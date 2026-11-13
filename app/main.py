"""FastAPI application: upload, monitor, cancel and download demux jobs."""

from __future__ import annotations

import os
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates

from .errors import DemuxError
from .jobs import (
    MAX_PARALLEL_JOBS,
    JobManager,
)

BASE_DIR = Path(__file__).resolve().parent
JOBS_ROOT = Path(os.environ.get("DEMUX_JOBS_ROOT", "/data/jobs"))
LANES = ("R1", "R2", "I1", "I2")
SHEET_FIELD = "sample_sheet"

@asynccontextmanager
async def lifespan(_app: FastAPI):
    # Task state lives in memory only; start every run with a clean disk.
    manager.wipe()
    yield


app = FastAPI(title="FASTQ Demultiplexer", version="1.0.0", lifespan=lifespan)
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))
manager = JobManager(JOBS_ROOT)


# ---------------------------------------------------------------- helpers


async def _stream_to_file(upload: UploadFile, dest: Path) -> None:
    """Save an upload chunk by chunk so whole files are never buffered."""
    try:
        with open(dest, "wb") as fh:
            while True:
                chunk = await upload.read(1024 * 1024)
                if not chunk:
                    break
                fh.write(chunk)
    except OSError as exc:
        raise HTTPException(status_code=500, detail=f"failed to store upload: {exc}")
    finally:
        await upload.close()


# ----------------------------------------------------------------- routes


@app.get("/", response_class=HTMLResponse)
async def index(request: Request) -> HTMLResponse:
    return templates.TemplateResponse(
        request,
        "index.html",
        {"max_parallel": MAX_PARALLEL_JOBS},
    )


@app.get("/api/health")
async def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/api/jobs", status_code=201)
async def create_job(
    R1: UploadFile = File(...),
    R2: UploadFile = File(...),
    I1: UploadFile = File(...),
    I2: UploadFile = File(...),
    sample_sheet: UploadFile = File(...),
) -> JSONResponse:
    job = manager.try_create()
    if job is None:
        raise HTTPException(
            status_code=409,
            detail=(
                f"already running the maximum of {MAX_PARALLEL_JOBS} batches; "
                "wait for one to finish or cancel it and retry"
            ),
        )

    uploads = {
        "R1": (R1, "R1.fastq.gz"),
        "R2": (R2, "R2.fastq.gz"),
        "I1": (I1, "I1.fastq.gz"),
        "I2": (I2, "I2.fastq.gz"),
        SHEET_FIELD: (sample_sheet, "samples.csv"),
    }
    upload_names: dict[str, str] = {}
    try:
        for _key, (upload, target_name) in uploads.items():
            if not upload.filename:
                raise HTTPException(status_code=400, detail=f"missing filename for {target_name}")
            upload_names[target_name] = upload.filename
            await _stream_to_file(upload, job.upload_dir / target_name)
    except HTTPException:
        manager.discard(job)
        raise
    except Exception as exc:  # noqa: BLE001
        manager.discard(job)
        raise HTTPException(status_code=400, detail=f"upload failed: {exc}") from exc

    manager.launch(job, upload_names)
    return JSONResponse(status_code=201, content=job.snapshot())


@app.get("/api/jobs")
async def list_jobs() -> dict[str, list[dict]]:
    return {"jobs": manager.list_jobs()}


@app.get("/api/jobs/{job_id}")
async def get_job(job_id: str) -> dict:
    job = manager.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found (state is memory-only)")
    return job.snapshot()


@app.post("/api/jobs/{job_id}/cancel")
async def cancel_job(job_id: str) -> dict:
    job = manager.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found (state is memory-only)")
    try:
        job = manager.cancel(job_id)
    except TimeoutError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    assert job is not None
    return job.snapshot()


@app.get("/api/jobs/{job_id}/download")
async def download_job(job_id: str) -> FileResponse:
    job = manager.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found (state is memory-only)")
    if job.status != "completed" or job.zip_path is None or not job.zip_path.exists():
        raise HTTPException(
            status_code=404,
            detail=(
                "no download available: only completed jobs publish a ZIP "
                f"(current status: {job.status})"
            ),
        )
    return FileResponse(
        job.zip_path,
        media_type="application/zip",
        filename=f"demux_{job_id}.zip",
    )


@app.exception_handler(DemuxError)
async def demux_error_handler(_request: Request, exc: DemuxError) -> JSONResponse:
    # Pipeline errors are recorded on the job; this is a defensive fallback.
    return JSONResponse(status_code=400, content={"detail": str(exc)})
