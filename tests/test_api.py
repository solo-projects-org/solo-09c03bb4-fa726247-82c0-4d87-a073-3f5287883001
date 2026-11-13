"""End-to-end tests through the FastAPI app and on-disk pipeline."""

from __future__ import annotations

import gzip
import io
import os
import time
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.fastq_io import open_fastq
from app.jobs import MAX_PARALLEL_JOBS
from app.main import app
from app.pipeline import Progress, run_demux

SAMPLES_DIR = Path(__file__).resolve().parent.parent / "samples"


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    jobs_root = tmp_path_factory.mktemp("jobs")
    os.environ["DEMUX_JOBS_ROOT"] = str(jobs_root)
    # Reload the module-level manager path.
    from app import main as main_module

    main_module.JOBS_ROOT = jobs_root
    main_module.manager = main_module.JobManager(jobs_root)
    main_module.manager.wipe()
    with TestClient(app) as c:
        yield c


def _files_for(dataset: str) -> dict[str, bytes]:
    base = SAMPLES_DIR / dataset
    return {
        "R1": (base / "R1.fastq.gz").read_bytes(),
        "R2": (base / "R2.fastq.gz").read_bytes(),
        "I1": (base / "I1.fastq.gz").read_bytes(),
        "I2": (base / "I2.fastq.gz").read_bytes(),
        "sample_sheet": (base / "samples.csv").read_bytes(),
    }


def _post_job(client: TestClient, dataset: str):
    files = _files_for(dataset)
    files_param = []
    for key in ("R1", "R2", "I1", "I2", "sample_sheet"):
        filename = "samples.csv" if key == "sample_sheet" else f"{key}.fastq.gz"
        files_param.append((key, (filename, files[key], "application/gzip")))
    return client.post("/api/jobs", files=files_param)


def _wait_terminal(client: TestClient, job_id: str, timeout: float = 20.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        snap = client.get(f"/api/jobs/{job_id}").json()
        if snap["status"] in ("completed", "failed", "cancelled"):
            return snap
        time.sleep(0.05)
    raise AssertionError("job did not finish in time")


# ------------------------------------------------------------------- happy


def test_normal_batch_completes_and_downloads_zip(client):
    resp = _post_job(client, "normal")
    assert resp.status_code == 201, resp.text
    job_id = resp.json()["id"]
    snap = _wait_terminal(client, job_id)
    assert snap["status"] == "completed", snap
    assert snap["total_reads"] == 6
    assert snap["final_counts"]["undetermined"] == 0
    assert sum(snap["final_counts"].values()) == 6

    dl = client.get(f"/api/jobs/{job_id}/download")
    assert dl.status_code == 200
    zf = zipfile.ZipFile(io.BytesIO(dl.content))
    names = set(zf.namelist())
    assert "assignment_counts.csv" in names
    for sample in ("sampleA", "sampleB", "sampleC", "undetermined"):
        assert f"{sample}_R1.fastq.gz" in names
        assert f"{sample}_R2.fastq.gz" in names

    counts = dict(
        line.split(",")
        for line in zf.read("assignment_counts.csv").decode().strip().splitlines()[1:]
    )
    # normal dataset: A gets 3, B gets 2, C gets 1
    assert counts == {
        "sampleA": "3",
        "sampleB": "2",
        "sampleC": "1",
        "undetermined": "0",
    }

    # R1/R2 original bytes and order are preserved inside a sample group.
    with gzip.open(
        io.BytesIO(zf.read("sampleA_R1.fastq.gz")), "rb"
    ) as fh:
        r1_text = fh.read().decode()
    names_in_r1 = [
        line.split(" ")[0].lstrip("@")
        for line in r1_text.splitlines()
        if line.startswith("@")
    ]
    assert names_in_r1 == ["SEQ-0001/1", "SEQ-0002/1", "SEQ-0006/1"]


def test_ambiguous_batch_routes_tie_and_no_candidate_to_undetermined(client):
    resp = _post_job(client, "ambiguous")
    assert resp.status_code == 201
    snap = _wait_terminal(client, resp.json()["id"])
    assert snap["status"] == "completed", snap
    assert snap["total_reads"] == 5
    # 3 unique samples + no-candidate + tie
    assert snap["final_counts"]["undetermined"] == 2
    assert snap["final_counts"]["sampleA"] == 1
    assert snap["final_counts"]["sampleB"] == 1
    assert snap["final_counts"]["sampleC"] == 1


def test_corrupt_gzip_fails_whole_batch_and_blocks_download(client):
    resp = _post_job(client, "corrupt")
    assert resp.status_code == 201
    snap = _wait_terminal(client, resp.json()["id"])
    assert snap["status"] == "failed"
    assert "corrupt" in snap["error"].lower() or "gzip" in snap["error"].lower()
    dl = client.get(f"/api/jobs/{resp.json()['id']}/download")
    assert dl.status_code == 404


# ------------------------------------------------------------- validation


def _make_fastq(records: list[tuple[str, str]]) -> bytes:
    buf = io.BytesIO()
    with gzip.open(buf, "wb") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode())
    return buf.getvalue()


def _post_raw(client, files: dict[str, bytes], sheet: bytes):
    payload = []
    for key in ("R1", "R2", "I1", "I2"):
        payload.append((key, (f"{key}.fastq.gz", files[key], "application/gzip")))
    payload.append(("sample_sheet", ("samples.csv", sheet, "text/csv")))
    return client.post("/api/jobs", files=payload)


def _sheet() -> bytes:
    return b"sample_id,i7,i5\nA,ACGT,ACGT\nB,TGCA,TGCA\n"


def test_mismatched_record_counts_fail_batch(client):
    good = _make_fastq([("r1", "ACGT"), ("r2", "ACGT")])
    short = _make_fastq([("r1", "ACGT")])
    resp = _post_raw(
        client,
        {"R1": good, "R2": good, "I1": short, "I2": good},
        _sheet(),
    )
    snap = _wait_terminal(client, resp.json()["id"])
    assert snap["status"] == "failed"
    assert "different numbers of records" in snap["error"]


def test_misordered_read_names_fail_batch(client):
    r1 = _make_fastq([("r1", "ACGT"), ("r2", "ACGT")])
    r2 = _make_fastq([("r2", "ACGT"), ("r1", "ACGT")])
    resp = _post_raw(
        client, {"R1": r1, "R2": r2, "I1": r1, "I2": r1}, _sheet()
    )
    snap = _wait_terminal(client, resp.json()["id"])
    assert snap["status"] == "failed"
    assert "not aligned" in snap["error"]


def test_sequence_quality_length_mismatch_fails_batch(client):
    buf = io.BytesIO()
    with gzip.open(buf, "wb") as fh:
        fh.write(b"@r1\nACGT\n+\nIII\n")  # qual one char short
    bad = buf.getvalue()
    good = _make_fastq([("r1", "ACGT")])
    resp = _post_raw(
        client, {"R1": bad, "R2": good, "I1": good, "I2": good}, _sheet()
    )
    snap = _wait_terminal(client, resp.json()["id"])
    assert snap["status"] == "failed"
    assert "sequence length" in snap["error"]


def test_garbage_gzip_fails_batch(client):
    good = _make_fastq([("r1", "ACGT")])
    garbage = b"this is not gzip at all"
    resp = _post_raw(
        client, {"R1": good, "R2": garbage, "I1": good, "I2": good}, _sheet()
    )
    snap = _wait_terminal(client, resp.json()["id"])
    assert snap["status"] == "failed"
    assert "gzip" in snap["error"].lower()


def test_bad_sample_sheet_fails_batch(client):
    good = _make_fastq([("r1", "ACGT")])
    resp = _post_raw(
        client,
        {"R1": good, "R2": good, "I1": good, "I2": good},
        b"id,i7,i5\nA,ACGT,ACGT\n",
    )
    snap = _wait_terminal(client, resp.json()["id"])
    assert snap["status"] == "failed"
    assert "header" in snap["error"]


# ------------------------------------------------------------ concurrency


def test_parallel_cap_and_cancel_frees_slot(client, monkeypatch):
    import threading

    import app.jobs as jobs_module
    from app import main as main_module

    gate = threading.Event()
    entered = threading.Event()
    real_run = jobs_module.run_demux

    def slow_run(upload_dir, work_dir, progress, is_cancelled):
        entered.set()
        # Honour cancellation promptly; otherwise park until the gate opens
        # so the two slots stay occupied for the cap assertion.
        while not gate.wait(timeout=0.02):
            if is_cancelled():
                from app.errors import DemuxError

                raise DemuxError("job cancelled by user")
        return real_run(upload_dir, work_dir, progress, is_cancelled)

    monkeypatch.setattr(jobs_module, "run_demux", slow_run)
    # jobs._run resolves run_demux as a module global, so patching
    # app.jobs is sufficient; main.py never imports the name directly.
    try:
        first = _post_job(client, "normal")
        assert entered.wait(timeout=5)
        entered.clear()
        second = _post_job(client, "normal")
        assert entered.wait(timeout=5)
        assert first.status_code == 201 and second.status_code == 201

        # Both slots occupied -> a third batch is rejected, not queued.
        third = _post_job(client, "normal")
        assert third.status_code == 409
        assert "maximum" in third.json()["detail"]

        # Cancel the first batch: the endpoint returns only after the worker
        # has actually stopped and released its slot.
        cancel_resp = client.post(f"/api/jobs/{first.json()['id']}/cancel")
        assert cancel_resp.status_code == 200
        snap = client.get(f"/api/jobs/{first.json()['id']}").json()
        assert snap["status"] == "cancelled"
        assert client.get(
            f"/api/jobs/{first.json()['id']}/download"
        ).status_code == 404

        # A new batch can take the freed slot immediately, while the second
        # batch continues unaffected.
        replacement = _post_job(client, "normal")
        assert replacement.status_code == 201

        gate.set()
        assert (
            _wait_terminal(client, second.json()["id"])["status"] == "completed"
        )
        assert (
            _wait_terminal(client, replacement.json()["id"])["status"] == "completed"
        )
    finally:
        monkeypatch.undo()


def test_memory_only_state_404_for_unknown_job(client):
    assert client.get("/api/jobs/deadbeef").status_code == 404
    assert client.post("/api/jobs/deadbeef/cancel").status_code == 404


def test_index_page_renders(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert "text/html" in resp.headers["content-type"]
    assert "jobForm" in resp.text
    assert "FASTQ" in resp.text


# --------------------------------------------------------------- pipeline


def test_pipeline_preserves_original_bytes_and_order(tmp_path):
    upload = tmp_path / "up"
    work = tmp_path / "work"
    upload.mkdir()
    records = [("r1", "ACGTACGT"), ("r2", "TTTTCCCC"), ("r3", "GGGGAAAA")]
    for lane in ("R1", "R2", "I1", "I2"):
        (upload / f"{lane}.fastq.gz").write_bytes(
            _make_fastq([(n, s) for n, s in records])
        )
    (upload / "samples.csv").write_text(
        "sample_id,i7,i5\nA,ACGTACGT,ACGTACGT\nB,TTTTCCCC,TTTTCCCC\n",
        encoding="utf-8",
    )
    zip_path, assigned, unassigned = run_demux(
        str(upload), str(work), Progress(), lambda: False
    )
    assert unassigned == 1
    assert assigned == {"A": 1, "B": 1}
    with zipfile.ZipFile(zip_path) as zf:
        with gzip.open(io.BytesIO(zf.read("A_R1.fastq.gz")), "rb") as fh:
            assert fh.read().startswith(b"@r1\nACGTACGT\n+\nIIIIIIII\n")


def test_reader_detects_truncated_record(tmp_path):
    buf = io.BytesIO()
    with gzip.open(buf, "wb") as fh:
        fh.write(b"@r1\nACGT\n+\nIIII\n@r2\nAC\n")  # second record cut short
    path = tmp_path / "x.fastq.gz"
    path.write_bytes(buf.getvalue())
    reader = open_fastq(str(path), "X")
    with pytest.raises(Exception):
        with reader:
            list(reader)
