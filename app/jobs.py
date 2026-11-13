"""In-memory job registry with a hard cap of two running batches."""

from __future__ import annotations

import threading
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .errors import DemuxError
from .pipeline import Progress, run_demux

# Lifecycle states.  Everything is held in memory only.
RUNNING = "running"
COMPLETED = "completed"
FAILED = "failed"
CANCELLED = "cancelled"
TERMINAL_STATES = frozenset((COMPLETED, FAILED, CANCELLED))

MAX_PARALLEL_JOBS = 2
# Safety bound for the cancel endpoint: a record cycle is normally
# milliseconds, so this is far beyond any realistic shutdown.
CANCEL_WAIT_TIMEOUT = 300.0


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


@dataclass
class Job:
    id: str
    upload_dir: Path
    work_dir: Path
    status: str = RUNNING
    error: str | None = None
    progress: Progress = field(default_factory=Progress)
    final_counts: dict[str, int] | None = None
    unassigned_count: int = 0
    zip_path: Path | None = None
    created_at: datetime = field(default_factory=_utc_now)
    finished_at: datetime | None = None
    cancel_event: threading.Event = field(default_factory=threading.Event)
    terminated: threading.Event = field(default_factory=threading.Event)
    thread: threading.Thread | None = None
    upload_names: dict[str, str] = field(default_factory=dict)

    def snapshot(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "status": self.status,
            "error": self.error,
            "total_reads": self.progress.total_reads,
            "assigned": dict(self.progress.assigned),
            "unassigned": self.progress.unassigned,
            "final_counts": dict(self.final_counts) if self.final_counts else None,
            "created_at": self.created_at.isoformat(),
            "finished_at": self.finished_at.isoformat() if self.finished_at else None,
            "upload_names": dict(self.upload_names),
            "download_available": self.status == COMPLETED and self.zip_path is not None
            and self.zip_path.exists(),
        }


class JobManager:
    def __init__(self, jobs_root: Path) -> None:
        self.jobs_root = jobs_root
        self._semaphore = threading.BoundedSemaphore(MAX_PARALLEL_JOBS)
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()

    # ---- lifecycle -----------------------------------------------------

    def acquire_slot_or_none(self) -> bool:
        """Non-blocking attempt to occupy one of the two batch slots."""
        return self._semaphore.acquire(blocking=False)

    def create_job(self) -> Job:
        """Allocate a job.  Caller must hold a slot (see try_create)."""
        job_id = uuid.uuid4().hex
        base = self.jobs_root / job_id
        upload_dir = base / "uploads"
        work_dir = base / "work"
        upload_dir.mkdir(parents=True, exist_ok=False)
        work_dir.mkdir(parents=True, exist_ok=False)
        job = Job(id=job_id, upload_dir=upload_dir, work_dir=work_dir)
        with self._lock:
            self._jobs[job_id] = job
        return job

    def try_create(self) -> Job | None:
        """Atomically take a slot and create a job, or reject the batch."""
        if not self.acquire_slot_or_none():
            return None
        try:
            return self.create_job()
        except Exception:
            self._semaphore.release()
            raise

    def discard(self, job: Job) -> None:
        """Release resources for a job whose uploads never started."""
        with self._lock:
            self._jobs.pop(job.id, None)
        self._remove_tree(job.work_dir.parent)
        self._semaphore.release()

    def launch(self, job: Job, upload_names: dict[str, str]) -> None:
        job.upload_names = upload_names
        thread = threading.Thread(
            target=self._run, args=(job,), name=f"demux-{job.id}", daemon=True
        )
        job.thread = thread
        thread.start()

    def _run(self, job: Job) -> None:
        try:
            zip_path, _assigned, unassigned = run_demux(
                str(job.upload_dir),
                str(job.work_dir),
                job.progress,
                job.cancel_event.is_set,
            )
            job.zip_path = Path(zip_path)
            job.unassigned_count = unassigned
            job.final_counts = dict(job.progress.assigned)
            job.final_counts["undetermined"] = unassigned
            job.status = COMPLETED
        except Exception as exc:  # noqa: BLE001 - every failure becomes a job state
            message = str(exc) if isinstance(exc, DemuxError) else f"unexpected error: {exc}"
            job.error = message
            if job.cancel_event.is_set():
                job.status = CANCELLED
            else:
                job.status = FAILED
            # Failed / cancelled batches offer no download and leave no
            # files behind (neither outputs nor their uploads).
            self._remove_tree(job.work_dir.parent)
        finally:
            job.finished_at = _utc_now()
            # Release the parallelism slot *before* signalling termination,
            # so a cancelled batch's slot is free the moment cancel returns.
            self._semaphore.release()
            job.terminated.set()

    def cancel(self, job_id: str) -> Job | None:
        """Request cancellation and wait until the worker has really stopped."""
        with self._lock:
            job = self._jobs.get(job_id)
        if job is None:
            return None
        if job.status in TERMINAL_STATES:
            return job
        job.cancel_event.set()
        if not job.terminated.wait(timeout=CANCEL_WAIT_TIMEOUT):
            raise TimeoutError(f"job {job_id} did not stop within {CANCEL_WAIT_TIMEOUT}s")
        return job

    # ---- queries -------------------------------------------------------

    def get(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(job_id)

    def list_jobs(self) -> list[dict[str, Any]]:
        with self._lock:
            jobs = list(self._jobs.values())
        jobs.sort(key=lambda j: j.created_at, reverse=True)
        return [j.snapshot() for j in jobs]

    # --- helpers --------------------------------------------------------

    @staticmethod
    def _remove_tree(path: Path) -> None:
        if not path.exists():
            return
        import shutil

        shutil.rmtree(path, ignore_errors=True)

    def wipe(self) -> None:
        """Remove all on-disk job data at startup (state is memory-only)."""
        self._remove_tree(self.jobs_root)
        self.jobs_root.mkdir(parents=True, exist_ok=True)
