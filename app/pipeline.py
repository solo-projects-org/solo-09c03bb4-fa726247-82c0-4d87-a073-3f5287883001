"""Streaming demultiplex pipeline.

The four gzip lanes are opened together and walked record by record.  For
each position the read names must agree (after whitespace / ``/1``,``/2``
normalisation), the I1/I2 index sequences decide the destination sample,
and the original R1/R2 bytes are appended to that sample's gzip outputs.
Nothing is buffered beyond the current record.
"""

from __future__ import annotations

import gzip
import os
import zipfile
from dataclasses import dataclass, field
from typing import Callable

from .demux import UNASSIGNED, assign_read, parse_sample_sheet
from .errors import BatchMismatchError, DemuxError
from .fastq_io import FastqReader, open_fastq

LANE_LABELS = ("R1", "R2", "I1", "I2")


@dataclass
class Progress:
    total_reads: int = 0
    assigned: dict[str, int] = field(default_factory=dict)
    unassigned: int = 0
    last_error: str | None = None


class _GzWriters:
    """One gzip writer per output group, opened up front."""

    def __init__(self, output_dir: str, groups: list[str]) -> None:
        self._paths: dict[str, tuple[str, str]] = {}
        self._r1: dict[str, gzip.GzipFile] = {}
        self._r2: dict[str, gzip.GzipFile] = {}
        for group in groups:
            r1_path = os.path.join(output_dir, f"{group}_R1.fastq.gz")
            r2_path = os.path.join(output_dir, f"{group}_R2.fastq.gz")
            self._paths[group] = (r1_path, r2_path)
            self._r1[group] = gzip.open(r1_path, "wb", compresslevel=6)
            self._r2[group] = gzip.open(r2_path, "wb", compresslevel=6)

    def write(self, group: str, r1: bytes, r2: bytes) -> None:
        self._r1[group].write(r1)
        self._r2[group].write(r2)

    def close(self) -> None:
        for fh in (*self._r1.values(), *self._r2.values()):
            fh.close()


def run_demux(
    upload_dir: str,
    work_dir: str,
    progress: Progress,
    is_cancelled: Callable[[], bool],
) -> tuple[str, dict[str, int], int]:
    """Run the full pipeline.

    Returns ``(zip_path, counts_by_sample, unassigned_count)`` on success.
    Raises :class:`DemuxError` on any validation / IO failure.
    """
    samples = parse_sample_sheet(os.path.join(upload_dir, "samples.csv"))
    groups = [s.sample_id for s in samples] + [UNASSIGNED]
    progress.assigned = {s.sample_id: 0 for s in samples}

    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    writers = _GzWriters(output_dir, groups)

    lane_paths = {
        lane: os.path.join(upload_dir, f"{lane}.fastq.gz") for lane in LANE_LABELS
    }
    readers: dict[str, FastqReader] = {
        lane: open_fastq(path, lane) for lane, path in lane_paths.items()
    }

    record_index = 0
    try:
        with readers["R1"] as r1_reader, readers["R2"] as r2_reader, \
                readers["I1"] as i1_reader, readers["I2"] as i2_reader:
            while True:
                if is_cancelled():
                    raise DemuxError("job cancelled by user")

                r1 = next(r1_reader, None)
                r2 = next(r2_reader, None)
                i1 = next(i1_reader, None)
                i2 = next(i2_reader, None)
                ended = [rec is None for rec in (r1, r2, i1, i2)]
                if any(ended):
                    if not all(ended):
                        counts = {
                            lane: reader.records_read
                            for lane, reader in readers.items()
                        }
                        raise BatchMismatchError(
                            "lanes contain different numbers of records: "
                            + ", ".join(f"{lane}={n}" for lane, n in counts.items())
                        )
                    break  # all four ended together: clean EOF

                record_index += 1

                # Read order / identity alignment across the four lanes.
                names = {"R1": r1.name, "R2": r2.name, "I1": i1.name, "I2": i2.name}
                unique_names = set(names.values())
                if len(unique_names) != 1:
                    raise BatchMismatchError(
                        f"record #{record_index}: read names are not aligned "
                        f"across lanes: {names}"
                    )

                assignment = assign_read(i1.sequence, i2.sequence, samples)
                group = assignment.sample_id or UNASSIGNED
                if assignment.sample_id is None:
                    progress.unassigned += 1
                else:
                    progress.assigned[assignment.sample_id] += 1

                writers.write(group, r1.raw_bytes(), r2.raw_bytes())
                progress.total_reads = record_index
    finally:
        writers.close()

    if is_cancelled():
        raise DemuxError("job cancelled by user")

    counts = dict(progress.assigned)
    counts[UNASSIGNED] = progress.unassigned

    zip_path = os.path.join(work_dir, "result.zip")
    _build_zip(zip_path, output_dir, groups, counts)
    return zip_path, progress.assigned, progress.unassigned


def _build_zip(
    zip_path: str, output_dir: str, groups: list[str], counts: dict[str, int]
) -> None:
    counts_path = os.path.join(output_dir, "assignment_counts.csv")
    with open(counts_path, "w", encoding="utf-8", newline="") as fh:
        fh.write("sample_id,read_count\n")
        for group in groups:
            fh.write(f"{group},{counts.get(group, 0)}\n")

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_STORED) as zf:
        for group in groups:
            for suffix in ("R1", "R2"):
                name = f"{group}_{suffix}.fastq.gz"
                zf.write(os.path.join(output_dir, name), arcname=name)
        zf.write(counts_path, arcname="assignment_counts.csv")
