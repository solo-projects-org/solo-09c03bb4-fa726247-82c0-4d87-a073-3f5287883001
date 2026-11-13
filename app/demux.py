"""Sample-sheet parsing and index-based read assignment."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass

from .errors import SampleSheetError

UNASSIGNED = "undetermined"
_ALLOWED_INDEX_CHARS = set("ACGTNacgtn")
# sample_id is used to build output file and ZIP entry names, so restrict
# it to a safe character set and forbid any path-separator semantics.
_SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,62}$")


@dataclass(slots=True)
class Sample:
    sample_id: str
    i7: str
    i5: str


def parse_sample_sheet(csv_path: str) -> list[Sample]:
    """Parse and validate the three-column sample sheet.

    Required header (order matters): ``sample_id,i7,i5``.
    """
    samples: list[Sample] = []
    seen_ids: set[str] = set()
    seen_pairs: set[tuple[str, str]] = set()
    index_len: int | None = None

    try:
        with open(csv_path, "r", encoding="utf-8-sig", newline="") as fh:
            reader = csv.reader(fh)
            try:
                header = next(reader)
            except StopIteration:
                raise SampleSheetError("sample sheet CSV is empty")
            if [c.strip().lower() for c in header] != ["sample_id", "i7", "i5"]:
                raise SampleSheetError(
                    "sample sheet header must be exactly: sample_id,i7,i5 "
                    f"(got {header!r})"
                )

            for row_no, row in enumerate(reader, start=2):
                if not row or all(not c.strip() for c in row):
                    continue  # tolerate blank trailing lines
                if len(row) != 3:
                    raise SampleSheetError(
                        f"row {row_no}: expected 3 columns, got {len(row)}"
                    )
                sample_id, i7, i5 = (c.strip() for c in row)
                if not sample_id:
                    raise SampleSheetError(f"row {row_no}: empty sample_id")
                if not _SAFE_SAMPLE_ID.fullmatch(sample_id) or ".." in sample_id:
                    raise SampleSheetError(
                        f"row {row_no}: unsafe sample_id {sample_id!r}; use 1-63 "
                        "letters/digits/._- with no path separators or '..'"
                    )
                if sample_id == UNASSIGNED:
                    raise SampleSheetError(
                        f"row {row_no}: sample_id {UNASSIGNED!r} is reserved"
                    )
                if sample_id in seen_ids:
                    raise SampleSheetError(f"row {row_no}: duplicate sample_id {sample_id!r}")
                if not i7 or not i5:
                    raise SampleSheetError(f"row {row_no}: empty index for {sample_id!r}")
                if not set(i7) <= _ALLOWED_INDEX_CHARS or not set(i5) <= _ALLOWED_INDEX_CHARS:
                    raise SampleSheetError(
                        f"row {row_no}: indices may only contain A/C/G/T/N, got "
                        f"i7={i7!r}, i5={i5!r}"
                    )
                i7, i5 = i7.upper(), i5.upper()
                if index_len is None:
                    index_len = len(i7)
                    if len(i5) != index_len:
                        raise SampleSheetError(
                            f"row {row_no}: i7 length {len(i7)} differs from i5 "
                            f"length {len(i5)}; all indices must be equal length"
                        )
                elif len(i7) != index_len or len(i5) != index_len:
                    raise SampleSheetError(
                        f"row {row_no}: index length {max(len(i7), len(i5))} differs "
                        f"from the sheet's index length {index_len}; all indices "
                        "must be equal length"
                    )
                pair = (i7, i5)
                if pair in seen_pairs:
                    raise SampleSheetError(
                        f"row {row_no}: duplicate i7+i5 combination {pair}"
                    )
                seen_ids.add(sample_id)
                seen_pairs.add(pair)
                samples.append(Sample(sample_id=sample_id, i7=i7, i5=i5))
    except UnicodeDecodeError as exc:
        raise SampleSheetError(f"sample sheet is not valid UTF-8 ({exc})") from exc
    except OSError as exc:
        raise SampleSheetError(f"cannot read sample sheet ({exc})") from exc

    if not samples:
        raise SampleSheetError("sample sheet contains no sample rows")
    return samples


def hamming_distance(a: str, b: str) -> int:
    """Position-wise mismatch count; different lengths cannot match."""
    if len(a) != len(b):
        return len(a) + len(b) + 1
    return sum(1 for x, y in zip(a, b) if x != y)


@dataclass(slots=True)
class Assignment:
    sample_id: str | None  # None => undetermined
    total_distance: int
    reason: str  # "unique" | "no_candidate" | "tie" | "length_mismatch"


def assign_read(
    i1_seq: str, i2_seq: str, samples: list[Sample]
) -> Assignment:
    """Assign one read pair to a sample.

    Each index end may differ from the sample's index by at most one
    character.  Among all candidates the sample with the smallest combined
    distance wins, but only when that best distance is held by exactly one
    sample.  No candidate or a tie sends the read to the undetermined group.
    """
    candidates: list[tuple[int, Sample]] = []
    length_ok = True
    for sample in samples:
        d7 = hamming_distance(i1_seq, sample.i7)
        d5 = hamming_distance(i2_seq, sample.i5)
        if d7 == len(i1_seq) + len(sample.i7) + 1 or d5 == len(i2_seq) + len(sample.i5) + 1:
            length_ok = False
            continue
        if d7 <= 1 and d5 <= 1:
            candidates.append((d7 + d5, sample))

    if not candidates:
        reason = "length_mismatch" if not length_ok else "no_candidate"
        return Assignment(sample_id=None, total_distance=-1, reason=reason)

    best = min(candidates, key=lambda item: item[0])[0]
    winners = [s for dist, s in candidates if dist == best]
    if len(winners) > 1:
        return Assignment(
            sample_id=None,
            total_distance=best,
            reason="tie",
        )
    return Assignment(
        sample_id=winners[0].sample_id,
        total_distance=best,
        reason="unique",
    )
