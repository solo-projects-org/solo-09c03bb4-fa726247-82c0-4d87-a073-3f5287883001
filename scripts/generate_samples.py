"""Generate the bundled demo datasets for the web UI.

Three datasets are produced, each with the four lanes
(R1/R2/I1/I2 .fastq.gz) and a samples.csv:

* normal    - every read matches a unique sample (perfect or 1 mismatch)
* ambiguous - one read with no candidate and one read whose best total
              distance is shared by two samples (a tie); both go to
              undetermined
* corrupt   - same content as normal but the R2 gzip stream is truncated
              at its tail, so the CRC / ISIZE check fails and the whole
              batch fails.

The sample sheet is chosen so that sampleA and sampleB differ at exactly
two positions on both index ends; flipping one of those positions yields a
read that sits at total distance 2 from each of them -> a genuine tie.

Run:  python -I scripts/generate_samples.py [output_root]
"""

from __future__ import annotations

import gzip
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = os.path.join(HERE, "..", "samples")

SAMPLE_SHEET = [
    ("sampleA", "AACGTGAT", "AACGTGAT"),
    ("sampleB", "ATCGTGAA", "ATCGTGAA"),  # 2 flips vs A on each end
    ("sampleC", "GCTACGAA", "TTACGCAG"),
]


def _record(name: str, seq: str) -> bytes:
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode("ascii")


def _write_lane(path: str, records: list[bytes]) -> None:
    with gzip.open(path, "wb", compresslevel=6) as fh:
        for rec in records:
            fh.write(rec)


def _read_names(n: int) -> list[str]:
    # Whitespace suffixes exercise the "split at first whitespace" rule.
    return [f"SEQ-{i:04d} run=L1" for i in range(1, n + 1)]


def _build_lanes(
    index_pairs: list[tuple[str, str]], names: list[str], out_dir: str
) -> None:
    r1, r2, i1, i2 = [], [], [], []
    for idx, (name, (i7, i5)) in enumerate(zip(names, index_pairs), start=1):
        base = name.split(" ")[0], name.split(" ", 1)[1]
        # R1/R2 carry /1 and /2 suffixes to confirm normalisation.
        r1.append(_record(f"{base[0]}/1 {base[1]}", f"R1SEQ{idx:04d}AAAA"))
        r2.append(_record(f"{base[0]}/2 {base[1]}", f"R2SEQ{idx:04d}TTTT"))
        i1.append(_record(name, i7))
        i2.append(_record(name, i5))
    _write_lane(os.path.join(out_dir, "R1.fastq.gz"), r1)
    _write_lane(os.path.join(out_dir, "R2.fastq.gz"), r2)
    _write_lane(os.path.join(out_dir, "I1.fastq.gz"), i1)
    _write_lane(os.path.join(out_dir, "I2.fastq.gz"), i2)


def _write_sheet(out_dir: str) -> None:
    with open(os.path.join(out_dir, "samples.csv"), "w", encoding="utf-8") as fh:
        fh.write("sample_id,i7,i5\n")
        for row in SAMPLE_SHEET:
            fh.write(",".join(row) + "\n")


def normal_dataset(out_dir: str) -> list[tuple[str, str]]:
    a7, a5 = SAMPLE_SHEET[0][1], SAMPLE_SHEET[0][2]
    b7, b5 = SAMPLE_SHEET[1][1], SAMPLE_SHEET[1][2]
    c7, c5 = SAMPLE_SHEET[2][1], SAMPLE_SHEET[2][2]
    pairs = [
        (a7, a5),                       # perfect -> sampleA
        ("AACGTGTT", a5),               # 1 mismatch on i7 -> sampleA
        (b7, "ATCGTGAC"),               # 1 mismatch on i5 -> sampleB
        ("ATCGTGAG", "ATCGTGAG"),       # 1+1 mismatch -> sampleB (total 2, unique)
        (c7, c5),                       # perfect -> sampleC
        (a7, a5),                       # perfect -> sampleA
    ]
    _build_lanes(pairs, _read_names(len(pairs)), out_dir)
    _write_sheet(out_dir)
    return pairs


def ambiguous_dataset(out_dir: str) -> None:
    a7, a5 = SAMPLE_SHEET[0][1], SAMPLE_SHEET[0][2]
    b7, b5 = SAMPLE_SHEET[1][1], SAMPLE_SHEET[1][2]
    # A and B differ at positions 1 and 7 on each end.  A read that keeps
    # A's base at position 1 and B's base at position 7 (AACGTGAA) is
    # distance 1 from each on each end -> both pass the <=1-per-end gate
    # and share total distance 2: a genuine tie -> undetermined.
    tie_pair_seq = "AACGTGAA"
    pairs = [
        (a7, a5),                       # unique -> sampleA
        (b7, b5),                       # unique -> sampleB
        (SAMPLE_SHEET[2][1], SAMPLE_SHEET[2][2]),  # unique -> sampleC
        ("GGGGGGGG", "CCCCCCCC"),       # >1 away from every sample -> no candidate
        (tie_pair_seq, tie_pair_seq),   # tied best (2 vs 2) -> undetermined
    ]
    _build_lanes(pairs, _read_names(len(pairs)), out_dir)
    _write_sheet(out_dir)


def corrupt_dataset(out_dir: str, normal_dir: str) -> None:
    # Copy normal, then cut the R2 gzip stream nine bytes short, removing
    # its 8-byte trailer (CRC32 + ISIZE) and one byte of deflate data.
    for name in ("R1.fastq.gz", "I1.fastq.gz", "I2.fastq.gz", "samples.csv"):
        shutil.copyfile(os.path.join(normal_dir, name), os.path.join(out_dir, name))
    src = os.path.join(normal_dir, "R2.fastq.gz")
    dst = os.path.join(out_dir, "R2.fastq.gz")
    with open(src, "rb") as fh:
        data = fh.read()
    with open(dst, "wb") as fh:
        fh.write(data[:-9])


def main(root: str) -> None:
    if os.path.isdir(root):
        shutil.rmtree(root)
    normal_dir = os.path.join(root, "normal")
    ambiguous_dir = os.path.join(root, "ambiguous")
    corrupt_dir = os.path.join(root, "corrupt")
    for d in (normal_dir, ambiguous_dir, corrupt_dir):
        os.makedirs(d, exist_ok=True)
    normal_dataset(normal_dir)
    ambiguous_dataset(ambiguous_dir)
    corrupt_dataset(corrupt_dir, normal_dir)
    print(f"sample datasets written under {root}")


if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ROOT
    main(os.path.abspath(root))
