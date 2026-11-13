"""Streaming gzip FASTQ reader.

Records are read one at a time so a whole batch never has to be held in
memory.  Any gzip-level failure (bad magic, truncated stream, failed CRC /
size check) is converted into :class:`GzipCorruptError`, and structural
FASTQ problems into :class:`FastqFormatError`.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass

from .errors import FastqFormatError, GzipCorruptError


@dataclass(slots=True)
class FastqRecord:
    """One FASTQ record.

    ``raw_lines`` keeps the exact four lines (including newline characters)
    as read from disk, so output files preserve the original bytes.
    """

    name: str
    sequence: str
    plus: str
    quality: str
    raw_lines: list[bytes]

    def raw_bytes(self) -> bytes:
        return b"".join(self.raw_lines)


def _normalise_name(header_line: bytes) -> str:
    """Return the read identifier used to align the four lanes.

    The header is split at the first whitespace, and a trailing ``/1`` or
    ``/2`` on the identifier is removed.
    """
    header = header_line.strip()
    token = header.split(b" ", 1)[0].split(b"\t", 1)[0]
    if token.endswith((b"/1", b"/2")):
        token = token[:-2]
    return token.decode("ascii", errors="replace")


def _readline(gz: gzip.GzipFile, label: str) -> bytes:
    try:
        return gz.readline()
    except (OSError, EOFError, gzip.BadGzipFile) as exc:  # pragma: no cover - defensive
        raise GzipCorruptError(f"{label}: gzip stream is corrupt ({exc})") from exc
    except Exception as exc:  # zlib.error subclasses Exception
        raise GzipCorruptError(f"{label}: gzip stream is corrupt ({exc})") from exc


class FastqReader:
    """Yield :class:`FastqRecord` objects from a gzipped FASTQ file."""

    def __init__(self, path: str, label: str) -> None:
        self.path = path
        self.label = label
        self._gz: gzip.GzipFile | None = None
        self.records_read = 0

    def __enter__(self) -> "FastqReader":
        try:
            # gzip.GzipFile reads the header lazily, so a missing/bad magic
            # surfaces on the first readline() rather than here.
            self._gz = gzip.open(self.path, "rb")
        except (OSError, gzip.BadGzipFile) as exc:
            raise GzipCorruptError(
                f"{self.label}: cannot open gzip stream ({exc})"
            ) from exc
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        assert self._gz is not None
        try:
            self._gz.close()
        except (OSError, EOFError) as close_exc:
            # Closing verifies the trailing CRC and ISIZE.
            if exc_type is None:
                raise GzipCorruptError(
                    f"{self.label}: gzip stream is corrupt ({close_exc})"
                ) from close_exc

    def __iter__(self) -> "FastqReader":
        return self

    def __next__(self) -> FastqRecord:
        assert self._gz is not None
        lines: list[bytes] = []
        for i in range(4):
            line = _readline(self._gz, self.label)
            if not line:
                if i == 0:
                    # Clean end of file.
                    raise StopIteration
                raise FastqFormatError(
                    f"{self.label}: record #{self.records_read + 1} ends after "
                    f"{i} line(s), expected 4 lines"
                )
            if not line.endswith(b"\n"):
                raise FastqFormatError(
                    f"{self.label}: record #{self.records_read + 1} line {i + 1} "
                    "is missing its newline (truncated file?)"
                )
            lines.append(line)

        header, seq, plus, qual = lines
        if not header.startswith(b"@"):
            raise FastqFormatError(
                f"{self.label}: record #{self.records_read + 1} header must start "
                f"with '@', got {header[:20]!r}"
            )
        if not plus.startswith(b"+"):
            raise FastqFormatError(
                f"{self.label}: record #{self.records_read + 1} third line must "
                f"start with '+', got {plus[:20]!r}"
            )
        sequence = seq.strip()
        quality = qual.strip()
        if not sequence:
            raise FastqFormatError(
                f"{self.label}: record #{self.records_read + 1} has an empty sequence"
            )
        if len(sequence) != len(quality):
            raise FastqFormatError(
                f"{self.label}: record #{self.records_read + 1} sequence length "
                f"({len(sequence)}) does not match quality length ({len(quality)})"
            )

        self.records_read += 1
        return FastqRecord(
            name=_normalise_name(header),
            sequence=sequence.decode("ascii", errors="replace"),
            plus=plus.decode("ascii", errors="replace"),
            quality=quality.decode("ascii", errors="replace"),
            raw_lines=lines,
        )


def open_fastq(path: str, label: str) -> FastqReader:
    return FastqReader(path, label)
