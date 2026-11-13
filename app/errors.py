"""Domain-level errors used across the demultiplex pipeline."""


class DemuxError(Exception):
    """Base class for expected, user-facing pipeline failures."""


class GzipCorruptError(DemuxError):
    """Raised when a gzip stream cannot be decoded or fails its checksum."""


class FastqFormatError(DemuxError):
    """Raised when a FASTQ record is structurally invalid."""


class SampleSheetError(DemuxError):
    """Raised when the uploaded CSV sample sheet is invalid."""


class BatchMismatchError(DemuxError):
    """Raised when the four lanes disagree on record count or read order."""
