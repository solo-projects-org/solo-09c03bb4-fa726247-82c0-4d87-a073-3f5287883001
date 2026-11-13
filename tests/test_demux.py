"""Tests for index matching and sample-sheet validation."""

import pytest

from app.demux import Sample, assign_read, hamming_distance, parse_sample_sheet
from app.errors import SampleSheetError

SAMPLES = [
    Sample("sampleA", "AACGTGAT", "AACGTGAT"),
    Sample("sampleB", "ATCGTGAA", "ATCGTGAA"),
    Sample("sampleC", "GCTACGAA", "TTACGCAG"),
]


def test_hamming():
    assert hamming_distance("ACGT", "ACGT") == 0
    assert hamming_distance("ACGT", "ACGA") == 1
    assert hamming_distance("AC", "ACGT") > 2  # length mismatch never matches


def test_perfect_assignment():
    res = assign_read("AACGTGAT", "AACGTGAT", SAMPLES)
    assert res.sample_id == "sampleA"
    assert res.total_distance == 0
    assert res.reason == "unique"


def test_one_mismatch_per_end_is_allowed():
    res = assign_read("AACGTGTT", "AACGTGAC", SAMPLES)
    assert res.sample_id == "sampleA"
    assert res.total_distance == 2


def test_two_mismatches_on_one_end_rejected():
    # Two flips on i7 alone disqualifies every sample on that end, while
    # i5 still matches sampleA exactly -> no candidate survives the gate.
    res = assign_read("AAAGTGAA", "AACGTGAT", SAMPLES)
    # i7 AAAGTGAA: vs A AACGTGAT -> pos2,3,7 = 3; vs B ATCGTGAA -> 2;
    # vs C GCTACGAA -> >=3.  Nothing is within 1 on i7, so no candidate.
    assert res.sample_id is None
    assert res.reason == "no_candidate"


def test_demo_sheet_tie_goes_to_undetermined():
    # Demo sheet: A=AACGTGAT, B=ATCGTGAA (flips at positions 1 and 7 on
    # each end).  Read AACGTGAA keeps A at position 1 and B at position 7,
    # so it is distance (1,1) -> total 2 from both A and B: a real tie.
    res = assign_read("AACGTGAA", "AACGTGAA", SAMPLES)
    assert res.sample_id is None
    assert res.reason == "tie"
    assert res.total_distance == 2
    # SampleC is far away and never enters the candidate set.
    assert all(
        hamming_distance("AACGTGAA", s.i7) > 1 or hamming_distance("AACGTGAA", s.i5) > 1
        for s in SAMPLES if s.sample_id == "sampleC"
    )


def test_explicit_tie_with_custom_sheet():
    # Two samples two flips apart on both ends; intermediate read is
    # distance (1,1) from each -> tied total distance 2.
    samples = [
        Sample("x", "AAAAAAAA", "AAAAAAAA"),
        Sample("y", "AATTAAAA", "AATTAAAA"),
    ]
    res = assign_read("AATAAAAA", "AATAAAAA", samples)
    # vs x: pos2 T (1) per end -> total 2
    # vs y: pos3 T->A (1) per end -> total 2
    assert res.sample_id is None
    assert res.reason == "tie"
    assert res.total_distance == 2


def test_length_mismatch_is_unassigned_not_fatal():
    res = assign_read("ACG", "AACGTGAT", SAMPLES)
    assert res.sample_id is None
    assert res.reason == "length_mismatch"


def test_sample_sheet_ok(tmp_path):
    csv_path = tmp_path / "s.csv"
    csv_path.write_text("sample_id,i7,i5\nA,ACGT,ACGT\nB,TGCA,TGCA\n", encoding="utf-8")
    samples = parse_sample_sheet(str(csv_path))
    assert [s.sample_id for s in samples] == ["A", "B"]


def test_sample_sheet_rejects_bad_header(tmp_path):
    csv_path = tmp_path / "s.csv"
    csv_path.write_text("id,i7,i5\nA,ACGT,ACGT\n", encoding="utf-8")
    with pytest.raises(SampleSheetError):
        parse_sample_sheet(str(csv_path))


def test_sample_sheet_rejects_duplicate_combo(tmp_path):
    csv_path = tmp_path / "s.csv"
    csv_path.write_text(
        "sample_id,i7,i5\nA,ACGT,ACGT\nB,ACGT,ACGT\n", encoding="utf-8"
    )
    with pytest.raises(SampleSheetError):
        parse_sample_sheet(str(csv_path))


def test_sample_sheet_rejects_unequal_lengths(tmp_path):
    csv_path = tmp_path / "s.csv"
    csv_path.write_text("sample_id,i7,i5\nA,ACGT,AC\nB,TTGG,TTGG\n", encoding="utf-8")
    with pytest.raises(SampleSheetError):
        parse_sample_sheet(str(csv_path))


def test_sample_sheet_rejects_bad_chars(tmp_path):
    csv_path = tmp_path / "s.csv"
    csv_path.write_text("sample_id,i7,i5\nA,ACGX,ACGT\n", encoding="utf-8")
    with pytest.raises(SampleSheetError):
        parse_sample_sheet(str(csv_path))


def test_sample_sheet_rejects_path_traversal_id(tmp_path):
    csv_path = tmp_path / "s.csv"
    csv_path.write_text(
        "sample_id,i7,i5\n../escape,ACGT,ACGT\nok,TGCA,TGCA\n", encoding="utf-8"
    )
    with pytest.raises(SampleSheetError, match="unsafe sample_id"):
        parse_sample_sheet(str(csv_path))


def test_sample_sheet_rejects_reserved_undetermined_id(tmp_path):
    csv_path = tmp_path / "s.csv"
    csv_path.write_text(
        "sample_id,i7,i5\nundetermined,ACGT,ACGT\n", encoding="utf-8"
    )
    with pytest.raises(SampleSheetError, match="reserved"):
        parse_sample_sheet(str(csv_path))
