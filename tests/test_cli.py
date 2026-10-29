"""Tests for the command-line summary."""

from pathlib import Path

from protein_explorer.__main__ import main

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_FASTA = PROJECT_ROOT / "data" / "sample.fasta"
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_summarize_sample_lists_every_record(capsys):
    assert main(["summarize", str(SAMPLE_FASTA)]) == 0
    out = capsys.readouterr().out
    assert "sp|P01308|INS_HUMAN" in out
    assert "sp|P02769|ALBU_BOVIN" in out
    assert "Records: 8  Valid: 8  Invalid: 0" in out
    assert "Length range: 110-607 residues" in out


def test_summarize_reports_invalid_sequence_without_failing(tmp_path, capsys):
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(">good\nMKTAYI\n>gapped\nMK-TAYI\n")
    assert main(["summarize", str(fasta)]) == 0
    out = capsys.readouterr().out
    assert "NO" in out
    assert "Records: 2  Valid: 1  Invalid: 1" in out


def test_summarize_missing_file_returns_error(tmp_path, capsys):
    assert main(["summarize", str(tmp_path / "missing.fasta")]) == 1
    assert "file not found" in capsys.readouterr().err


def test_summarize_malformed_file_returns_error(capsys):
    assert main(["summarize", str(FIXTURES / "orphan_sequence.fasta")]) == 1
    err = capsys.readouterr().err
    assert "not valid FASTA" in err
    assert "Line 1" in err
