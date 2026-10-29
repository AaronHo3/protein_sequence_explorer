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


def test_analyze_writes_all_outputs(tmp_path, capsys):
    output = tmp_path / "results"
    assert main(["analyze", str(SAMPLE_FASTA), "--output", str(output)]) == 0
    assert {p.name for p in output.iterdir()} == {
        "features.csv", "composition.csv", "summary.csv", "run_info.json",
    }
    assert "Analyzed 8 records (8 valid, 0 invalid)" in capsys.readouterr().out


def test_analyze_ph_option_names_charge_column(tmp_path):
    output = tmp_path / "results"
    assert main(["analyze", str(SAMPLE_FASTA), "--output", str(output), "--ph", "7.4"]) == 0
    header = (output / "features.csv").read_text().splitlines()[0]
    assert "net_charge_ph7.4" in header


def test_analyze_keeps_invalid_records(tmp_path, capsys):
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(">good\nMKTAYI\n>gapped\nMK-TAYI\n")
    output = tmp_path / "results"
    assert main(["analyze", str(fasta), "--output", str(output)]) == 0
    assert "(1 valid, 1 invalid)" in capsys.readouterr().out
    assert len((output / "features.csv").read_text().splitlines()) == 3  # header + 2 rows


def test_analyze_malformed_file_returns_error_and_writes_nothing(tmp_path, capsys):
    output = tmp_path / "results"
    assert main(["analyze", str(FIXTURES / "empty_record.fasta"), "--output", str(output)]) == 1
    assert "not valid FASTA" in capsys.readouterr().err
    assert not output.exists()


def test_plot_writes_figures(tmp_path, capsys):
    output = tmp_path / "figures"
    assert main(["plot", str(SAMPLE_FASTA), "--output", str(output), "--format", "svg"]) == 0
    assert "Wrote 20 figures" in capsys.readouterr().out
    assert (output / "mw_vs_length.svg").exists()
    assert (output / "hydropathy" / "INS_HUMAN.svg").exists()


def test_plot_rejects_even_window(tmp_path, capsys):
    assert main(["plot", str(SAMPLE_FASTA), "--output", str(tmp_path), "--window", "8"]) == 1
    assert "odd" in capsys.readouterr().err
