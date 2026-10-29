# Figure exports: a ZIP of the publication-resolution static figures (matplotlib, 300 dpi).
# Interactive charts export individually through their toolbar (theme.plotly_config).

import io
import tempfile
import zipfile
from pathlib import Path

from protein_explorer.fasta import FastaRecord
from protein_explorer.plots import write_all_figures

def static_figure_bundle(records: list[FastaRecord], window: int = 9, fmt: str = "png") -> bytes:
    """All static figures (dataset-level plus per-protein) as a ZIP archive in memory."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        paths = write_all_figures(records, root, window=window, fmt=fmt)
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in paths:
                archive.write(path, path.relative_to(root).as_posix())
    return buffer.getvalue()
