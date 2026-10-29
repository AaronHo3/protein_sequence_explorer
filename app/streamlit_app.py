# Entry point: uv run streamlit run app/streamlit_app.py
# All interface code lives in src/protein_explorer/ui; pages are app/views/*.py.

from protein_explorer.ui.app import run

run()
