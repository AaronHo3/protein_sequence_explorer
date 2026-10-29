# Streamlit wrapper for the interactive sequence viewer (components v2). The browser half is
# ui/assets/sequence_viewer.{js,css}; the payload comes from visualization.sequence_viewer.

from collections.abc import Callable
from importlib.resources import files

import streamlit as st

from protein_explorer.visualization.sequence_viewer import ResidueRange

_ASSETS = files("protein_explorer.ui") / "assets"

_COMPONENT = st.components.v2.component(
    "protein_explorer_sequence_viewer",
    js=(_ASSETS / "sequence_viewer.js").read_text(encoding="utf-8"),
    css=(_ASSETS / "sequence_viewer.css").read_text(encoding="utf-8"),
)

def parse_selection_event(event: object) -> tuple[bool, ResidueRange | None]:
    """Interpret a viewer trigger value. Returns (is_event, selection): a clear event gives
    (True, None); a missing or malformed value gives (False, None)."""
    if not isinstance(event, dict):
        return False, None
    if event.get("clear"):
        return True, None
    try:
        return True, ResidueRange(int(event["start"]), int(event["end"]))
    except (KeyError, TypeError, ValueError):
        return False, None

def sequence_viewer(
    data: dict, key: str, on_select: Callable[[ResidueRange | None], None]
) -> None:
    """Mount the viewer. on_select runs (before the rerun) when the user selects or clears."""
    def handle() -> None:
        result = st.session_state.get(key)
        event = getattr(result, "select", None)
        if event is None and isinstance(result, dict):
            event = result.get("select")
        is_event, selection = parse_selection_event(event)
        if is_event:
            on_select(selection)

    _COMPONENT(key=key, data=data, on_select_change=handle)
