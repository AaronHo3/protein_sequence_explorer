# Global CSS: layout width, spacing rhythm, and the shared card/header classes used by the
# HTML widgets. Colours come from the active palette so light and dark both work.

import streamlit as st

from protein_explorer.visualization.theme import ChartPalette

def inject_styles(palette: ChartPalette) -> None:
    st.html(f"""<style>
[data-testid="stMainBlockContainer"] {{ max-width: 1360px; padding-top: 3.2rem; padding-bottom: 4rem; }}
h1 {{ font-size: 1.85rem !important; letter-spacing: -0.01em; }}
h2 {{ font-size: 1.3rem !important; margin-top: 0.6rem !important; }}
h3 {{ font-size: 1.08rem !important; }}
.pe-page-intro {{ color: {palette.ink_secondary}; font-size: 0.98rem; max-width: 72ch; margin: -0.4rem 0 0.4rem; }}
.pe-cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; margin: 0.25rem 0 0.5rem; }}
.pe-card {{ background: {palette.raised}; border: 1px solid {palette.gridline}; border-radius: 8px;
            padding: 10px 14px 11px; box-shadow: 0 1px 2px rgba(0,0,0,0.04); min-width: 0; }}
.pe-card-label {{ font-size: 12.5px; color: {palette.ink_secondary}; line-height: 1.3; }}
.pe-card-label[title] {{ text-decoration: underline dotted {palette.baseline}; text-underline-offset: 3px; cursor: help; }}
.pe-card-value {{ font-size: 22px; font-weight: 600; color: {palette.ink}; font-variant-numeric: tabular-nums; margin-top: 2px; white-space: nowrap; }}
.pe-card-unit {{ font-size: 13px; font-weight: 400; color: {palette.muted}; margin-left: 3px; }}
.pe-card-detail {{ font-size: 12px; color: {palette.muted}; margin-top: 1px; line-height: 1.3; }}
.pe-protein-title {{ display: flex; flex-wrap: wrap; align-items: baseline; gap: 4px 14px; }}
.pe-protein-title h1 {{ margin: 0; padding: 0; }}
.pe-protein-sub {{ color: {palette.ink_secondary}; font-size: 0.98rem; }}
.pe-protein-sub em {{ color: {palette.ink_secondary}; }}
.pe-protein-sub a {{ color: {palette.accent}; }}
.pe-inspector {{ background: {palette.raised}; border: 1px solid {palette.gridline}; border-radius: 8px; padding: 14px 16px; }}
.pe-inspector-head {{ display: flex; align-items: center; gap: 12px; margin-bottom: 10px; }}
.pe-inspector-tile {{ width: 44px; height: 44px; border-radius: 8px; display: flex; align-items: center; justify-content: center;
                     font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 22px; font-weight: 600; flex: none; }}
.pe-inspector-name {{ font-weight: 600; font-size: 1.02rem; color: {palette.ink}; }}
.pe-inspector-codes {{ color: {palette.muted}; font-size: 13px; }}
.pe-inspector dl {{ display: grid; grid-template-columns: minmax(0, 1.15fr) minmax(0, 1fr); gap: 5px 12px; margin: 0; padding: 0; font-size: 13.5px; }}
.pe-inspector dt {{ color: {palette.ink_secondary}; margin: 0; padding: 0; }}
.pe-inspector dd {{ margin: 0; color: {palette.ink}; font-variant-numeric: tabular-nums; }}
.pe-inspector-seq {{ font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 12.5px; word-break: break-all;
                     color: {palette.ink_secondary}; margin-top: 10px; }}
.pe-empty {{ color: {palette.muted}; font-size: 13.5px; }}
.pe-status {{ display: inline-block; font-size: 12px; padding: 1px 8px; border-radius: 10px; border: 1px solid {palette.gridline}; color: {palette.ink_secondary}; }}
</style>""")

def page_header(title: str, intro: str = "") -> None:
    st.title(title)
    if intro:
        st.html(f'<p class="pe-page-intro">{intro}</p>')
