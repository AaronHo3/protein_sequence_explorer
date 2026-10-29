# Compact metric cards rendered as one HTML grid (denser and more consistent than one
# st.metric per column). All text is escaped.

from dataclasses import dataclass
from html import escape

import streamlit as st

@dataclass(frozen=True)
class Metric:
    label: str
    value: str
    unit: str = ""
    detail: str = ""
    help: str = ""  # shown as a tooltip on the label

def metric_cards_html(metrics: list[Metric]) -> str:
    cards = []
    for m in metrics:
        title = f' title="{escape(m.help)}"' if m.help else ""
        unit = f'<span class="pe-card-unit">{escape(m.unit)}</span>' if m.unit else ""
        detail = f'<div class="pe-card-detail">{escape(m.detail)}</div>' if m.detail else ""
        cards.append(
            f'<div class="pe-card"><div class="pe-card-label"{title}>{escape(m.label)}</div>'
            f'<div class="pe-card-value">{escape(m.value)}{unit}</div>{detail}</div>'
        )
    return f'<div class="pe-cards">{"".join(cards)}</div>'

def metric_cards(metrics: list[Metric]) -> None:
    st.html(metric_cards_html(metrics))
