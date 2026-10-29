// Interactive sequence viewer (Streamlit components v2). Renders the payload built by
// protein_explorer.visualization.sequence_viewer.build_viewer_data and reports residue
// selections to Python as a trigger value: {start, end} (1-based, inclusive) or {clear: true}.
//
// Interaction: click selects a residue, shift-click extends to a range, arrow keys move the
// selection (shift extends), Escape clears. Hover shows residue details.

const CELL_WIDTH = 17; // px, must match .pev-r width in the CSS
const CELL_GAP = 1;
const BLOCK_GAP = 10;
const NUMBER_COLUMN = 56;

function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, (c) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  }[c]));
}

function hexToRgb(hex) {
  const v = hex.replace("#", "");
  return [0, 2, 4].map((i) => parseInt(v.slice(i, i + 2), 16));
}

function trackColor(track, value) {
  if (value === null || value === undefined) return "transparent";
  const t = Math.max(-1, Math.min(1, value >= 0 ? value / track.high : -value / track.low));
  const from = hexToRgb(track.midColor);
  const to = hexToRgb(value >= 0 ? track.highColor : track.lowColor);
  const mixed = from.map((c, i) => Math.round(c + (to[i] - c) * Math.abs(t)));
  return `rgb(${mixed.join(",")})`;
}

function blocksPerLine(data, width) {
  if (data.lineLength) return Math.max(1, Math.round(data.lineLength / data.blockSize));
  const blockWidth = data.blockSize * (CELL_WIDTH + CELL_GAP) + BLOCK_GAP;
  return Math.max(1, Math.floor((width - NUMBER_COLUMN) / blockWidth));
}

function positionAt(data, index) {
  return data.positions ? data.positions[index] : index + 1;
}

function renderLines(data, perLine) {
  const seq = data.sequence;
  const lineSize = perLine * data.blockSize;
  const parts = [];
  for (let start = 0; start < seq.length; start += lineSize) {
    const end = Math.min(start + lineSize, seq.length);
    const firstPos = positionAt(data, start);
    const annotationRow = renderAnnotationRow(data, start, end);
    let cells = "";
    let tracks = data.tracks.map(() => "");
    for (let b = start; b < end; b += data.blockSize) {
      let blockCells = "";
      const trackBlocks = data.tracks.map(() => "");
      for (let i = b; i < Math.min(b + data.blockSize, end); i++) {
        const symbol = seq[i];
        const info = data.lookup[symbol];
        blockCells += `<span class="pev-r" data-i="${i}" style="background:${info.fill};color:${info.text}">${escapeHtml(symbol)}</span>`;
        data.tracks.forEach((track, t) => {
          trackBlocks[t] += `<span class="pev-t" data-i="${i}" data-track="${t}" style="background:${trackColor(track, track.values[i])}"></span>`;
        });
      }
      cells += `<span class="pev-block">${blockCells}</span>`;
      tracks = tracks.map((row, t) => row + `<span class="pev-block">${trackBlocks[t]}</span>`);
    }
    const trackRows = tracks
      .map((row) => `<div class="pev-row pev-track-row"><span class="pev-num"></span>${row}</div>`)
      .join("");
    parts.push(
      `<div class="pev-line">${annotationRow}` +
      `<div class="pev-row"><span class="pev-num">${firstPos ?? ""}</span>${cells}</div>` +
      `${trackRows}</div>`
    );
  }
  return parts.join("");
}

function renderAnnotationRow(data, start, end) {
  if (!data.annotations.length) return "";
  let row = "";
  for (let b = start; b < end; b += data.blockSize) {
    let block = "";
    for (let i = b; i < Math.min(b + data.blockSize, end); i++) {
      const pos = positionAt(data, i);
      const hit = data.annotations.find((a) => pos !== null && pos >= a.start && pos <= a.end);
      block += hit
        ? `<span class="pev-a" style="background:${hit.color}" title="${escapeHtml(hit.label)}"></span>`
        : `<span class="pev-a"></span>`;
    }
    row += `<span class="pev-block">${block}</span>`;
  }
  return `<div class="pev-row pev-annotation-row"><span class="pev-num"></span>${row}</div>`;
}

function renderLegend(data) {
  const items = data.legend
    .map((c) => `<span class="pev-key"><span class="pev-swatch" style="background:${c.color}"></span>${escapeHtml(c.label)}</span>`)
    .join("");
  const tracks = data.tracks
    .map((t) => `<span class="pev-key" title="${escapeHtml(t.description)}"><span class="pev-gradient" style="background:linear-gradient(90deg, ${t.lowColor}, ${t.midColor}, ${t.highColor})"></span>Strip: ${escapeHtml(t.name)} (${t.low} to +${t.high})</span>`)
    .join("");
  return `<div class="pev-legend">${items}${tracks}</div>`;
}

function selectionLabel(data, selection) {
  if (!selection) return "Click a residue to inspect it. Shift-click selects a range.";
  const [start, end] = selection;
  if (start === end) {
    const index = data.positions ? data.positions.indexOf(start) : start - 1;
    const info = data.lookup[data.sequence[index]] || {};
    return `Selected residue ${start} (${escapeHtml(info.three || "?")})`;
  }
  return `Selected residues ${start}–${end} (${end - start + 1} residues)`;
}

function tooltipHtml(data, index) {
  const symbol = data.sequence[index];
  const info = data.lookup[symbol];
  const pos = positionAt(data, index);
  const rows = [
    ["Position", pos ?? "gap"],
    ["Class", info.chemicalClass],
  ];
  if (info.category !== info.chemicalClass) rows.push(["Colour group", info.category]);
  rows.push(["Hydropathy", info.hydropathy === null ? "not defined" : `${info.hydropathy > 0 ? "+" : ""}${info.hydropathy} (residue)`]);
  data.tracks.forEach((t) => {
    const v = t.values[index];
    rows.push([t.name, v === null || v === undefined ? "no full window" : `${v > 0 ? "+" : ""}${v.toFixed(2)} (window mean)`]);
  });
  const body = rows.map(([k, v]) => `<span class="pev-tt-k">${escapeHtml(k)}</span><span>${escapeHtml(v)}</span>`).join("");
  return `<div class="pev-tt-title">${escapeHtml(info.name)} <span class="pev-tt-code">${escapeHtml(info.three)}, ${escapeHtml(symbol)}</span></div><div class="pev-tt-grid">${body}</div>`;
}

async function copyText(text) {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch (e) {
    const area = document.createElement("textarea");
    area.value = text;
    document.body.appendChild(area);
    area.select();
    const ok = document.execCommand("copy");
    area.remove();
    return ok;
  }
}

export default function (component) {
  const { data, parentElement, setTriggerValue } = component;
  if (!data || !data.sequence) return;

  let root = parentElement.querySelector(".pev");
  if (!root) {
    root = document.createElement("div");
    root.className = "pev";
    root.tabIndex = 0;
    root.setAttribute("role", "application");
    root.setAttribute("aria-label", "Protein sequence viewer");
    parentElement.appendChild(root);
  }
  root.style.setProperty("--pev-selection", data.colors.selection);
  root.style.setProperty("--pev-muted", data.colors.muted);
  root.style.setProperty("--pev-raised", data.colors.raised);

  const selection = data.selection;
  // The shift-click anchor lives on the element, not in the payload, so a shift-click made
  // before Python's rerun returns still extends from the residue clicked first.
  if (!selection) root.pevAnchor = null;
  else if (root.pevAnchor == null) root.pevAnchor = selection[0];

  function render() {
    const scrollTop = root.querySelector(".pev-scroll")?.scrollTop ?? 0;
    const perLine = blocksPerLine(data, root.clientWidth || 800);
    root.innerHTML =
      `<div class="pev-toolbar"><span class="pev-status">${selectionLabel(data, selection)}</span>` +
      `<span class="pev-actions">` +
      `<button type="button" data-action="copy">Copy sequence</button>` +
      (selection ? `<button type="button" data-action="copy-selection">Copy selection</button><button type="button" data-action="clear">Clear</button>` : "") +
      `</span></div>` +
      renderLegend(data) +
      `<div class="pev-scroll"><div class="pev-lines">${renderLines(data, perLine)}</div></div>` +
      `<div class="pev-tooltip" hidden></div>`;
    root.querySelector(".pev-scroll").scrollTop = scrollTop;
    markSelection();
  }

  function markSelection() {
    if (!selection) return;
    const [start, end] = selection;
    root.querySelectorAll(".pev-r").forEach((cell) => {
      const pos = positionAt(data, Number(cell.dataset.i));
      if (pos !== null && pos >= start && pos <= end) {
        cell.classList.add("pev-sel");
        if (pos === start) cell.classList.add("pev-sel-start");
        if (pos === end) cell.classList.add("pev-sel-end");
      }
    });
    const first = root.querySelector(".pev-sel");
    const scroller = root.querySelector(".pev-scroll");
    if (first && scroller) {
      const top = first.offsetTop - scroller.offsetTop;
      if (top < scroller.scrollTop || top > scroller.scrollTop + scroller.clientHeight - 30) {
        scroller.scrollTop = Math.max(0, top - 40);
      }
    }
  }

  function select(start, end, anchor = start) {
    root.pevAnchor = anchor;
    setTriggerValue("select", { start: Math.min(start, end), end: Math.max(start, end) });
  }

  function onClick(event) {
    const button = event.target.closest("button[data-action]");
    if (button) {
      const action = button.dataset.action;
      if (action === "clear") {
        root.pevAnchor = null;
        setTriggerValue("select", { clear: true });
      }
      if (action === "copy" || action === "copy-selection") {
        let text = data.sequence;
        if (action === "copy-selection" && selection) {
          const idx = (p) => (data.positions ? data.positions.indexOf(p) : p - 1);
          text = data.sequence.slice(idx(selection[0]), idx(selection[1]) + 1);
        }
        copyText(text.replaceAll("-", "")).then((ok) => {
          button.textContent = ok ? "Copied" : "Copy failed";
          setTimeout(() => render(), 1200);
        });
      }
      return;
    }
    const cell = event.target.closest(".pev-r, .pev-t");
    if (!cell) return;
    const pos = positionAt(data, Number(cell.dataset.i));
    if (pos === null) return;
    if (event.shiftKey && root.pevAnchor != null) select(root.pevAnchor, pos, root.pevAnchor);
    else select(pos, pos);
  }

  function onKey(event) {
    if (event.key === "Escape" && selection) {
      root.pevAnchor = null;
      setTriggerValue("select", { clear: true });
      return;
    }
    if (!selection || !["ArrowLeft", "ArrowRight"].includes(event.key)) return;
    event.preventDefault();
    const maxPos = data.positions ? Math.max(...data.positions.filter((p) => p !== null)) : data.sequence.length;
    const step = event.key === "ArrowRight" ? 1 : -1;
    if (event.shiftKey) {
      const anchor = root.pevAnchor ?? selection[0];
      const moving = selection[0] === anchor ? selection[1] : selection[0];
      select(anchor, Math.min(maxPos, Math.max(1, moving + step)), anchor);
    } else {
      const next = Math.min(maxPos, Math.max(1, selection[1] + step));
      select(next, next);
    }
  }

  function onOver(event) {
    const cell = event.target.closest(".pev-r, .pev-t");
    const tip = root.querySelector(".pev-tooltip");
    if (!cell || !tip) return;
    tip.innerHTML = tooltipHtml(data, Number(cell.dataset.i));
    tip.hidden = false;
    const rootBox = root.getBoundingClientRect();
    const box = cell.getBoundingClientRect();
    const left = Math.min(Math.max(0, box.left - rootBox.left - 8), rootBox.width - tip.offsetWidth);
    const above = box.top - rootBox.top - tip.offsetHeight - 6;
    tip.style.left = `${left}px`;
    tip.style.top = `${above >= 0 ? above : box.bottom - rootBox.top + 6}px`;
  }

  function onOut(event) {
    if (!event.relatedTarget || !root.contains(event.relatedTarget) || !event.relatedTarget.closest(".pev-r, .pev-t")) {
      const tip = root.querySelector(".pev-tooltip");
      if (tip) tip.hidden = true;
    }
  }

  render();
  root.addEventListener("click", onClick);
  root.addEventListener("keydown", onKey);
  root.addEventListener("mouseover", onOver);
  root.addEventListener("mouseout", onOut);

  let lastWidth = root.clientWidth;
  const observer = new ResizeObserver(() => {
    if (!data.lineLength && Math.abs(root.clientWidth - lastWidth) > 4) {
      lastWidth = root.clientWidth;
      render();
    }
  });
  observer.observe(root);

  return () => {
    observer.disconnect();
    root.removeEventListener("click", onClick);
    root.removeEventListener("keydown", onKey);
    root.removeEventListener("mouseover", onOver);
    root.removeEventListener("mouseout", onOut);
  };
}
