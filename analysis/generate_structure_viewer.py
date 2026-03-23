"""
analysis/generate_structure_viewer.py
--------------------------------------
Generates a self-contained HTML structure viewer for the 6 FoldSeek bars.

Embeds all 5 Boltz-2 diffusion models per bar (30 PDB files total).
Uses 3Dmol.js (CDN) for interactive 3D rendering.

Features:
  - 2×3 grid, one panel per bar
  - Color by pLDDT (B-factor) — AlphaFold-style gradient
  - Toggle: rainbow (N→C) / surface / stick
  - Per-panel model selector (0–4)
  - Global: spin, reset, style controls
  - Metrics: pLDDT, pTM, seq len, structural class, song

Output: outputs/structure_viewer/viewer_foldseek_bars.html

Usage:
  python analysis/generate_structure_viewer.py
"""

import csv
import json
from pathlib import Path

ROOT        = Path(__file__).resolve().parent.parent
ID_MAP      = ROOT / "data/boltz_id_map.json"
SNAPSHOT    = ROOT / "data/bar_index_snapshot.json"
SUMMARY     = ROOT / "outputs/boltz/boltz_summary.csv"
MODELS_CSV  = ROOT / "outputs/boltz/boltz_models.csv"
PRED_DIR    = ROOT / "outputs/boltz_outputs/boltz_results_boltz_inputs/predictions"
OUT_DIR     = ROOT / "outputs/structure_viewer"

BARS = ["bar_11", "bar_17", "bar_27", "bar_38", "bar_49", "bar_53"]
N_MODELS = 5

CLASS_BADGE = {
    "confident_protein_like":  ("Confident protein-like", "#27ae60"),
    "uncertain_protein_like":  ("Uncertain protein-like", "#f39c12"),
    "disordered":              ("Disordered",             "#95a5a6"),
}
BAR_ACCENT = {
    "bar_11": "#e74c3c",
    "bar_17": "#3498db",
    "bar_27": "#2ecc71",
    "bar_38": "#f39c12",
    "bar_49": "#9b59b6",
    "bar_53": "#1abc9c",
}


def load_meta():
    id_map   = json.load(open(ID_MAP))
    inv_map  = {v: k for k, v in id_map.items()}
    snap     = json.load(open(SNAPSHOT))
    summary  = {r["bar_id"]: r for r in csv.DictReader(open(SUMMARY))}
    models   = {}
    for row in csv.DictReader(open(MODELS_CSV)):
        bid = row["bar_id"]
        if bid not in models:
            models[bid] = {}
        models[bid][int(row["model"])] = {
            "plddt": float(row.get("plddt") or 0),
            "ptm":   float(row.get("ptm") or 0),
        }

    meta = {}
    for bar_id in BARS:
        r  = summary[bar_id]
        s  = snap[bar_id]
        meta[bar_id] = {
            "boltz_id":     inv_map[bar_id],
            "song":         s.get("genius_song_title", ""),
            "bar":          s.get("canonical_bar", ""),
            "iconicity":    float(s.get("aggregate_iconicity") or 0),
            "seq_len":      int(r.get("seq_len") or 0),
            "best_model":   int(r.get("best_model") or 0),
            "plddt_mean":   float(r.get("plddt_mean") or 0),
            "plddt_best":   float(r.get("plddt_best") or 0),
            "ptm_mean":     float(r.get("ptm_mean") or 0),
            "struct_class": r.get("structural_class", "disordered"),
            "models":       models.get(bar_id, {}),
        }
    return meta


def fix_boltz_chain(pdb_text: str) -> str:
    """Boltz writes 3-char chain IDs (e.g. 'b11') which shift coordinate columns
    by 2 relative to standard PDB fixed-width format, breaking 3Dmol's parser.
    Replace the 3-char chain with 'A' so cols 22-26 (resSeq) and 31-54 (xyz)
    are correctly positioned."""
    out = []
    for line in pdb_text.splitlines(keepends=True):
        if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 32:
            line = line[:21] + "A" + line[24:]
        out.append(line)
    return "".join(out)


def load_pdbs(meta):
    pdbs = {}
    for bar_id in BARS:
        boltz_id = meta[bar_id]["boltz_id"]
        pdbs[bar_id] = {}
        for mi in range(N_MODELS):
            path = PRED_DIR / boltz_id / f"{boltz_id}_model_{mi}.pdb"
            content = path.read_text(errors="replace")
            content = fix_boltz_chain(content)
            # Escape backticks and backslashes for JS template literal
            content = content.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")
            pdbs[bar_id][mi] = content
    return pdbs


def render_html(meta, pdbs) -> str:
    # --- Build JS data blobs ---
    meta_js_parts = []
    for bar_id in BARS:
        m = meta[bar_id]
        cls_label, cls_color = CLASS_BADGE.get(m["struct_class"], ("Unknown", "#888"))
        model_metrics = {
            str(i): {"plddt": m["models"].get(i, {}).get("plddt", 0),
                     "ptm":   m["models"].get(i, {}).get("ptm", 0)}
            for i in range(N_MODELS)
        }
        entry = f"""  "{bar_id}": {{
    song: {json.dumps(m["song"])},
    bar: {json.dumps(m["bar"][:80])},
    seqLen: {m["seq_len"]},
    bestModel: {m["best_model"]},
    plddtMean: {m["plddt_mean"]:.4f},
    plddtBest: {m["plddt_best"]:.4f},
    ptmMean: {m["ptm_mean"]:.4f},
    iconicity: {m["iconicity"]:.4f},
    structClass: {json.dumps(m["struct_class"])},
    classLabel: {json.dumps(cls_label)},
    classColor: {json.dumps(cls_color)},
    accent: {json.dumps(BAR_ACCENT[bar_id])},
    modelMetrics: {json.dumps(model_metrics)}
  }}"""
        meta_js_parts.append(entry)

    pdb_js_parts = []
    for bar_id in BARS:
        model_parts = ", ".join(
            f'{i}: `{pdbs[bar_id][i]}`' for i in range(N_MODELS)
        )
        pdb_js_parts.append(f'  "{bar_id}": {{ {model_parts} }}')

    meta_js = "const BAR_META = {\n" + ",\n".join(meta_js_parts) + "\n};"
    pdb_js  = "const BAR_PDBS = {\n" + ",\n".join(pdb_js_parts) + "\n};"

    bars_json = json.dumps(BARS)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Structure Viewer — FoldSeek Bars</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  :root {{
    --bg:      #0f0f1a;
    --surface: #1a1a2e;
    --card:    #16213e;
    --border:  #2a2a4a;
    --text:    #e8e8f0;
    --muted:   #8888aa;
    --radius:  10px;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    background: var(--bg); color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    min-height: 100vh;
  }}
  header {{
    padding: 18px 28px 14px;
    background: var(--surface);
    border-bottom: 1px solid var(--border);
    display: flex; align-items: center; justify-content: space-between;
    flex-wrap: wrap; gap: 12px;
  }}
  header h1 {{ font-size: 1.1rem; font-weight: 600; letter-spacing: 0.01em; }}
  header p  {{ font-size: 0.78rem; color: var(--muted); margin-top: 2px; }}
  .controls {{
    display: flex; gap: 8px; flex-wrap: wrap; align-items: center;
  }}
  .btn {{
    padding: 6px 14px; border-radius: 6px; border: 1px solid var(--border);
    background: var(--card); color: var(--text); font-size: 0.78rem;
    cursor: pointer; transition: all .15s; user-select: none;
  }}
  .btn:hover  {{ background: #2a2a4a; border-color: #4a4a8a; }}
  .btn.active {{ background: #2a2a6a; border-color: #6a6aaa; color: #aaddff; }}
  select.btn  {{ appearance: none; padding-right: 24px;
                background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%238888aa'/%3E%3C/svg%3E");
                background-repeat: no-repeat; background-position: right 8px center; }}
  .grid {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 16px;
    padding: 20px;
    max-width: 1400px;
    margin: 0 auto;
  }}
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    overflow: hidden;
    display: flex; flex-direction: column;
    transition: box-shadow .2s;
  }}
  .card:hover {{ box-shadow: 0 0 0 1px var(--border), 0 4px 24px rgba(0,0,0,.4); }}
  .card-header {{
    padding: 12px 14px 8px;
    border-bottom: 2px solid var(--border);
    display: flex; justify-content: space-between; align-items: flex-start;
  }}
  .card-header .bar-id {{
    font-size: 0.72rem; font-weight: 700; letter-spacing: 0.05em;
    text-transform: uppercase; color: var(--muted); margin-bottom: 2px;
  }}
  .card-header .song {{ font-size: 0.9rem; font-weight: 600; }}
  .badge {{
    font-size: 0.65rem; padding: 2px 8px; border-radius: 20px;
    font-weight: 600; letter-spacing: 0.03em; white-space: nowrap;
    margin-top: 2px;
  }}
  .viewer-wrap {{
    position: relative; width: 100%; height: 360px;
    background: #0a0a18; overflow: hidden;
  }}
  .viewer-wrap .viewer-inner {{
    width: 100%; height: 360px; position: relative; display: block;
  }}
  .card-footer {{
    padding: 10px 14px;
    display: flex; flex-direction: column; gap: 6px;
  }}
  .metrics {{
    display: flex; gap: 10px; flex-wrap: wrap;
  }}
  .metric {{
    display: flex; flex-direction: column;
    font-size: 0.7rem;
  }}
  .metric .label {{ color: var(--muted); margin-bottom: 1px; }}
  .metric .value {{ font-weight: 700; font-size: 0.82rem; }}
  .model-row {{
    display: flex; align-items: center; gap: 6px; margin-top: 2px;
  }}
  .model-row span {{ font-size: 0.7rem; color: var(--muted); }}
  .model-btns {{ display: flex; gap: 4px; }}
  .model-btn {{
    width: 26px; height: 22px; border-radius: 4px; border: 1px solid var(--border);
    background: var(--surface); color: var(--muted); font-size: 0.68rem;
    cursor: pointer; transition: all .12s; display: flex; align-items: center;
    justify-content: center; font-weight: 600;
  }}
  .model-btn:hover   {{ background: #2a2a4a; color: var(--text); }}
  .model-btn.active  {{ color: white; border-color: transparent; }}
  .lyric-snippet {{
    font-size: 0.68rem; color: var(--muted); font-style: italic;
    line-height: 1.4; margin-top: 2px;
    white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
  }}
  .legend {{
    display: flex; align-items: center; justify-content: center;
    gap: 20px; padding: 12px 20px; flex-wrap: wrap;
    background: var(--surface); border-top: 1px solid var(--border);
    font-size: 0.72rem; color: var(--muted);
  }}
  .legend-gradient {{
    display: flex; align-items: center; gap: 8px;
  }}
  .gradient-bar {{
    width: 120px; height: 10px; border-radius: 5px;
    background: linear-gradient(to right, #FF7D45, #FFDB13, #65CBF3, #0053D6);
  }}
  .loading {{
    position: absolute; inset: 0; display: flex; align-items: center;
    justify-content: center; font-size: 0.75rem; color: var(--muted);
    background: #0a0a18; z-index: 10; pointer-events: none;
    transition: opacity .3s;
  }}
  @media (max-width: 900px) {{
    .grid {{ grid-template-columns: repeat(2, 1fr); }}
  }}
  @media (max-width: 600px) {{
    .grid {{ grid-template-columns: 1fr; }}
  }}
</style>
</head>
<body>

<header>
  <div>
    <h1>Boltz-2 Structure Viewer &nbsp;·&nbsp; FoldSeek Bars (pTM &ge; 0.4)</h1>
    <p>6 bars &nbsp;·&nbsp; 5 diffusion models each &nbsp;·&nbsp; concordance mapping &nbsp;·&nbsp; Nicki Minaj iconicity corpus</p>
  </div>
  <div class="controls">
    <button class="btn active" id="btn-cartoon"  onclick="setStyle('cartoon')">Cartoon</button>
    <button class="btn"        id="btn-surface"  onclick="setStyle('surface')">Surface</button>
    <button class="btn"        id="btn-stick"    onclick="setStyle('stick')">Stick</button>
    <span style="color:var(--border)">|</span>
    <button class="btn active" id="btn-plddt"    onclick="setColor('plddt')">pLDDT</button>
    <button class="btn"        id="btn-rainbow"  onclick="setColor('rainbow')">Rainbow</button>
    <span style="color:var(--border)">|</span>
    <button class="btn"        id="btn-spin"     onclick="toggleSpin()">&#9654; Spin</button>
    <button class="btn"        id="btn-reset"    onclick="resetAll()">Reset views</button>
  </div>
</header>

<div class="grid" id="grid"></div>

<div class="legend">
  <div class="legend-gradient">
    <span>pLDDT (B-factor):</span>
    <div class="gradient-bar"></div>
    <span>Low (30) &rarr; High (70)</span>
  </div>
  <span>&nbsp;|&nbsp; 3Dmol.js &nbsp;&middot;&nbsp; Boltz-2 predictions &nbsp;&middot;&nbsp; Cα RMSD: bar_17=3.4Å, bar_53=9.7Å mean across models</span>
</div>

<script>
{meta_js}

{pdb_js}

const BARS = {bars_json};

// ---------- Global state ----------
let currentStyle = 'cartoon';
let currentColor = 'plddt';
let spinning = false;
const viewers = {{}};
const activeModels = {{}};

// ---------- Style helpers ----------
function applyStyle(barId) {{
  const v = viewers[barId];
  const m = activeModels[barId];
  if (!v) return;
  v.removeAllModels();
  v.addModel(BAR_PDBS[barId][m], 'pdb');
  const style = buildStyle();
  v.setStyle({{}}, style);
  if (currentStyle === 'surface') {{
    v.addSurface('SES', {{
      opacity: 0.85,
      colorscheme: currentColor === 'plddt'
        ? {{prop: 'b', gradient: 'roygb', min: 30, max: 70}}
        : {{color: 'spectrum'}}
    }});
  }}
  v.zoomTo();
  v.render();
}}

function buildStyle() {{
  const cs = currentColor === 'plddt'
    ? {{prop: 'b', gradient: 'roygb', min: 30, max: 70}}
    : 'spectrum';
  if (currentStyle === 'cartoon') {{
    return {{cartoon: {{colorscheme: cs}}}};
  }} else if (currentStyle === 'stick') {{
    return {{stick: {{colorscheme: cs, radius: 0.15}}}};
  }}
  return {{}};
}}

// ---------- Controls ----------
function setStyle(s) {{
  currentStyle = s;
  ['cartoon','surface','stick'].forEach(x => {{
    document.getElementById('btn-' + x).classList.toggle('active', x === s);
  }});
  BARS.forEach(applyStyle);
}}

function setColor(c) {{
  currentColor = c;
  ['plddt','rainbow'].forEach(x => {{
    document.getElementById('btn-' + x).classList.toggle('active', x === c);
  }});
  BARS.forEach(applyStyle);
}}

function toggleSpin() {{
  spinning = !spinning;
  const btn = document.getElementById('btn-spin');
  btn.classList.toggle('active', spinning);
  btn.textContent = spinning ? '⏹ Stop' : '▶ Spin';
  BARS.forEach(b => {{
    if (viewers[b]) spinning ? viewers[b].spin('y') : viewers[b].spin(false);
  }});
}}

function resetAll() {{
  BARS.forEach(b => {{
    if (viewers[b]) {{ viewers[b].zoomTo(); viewers[b].render(); }}
  }});
}}

// ---------- Per-panel model switch ----------
function switchModel(barId, modelN) {{
  activeModels[barId] = modelN;
  // Update button states
  for (let i = 0; i < 5; i++) {{
    const btn = document.getElementById(`mb-${{barId}}-${{i}}`);
    if (btn) btn.classList.toggle('active', i === modelN);
  }}
  // Update metrics
  const mm = BAR_META[barId].modelMetrics[String(modelN)];
  document.getElementById(`plddt-${{barId}}`).textContent = mm.plddt.toFixed(3);
  document.getElementById(`ptm-${{barId}}`).textContent   = mm.ptm.toFixed(3);
  // Reload model
  applyStyle(barId);
}}

// ---------- Build grid ----------
function buildGrid() {{
  const grid = document.getElementById('grid');
  BARS.forEach(barId => {{
    const m = BAR_META[barId];
    activeModels[barId] = m.bestModel;
    const initPlddt = (m.modelMetrics[String(m.bestModel)]?.plddt || m.plddtMean).toFixed(3);
    const initPtm   = (m.modelMetrics[String(m.bestModel)]?.ptm   || m.ptmMean).toFixed(3);
    const [classLabel, classColor] = [m.classLabel, m.classColor];

    const modelBtns = Array.from({{length: 5}}, (_, i) => {{
      const isActive = i === m.bestModel;
      return `<button id="mb-${{barId}}-${{i}}" class="model-btn${{isActive ? ' active' : ''}}"
        style="${{isActive ? `background:${{m.accent}};` : ''}}"
        onclick="switchModel('${{barId}}', ${{i}})">M${{i}}</button>`;
    }}).join('');

    grid.insertAdjacentHTML('beforeend', `
      <div class="card" style="border-top: 3px solid ${{m.accent}}">
        <div class="card-header" style="border-bottom-color: ${{m.accent}}33">
          <div>
            <div class="bar-id">${{barId}}</div>
            <div class="song">${{m.song}}</div>
          </div>
          <span class="badge" style="background:${{classColor}}22; color:${{classColor}}; border:1px solid ${{classColor}}44">
            ${{classLabel}}
          </span>
        </div>
        <div class="viewer-wrap">
          <div id="viewer-${{barId}}" class="viewer-inner"></div>
          <div class="loading" id="loading-${{barId}}">Loading&hellip;</div>
        </div>
        <div class="card-footer">
          <div class="metrics">
            <div class="metric">
              <span class="label">pLDDT (mean)</span>
              <span class="value" id="plddt-${{barId}}" style="color:${{m.accent}}">${{initPlddt}}</span>
            </div>
            <div class="metric">
              <span class="label">pTM (mean)</span>
              <span class="value" id="ptm-${{barId}}">${{initPtm}}</span>
            </div>
            <div class="metric">
              <span class="label">Seq len</span>
              <span class="value">${{m.seqLen}} AA</span>
            </div>
            <div class="metric">
              <span class="label">Iconicity</span>
              <span class="value">${{m.iconicity.toFixed(3)}}</span>
            </div>
          </div>
          <div class="model-row">
            <span>Diffusion model:</span>
            <div class="model-btns">${{modelBtns}}</div>
          </div>
          <div class="lyric-snippet">&ldquo;${{m.bar}}&rdquo;</div>
        </div>
      </div>
    `);
  }});
}}

// ---------- Init viewers ----------
function initViewer(barId) {{
  return new Promise(resolve => {{
    setTimeout(() => {{
      const loader = document.getElementById(`loading-${{barId}}`);
      try {{
        if (typeof $3Dmol === 'undefined') throw new Error('3Dmol.js not loaded');
        const el = document.getElementById(`viewer-${{barId}}`);
        // 3Dmol needs explicit pixel dimensions on the container element itself
        const wrap = el.parentElement;
        const w = wrap.clientWidth  || wrap.offsetWidth  || 400;
        const h = 360;
        el.style.width  = w + 'px';
        el.style.height = h + 'px';
        // Pass ID string — classic 3Dmol usage
        const viewer = $3Dmol.createViewer(`viewer-${{barId}}`, {{
          backgroundColor: '#0a0a18',
          antialias: true,
          width: w,
          height: h,
        }});
        if (!viewer) throw new Error(`createViewer returned null (w=${{w}} h=${{h}})`);
        viewers[barId] = viewer;
        const modelN = activeModels[barId];
        viewer.addModel(BAR_PDBS[barId][modelN], 'pdb');
        viewer.setStyle({{}}, buildStyle());
        viewer.zoomTo();
        viewer.render();
        if (loader) {{ loader.style.opacity = '0'; setTimeout(() => loader.remove(), 300); }}
      }} catch(e) {{
        console.error(barId, e);
        if (loader) {{
          loader.style.opacity = '1';
          loader.style.color = '#e74c3c';
          loader.style.fontSize = '0.65rem';
          loader.style.padding = '8px';
          loader.style.textAlign = 'center';
          const msg = (e && e.message) ? e.message : String(e);
          loader.textContent = '3Dmol error: ' + msg;
        }}
      }}
      resolve();
    }}, 300);
  }});
}}

async function initAll() {{
  buildGrid();
  // Give the grid one full render cycle before initialising viewers
  await new Promise(r => setTimeout(r, 100));
  for (const barId of BARS) {{
    await initViewer(barId);
  }}
}}

window.addEventListener('load', initAll);
</script>
</body>
</html>"""


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("[generate_structure_viewer.py] Loading metadata...")
    meta = load_meta()
    print("  Loading 30 PDB files (5 models × 6 bars)...")
    pdbs = load_pdbs(meta)
    print("  Rendering HTML...")
    html = render_html(meta, pdbs)
    out_path = OUT_DIR / "viewer_foldseek_bars.html"
    out_path.write_text(html, encoding="utf-8")
    size_kb = out_path.stat().st_size // 1024
    print(f"  Saved {out_path}  ({size_kb} KB)")
    print(f"\nOpen in browser: open {out_path}")


if __name__ == "__main__":
    main()
