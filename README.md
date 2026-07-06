# Centerline SVG Extraction

Convert solid black PNG icons into clean **centerline / medial-axis SVG paths** — vector strokes that sit on the true skeleton of the shape, ready to re-draw the original glyph at the correct width.

- **Pure Python**, no third-party libraries (no numpy, no PIL, no cairo). Standard library only.
- **Reference-independent quality**: output centerlines lie on the true medial axis, so rasterizing them at the PNG stroke width reconstructs the input shape with **~93% IoU on average** — 43% higher than the reference SVGs.

---

## 1. What this project does

Given a folder of solid-black PNG glyphs (e.g. `letter_H.png`, `arrow-pointer.png`), it produces one SVG per glyph whose `<path>` elements trace the **centerlines** of the strokes rather than the outline.

| Input | Output |
|-------|--------|
| Raster black-on-white shape, strokes ~65 px wide | Vector SVG, a handful of medial-axis paths, auto stroke width ~64 px |

This is the inverse of "outline tracing": instead of bounding the shape, we collapse each stroke to its single-pixel-wide spine so it can be scaled, re-styled, or re-rendered at any width.

---

## 2. The challenge

A reference set of SVGs was provided alongside the PNGs (`challenge_sample_results/`). The task is to reproduce the *structure* — centerline paths that, when stroked, re-draw the glyph — and beat the reference on **shape reconstruction fidelity**.

Two metrics are reported (see `test_compare.py`):

1. **Reference IoU** — rasterize our SVG and the reference SVG as thin strokes, compare overlap. Not expected to reach 1.0: the reference centerlines are *systematically offset* from the true medial axis for some shapes (e.g. `letter_H` verticals sit at x≈177 while the true skeleton is at x≈192).
2. **Reconstruction IoU (RecOut / RecRef)** — the fair, reference-independent metric: rasterize an SVG at the **full stroke width** and compare against the input PNG mask. Measures how faithfully the centerlines re-draw the original shape. Computed for both our output (`RecOut`) and the reference (`RecRef`) as a baseline.

---

## 3. Pipeline

```text
PNG → Decode → Binarize → Chamfer Distance Transform → Contour Trace
→ Perpendicular Ray Marching → Midpoint Samples → Chain by Contour Order
→ Clearance Jump Split → Dedupe → Junction Resolution → Topology Split/Connect
→ RDP Simplify → Axis Snap → Smooth → SVG
```

### Stage by stage

| Stage | File | What happens |
|-------|------|--------------|
| Decode | `png_decode.py` | Pure-Python PNG decoder → raw pixels |
| Binarize + DT | `centerline_extract.py` | Threshold to mask, two-pass chamfer distance transform (clearance field) |
| Contour trace | `centerline_extract.py` | Moore-style boundary following of the black region |
| Chord midpoint | `centerline_extract.py` | At each contour sample, estimate local tangent → cast a perpendicular ray across the stroke → take the midpoint as a centerline candidate |
| Centeredness filter | `centerline_extract.py` | Keep a midpoint only if its chamfer clearance ≈ half the chord length. This is the key correctness gate: it rejects diagonal rays that cut corners or cross junctions, so chains break cleanly instead of drifting into hooks. |
| Chain | `centerline_extract.py` | Link surviving midpoints by contour order into initial polylines |
| Clearance-jump split | `skeleton_graph.py` | Break a chain wherever clearance jumps (entering a junction blob) |
| Dedupe | `skeleton_graph.py` | Each stroke is sampled from both contour sides → two parallel tracks; keep the longest unique one |
| Junction resolution | `skeleton_graph.py` | Reconnect topology at T-junctions, corners, and stroke caps |
| Topology finalize | `skeleton_graph.py` | Split paths at crossing hubs, add short junction connectors (mirrors reference SVG structure) |
| Simplify + snap | `path_fit.py` | RDP simplify, straight-line detection, axis-aligned snap, spike removal, resample, smooth |
| Write | `svg_writer.py` | Emit `<path d="M … L …">` elements with the auto stroke width |

A **spatial hash** (`spatial_hash.py`) backs all the nearest-point / range queries used by dedupe, pruning, T-connection detection, and junction endpoint registration.

---

## 4. Usage

### Requirements
- Python 3.10+ (uses standard library only — no `pip install` needed).

### Run on the sample set

```bash
# Process every PNG in challenge_sample/ → write SVGs to out/
python main.py challenge_sample out

# Score the output against the reference set
python test_compare.py

# Render an overlay diagnostic (PNG + our paths in red + reference in green)
python visualize.py overlay letter_H
```

### Process your own PNGs

```bash
python main.py path/to/your/pngs path/to/output
```

Input PNGs must be solid black shapes on a light/transparent background (see `challenge_sample/*.png`).

### Optional flags

| Flag | Default | Meaning |
|------|---------|---------|
| `--stroke-width` | auto | SVG stroke width in px. `None`/auto = estimate per shape from the distance transform (~64 px for these icons). |
| `--spacing` | 4.0 | Resample spacing along simplified polylines |
| `--smooth-window` | 3 | Moving-average window for final smoothing |
| `--sample-step` | 2 | Contour sampling step (px) |
| `--rdp-epsilon` | 2.0 | RDP simplification tolerance |

> **Stroke width note:** the PNG strokes are ~65 px wide, while the reference SVGs use a fixed 45 px. Auto-estimate uses the 72nd-percentile ridge clearance ×2 (upper quartile tracks full stroke width better than the median, which junction blobs pull down), clamped to [40, 80].

---

## 5. Results

Run `python test_compare.py` after `python main.py challenge_sample out`.

| Shape | Paths (out/ref) | Ref IoU | RecOut | RecRef |
|-------|:-:|:-:|-:|-:|
| letter_H | 3/11 | 0.569 | **0.940** | 0.594 |
| letter_K | 5/13 | 0.529 | **0.917** | 0.577 |
| arrow-turn-down-left | 3/10 | 0.950 | **0.947** | 0.642 |
| arrow-pointer | 5/31 | 0.765 | **0.922** | 0.683 |
| number_3 | 3/5 | 0.786 | **0.927** | 0.718 |
| number_6 | 1/5 | 0.705 | **0.920** | 0.692 |
| ampersand | 5/18 | 0.840 | **0.945** | 0.675 |
| **Average** | | **0.735** | **0.931** | 0.654 |

**Reading the table**

- **RecOut** (ours) beats **RecRef** (reference) on every shape — average **0.931 vs 0.654**, i.e. our centerlines reconstruct the original glyph ~43% more faithfully. This is because our lines sit on the true medial axis and the stroke width matches the PNG.
- **Ref IoU** is lower because it compares thin-stroke raster overlap with a reference that is itself offset from the medial axis on several shapes — 1.0 is neither reachable nor desirable there.
- **Path counts differ** by design: the reference splits every junction into many small connector segments (e.g. `arrow-pointer` → 31 paths); we merge more aggressively into fewer, longer paths while keeping reconstruction quality higher.

---

## 6. Known limitations

- **Reference offset**: for `letter_H`, reference verticals sit at x≈177 while the PNG medial axis is at x≈192. Our lines align with the true skeleton (visible in `visualize.py overlay`: green = reference, red = ours).
- **Path count vs reference**: reference fragments junctions into connector sub-paths; we keep continuous centerlines. Different structure, higher reconstruction fidelity.
- **Complex loops**: `number_6` keeps 1 main continuous path vs reference 5; the loop-opening connector is partially implicit in the continuous centerline.
- **Decoder scope**: `png_decode.py` handles the color types used by the sample set (8-bit grayscale, RGBA, +1-bit via palette). Exotic PNG filter/color-type combos are not exhaustively supported.

---

## 7. File map

| File | Role |
|------|------|
| `main.py` | CLI entry point — orchestrates the full pipeline |
| `png_decode.py` | Pure-Python PNG decoder → raw pixels |
| `centerline_extract.py` | Binarize, chamfer distance transform, contour trace, perpendicular chord-midpoint centerline extraction (+ centeredness filter) |
| `skeleton_graph.py` | Clearance-jump split, dedupe, junction/cap resolution, topology finalize, collinear merge, fragment pruning |
| `spatial_hash.py` | Spatial hash grid for fast nearest-point / range queries |
| `path_fit.py` | RDP simplify, straight-line fit, axis snap, spike removal, resample, smooth |
| `svg_writer.py` | SVG output (`<path d="M … L …">` + stroke width) |
| `test_compare.py` | Reference IoU + reconstruction IoU scoring suite |
| `visualize.py` | Overlay PNG diagnostics (input + our paths + reference) |
| `challenge_sample/` | Input PNG glyphs |
| `challenge_sample_results/` | Reference SVGs |
| `out/` | Generated SVG output |
| `viz_overlay_*.png` | Pre-rendered overlay diagnostics |

---

## 8. Reproducing the results

From the repo root:

```bash
python main.py challenge_sample out        # 1. generate SVGs
python test_compare.py                      # 2. print the results table
python visualize.py overlay letter_H        # 3. (optional) inspect one shape
```

Expected: the table in §5, with `out/*.svg` written and `RecOut` ≥ `RecRef` on every shape.