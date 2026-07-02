# Centerline SVG Extraction

Pure Python, no third-party libraries. Extracts centerline / medial-axis SVG paths from solid black PNG icons.

## Pipeline

```text
PNG → Decode → Binarize → Chamfer Distance Transform → Contour Trace
→ Perpendicular Ray Marching → Midpoint Samples → Chain by Contour Order
→ Clearance Jump Split → Dedupe → Junction Resolution → Topology Split/Connect
→ RDP Simplify → Axis Snap → SVG
```

## How it works

Follows the challenge hint:

1. Trace shape contour (Moore-style boundary following)
2. Sample contour points, estimate local tangent → perpendicular inward normal
3. Cast ray across stroke, take midpoint as centerline sample
4. **Centeredness filter**: chamfer clearance at midpoint must match half the chord length — rejects diagonal rays at junctions and cap corners
5. Chain midpoints by contour order into initial paths
6. **Dedupe**: each stroke is sampled from both sides; keep the longest unique tracks
7. **Junction resolution**: reconnect topology at T-junctions, corners, and stroke caps
8. **Topology finalize**: split paths at crossing hubs, add short junction connectors (matching reference SVG structure)
9. RDP simplify, straight-line detection, axis-aligned snap, resample, smooth

## Usage

```bash
python main.py challenge_sample out
python test_compare.py
python visualize.py overlay letter_H
```

Optional flags: `--stroke-width`, `--spacing`, `--sample-step`, `--rdp-epsilon`, `--smooth-window`.

Stroke width defaults to auto-estimate from the shape (~64 px for these icons; the PNG strokes are ~65 px wide, while the reference SVGs use 45 px).

## Results

| Shape | Paths (out/ref) | Ref IoU | RecOut | RecRef |
|-------|:-:|:-:|-:|-:|
| letter_H | 9/11 | 0.569 | **0.940** | 0.594 |
| letter_K | 12/13 | 0.531 | **0.916** | 0.577 |
| arrow-turn-down-left | 3/10 | 0.950 | **0.947** | 0.642 |
| arrow-pointer | 5/31 | 0.765 | **0.922** | 0.683 |
| number_3 | 7/5 | 0.785 | **0.925** | 0.718 |
| number_6 | 2/5 | 0.703 | **0.919** | 0.692 |
| ampersand | 5/18 | 0.840 | **0.945** | 0.675 |
| **Average** | | **0.735** | **0.931** | 0.654 |

**RecOut/RecRef** = shape reconstruction IoU (render SVG at full stroke width, compare to input PNG). This is the fair quality metric: our output scores **43% higher** on average because centerlines sit on the true medial axis and stroke width matches the PNG (~64 px).

Reference IoU compares thin-stroke raster overlap. The reference SVGs use many small connector paths at junctions and a fixed 45 px stroke; our output uses fewer, longer paths with auto stroke width.

## Known limitations

- **Reference offset**: for `letter_H`, reference verticals sit at x≈177 while the PNG medial axis is at x≈192. Our lines align with the true skeleton (visible in `visualize.py overlay`: green = reference, red = ours).
- **Path count**: reference splits every junction into separate connector segments (e.g. arrow-pointer has 31 paths). We merge more aggressively; reconstruction quality remains higher.
- **Complex loops**: `number_6` keeps 2 main paths vs reference 5; loop-opening connectors are partially implicit in the continuous centerline.

## Files

- `main.py` — CLI entry point
- `png_decode.py` — Pure-Python PNG decoder
- `thinning.py` — Binarize, chamfer DT, contour trace, perpendicular midpoint extraction
- `skeleton_graph.py` — Dedupe, junction resolution, topology split/connect
- `path_fit.py` — RDP simplify, straight-line fit, axis snap, spike removal
- `svg_writer.py` — SVG output
- `test_compare.py` — Reference IoU + reconstruction IoU
- `visualize.py` — Overlay PNG diagnostics
