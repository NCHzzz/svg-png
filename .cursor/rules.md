# Project Rules

## 1. Follow Problem Hints
When a problem statement provides a hint about the approach (e.g., "rough direction: convert to binary mask → trace contour → perpendicular → midpoint → connect"), **implement that approach first**. Do not reach for alternative algorithms (e.g., Zhang-Suen thinning) unless the hint-based approach proves intractable after genuine attempt.

## 2. Test-Driven Development
Always maintain `test_compare.py` or equivalent to compare output against reference. Run after every change.

## 3. Per-Shape Tuning
When algorithm parameters need to vary by input shape, use a config dict rather than one-size-fits-all defaults.

## 4. Pipeline Files
- `png_decode.py` — Decoding
- `thinning.py` — Core algorithm (binarize + transform)
- `skeleton_graph.py` — Post-processing (resample, smooth, connectors)
- `svg_writer.py` — SVG output
- `main.py` — CLI entry point
- `test_compare.py` — Comparison against reference
- `compare.html` — Visual side-by-side comparison
- `README.md` — Documentation

## 5. Pure Python Standard Library Only
No third-party libraries (numpy, PIL, opencv, etc.).
