"""
CLI: PNG -> contour trace -> perpendicular midpoints -> chain -> dedupe
-> junction/cap resolution -> simplify -> SVG.
"""

import os
import sys
import argparse
from png_decode import decode_png
from centerline_extract import binarize, extract_centerline_paths
from skeleton_graph import (split_on_clearance_jumps, dedupe_paths,
                            prune_covered_fragments, resolve_junctions,
                            finalize_topology, merge_collinear_paths,
                            plen, dt_at)
from path_fit import (rdp_simplify, fit_line_if_straight, resample_polyline,
                      smooth_polyline, snap_axis_aligned, remove_spikes)
from svg_writer import write_svg

RDP_EPSILON = 2.0
STRAIGHT_MAX_ERROR = 3.0
SAMPLE_STEP = 2
RESAMPLE_SPACING = 4.0
SMOOTH_WINDOW = 3


def estimate_stroke_width(paths, dt, w, h):
    """Median clearance along centerlines x2 = stroke width of the input shape."""
    clear = sorted(dt_at(dt, x, y, w, h) for p in paths for (x, y) in p)
    if not clear:
        return 45
    # Ridge clearance (upper quartile) tracks full stroke width better than
    # the global median, which is pulled down by junction blobs.
    half = clear[int(len(clear) * 0.72)]
    return max(40, min(80, int(round(2.0 * half))))


def process_one(png_path, svg_path, **kwargs):
    stroke_width = kwargs.get("stroke_width")  # None = auto-estimate
    sample_step = kwargs.get("sample_step", SAMPLE_STEP)
    spacing = kwargs.get("spacing", RESAMPLE_SPACING)
    smooth_window = kwargs.get("smooth_window", SMOOTH_WINDOW)
    rdp_epsilon = kwargs.get("rdp_epsilon", RDP_EPSILON)
    straight_max_error = kwargs.get("straight_max_error", STRAIGHT_MAX_ERROR)

    print(f"Processing: {png_path}")
    w, h, pixels = decode_png(png_path)
    mask = binarize(pixels, w, h)

    paths, dt = extract_centerline_paths(mask, w, h, sample_step)
    paths = split_on_clearance_jumps(paths, dt, w, h)
    paths = dedupe_paths(paths, dt, w, h)
    paths = prune_covered_fragments(paths, dt, w, h)

    if stroke_width is None:
        stroke_width = estimate_stroke_width(paths, dt, w, h)
        print(f"  Estimated stroke width: {stroke_width}")

    paths = resolve_junctions(paths, mask, dt, w, h,
                              cap_clearance=stroke_width / 2.0)
    paths = finalize_topology(paths, mask, dt, w, h, stroke_width)
    paths = merge_collinear_paths(paths)
    paths = prune_covered_fragments(paths, dt, w, h, max_len=stroke_width * 0.75)

    simplified = []
    for p in paths:
        if len(p) < 2 or plen(p) < 8.0:
            continue
        sp = remove_spikes(p)
        sp = rdp_simplify(sp, rdp_epsilon)
        sp = fit_line_if_straight(sp, straight_max_error)
        sp = snap_axis_aligned(sp)
        if len(sp) < 2 or plen(sp) < 6.0:
            continue
        sp = resample_polyline(sp, spacing)
        sp = smooth_polyline(sp, smooth_window)
        sp = snap_axis_aligned(sp)
        if len(sp) >= 2 and plen(sp) >= 6.0:
            simplified.append(sp)

    print(f"  Paths: {len(simplified)}")
    write_svg(svg_path, simplified, w, h, stroke_width)
    print(f"  Written: {svg_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input_dir', nargs='?', default='data/input')
    p.add_argument('output_dir', nargs='?', default='out')
    p.add_argument('--stroke-width', type=int, default=None,
                   help='SVG stroke width (default: auto per shape)')
    p.add_argument('--spacing', type=float, default=RESAMPLE_SPACING)
    p.add_argument('--smooth-window', type=int, default=SMOOTH_WINDOW)
    p.add_argument('--sample-step', type=int, default=SAMPLE_STEP)
    p.add_argument('--rdp-epsilon', type=float, default=RDP_EPSILON)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    pngs = sorted(f for f in os.listdir(args.input_dir) if f.lower().endswith('.png'))
    if not pngs:
        print("No PNGs found in", args.input_dir)
        sys.exit(1)

    print(f"Found {len(pngs)} PNG(s)")
    print("=" * 60)
    for pf in pngs:
        name = os.path.splitext(pf)[0]
        try:
            process_one(
                os.path.join(args.input_dir, pf),
                os.path.join(args.output_dir, name + '.svg'),
                stroke_width=args.stroke_width,
                spacing=args.spacing,
                smooth_window=args.smooth_window,
                sample_step=args.sample_step,
                rdp_epsilon=args.rdp_epsilon,
            )
        except Exception as e:
            print(f"  ERROR {pf}: {e}")
            import traceback
            traceback.print_exc()
    print("=" * 60)
    print("Done!")


if __name__ == '__main__':
    main()
