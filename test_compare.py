"""
Test suite with two kinds of metrics, pure Python:

1. Reference similarity: rasterize our SVG and the reference SVG as thin
   strokes and compare (IoU/Dice). NOTE: the reference centerlines are
   systematically offset from the true medial axis of the PNGs for some
   shapes (e.g. letter_H verticals), so 1.0 is not reachable nor desirable.
2. Shape reconstruction (reference-independent): rasterize the SVG at the
   full stroke width and compare against the input PNG mask. This measures
   how well the centerlines re-draw the original shape and is the fairer
   quality signal. It is also computed for the reference SVGs as a baseline.
"""

import os
import sys
import math
import xml.etree.ElementTree as ET


def parse_svg_paths(svg_path):
    """Extract all path data strings from an SVG file."""
    try:
        tree = ET.parse(svg_path)
        root = tree.getroot()
        ns = {'svg': 'http://www.w3.org/2000/svg'}
        paths = []
        for path_elem in root.findall('.//{http://www.w3.org/2000/svg}path'):
            d = path_elem.get('d', '')
            paths.append(d)
        # Also try without namespace
        if not paths:
            for path_elem in root.findall('.//path'):
                d = path_elem.get('d', '')
                paths.append(d)
        return paths
    except Exception as e:
        print(f"  Parse error: {e}")
        return []


def path_to_subpaths(d_string, step=2.0):
    """
    Convert SVG path d-string to a list of polylines (one per M command).
    Handles M (move) and L (line) commands. Keeping subpaths separate is
    important: rasterizing a flat point list would draw phantom segments
    between the end of one subpath and the start of the next.
    """
    subpaths = []
    points = []
    # Tokenize: split on command letters, keep letters as tokens
    d_string = d_string.replace(',', ' ').strip()
    tokens = []
    current = ''
    for ch in d_string:
        if ch.isalpha():
            if current:
                tokens.append(current.strip())
                current = ''
            tokens.append(ch)
        else:
            current += ch
    if current:
        tokens.append(current.strip())

    x, y = 0.0, 0.0
    i = 0
    while i < len(tokens):
        cmd = tokens[i]
        i += 1
        if cmd in ('M', 'L'):
            # Collect all numbers until next command
            nums = []
            while i < len(tokens) and (tokens[i][0].isdigit() or tokens[i][0] in '-.'):
                # Split the token by whitespace (it may contain multiple numbers)
                for num_str in tokens[i].split():
                    if num_str:
                        nums.append(float(num_str))
                i += 1
            for j in range(0, len(nums), 2):
                if j + 1 >= len(nums):
                    break
                x, y = nums[j], nums[j+1]
                if cmd == 'M' and j == 0:
                    if len(points) >= 2:
                        subpaths.append(points)
                    points = [(x, y)]
                elif not points:
                    points.append((x, y))
                else:
                    # Line: interpolate from previous point
                    px, py = points[-1]
                    dist = math.sqrt((x-px)**2 + (y-py)**2)
                    if dist > step:
                        segs = max(1, int(dist / step))
                        for s in range(1, segs):
                            t = s / segs
                            points.append((px + t*(x-px), py + t*(y-py)))
                    points.append((x, y))
        else:
            # Skip unknown commands and their numbers
            while i < len(tokens) and (tokens[i][0].isdigit() or tokens[i][0] in '-.'):
                i += 1

    if len(points) >= 1:
        subpaths.append(points)
    return subpaths


def path_to_points(d_string, step=2.0):
    """Flattened point list (for endpoint sampling / rendering overlays)."""
    pts = []
    for sp in path_to_subpaths(d_string, step):
        pts.extend(sp)
    return pts


def rasterize_polylines(polylines, width=256, height=256, stroke_width=11):
    """
    Rasterize a list of polylines to a binary grid.
    Uses Bresenham-like line drawing with stroke width.
    """
    grid = bytearray(width * height)
    for poly in polylines:
        _draw_polyline(grid, poly, width, height, stroke_width)
    return grid


def _draw_polyline(grid, points, width, height, stroke_width):
    sw = stroke_width

    if len(points) == 1:
        points = [points[0], points[0]]

    for i in range(1, len(points)):
        x0, y0 = points[i-1]
        x1, y1 = points[i]

        # Scale from 1024 to grid size
        sx0 = x0 * width / 1024.0
        sy0 = y0 * height / 1024.0
        sx1 = x1 * width / 1024.0
        sy1 = y1 * height / 1024.0

        # Draw thick line (simplified: fill rect around each segment point)
        dist = math.sqrt((sx1-sx0)**2 + (sy1-sy0)**2)
        steps = max(1, int(dist))
        for s in range(steps + 1):
            t = s / max(steps, 1)
            cx = sx0 + t * (sx1 - sx0)
            cy = sy0 + t * (sy1 - sy0)
            # Draw circle of radius sw/2
            r = sw / 2.0
            for dy in range(-int(r)-1, int(r)+2):
                for dx in range(-int(r)-1, int(r)+2):
                    if dx*dx + dy*dy <= r*r:
                        gx = int(cx + dx)
                        gy = int(cy + dy)
                        if 0 <= gx < width and 0 <= gy < height:
                            grid[gy * width + gx] = 1

    return grid


def compute_similarity(grid1, grid2):
    """Compute IoU and Dice score between two binary grids."""
    total = len(grid1)
    intersection = sum(1 for i in range(total) if grid1[i] and grid2[i])
    union = sum(1 for i in range(total) if grid1[i] or grid2[i])
    both = sum(grid1) + sum(grid2)

    if union == 0:
        iou = 1.0
    else:
        iou = intersection / union

    if both == 0:
        dice = 1.0
    else:
        dice = 2 * intersection / both

    return iou, dice


def _mask_grid(name, size=256):
    """Downsample the input PNG's binary mask to the metric grid."""
    from png_decode import decode_png
    from thinning import binarize
    w, h, px = decode_png(f"challenge_sample/{name}.png")
    mask = binarize(px, w, h)
    grid = bytearray(size * size)
    for y in range(size):
        row = mask[int(y * h / size)]
        base = y * size
        for x in range(size):
            if row[int(x * w / size)]:
                grid[base + x] = 1
    return grid


def _svg_stroke_width(svg_path, default=45.0):
    try:
        tree = ET.parse(svg_path)
        for el in tree.getroot().iter():
            sw = el.get('stroke-width')
            if sw:
                return float(sw)
    except Exception:
        pass
    return default


def _svg_full_stroke_grid(svg_path, size=256):
    stroke = _svg_stroke_width(svg_path)
    polys = []
    for d in parse_svg_paths(svg_path):
        polys.extend(path_to_subpaths(d, step=1.0))
    return rasterize_polylines(polys, size, size, stroke * size / 1024.0)


def reconstruction_iou(name, svg_path):
    """Reference-independent: how well does the SVG, drawn at full stroke
    width, reproduce the input shape?"""
    iou, _ = compute_similarity(_mask_grid(name), _svg_full_stroke_grid(svg_path))
    return iou


def test_one_shape(name):
    """Compare one shape's output against reference."""
    ref_path = f"challenge_sample_results/{name}.svg"
    out_path = f"out/{name}.svg"

    if not os.path.exists(ref_path):
        return {"error": f"Reference not found: {ref_path}"}
    if not os.path.exists(out_path):
        return {"error": f"Output not found: {out_path}"}

    ref_paths = parse_svg_paths(ref_path)
    out_paths = parse_svg_paths(out_path)

    # Path counts
    ref_count = len(ref_paths)
    out_count = len(out_paths)

    # Sample points per subpath (no phantom segments between paths)
    ref_polys = []
    for d in ref_paths:
        ref_polys.extend(path_to_subpaths(d, step=1.0))
    out_polys = []
    for d in out_paths:
        out_polys.extend(path_to_subpaths(d, step=1.0))
    ref_points = [p for poly in ref_polys for p in poly]
    out_points = [p for poly in out_polys for p in poly]

    # Rasterize
    ref_grid = rasterize_polylines(ref_polys, width=256, height=256, stroke_width=11)
    out_grid = rasterize_polylines(out_polys, width=256, height=256, stroke_width=11)

    iou, dice = compute_similarity(ref_grid, out_grid)

    # Endpoint analysis
    def get_endpoints(paths_list):
        eps = []
        for d in paths_list:
            pts = path_to_points(d, step=10.0)
            if pts:
                eps.append(pts[0])
                eps.append(pts[-1])
        return eps

    ref_eps = get_endpoints(ref_paths)
    out_eps = get_endpoints(out_paths)

    # Count how many output endpoints are within 10px of reference endpoints
    matched = 0
    for (ox, oy) in out_eps:
        for (rx, ry) in ref_eps:
            if math.sqrt((ox-rx)**2 + (oy-ry)**2) <= 15:
                matched += 1
                break

    endpoint_score = matched / max(len(out_eps), 1) if out_eps else 1.0

    return {
        "ref_paths": ref_count,
        "out_paths": out_count,
        "path_diff": out_count - ref_count,
        "ref_points": len(ref_points),
        "out_points": len(out_points),
        "iou": round(iou, 4),
        "dice": round(dice, 4),
        "endpoint_match": round(endpoint_score, 4),
        "recon_out": round(reconstruction_iou(name, out_path), 4),
        "recon_ref": round(reconstruction_iou(name, ref_path), 4),
    }


def run_all_tests():
    shapes = [
        "letter_H", "letter_K", "arrow-turn-down-left",
        "arrow-pointer", "number_3", "number_6", "ampersand",
    ]

    print("=" * 96)
    print("CENTERLINE SVG COMPARISON TEST")
    print("(RecOut/RecRef = shape reconstruction IoU of our output / of the reference)")
    print("=" * 96)
    print(f"{'Shape':<25} {'Ref':>5} {'Out':>5} {'Diff':>5} {'IoU':>8} {'Dice':>8} "
          f"{'EndPt':>8} {'RecOut':>8} {'RecRef':>8}")
    print("-" * 96)

    results = {}
    total_iou = 0
    total_dice = 0
    total_ro = 0
    total_rr = 0
    n = 0

    for shape in shapes:
        r = test_one_shape(shape)
        results[shape] = r
        if "error" in r:
            print(f"{shape:<25} ERROR: {r['error']}")
        else:
            total_iou += r['iou']
            total_dice += r['dice']
            total_ro += r['recon_out']
            total_rr += r['recon_ref']
            n += 1
            flag = "OK" if r['recon_out'] >= r['recon_ref'] - 0.02 else "??"
            print(f"{shape:<25} {r['ref_paths']:>5} {r['out_paths']:>5} "
                  f"{r['path_diff']:>+5} {r['iou']:>8.3f} {r['dice']:>8.3f} "
                  f"{r['endpoint_match']:>8.3f} {r['recon_out']:>8.3f} "
                  f"{r['recon_ref']:>8.3f} {flag}")

    print("-" * 96)
    if n > 0:
        print(f"{'AVERAGE':<25} {'':>5} {'':>5} {'':>5} "
              f"{total_iou/n:>8.3f} {total_dice/n:>8.3f} {'':>8} "
              f"{total_ro/n:>8.3f} {total_rr/n:>8.3f}")
    print("=" * 96)

    # Per-shape details
    print("\nDETAILS:")
    for shape in shapes:
        r = results.get(shape, {})
        if "error" not in r:
            print(f"  {shape}: ref={r['ref_paths']}paths/{r['ref_points']}pts, "
                  f"out={r['out_paths']}paths/{r['out_points']}pts, "
                  f"IoU={r['iou']:.3f}, endPt={r['endpoint_match']:.3f}")

    return results


if __name__ == '__main__':
    run_all_tests()
