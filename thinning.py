"""
Core centerline extraction, following the challenge hint:
binary mask -> contour trace -> perpendicular rays across the stroke
-> midpoints -> chain midpoints by contour order.

The key correctness filter is "centeredness": on the true medial axis the
clearance of the midpoint (distance transform value) equals half the chord
length of the perpendicular ray. Rays that cross junctions or cut corners
diagonally violate this and are rejected, so chains break cleanly at
junctions/caps instead of producing hooks and drift.
"""

import math


def binarize(pixels, width, height, threshold=128):
    mask = [bytearray(width) for _ in range(height)]
    for y in range(height):
        row = pixels[y]
        mrow = mask[y]
        for x in range(width):
            mrow[x] = 1 if row[x] < threshold else 0
    return mask


def chamfer_dt(mask, width, height):
    """Two-pass 3-4 chamfer distance transform. Values are in units of
    1/3 pixel (divide by 3 to get approximate euclidean pixels)."""
    INF = float('inf')
    w, h = width, height
    dt = [[INF] * w for _ in range(h)]
    for y in range(h):
        for x in range(w):
            if mask[y][x] == 0:
                dt[y][x] = 0.0
    for y in range(h):
        for x in range(w):
            if mask[y][x] == 0:
                continue
            best = dt[y][x]
            if y > 0:
                if x > 0 and dt[y - 1][x - 1] + 4 < best:
                    best = dt[y - 1][x - 1] + 4
                if dt[y - 1][x] + 3 < best:
                    best = dt[y - 1][x] + 3
                if x + 1 < w and dt[y - 1][x + 1] + 4 < best:
                    best = dt[y - 1][x + 1] + 4
            if x > 0 and dt[y][x - 1] + 3 < best:
                best = dt[y][x - 1] + 3
            dt[y][x] = best
    for y in range(h - 1, -1, -1):
        for x in range(w - 1, -1, -1):
            if mask[y][x] == 0:
                continue
            best = dt[y][x]
            if y + 1 < h:
                if dt[y + 1][x] + 3 < best:
                    best = dt[y + 1][x] + 3
                if x + 1 < w and dt[y + 1][x + 1] + 4 < best:
                    best = dt[y + 1][x + 1] + 4
                if x > 0 and dt[y + 1][x - 1] + 4 < best:
                    best = dt[y + 1][x - 1] + 4
            if x + 1 < w and dt[y][x + 1] + 3 < best:
                best = dt[y][x + 1] + 3
            dt[y][x] = best
    return dt


def _is_boundary(mask, x, y, w, h):
    if mask[y][x] == 0:
        return False
    for dx, dy in [(0, -1), (1, 0), (0, 1), (-1, 0)]:
        nx, ny = x + dx, y + dy
        if nx < 0 or nx >= w or ny < 0 or ny >= h:
            return True
        if mask[ny][nx] == 0:
            return True
    return False


def trace_all_contours(mask, width, height):
    """Moore-style boundary following. Returns all contours (outer + holes)
    as ordered pixel lists."""
    w, h = width, height
    visited = bytearray(w * h)
    dirs = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]
    contours = []
    for y0 in range(h):
        for x0 in range(w):
            if not _is_boundary(mask, x0, y0, w, h) or visited[y0 * w + x0]:
                continue
            start = (x0, y0)
            contour = [start]
            visited[start[1] * w + start[0]] = 1
            cx, cy = start
            px, py = -1, -1
            ss = 0
            for _ in range(w * h):
                found = False
                for d in range(8):
                    di = (ss + d) % 8
                    dx, dy = dirs[di]
                    nx, ny = cx + dx, cy + dy
                    if 0 <= nx < w and 0 <= ny < h and _is_boundary(mask, nx, ny, w, h):
                        if (nx, ny) == (px, py):
                            continue
                        if (nx, ny) == start:
                            if len(contour) > 2:
                                found = True
                                break
                            continue
                        contour.append((nx, ny))
                        visited[ny * w + nx] = 1
                        px, py = cx, cy
                        cx, cy = nx, ny
                        ss = (di + 5) % 8
                        found = True
                        break
                if not found or (cx, cy) == start:
                    break
            if len(contour) >= 10:
                contours.append(contour)
    return contours


def _tangent(contour, i, window):
    n = len(contour)
    im = (i - window) % n
    ip = (i + window) % n
    dx = contour[ip][0] - contour[im][0]
    dy = contour[ip][1] - contour[im][1]
    mag = math.hypot(dx, dy)
    return (dx / mag, dy / mag) if mag > 1e-6 else (0.0, 0.0)


def _inward_normal(mask, x, y, tx, ty, w, h):
    """Of the two perpendiculars to the tangent, pick the one pointing
    into the shape (deeper probe stays inside longer)."""
    n1x, n1y = -ty, tx
    n2x, n2y = ty, -tx

    def _steps(nx, ny):
        fx, fy = float(x) + nx * 2, float(y) + ny * 2
        for s in range(15):
            xi, yi = int(fx), int(fy)
            if xi < 0 or xi >= w or yi < 0 or yi >= h:
                return s
            if mask[yi][xi] == 0:
                return s
            fx += nx
            fy += ny
        return 15

    return (n1x, n1y) if _steps(n1x, n1y) >= _steps(n2x, n2y) else (n2x, n2y)


def extract_centerline_paths(mask, width, height, sample_step=2):
    """Hint pipeline: contour -> perpendicular chord -> midpoint -> chain.

    Returns (paths, dt). Each path is a polyline [(x, y), ...] in pixel
    coordinates; closed loops have path[0] == path[-1]. Each stroke will
    typically appear TWICE (once per side) - see skeleton_graph.dedupe_paths.
    """
    w, h = width, height
    dt = chamfer_dt(mask, w, h)
    contours = trace_all_contours(mask, w, h)
    print(f"  Contours: {len(contours)} ({sum(len(c) for c in contours)} pts)")

    max_igap = sample_step * 8
    max_sgap = 25.0
    paths = []

    for contour in contours:
        n = len(contour)
        window = max(4, min(12, n // 100))
        samples = []  # (contour_index, mx, my)

        for idx in range(0, n, sample_step):
            x, y = contour[idx]
            tx, ty = _tangent(contour, idx, window)
            if tx == 0 and ty == 0:
                continue
            nx, ny = _inward_normal(mask, x, y, tx, ty, w, h)

            # March the perpendicular chord across the stroke.
            fx, fy = float(x), float(y)
            hit = False
            for _ in range(300):
                fx += nx
                fy += ny
                xi, yi = int(fx), int(fy)
                if xi < 0 or xi >= w or yi < 0 or yi >= h:
                    break
                if mask[yi][xi] == 0:
                    hit = True
                    break
            if not hit:
                continue

            ox, oy = fx - nx, fy - ny  # last point still inside
            mx = (x + ox) / 2.0
            my = (y + oy) / 2.0
            chord = math.hypot(ox - x, oy - y)
            if chord < 4.0:
                continue
            mxi, myi = int(mx), int(my)
            if not (0 <= mxi < w and 0 <= myi < h and mask[myi][mxi]):
                continue

            # Centeredness filter: clearance(midpoint) must match half chord.
            # Strict tolerance rejects diagonal chords at cap corners and
            # chords reaching into junction blobs; each stroke is sampled
            # from both sides, so losing one side's samples is harmless.
            half = chord / 2.0
            clearance = dt[myi][mxi] / 3.0
            if abs(half - clearance) > max(2.5, 0.15 * half):
                continue

            # Refine: the chord crosses the medial axis where clearance
            # peaks. The plain midpoint is biased on curved strokes (the
            # chord is perpendicular to one side only); slide to the max.
            best_c = clearance
            bx, by = mx, my
            steps = int(chord)
            for s in range(steps + 1):
                t = s / max(steps, 1)
                px_ = x + t * (ox - x)
                py_ = y + t * (oy - y)
                pxi, pyi = int(px_), int(py_)
                if 0 <= pxi < w and 0 <= pyi < h:
                    c = dt[pyi][pxi] / 3.0
                    if c > best_c:
                        best_c = c
                        bx, by = px_, py_

            samples.append((idx, bx, by, best_c))

        if len(samples) < 3:
            continue

        # Chain consecutive samples; break on contour-index or spatial gaps.
        chains = [[samples[0]]]
        for k in range(1, len(samples)):
            pidx, px_, py_, _ = chains[-1][-1]
            cidx, cx_, cy_, _ = samples[k]
            if (cidx - pidx) > max_igap or math.hypot(cx_ - px_, cy_ - py_) > max_sgap:
                chains.append([samples[k]])
            else:
                chains[-1].append(samples[k])

        # The contour is a closed loop: the last chain may continue into the
        # first one across the index wrap-around.
        def _wrap_connects(tail, head):
            igap = head[0][0] + n - tail[-1][0]
            sgap = math.hypot(head[0][1] - tail[-1][1], head[0][2] - tail[-1][2])
            return igap <= max_igap and sgap <= max_sgap

        closed = False
        if len(chains) > 1 and _wrap_connects(chains[-1], chains[0]):
            chains[0] = chains.pop() + chains[0]
        elif len(chains) == 1 and _wrap_connects(chains[0], chains[0]):
            closed = True

        for ch in chains:
            if len(ch) < 3:
                continue
            pts = [(mx, my) for _, mx, my, _ in ch]
            if closed:
                pts.append(pts[0])
            paths.append(pts)

    print(f"  Raw midpoint tracks: {len(paths)}")
    return paths, dt
