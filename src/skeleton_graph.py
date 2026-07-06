"""
Post-processing of raw midpoint tracks:

1. dedupe_paths     - every stroke is sampled from both of its sides, giving
                      two nearly identical midpoint tracks; keep one.
2. resolve_junctions - reconnect the topology that the centeredness filter
                      intentionally broke:
                      a) endpoint clusters: 2 endpoints -> corner join at the
                         intersection of their tangents (elbows, sharp tips);
                         3+ endpoints -> shared junction point (Y/T hubs);
                      b) T-connections: walk a free endpoint forward until it
                         lands on a path passing by (e.g. H crossbar -> stem);
                      c) cap extension: push remaining free ends into the
                         stroke tips until clearance drops to the drawn
                         stroke's cap radius, so round caps fill the tip.
"""

import math
from collections import defaultdict
from spatial_hash import PointGrid

# ── Pixel-distance thresholds (tuned for ~256–1024 px icons) ────────────

# dedupe
_DEDUPE_MIN_LEN = 12.0          # min arc length to keep a segment
_DEDUPE_MAX_RADIUS = 13.0       # max radius for duplicate-point query
_DEDUPE_MIN_RADIUS = 5.0        # min radius for duplicate-point query

# prune
_PRUNE_MAX_LEN = 60.0           # paths shorter than this may be pruned
_PRUNE_COVERAGE_MIN_R = 15.0    # min radius when testing coverage

# junction clustering
_CLOSED_COINCIDE = 3.0          # max start→end gap to declare a loop
_END_DIR_SPAN = 18.0            # arc length for tangent direction estimation
_CLUSTER_GAP_BASE = 28.0        # base max gap for endpoint union-find
_CORNER_TMAX_BASE = 40.0        # base forward limit for corner ray intersection
_HUB_TMAX_BASE = 40.0           # base forward limit for hub ray intersection

# cap walk (T-connection / cap extension)
_CAP_CLEARANCE_DEFAULT = 22.5   # default half-stroke clearance
_T_TOUCH_RADIUS = 2.5           # if endpoint near another path, skip walk
_T_HIT_RADIUS = 4.0             # distance to register T-connection hit
_CAP_EXTEND_MIN = 2.0           # min walk distance to apply cap extension

# topology: crossing detection, connectors, bridges
_MERGE_POINT_RADIUS = 18.0      # radius for clustering hit points
_PATH_CROSSING_RADIUS = 22.0    # default radius for hub/endpoint detection
_CONNECTOR_MAX_GAP = 42.0       # max gap for junction connector segs
_CONNECTOR_MIN_GAP = 1.5        # min gap for junction connector segs
_BRIDGE_MIN_GAP = 3.0           # min gap for endpoint bridge
_BRIDGE_MAX_GAP = 22.0          # max gap for endpoint bridge

# merge collinear
_MERGE_MAX_GAP = 28.0           # max gap between collinear endpoints


def plen(path):
    """Total arc length of a polyline."""
    return sum(math.hypot(path[i][0] - path[i - 1][0], path[i][1] - path[i - 1][1])
               for i in range(1, len(path)))


def dt_at(dt, x, y, w, h):
    """Distance-transform value at (x,y) in pixel units (chamfer ÷ 3)."""
    xi, yi = int(x), int(y)
    if 0 <= xi < w and 0 <= yi < h:
        return dt[yi][xi] / 3.0
    return 0.0


def _inside(mask, x, y, w, h):
    """True if (x,y) is inside the foreground region."""
    xi, yi = int(x), int(y)
    return 0 <= xi < w and 0 <= yi < h and mask[yi][xi] == 1


def _segment_inside_mask(mask, w, h, ax, ay, bx, by):
    """True if every sampled point on segment a→b is inside foreground."""
    steps = max(1, int(math.hypot(bx - ax, by - ay)))
    for s in range(steps + 1):
        t = s / steps
        if not _inside(mask, ax + t * (bx - ax), ay + t * (by - ay), w, h):
            return False
    return True


def split_on_clearance_jumps(paths, dt, w, h, up=1.30, down=0.72):
    """Near a junction the perpendicular chord escapes through the junction
    opening: the midpoint drifts into the junction blob and its clearance
    jumps well above the stroke's typical half-width. Split chains at such
    outliers so junction resolution can rebuild the topology cleanly."""
    result = []
    for path in paths:
        clear = [dt_at(dt, x, y, w, h) for (x, y) in path]
        cur = []
        window = []
        for i, pt in enumerate(path):
            if window:
                med = sorted(window)[len(window) // 2]
                if clear[i] > med * up or clear[i] < med * down:
                    if len(cur) >= 3:
                        result.append(cur)
                    cur = []
                    window = []
                    continue
            cur.append(pt)
            window.append(clear[i])
            if len(window) > 15:
                window.pop(0)
        if len(cur) >= 3:
            result.append(cur)
    return result


def dedupe_paths(paths, dt, w, h, min_len=_DEDUPE_MIN_LEN):
    """Every stroke yields two midpoint tracks (one per contour side), and a
    track can additionally fold back on itself where the contour wraps a
    rounded cap. Keep longest tracks first; drop points that duplicate
    already-kept geometry, or the track's own geometry from a distant
    position along the same track (fold-back)."""
    order = sorted(range(len(paths)), key=lambda i: -plen(paths[i]))
    kept = PointGrid()
    result = []

    for pi in order:
        path = paths[pi]
        own = PointGrid()
        arc = 0.0
        dup = []
        for k, (x, y) in enumerate(path):
            if k > 0:
                arc += math.hypot(x - path[k - 1][0], y - path[k - 1][1])
            clearance = dt_at(dt, x, y, w, h)
            r = min(_DEDUPE_MAX_RADIUS, max(_DEDUPE_MIN_RADIUS, 0.4 * clearance))
            is_dup = kept.near(x, y, r) is not None
            if not is_dup:
                # Self fold-back: same location reached much earlier in arc
                # (track wrapped a cap and doubled back along the far side).
                guard = 2.0 * r + 6.0
                cur_arc = arc
                if own.any_match(x, y, r, lambda a: cur_arc - a > guard):
                    is_dup = True
            dup.append(is_dup)
            own.add(x, y, arc)

        # Split into runs of unique points.
        runs = []
        s = None
        for i, f in enumerate(dup):
            if not f and s is None:
                s = i
            elif f and s is not None:
                runs.append((s, i))
                s = None
        if s is not None:
            runs.append((s, len(path)))

        for (a, b) in runs:
            seg = path[a:b]
            if len(seg) >= 3 and plen(seg) >= min_len:
                for (x, y) in seg:
                    kept.add(x, y)
                result.append(seg)

    print(f"  After dedupe: {len(result)} paths")
    return result


def prune_covered_fragments(paths, dt, w, h, max_len=_PRUNE_MAX_LEN):
    """Drop short leftover tracks (cap-wrap remnants, junction-blob branch
    stubs) that lie entirely within the stroke area already covered by
    longer paths."""
    if len(paths) < 2:
        return paths
    lengths = [plen(p) for p in paths]
    grids = []
    for p in paths:
        g = PointGrid()
        for (x, y) in p:
            g.add(x, y)
        grids.append(g)

    keep = [True] * len(paths)
    for pi, path in enumerate(paths):
        if lengths[pi] >= max_len:
            continue
        covered = True
        for (x, y) in path:
            r = max(1.1 * dt_at(dt, x, y, w, h), _PRUNE_COVERAGE_MIN_R)
            hit = False
            for pj in range(len(paths)):
                if pj == pi or not keep[pj] or lengths[pj] <= lengths[pi]:
                    continue
                if grids[pj].near(x, y, r) is not None:
                    hit = True
                    break
            if not hit:
                covered = False
                break
        if covered:
            keep[pi] = False

    pruned = sum(1 for k in keep if not k)
    if pruned:
        print(f"  Pruned {pruned} covered fragments")
    return [p for i, p in enumerate(paths) if keep[i]]


def _is_closed(path):
    """True if path's first and last points are nearly coincident (≈loop)."""
    return len(path) >= 4 and math.hypot(path[0][0] - path[-1][0],
                                         path[0][1] - path[-1][1]) < _CLOSED_COINCIDE


def _end_dir(path, side, span=_END_DIR_SPAN):
    """Unit direction pointing OUT of the path at the given end, measured
    over `span` pixels of arc length (robust to sub-pixel point noise)."""
    if side == 0:
        bx, by = path[0]
        idx = range(1, len(path))
    else:
        bx, by = path[-1]
        idx = range(len(path) - 2, -1, -1)
    ax, ay = bx, by
    arc = 0.0
    px, py = bx, by
    for k in idx:
        x, y = path[k]
        arc += math.hypot(x - px, y - py)
        px, py = x, y
        ax, ay = x, y
        if arc >= span:
            break
    dx, dy = bx - ax, by - ay
    m = math.hypot(dx, dy)
    return (dx / m, dy / m) if m > 1e-9 else (0.0, 0.0)


def _append_end(path, side, pt):
    """Append pt to start (side=0) or end (side=1) of path in-place."""
    if side == 0:
        path.insert(0, pt)
    else:
        path.append(pt)


def _ray_intersection(p1, d1, p2, d2):
    """Intersection of ray p1+t*d1 and ray p2+t*d2. Returns (t1, t2, x, y) or None."""
    det = d1[0] * (-d2[1]) - (-d2[0]) * d1[1]
    if abs(det) < 1e-6:
        return None
    rx, ry = p2[0] - p1[0], p2[1] - p1[1]
    t1 = (rx * (-d2[1]) - (-d2[0]) * ry) / det
    t2 = (d1[0] * ry - rx * d1[1]) / det
    return (t1, t2, p1[0] + t1 * d1[0], p1[1] + t1 * d1[1])


def resolve_junctions(paths, mask, dt, w, h, cap_clearance=_CAP_CLEARANCE_DEFAULT):
    """Reconnect topology at junctions, corners and caps. Mutates copies."""
    paths = [list(p) for p in paths]

    # ---- collect free endpoints -------------------------------------------
    endpoints = []  # (path_idx, side, x, y)
    for pi, p in enumerate(paths):
        if len(p) < 2 or _is_closed(p):
            continue
        endpoints.append([pi, 0, p[0][0], p[0][1]])
        endpoints.append([pi, 1, p[-1][0], p[-1][1]])

    resolved = set()  # endpoint list indices already handled

    # ---- pass A: endpoint clusters (corners and junction hubs) ------------
    m = len(endpoints)
    parent = list(range(m))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(m):
        for j in range(i + 1, m):
            xi, yi = endpoints[i][2], endpoints[i][3]
            xj, yj = endpoints[j][2], endpoints[j][3]
            di = dt_at(dt, xi, yi, w, h)
            dj = dt_at(dt, xj, yj, w, h)
            gap = math.hypot(xi - xj, yi - yj)
            if gap < max(_CLUSTER_GAP_BASE, 1.1 * (di + dj)):
                if _segment_inside_mask(mask, w, h, xi, yi, xj, yj):
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        parent[ra] = rb

    clusters = defaultdict(list)
    for i in range(m):
        clusters[find(i)].append(i)

    def _facing_score(i, j):
        """How well two endpoints continue into each other: both must point
        toward the other, and their directions must be anti-parallel."""
        pi, si, xi, yi = endpoints[i]
        pj, sj, xj, yj = endpoints[j]
        d1 = _end_dir(paths[pi], si)
        d2 = _end_dir(paths[pj], sj)
        gap = math.hypot(xj - xi, yj - yi)
        if gap < 1e-9:
            return 1.0
        ux, uy = (xj - xi) / gap, (yj - yi) / gap
        toward1 = d1[0] * ux + d1[1] * uy
        toward2 = -(d2[0] * ux + d2[1] * uy)
        anti = -(d1[0] * d2[0] + d1[1] * d2[1])
        if toward1 < 0.4 or toward2 < 0.4:
            return -1.0
        return min(anti, toward1, toward2)

    def _corner_join(i, j):
        """Join two endpoints at the intersection of their tangents (falls
        back to their midpoint), e.g. elbows, sharp tips, loop closures."""
        pi, si, xi, yi = endpoints[i]
        pj, sj, xj, yj = endpoints[j]
        d1 = _end_dir(paths[pi], si)
        d2 = _end_dir(paths[pj], sj)
        gap = math.hypot(xi - xj, yi - yj)
        corner = None
        hit = _ray_intersection((xi, yi), d1, (xj, yj), d2)
        if hit is not None:
            t1, t2, cx, cy = hit
            tmax = 2.5 * gap + _CORNER_TMAX_BASE
            if -4.0 <= t1 <= tmax and -4.0 <= t2 <= tmax and \
                    _inside(mask, cx, cy, w, h) and \
                    _segment_inside_mask(mask, w, h, xi, yi, cx, cy) and \
                    _segment_inside_mask(mask, w, h, xj, yj, cx, cy):
                corner = (cx, cy)
        if corner is None:
            corner = ((xi + xj) / 2.0, (yi + yj) / 2.0)
        _append_end(paths[pi], si, corner)
        _append_end(paths[pj], sj, corner)
        resolved.add(i)
        resolved.add(j)

    def _hub_join(members):
        """Connect endpoints to a shared junction point. The hub is placed
        where the endpoint tangent rays intersect (that is where the branch
        centerlines actually meet, e.g. deep inside an arrowhead vertex);
        falls back to the endpoint centroid. Endpoints that do not point
        toward the hub are left for pass B."""
        pts = []
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                ia, ib = members[a], members[b]
                pa, sa, xa, ya = endpoints[ia]
                pb, sb, xb, yb = endpoints[ib]
                da = _end_dir(paths[pa], sa)
                db = _end_dir(paths[pb], sb)
                gap = math.hypot(xb - xa, yb - ya)
                hit = _ray_intersection((xa, ya), da, (xb, yb), db)
                if hit is None:
                    continue
                t1, t2, ix, iy = hit
                tmax = 2.0 * gap + _HUB_TMAX_BASE
                if -4.0 <= t1 <= tmax and -4.0 <= t2 <= tmax and \
                        _inside(mask, ix, iy, w, h):
                    pts.append((ix, iy))
        if pts:
            cx = sum(p[0] for p in pts) / len(pts)
            cy = sum(p[1] for p in pts) / len(pts)
        else:
            cx = sum(endpoints[i][2] for i in members) / len(members)
            cy = sum(endpoints[i][3] for i in members) / len(members)
        joined = 0
        for i in members:
            pi, si, xi, yi = endpoints[i]
            gx, gy = cx - xi, cy - yi
            gm = math.hypot(gx, gy)
            if gm > 1e-9:
                d = _end_dir(paths[pi], si)
                if (d[0] * gx + d[1] * gy) / gm < 0.0:
                    continue  # points away from the hub
            if _segment_inside_mask(mask, w, h, xi, yi, cx, cy):
                _append_end(paths[pi], si, (cx, cy))
                resolved.add(i)
                joined += 1
        return joined

    corner_joins = 0
    through_joins = 0
    hub_joins = 0
    for members in clusters.values():
        if len(members) < 2:
            continue
        if len(members) == 2:
            _corner_join(members[0], members[1])
            corner_joins += 1
            continue

        # 3+ endpoints: a stroke passing straight through the junction was
        # split into two collinear pieces. Reconnect the best facing pairs
        # through the junction, then join remaining branch ends to a hub.
        pending = list(members)
        while len(pending) >= 2:
            best = None
            for a in range(len(pending)):
                for b in range(a + 1, len(pending)):
                    s = _facing_score(pending[a], pending[b])
                    if s > 0.75 and (best is None or s > best[0]):
                        best = (s, pending[a], pending[b])
            if best is None:
                break
            _, i, j = best
            pi, si, xi, yi = endpoints[i]
            pj, sj, xj, yj = endpoints[j]
            mid = ((xi + xj) / 2.0, (yi + yj) / 2.0)
            _append_end(paths[pi], si, mid)
            _append_end(paths[pj], sj, mid)
            resolved.add(i)
            resolved.add(j)
            pending.remove(i)
            pending.remove(j)
            through_joins += 1
        # 1-2 leftovers are usually branch ends beside a through-stroke;
        # the T-connect walk (pass B) attaches them more faithfully than
        # a hub point would.
        if len(pending) >= 3:
            if _hub_join(pending):
                hub_joins += 1

    # ---- pass B: T-connections and cap extension --------------------------
    # Registry of all current path points for T-hit tests.
    registry = PointGrid()
    for pi, p in enumerate(paths):
        for (x, y) in p:
            registry.add(x, y, pi)

    t_connects = 0
    cap_extends = 0
    for ei, (pi, si, x, y) in enumerate(endpoints):
        if ei in resolved:
            continue
        dx, dy = _end_dir(paths[pi], si)
        if dx == 0 and dy == 0:
            continue
        d_end = dt_at(dt, x, y, w, h)
        if d_end < 0.85 * cap_clearance:
            # Tapering terminal (stroke narrows toward the tip): follow it
            # deeper so the drawn cap reaches the tip of the shape.
            stop_dt = max(0.55 * d_end, 6.0)
        else:
            # Flat/round cap of a full-width stroke: stop half a drawn
            # stroke-width from the boundary so the round cap lands flush.
            stop_dt = min(cap_clearance, 0.9 * d_end)
        max_walk = min(100.0, 4.0 * max(d_end, 10.0))

        # Already touching another path?
        if registry.near(x, y, _T_TOUCH_RADIUS, exclude_tag=pi) is not None:
            continue

        cx, cy = x, y
        cap_trail = []  # ridge positions with clearance >= stop_dt
        walked = 0.0
        connected = False
        while walked < max_walk:
            # Steer along the clearance ridge: of straight-ahead and two
            # slightly turned steps, take the one with max clearance. This
            # lets the walk follow curved tapering terminals.
            nx, ny = -dy, dx
            best = None
            for off in (0.0, 0.7, -0.7):
                px_ = cx + dx + nx * off
                py_ = cy + dy + ny * off
                c = dt_at(dt, px_, py_, w, h)
                if best is None or c > best[0]:
                    best = (c, px_, py_)
            step_c, sx_, sy_ = best
            m = math.hypot(sx_ - cx, sy_ - cy)
            dx, dy = (sx_ - cx) / m, (sy_ - cy) / m
            cx, cy = sx_, sy_
            walked += m
            if not _inside(mask, cx, cy, w, h):
                break
            hit = registry.near(cx, cy, _T_HIT_RADIUS, exclude_tag=pi)
            if hit is not None:
                _append_end(paths[pi], si, (cx, cy))
                registry.add(cx, cy, pi)
                t_connects += 1
                connected = True
                break
            if step_c >= stop_dt:
                cap_trail.append((cx, cy))
            else:
                break
        if not connected and cap_trail and \
                math.hypot(cap_trail[-1][0] - x, cap_trail[-1][1] - y) >= _CAP_EXTEND_MIN:
            for pt in cap_trail:
                _append_end(paths[pi], si, pt)
                registry.add(pt[0], pt[1], pi)
            cap_extends += 1

    print(f"  Junctions: {corner_joins} corner joins, {through_joins} through joins, "
          f"{hub_joins} hubs, {t_connects} T-connects, {cap_extends} cap extensions")
    return paths


def _segments(poly):
    """Iterate over consecutive point pairs (edges) of a polyline."""
    for i in range(1, len(poly)):
        yield poly[i - 1], poly[i]


def _seg_intersect(a1, a2, b1, b2):
    """Intersection point of segments a1→a2 and b1→b2, or None if parallel/non-intersecting."""
    x1, y1 = a1
    x2, y2 = a2
    x3, y3 = b1
    x4, y4 = b2
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(den) < 1e-9:
        return None
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / den
    if -0.02 <= t <= 1.02 and -0.02 <= u <= 1.02:
        return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))
    return None


def _merge_points(points, radius=_MERGE_POINT_RADIUS):
    """Cluster nearby points; return centroids."""
    if not points:
        return []
    used = [False] * len(points)
    hubs = []
    for i, (x, y) in enumerate(points):
        if used[i]:
            continue
        group = [(x, y)]
        used[i] = True
        for j in range(i + 1, len(points)):
            if used[j]:
                continue
            if math.hypot(points[j][0] - x, points[j][1] - y) <= radius:
                group.append(points[j])
                used[j] = True
        hubs.append((sum(p[0] for p in group) / len(group),
                     sum(p[1] for p in group) / len(group)))
    return hubs


def _path_crossing_hubs(paths, merge_radius=_PATH_CROSSING_RADIUS):
    """Junction points from segment crossings and endpoint clusters."""
    hits = []
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            for a1, a2 in _segments(paths[i]):
                for b1, b2 in _segments(paths[j]):
                    pt = _seg_intersect(a1, a2, b1, b2)
                    if pt is not None:
                        hits.append(pt)

    # Endpoints from different paths that meet at a junction often cluster
    # without sharing an exact segment intersection (caps stop short).
    endpoints = []
    for pi, path in enumerate(paths):
        if len(path) < 2:
            continue
        endpoints.append((path[0][0], path[0][1], pi))
        endpoints.append((path[-1][0], path[-1][1], pi))
    for i in range(len(endpoints)):
        xi, yi, pi = endpoints[i]
        group = [(xi, yi)]
        paths_in = {pi}
        for j in range(i + 1, len(endpoints)):
            xj, yj, pj = endpoints[j]
            if pi == pj:
                continue
            if math.hypot(xj - xi, yj - yi) <= merge_radius * 1.35:
                group.append((xj, yj))
                paths_in.add(pj)
        if len(paths_in) >= 2 and len(group) >= 2:
            hits.append((sum(p[0] for p in group) / len(group),
                         sum(p[1] for p in group) / len(group)))

    # Interior T-hits: a path endpoint lands near another path's interior.
    for pi, path in enumerate(paths):
        if len(path) < 2:
            continue
        for end in (path[0], path[-1]):
            ex, ey = end
            for pj, other in enumerate(paths):
                if pi == pj or len(other) < 2:
                    continue
                for k in range(1, len(other) - 1):
                    ox, oy = other[k]
                    if math.hypot(ox - ex, oy - ey) <= merge_radius:
                        hits.append((ox, oy))
                        break

    return _merge_points(hits, merge_radius)


def _nearest_index(path, hub, radius=_PATH_CROSSING_RADIUS):
    """Index of the path point closest to hub within radius, or None."""
    best = None
    for i, (x, y) in enumerate(path):
        d = math.hypot(x - hub[0], y - hub[1])
        if d <= radius and (best is None or d < best[0]):
            best = (d, i)
    return best[1] if best else None


def _split_path_at(path, idx, hub):
    """Split path at idx, projecting the cut to the hub."""
    if idx <= 0 or idx >= len(path) - 1:
        return None
    left = path[:idx + 1] + [hub]
    right = [hub] + path[idx + 1:]
    if len(left) < 2 or len(right) < 2:
        return None
    if plen(left) < 8.0 or plen(right) < 8.0:
        return None
    return left, right


def split_at_crossings(paths, hubs, cut_radius=_PATH_CROSSING_RADIUS, min_len=8.0):
    """Split paths at crossing hubs so each arm stops at the junction."""
    if not hubs:
        return paths
    result = list(paths)
    changed = True
    while changed:
        changed = False
        next_paths = []
        for path in result:
            split_idx = None
            split_hub = None
            for hub in hubs:
                idx = _nearest_index(path, hub, cut_radius)
                if idx is None:
                    continue
                if 0 < idx < len(path) - 1:
                    if split_idx is None or abs(idx - len(path) // 2) < abs(split_idx - len(path) // 2):
                        split_idx = idx
                        split_hub = hub
            if split_idx is not None:
                parts = _split_path_at(path, split_idx, split_hub)
                if parts:
                    next_paths.extend(parts)
                    changed = True
                    continue
            if len(path) >= 2 and plen(path) >= min_len:
                next_paths.append(path)
        result = next_paths
    return result


def add_junction_connectors(paths, hubs, max_gap=_CONNECTOR_MAX_GAP, min_gap=_CONNECTOR_MIN_GAP):
    """Short connector paths from segment ends to shared junction hubs."""
    if not hubs:
        return paths
    connectors = []
    seen = set()
    endpoints = []
    for path in paths:
        if len(path) < 2:
            continue
        endpoints.append(path[0])
        endpoints.append(path[-1])

    for hub in hubs:
        hx, hy = hub
        for (ex, ey) in endpoints:
            gap = math.hypot(ex - hx, ey - hy)
            if gap < min_gap or gap > max_gap:
                continue
            key = (round(ex, 1), round(ey, 1), round(hx, 1), round(hy, 1))
            rev = (round(hx, 1), round(hy, 1), round(ex, 1), round(ey, 1))
            if key in seen or rev in seen:
                continue
            seen.add(key)
            connectors.append([(ex, ey), (hx, hy)])
    return paths + connectors


def bridge_endpoint_gaps(paths, mask, w, h, hubs=(), min_gap=_BRIDGE_MIN_GAP, max_gap=_BRIDGE_MAX_GAP):
    """Connect nearby free endpoints (loop openings, stem gaps)."""
    endpoints = []
    for pi, path in enumerate(paths):
        if len(path) < 2:
            continue
        endpoints.append((pi, 0, path[0]))
        endpoints.append((pi, 1, path[-1]))

    def _near_hub(pt, radius=14.0):
        for hx, hy in hubs:
            if math.hypot(pt[0] - hx, pt[1] - hy) <= radius:
                return True
        return False

    used = set()
    connectors = []
    for i in range(len(endpoints)):
        if i in used:
            continue
        pi, si, (x1, y1) = endpoints[i]
        if _near_hub((x1, y1)):
            continue
        best = None
        for j in range(len(endpoints)):
            if i == j or j in used:
                continue
            pj, sj, (x2, y2) = endpoints[j]
            if pi == pj or _near_hub((x2, y2)):
                continue
            gap = math.hypot(x2 - x1, y2 - y1)
            if gap < min_gap or gap > max_gap:
                continue
            if not _segment_inside_mask(mask, w, h, x1, y1, x2, y2):
                continue
            d1 = _end_dir(paths[pi], si, span=min(12.0, gap))
            d2 = _end_dir(paths[pj], sj, span=min(12.0, gap))
            ux, uy = (x2 - x1) / gap, (y2 - y1) / gap
            score = min(d1[0] * ux + d1[1] * uy, -(d2[0] * ux + d2[1] * uy))
            if score < 0.35:
                continue
            if best is None or gap < best[0]:
                best = (gap, j, x2, y2)
        if best is not None:
            _, j, x2, y2 = best
            connectors.append([(x1, y1), (x2, y2)])
            used.add(i)
            used.add(j)
    return paths + connectors


def finalize_topology(paths, mask, dt, w, h, stroke_width=45):
    """Split arms at crossings and add junction connectors like the reference."""
    hubs = _path_crossing_hubs(paths)
    split = split_at_crossings(paths, hubs) if hubs else paths
    connected = add_junction_connectors(split, hubs, max_gap=stroke_width * 0.85)
    connected = bridge_endpoint_gaps(connected, mask, w, h, hubs,
                                     max_gap=stroke_width * 0.38)
    print(f"  Topology: {len(hubs)} hubs, {len(connected)} paths "
          f"({len(connected) - len(split)} connectors)")
    return connected

def merge_collinear_paths(paths, cos_threshold=0.90, max_gap=_MERGE_MAX_GAP):
    """Merge collinear pieces split by junction processing.

    After split_at_crossings cuts stems at junction hubs, the two halves
    of each stem are collinear and their endpoints meet near (or at) the
    hub. add_junction_connectors may add short stubs. This phase merges
    collinear pieces whose endpoint tangent directions are nearly parallel
    (|dot| > cos_threshold), skipping perpendicular arms (e.g. vertical
    stem into horizontal crossbar at 90 deg) and short connector stubs
    that are not collinear with the main stem direction.
    """
    if len(paths) < 2:
        return paths

    result_paths = [list(p) for p in paths]
    active = [True] * len(paths)

    changed = True
    while changed:
        changed = False
        best_val = -1.0
        best = None  # (i, j, reverse_j, i_first)

        for i in range(len(paths)):
            if not active[i] or len(result_paths[i]) < 2:
                continue
            pi = result_paths[i]
            for j in range(len(paths)):
                if i == j or not active[j] or len(result_paths[j]) < 2:
                    continue
                pj = result_paths[j]

                for si in (0, 1):
                    if si == 0:
                        xi, yi = pi[0]
                    else:
                        xi, yi = pi[-1]
                    di = _end_dir(pi, si)

                    for sj in (0, 1):
                        if sj == 0:
                            xj, yj = pj[0]
                        else:
                            xj, yj = pj[-1]
                        gap = math.hypot(xj - xi, yj - yi)
                        if gap > max_gap:
                            continue

                        dj = _end_dir(pj, sj)
                        dot = di[0] * dj[0] + di[1] * dj[1]
                        abs_dot = abs(dot)
                        if abs_dot <= cos_threshold:
                            continue

                        score = abs_dot - gap / max_gap
                        if score > best_val:
                            best_val = score
                            dot_pos = dot > 0

                            if not dot_pos:
                                # anti-parallel outward: natural connection
                                if si == 1 and sj == 0:
                                    best = (i, j, False, True)
                                elif si == 0 and sj == 1:
                                    best = (i, j, False, False)
                            else:
                                # parallel outward: one path needs reversing
                                if si == 1 and sj == 1:
                                    best = (i, j, True, True)
                                elif si == 0 and sj == 0:
                                    best = (i, j, True, False)

        if best is None:
            break

        i, j, reverse_j, i_first = best
        pj = result_paths[j]
        if reverse_j:
            pj = list(reversed(pj))

        if i_first:
            dup = math.hypot(result_paths[i][-1][0] - pj[0][0],
                             result_paths[i][-1][1] - pj[0][1]) < 1.0
            merged = result_paths[i] + (pj[1:] if dup else pj)
        else:
            dup = math.hypot(pj[-1][0] - result_paths[i][0][0],
                             pj[-1][1] - result_paths[i][0][1]) < 1.0
            merged = pj + (result_paths[i][1:] if dup else result_paths[i])

        result_paths[i] = merged
        active[j] = False
        changed = True

    merged = [result_paths[i] for i in range(len(paths))
              if active[i] and len(result_paths[i]) >= 2]
    removed = len(paths) - len(merged)
    if removed:
        print(f"  Merge collinear: {removed} -> {len(merged)} paths")
    return merged
