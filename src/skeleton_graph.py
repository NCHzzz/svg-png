"""
Post-processing of raw midpoint tracks into a connected skeleton graph.

Pipeline (5 steps):
1. split_on_clearance_jumps — near a junction the perpendicular chord escapes
   through the opening, causing a clearance spike; split chains at these
   outliers so junction resolution can rebuild topology cleanly.
2. dedupe_paths            — every stroke is sampled from both of its sides,
   giving two nearly identical midpoint tracks; keep one, drop fold-backs.
3. prune_covered_fragments — drop short leftover fragments (cap-wrap remnants,
   junction-blob stubs) already covered by longer paths.
4. resolve_junctions       — reconnect the topology that the centeredness
   filter intentionally broke:
   a) endpoint clusters: 2 endpoints -> corner join at the intersection of
      their tangents (elbows, sharp tips); 3+ endpoints -> shared junction
      point (Y/T hubs), with through-joins for collinear pairs;
   b) T-connections: walk a free endpoint forward along the clearance ridge
      until it lands on a path passing by (e.g. H crossbar -> stem);
   c) cap extension: push remaining free ends into the stroke tips until
      clearance drops to the drawn stroke's cap radius, so round caps fill
      the tip.
5. finalize_topology        — split arms at crossing hubs, add short junction
   connectors, bridge nearby endpoint gaps, then merge collinear pieces that
   were cut apart by the split.

Hậu xử lý các đường trung tuyến thô thành đồ thị khung xương liên thông.

Quy trình (5 bước):
1. split_on_clearance_jumps — gần điểm nối, dây cung vuông góc thoát ra
   qua khe hở, gây đột biến clearance; tách chuỗi tại các điểm ngoại lai
   này để bước nối điểm tái tạo cấu trúc liên kết sạch sẽ.
2. dedupe_paths — mỗi nét vẽ được lấy mẫu từ cả hai cạnh, tạo ra hai
   đường trung tuyến gần như giống hệt; giữ một, bỏ đường gập ngược.
3. prune_covered_fragments — bỏ đoạn thừa ngắn (tàn dư bọc đầu, đoạn cụt
   vùng nối) đã bị đường dài hơn che phủ.
4. resolve_junctions — tái kết nối cấu trúc liên kết mà bộ lọc định tâm
   đã cố ý phá vỡ:
   a) cụm điểm đầu: 2 điểm đầu -> nối góc tại giao điểm tiếp tuyến
      (khuỷu, mũi nhọn); 3+ điểm đầu -> điểm nối chung (ngã rẽ Y/T),
      với nối xuyên cho cặp thẳng hàng;
   b) kết nối chữ T: đi dọc đường đỉnh clearance từ điểm đầu tự do
      đến khi chạm đường đi ngang (vd: thanh ngang chữ H -> thân);
   c) mở rộng đầu: đẩy điểm đầu tự do còn lại vào mũi nét vẽ đến khi
      clearance giảm xuống bán kính đầu nét vẽ, để đầu tròn lấp đầy mũi.
5. finalize_topology — tách nhánh tại điểm giao, thêm đoạn nối ngắn,
   bắc cầu khoảng trống gần điểm đầu, rồi hợp nhất đoạn thẳng hàng
   đã bị cắt rời bởi bước tách.
"""



import math
from collections import defaultdict
from spatial_hash import PointGrid

# ── Pixel-distance thresholds (tuned for ~256–1024 px icons) ────────────
# ── Ngưỡng khoảng cách pixel (tinh chỉnh cho icon ~256–1024 px) ──────────

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
    """tính tổng chiều dài của tất cả đoạn nhỏ nối các điểm của một đường polyline ( Polyline là một đường gồm nhiều điểm nối lại với nhau )."""
    return sum(math.hypot(path[i][0] - path[i - 1][0], path[i][1] - path[i - 1][1])
               for i in range(1, len(path)))


def dt_at(dt, x, y, w, h):
    """Distance-transform value at (x,y) in pixel units (chamfer ÷ 3)."""
    """Giá trị biến đổi khoảng cách tại (x,y) theo đơn vị pixel (chamfer ÷ 3)."""
    xi, yi = int(x), int(y)
    if 0 <= xi < w and 0 <= yi < h:
        return dt[yi][xi] / 3.0
    return 0.0


def _inside(mask, x, y, w, h):
    """True if (x,y) is inside the foreground region."""
    """True nếu (x,y) nằm trong mask ( vùng đen )."""
    xi, yi = int(x), int(y)
    return 0 <= xi < w and 0 <= yi < h and mask[yi][xi] == 1


def _segment_inside_mask(mask, w, h, ax, ay, bx, by):
    """True if every sampled point on segment a→b is inside foreground."""
    """True nếu mọi điểm lấy mẫu trên đoạn a→b nằm trong mask."""
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
    """Khi đường trung tuyến đi gần junction, ví dụ chỗ giao của chữ H, K, T, chord có thể bị lệch vào vùng giao nhau lớn. 
    Khi đó điểm midpoint không còn nằm trên stroke bình thường nữa, và giá trị distance transform/clearance bị nhảy vọt. 
    Hàm này phát hiện các điểm bất thường đó và cắt path tại đó.  
    """
    result = []
    for path in paths:
        clear = [dt_at(dt, x, y, w, h) for (x, y) in path]
        cur = []
        window = []
        for i, pt in enumerate(path): # i là index, pt là (x,y) của điểm trong path
            if window:
                med = sorted(window)[len(window) // 2]
                if clear[i] > med * up or clear[i] < med * down: # Nếu giá trị clearance hiện tại vượt quá ngưỡng up hoặc nhỏ hơn ngưỡng down so với median của window, thì coi như điểm này là bất thường
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
    """Mỗi nét vẽ tạo ra hai đường trung tuyến (mỗi cạnh contour một đường),
    và một đường có thể gập ngược lại chính nó tại nơi contour bọc đầu tròn.
    Giữ đường dài nhất trước; bỏ điểm trùng với hình học đã giữ, hoặc trùng
    với chính đường đó từ vị trí xa dọc theo cùng đường (gập ngược)."""
    order = sorted(range(len(paths)), key=lambda i: -plen(paths[i]))
    kept = PointGrid()
    result = []

    # Duyệt từng path theo thứ tự từ dài nhất đến ngắn nhất
    for pi in order:
        """
        path = đường hiện tại đang xét
        own = PointGrid lưu các điểm của chính path này
        arc = chiều dài đã đi được dọc theo path
        dup = danh sách đánh dấu điểm nào bị trùng
        """ 
        path = paths[pi]
        own = PointGrid()
        arc = 0.0
        dup = []
        for k, (x, y) in enumerate(path):
            if k > 0:
                arc += math.hypot(x - path[k - 1][0], y - path[k - 1][1]) # Tính chiều dài đoạn thẳng từ điểm trước đến điểm hiện tại dọc theo path
            clearance = dt_at(dt, x, y, w, h)
            r = min(_DEDUPE_MAX_RADIUS, max(_DEDUPE_MIN_RADIUS, 0.4 * clearance))
            is_dup = kept.near(x, y, r) is not None
            if not is_dup:
                # Self fold-back: same location reached much earlier in arc
                # (track wrapped a cap and doubled back along the far side).
                # Điểm hiện tại có gần điểm cũ của chính path này không?
                # Và điểm cũ đó có cách khá xa theo chiều dài path không?
                guard = 2.0 * r + 6.0
                cur_arc = arc
                if own.any_match(x, y, r, lambda a: cur_arc - a > guard):
                    is_dup = True
            dup.append(is_dup)
            own.add(x, y, arc)

        # Split into unique_ranges of unique points.
        # Tách thành các đoạn chạy gồm điểm duy nhất.
        unique_ranges = []
        start = None
        for i, is_duplicate in enumerate(dup): # duyệt qua danh sách đánh dấu điểm trùng
            if not is_duplicate and start is None:  # Nếu điểm hiện tại không trùng và chưa bắt đầu một đoạn duy nhất nào
                start = i
            elif is_duplicate and start is not None: # Nếu điểm hiện tại trùng và đã bắt đầu một đoạn duy nhất, kết thúc đoạn đó
                unique_ranges.append((start, i))
                start = None
        if start is not None:
            unique_ranges.append((start, len(path)))

        for (a, b) in unique_ranges:
            seg = path[a:b]
            if len(seg) >= 3 and plen(seg) >= min_len:
                for (x, y) in seg:
                    kept.add(x, y)
                result.append(seg) # Thêm đoạn duy nhất vào kết quả

    print(f"  After dedupe: {len(result)} paths")
    return result # trả về danh sách các đoạn path đã được loại bỏ các điểm trùng lặp.


def prune_covered_fragments(paths, dt, w, h, max_len=_PRUNE_MAX_LEN):
    """Drop short leftover tracks (cap-wrap remnants, junction-blob branch
    stubs) that lie entirely within the stroke area already covered by
    longer paths."""
    """Bỏ các đường thừa ngắn (tàn dư bọc đầu, đoạn cụt nhánh vùng nối)
    nằm hoàn toàn trong vùng nét vẽ đã bị đường dài hơn che phủ."""
    if len(paths) < 2:
        return paths
    lengths = [plen(p) for p in paths]
    grids = []
    for p in paths:
        g = PointGrid()
        for (x, y) in p:
            g.add(x, y)
        grids.append(g) # Tạo PointGrid cho từng path

    keep = [True] * len(paths) # mặc định ban đầu giữ tất cả các path
    for pi, path in enumerate(paths): # hàm duyệt từng path để kiểm tra xem có bị ngắn và bị che phủ bởi các path khác không
        if lengths[pi] >= max_len: # nếu path hiện tại dài hơn max_len thì bỏ qua, không cần kiểm tra
            continue
        covered = True
        for (x, y) in path:
            r = max(1.1 * dt_at(dt, x, y, w, h), _PRUNE_COVERAGE_MIN_R)
            hit = False
            for pj in range(len(paths)):
                if pj == pi or not keep[pj] or lengths[pj] <= lengths[pi]: # Không so với chính nó, path đã bị loại bỏ, hoặc path ngắn hơn hoặc bằng path hiện tại
                    continue
                if grids[pj].near(x, y, r) is not None: # Trong những path dài hơn, có path nào có điểm gần (x,y) trong bán kính r không?
                    hit = True
                    break
            if not hit: # chỉ cần 1 điểm không bị che phủ bởi các path dài hơn là đủ để giữ path hiện tại
                covered = False
                break
        if covered:
            keep[pi] = False

    pruned = sum(1 for k in keep if not k) # đếm số lượng path bị loại bỏ
    if pruned:
        print(f"  Pruned {pruned} covered fragments")
    return [p for i, p in enumerate(paths) if keep[i]] # trả về danh sách các path còn lại sau khi loại bỏ các path bị che phủ.


def _is_closed(path):
    """True if path's first and last points are nearly coincident (≈loop)."""
    """True nếu điểm đầu và điểm cuối của đường gần trùng nhau (≈vòng kín)."""
    return len(path) >= 4 and math.hypot(path[0][0] - path[-1][0],
                                         path[0][1] - path[-1][1]) < _CLOSED_COINCIDE # nếu đường có ít nhất 4 điểm và khoảng cách giữa điểm đầu và điểm cuối nhỏ hơn _CLOSED_COINCIDE thì coi như đường là vòng kín.


def _end_dir(path, side, span=_END_DIR_SPAN):
    """Unit direction pointing OUT of the path at the given end, measured
    over `span` pixels of arc length (robust to sub-pixel point noise)."""
    """Hướng đơn vị chỉ RA NGOÀI đường tại đầu đã cho, đo trên `span` pixel
    chiều dài cung (chống nhiễu điểm dưới pixel).
    path = danh sách điểm centerline
    side = 0 nghĩa là lấy hướng ở đầu path
    side = 1 nghĩa là lấy hướng ở cuối path
    span = nhìn sâu vào path khoảng bao nhiêu pixel để tính hướng
    """
    if side == 0:
        bx, by = path[0]
        idx = range(1, len(path)) # đi từ đầu đến cuối
    else:
        bx, by = path[-1]
        idx = range(len(path) - 2, -1, -1) # đi từ cuối đến đầu
    inside_x, inside_y = bx, by # điểm mới nhất đã đi tới trong path - dùng cho việc tính direction sau vòng lặp
    arc = 0.0
    prev_x, prev_y = bx, by  # điểm mới để vòng sau tính khoảng cách - dùng cho việc cộng arc trong vòng lặp
    for k in idx:
        x, y = path[k]
        arc += math.hypot(x - prev_x, y - prev_y)
        prev_x, prev_y = x, y
        inside_x, inside_y = inside_x, inside_y
        if arc >= span:
            break
    dx, dy = bx - inside_x, by - inside_y
    m = math.hypot(dx, dy)
    return (dx / m, dy / m) if m > 1e-9 else (0.0, 0.0) # trả về vector direction đơn vị từ inside point ra ngoài điểm đầu hoặc cuối của path chuẩn hóa (0,1)


def _append_end(path, side, pt):
    """Append pt to start (side=0) or end (side=1) of path in-place."""
    """Thêm pt vào đầu (side=0) hoặc cuối (side=1) của đường tại chỗ."""
    if side == 0:
        path.insert(0, pt)
    else:
        path.append(pt)


def _ray_intersection(p1, d1, p2, d2):
    """Intersection of ray p1+t*d1 and ray p2+t*d2. Returns (t1, t2, x, y) or None."""
    """Giao điểm của tia p1+t*d1 và tia p2+t*d2. Trả về (t1, t2, x, y) hoặc None."""
    """
    p1 = điểm bắt đầu của tia 1
    d1 = hướng của tia 1
    p2 = điểm bắt đầu của tia 2
    d2 = hướng của tia 2
    tìm t1, t2 sao cho p1 + t1*d1 = p2 + t2*d2 ( công thức toán học giao điểm của hai tia )
    """
    
    det = d1[0] * (-d2[1]) - (-d2[0]) * d1[1] # 3 dòng này kiểm tra xem có song song không, nếu song song thì không có giao điểm
    if abs(det) < 1e-6: 
        return None
    rx, ry = p2[0] - p1[0], p2[1] - p1[1]
    t1 = (rx * (-d2[1]) - (-d2[0]) * ry) / det
    t2 = (d1[0] * ry - rx * d1[1]) / det
    return (t1, t2, p1[0] + t1 * d1[0], p1[1] + t1 * d1[1]) # trả về t1, t2 và tọa độ giao điểm (x,y) nếu có giao điểm, ngược lại trả về None nếu hai tia song song.


def resolve_junctions(paths, mask, dt, w, h, cap_clearance=_CAP_CLEARANCE_DEFAULT):
    """Reconnect topology at junctions, corners and caps. Mutates copies."""
    """Tái kết nối cấu trúc liên kết tại điểm nối, góc và đầu nét. Biến đổi bản sao."""
    """
    Input:
    paths = nhiều đường centerline đã được lọc bớt trùng
    mask  = ảnh nhị phân, pixel đen là 1
    dt    = distance transform / clearance map
    w, h  = kích thước ảnh
    cap_clearance = clearance mặc định của nửa stroke
    Output:
    paths = các đường centerline đã được nối lại ở góc, junction, T-shape, cap
    """
    """
    resolve_junctions() làm 2 lượt:
    Pass A:
    - Gom endpoint gần nhau thành cluster.
    - 2 endpoint → nối góc.
    - 3+ endpoint → nối cặp thẳng hàng, hoặc nối vào hub.

    Pass B:
    - Với endpoint còn dư:
    - đi tiếp theo hướng centerline
    - nếu gặp path khác → T-connect
    - nếu không gặp nhưng còn trong cap → kéo dài đầu nét
    """
    paths = [list(p) for p in paths]

    # ---- collect free endpoints -------------------------------------------
    # ---- thu thập điểm đầu tự do ------------------------------------------
    # Mỗi path có 2 đầu thu thập chúng. Ví dụ: path = [(10, 0), (20, 0), (30, 0)]
    # [pi, 0, 10, 0]  # đầu path
    # [pi, 1, 30, 0]  # cuối path
    endpoints = []  # (path_idx, side, x, y)
    for pi, p in enumerate(paths):
        if len(p) < 2 or _is_closed(p):
            continue
        endpoints.append([pi, 0, p[0][0], p[0][1]])
        endpoints.append([pi, 1, p[-1][0], p[-1][1]])

    resolved = set()  # endpoint list indices already handled

    # ---- pass A: endpoint clusters (corners and junction hubs) ------------
    # ---- lượt A:  nối endpoint gần nhau thành góc / junction / hub -----------------
    m = len(endpoints)
    parent = list(range(m))

    # dùng để tìm gốc cụm của endpoint a
    # đi ngược theo parent để xem endpoint a thuộc nhóm nào.
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
            if gap < max(_CLUSTER_GAP_BASE, 1.1 * (di + dj)): # Nếu Nếu hai endpoint đủ gần nhau→ có thể chúng thuộc cùng một junction
                if _segment_inside_mask(mask, w, h, xi, yi, xj, yj):
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        parent[ra] = rb

    # Tạo các cụm endpoint từ kết quả union-find
    clusters = defaultdict(list)
    for i in range(m):
        clusters[find(i)].append(i)

    def _facing_score(i, j):
        """How well two endpoints continue into each other: both must point
        toward the other, and their directions must be anti-parallel."""
        """Mức độ hai điểm đầu tiếp nối vào nhau: cả hai phải hướng về nhau,
        và hướng của chúng phải ngược chiều (anti-parallel)."""
        pi, si, xi, yi = endpoints[i]
        pj, sj, xj, yj = endpoints[j]
        d1 = _end_dir(paths[pi], si)
        d2 = _end_dir(paths[pj], sj)
        gap = math.hypot(xj - xi, yj - yi)
        if gap < 1e-9:
            return 1.0
        ux, uy = (xj - xi) / gap, (yj - yi) / gap
        toward1 = d1[0] * ux + d1[1] * uy # endpoint i có hướng về endpoint j không
        toward2 = -(d2[0] * ux + d2[1] * uy) # endpoint j có hướng về endpoint i không
        anti = -(d1[0] * d2[0] + d1[1] * d2[1]) # hai hướng có ngược chiều nhau không
        if toward1 < 0.4 or toward2 < 0.4:
            return -1.0
        return min(anti, toward1, toward2) # lấy điểm yếu nhất làm điểm tổng. => mặc định cả 3 đều thỏa

    def _corner_join(i, j):
        """Join two endpoints at the intersection of their tangents (falls
        back to their midpoint), e.g. elbows, sharp tips, loop closures."""
        """Nối hai điểm đầu tại giao điểm tiếp tuyến của chúng (dùng điểm
        giữa nếu không tìm được), vd: khuỷu, mũi nhọn, đóng vòng."""
        pi, si, xi, yi = endpoints[i]
        pj, sj, xj, yj = endpoints[j]
        # Lấy hướng của 2 endpoint
        d1 = _end_dir(paths[pi], si)
        d2 = _end_dir(paths[pj], sj)
        gap = math.hypot(xi - xj, yi - yj)
        corner = None
        # Tìm giao điểm của 2 tia đi ra từ 2 endpoint
        hit = _ray_intersection((xi, yi), d1, (xj, yj), d2)
        if hit is not None:
            t1, t2, cx, cy = hit
            tmax = 2.5 * gap + _CORNER_TMAX_BASE
            """
            1. Giao điểm không nằm quá xa.
            2. Giao điểm nằm trong vùng đen của icon.
            3. Đoạn từ endpoint i đến corner nằm trong vùng đen.
            4. Đoạn từ endpoint j đến corner nằm trong vùng đen.
            """
            if -4.0 <= t1 <= tmax and -4.0 <= t2 <= tmax and \
                    _inside(mask, cx, cy, w, h) and \
                    _segment_inside_mask(mask, w, h, xi, yi, cx, cy) and \
                    _segment_inside_mask(mask, w, h, xj, yj, cx, cy):
                corner = (cx, cy)
        if corner is None:
            corner = ((xi + xj) / 2.0, (yi + yj) / 2.0) # Nếu không tìm được giao điểm, lấy điểm giữa 2 endpoint làm corner
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
        """Nối các điểm đầu vào một điểm nối chung. Trung tâm được đặt tại
        nơi các tia tiếp tuyến giao nhau (đó là nơi các đường tâm nhánh
        thực sự gặp nhau, vd: sâu trong đỉnh mũi tên); dùng trọng tâm điểm
        đầu nếu không tìm được. Điểm đầu không hướng về trung tâm được để
        lại cho lượt B."""
        
        """
        Hàm resolve_junctions() là hàm nối lại các centerline path đang bị đứt ở junction, góc, chữ T, đầu nét/cap. 
        Đây là bước sau khi đã có nhiều path thô nhưng chúng chưa nối topology đúng.
        Input:
        paths = nhiều đường centerline đã được lọc bớt trùng
        Output:
        paths = các đường centerline đã được nối lại ở góc, junction, T-shape, cap
        """
        
        # Tìm giao điểm của tất cả các endpoint trong cụm members. 
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
        if pts: # nếu nhiều giao điểm, lay trung bình của chúng làm hub point
            cx = sum(p[0] for p in pts) / len(pts)
            cy = sum(p[1] for p in pts) / len(pts)
        else: # nếu không tìm được giao điểm, dùng trung điểm ( centeroid endpoints ) làm hub point
            cx = sum(endpoints[i][2] for i in members) / len(members)
            cy = sum(endpoints[i][3] for i in members) / len(members)
        joined = 0
        for i in members:
            pi, si, xi, yi = endpoints[i]
            gx, gy = cx - xi, cy - yi # vector từ endpoint đến hub point
            gm = math.hypot(gx, gy)
            if gm > 1e-9:
                d = _end_dir(paths[pi], si) # hướng của endpoint
                if (d[0] * gx + d[1] * gy) / gm < 0.0: # Nếu hướng endpoint không hướng về hub point, bỏ qua endpoint này.
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
        if len(members) < 2: # nếu cụm chỉ có 1 endpoint, không cần nối gì cả
            continue
        if len(members) == 2: # nếu cụm có 2 endpoint, nối chúng lại bằng góc
            _corner_join(members[0], members[1])
            corner_joins += 1
            continue

        # 3+ endpoints: a stroke passing straight through the junction was
        # split into two collinear pieces. Reconnect the best facing pairs
        # through the junction, then join remaining branch ends to a hub.
        pending = list(members)
        while len(pending) >= 2:
            best = None
            for a in range(len(pending)): # duyệt từng endpoint trong cụm
                for b in range(a + 1, len(pending)):
                    s = _facing_score(pending[a], pending[b]) # tính điểm số mức độ 2 endpoint có hướng về nhau và ngược chiều nhau
                    if s > 0.75 and (best is None or s > best[0]): # nếu điểm số > 0.75 và tốt hơn điểm số tốt nhất trước đó, lưu lại cặp endpoint này
                        best = (s, pending[a], pending[b]) # lưu lại điểm số và index của 2 endpoint
            if best is None:
                break
            _, i, j = best
            pi, si, xi, yi = endpoints[i]
            pj, sj, xj, yj = endpoints[j]
            mid = ((xi + xj) / 2.0, (yi + yj) / 2.0) # lấy trung điểm của 2 endpoint làm điểm nối
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
        # Nếu còn nhiều endpoint dư thì nối hub
        if len(pending) >= 3:
            if _hub_join(pending):
                hub_joins += 1

    # ---- pass B: T-connections and cap extension --------------------------
    # ---- lượt B: xử lý endpoint còn dư, nối T-shape hoặc kéo dài cap -------------------------
    # Registry of all current path points for T-hit tests.
    # Sổ đăng ký tất cả điểm đường hiện tại cho kiểm tra va chạm chữ T.
    registry = PointGrid() # tất cả path hiện tại để kiểm tra va chạm T-connect
    for pi, p in enumerate(paths):
        for (x, y) in p:
            registry.add(x, y, pi)

    t_connects = 0
    cap_extends = 0
    for ei, (pi, si, x, y) in enumerate(endpoints):
        if ei in resolved:
            continue # nếu endpoint đã được xử lý ở lượt A thì bỏ qua, B xử lý endpoint còn tự do / còn dư
        dx, dy = _end_dir(paths[pi], si)
        if dx == 0 and dy == 0: # nếu hướng endpoint không xác định được (ví dụ path chỉ có 1 điểm) thì bỏ qua endpoint này
            continue
        d_end = dt_at(dt, x, y, w, h) # lấy giá trị distance transform tại endpoint
        if d_end < 0.85 * cap_clearance: 
            # nằm ở vùng hẹp hơn stroke bình thường. 
            # ( hiện tại thì witđh tầm 45px 
            # thì clearance tầm 22px, nếu clearance < 0.85 * 22px 
            # thì coi như đầu nét vẽ đang thuôn nhọn )
            # cần đi sâu hơn để đầu nét vẽ chạm tới mũi của hình.
            # Tapering terminal (stroke narrows toward the tip): follow it
            # deeper so the drawn cap reaches the tip of the shape.
            stop_dt = max(0.55 * d_end, 6.0)
        else:
            # Flat/round cap of a full-width stroke: stop half a drawn
            # stroke-width from the boundary so the round cap lands flush.
            # Đầu phẳng/tròn của nét vẽ đầy đủ: dừng cách biên nửa
            # chiều rộng nét vẽ để đầu tròn nằm khít.
            stop_dt = min(cap_clearance, 0.9 * d_end)
        max_walk = min(100.0, 4.0 * max(d_end, 10.0)) # giới hạn chiều dài đi bộ tối đa để tìm T-connect hoặc cap extension, tránh đi quá xa khỏi endpoint

        # Already touching another path?
        # Đã chạm đường khác chưa?
        if registry.near(x, y, _T_TOUCH_RADIUS, exclude_tag=pi) is not None:
            continue
        """
                | Biến gốc    | Tên dễ hiểu hơn            | Ý nghĩa                               |
        | ----------- | -------------------------- | ------------------------------------- |
        | `cx, cy`    | `cur_x, cur_y`             | vị trí hiện tại khi đang walk         |
        | `dx, dy`    | `dir_x, dir_y`             | hướng đang đi                         |
        | `nx, ny`    | `normal_x, normal_y`       | hướng vuông góc để thử lệch trái/phải |
        | `off`       | `side_offset`              | độ lệch trái/phải                     |
        | `px_, py_`  | `candidate_x, candidate_y` | điểm thử tiếp theo                    |
        | `c`         | `candidate_clearance`      | clearance tại điểm thử                |
        | `best`      | `best_candidate`           | điểm thử tốt nhất                     |
        | `step_c`    | `next_clearance`           | clearance của điểm được chọn          |
        | `sx_, sy_`  | `next_x, next_y`           | điểm tiếp theo được chọn              |
        | `m`         | `step_len`                 | độ dài bước vừa đi                    |
        | `walked`    | `walked_dist`              | tổng quãng đường đã đi                |
        | `cap_trail` | `extension_points`         | các điểm hợp lệ để kéo dài cap        |
        | `connected` | `connected_to_other_path`  | đã nối được vào path khác chưa        |
        | `hit`       | `nearby_path`              | có path khác gần điểm hiện tại không  |

        """
        cx, cy = x, y # current position, start at endpoint
        cap_trail = []  # ridge positions (danh sách điểm đi thêm hợp lệ ) with clearance >= stop_dt  (vị trí đỉnh clearance >= stop_dt)
        walked = 0.0
        connected = False
        while walked < max_walk:
            # Steer along the clearance ridge: of straight-ahead and two
            # slightly turned steps, take the one with max clearance. This
            # lets the walk follow curved tapering terminals.
            # Lái dọc theo đường đỉnh clearance: trong số bước thẳng và hai
            # bước hơi rẽ, chọn bước có clearance cao nhất. Việc này cho phép
            # đường đi bám theo đầu thuôn cong.
            nx, ny = -dy, dx
            best = None
            for off in (0.0, 0.7, -0.7): # thử 3 hướng: thẳng, lệch trái, lệch phải
                px_ = cx + dx + nx * off
                py_ = cy + dy + ny * off
                c = dt_at(dt, px_, py_, w, h)
                if best is None or c > best[0]: # chọn điểm thử có clearance cao nhất ( centerline thường nằm ở nơi xa biên nhất, )
                    best = (c, px_, py_)
            step_c, sx_, sy_ = best # step_c = clearance của điểm thử tốt nhất, sx_, sy_ = tọa độ điểm thử tốt nhất
            m = math.hypot(sx_ - cx, sy_ - cy) # quãng đường từ điểm cũ tới điểm mới
            dx, dy = (sx_ - cx) / m, (sy_ - cy) / m # hướng đi từ điểm cũ tới điểm mới
            cx, cy = sx_, sy_ # cập nhật vị trí hiện tại
            walked += m # cộng quãng đường đã đi
            if not _inside(mask, cx, cy, w, h):
                break
            hit = registry.near(cx, cy, _T_HIT_RADIUS, exclude_tag=pi) # kiểm tra xem có path khác gần điểm hiện tại không, nếu có thì nối vào path đó
            if hit is not None:
                _append_end(paths[pi], si, (cx, cy)) # nối endpoint hiện tại vào điểm chạm với path khác
                registry.add(cx, cy, pi) # thêm điểm chạm vào registry để các endpoint khác không nối vào điểm này nữa
                t_connects += 1
                connected = True
                break
            if step_c >= stop_dt: #  Nếu clearance tại điểm thử >= stop_dt thì thêm điểm thử vào cap_trail để kéo dài cap
                cap_trail.append((cx, cy))
            else:
                break
        """
        - Không nối được vào path khác
        - Nhưng có đi thêm được vài điểm hợp lệ
        - Và đoạn đi thêm đủ dài
        """
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
    """
    Duyệt qua các cặp điểm liên tiếp (cạnh) của polyline. 
    và trả về các cặp điểm (x1, y1), (x2, y2), ...
    """
    for i in range(1, len(poly)):
        yield poly[i - 1], poly[i]


def _seg_intersect(a1, a2, b1, b2):
    """Intersection point of segments a1→a2 and b1→b2, or None if parallel/non-intersecting."""
    """Giao điểm của đoạn a1→a2 và b1→b2, hoặc None nếu song song/không giao.
    Hàm này kiểm tra hai đoạn thẳng có cắt nhau không
    """
    x1, y1 = a1
    x2, y2 = a2
    x3, y3 = b1
    x4, y4 = b2
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4) # tính định thức để kiểm tra hai đoạn có song song không
    if abs(den) < 1e-9:
        return None
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den # tính tham số t để tìm giao điểm trên đoạn a1→a2
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / den # tính tham số u để tìm giao điểm trên đoạn b1→b2
    if -0.02 <= t <= 1.02 and -0.02 <= u <= 1.02: # nếu t và u nằm trong khoảng [0, 1] thì hai đoạn cắt nhau, nới rộng tránh sai số float
        return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))
    return None


def _merge_points(points, radius=_MERGE_POINT_RADIUS):
    """Cluster nearby points; return centroids."""
    """Gom cụm các điểm gần nhau; trả về trọng tâm."""
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
    """Điểm nối từ giao điểm đoạn và cụm điểm đầu."""
    """
    Input: nhiều centerline path
    Output: danh sách các điểm hub, 
    tức là nơi các path giao nhau / gần chạm nhau / tạo chữ T
    Hàm này tìm hub bằng 3 cách
    1. Tìm giao điểm thật giữa các đoạn path.
    2. Tìm các endpoint của nhiều path nằm gần nhau.
    3. Tìm endpoint của path này nằm gần thân path khác, tức T-junction.
    """
    
    """
    Duyệt từng cặp path khác nhau.
    Tách mỗi path thành các đoạn nhỏ.
    Kiểm tra từng đoạn của path A có cắt từng đoạn của path B không.
    Nếu cắt nhau → lưu giao điểm vào hits.
    """
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
    # Điểm đầu từ các đường khác nhau gặp tại điểm nối thường tụ lại
    # mà không chia sẻ giao điểm đoạn chính xác (đầu nét dừng sớm).
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
                paths_in.add(pj) # Lấy một endpoint->Tìm endpoint khác từ path khác nằm gần nó. -> Nếu gần → gom vào group."""
        if len(paths_in) >= 2 and len(group) >= 2:
            hits.append((sum(p[0] for p in group) / len(group),
                         sum(p[1] for p in group) / len(group)))  #Nếu group có 2 path khác nhau và có ít nhất 2 endpoint → tính trọng tâm của group làm hub point.


    # Interior T-hits: a path endpoint lands near another path's interior.
    # Va chạm chữ T nội bộ: điểm đầu đường nằm gần phần thân của đường khác.
    for pi, path in enumerate(paths):
        if len(path) < 2:
            continue
        for end in (path[0], path[-1]): # chỉ xét endpoint đầu và cuối của path
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
    """tìm điểm trên path gần hub nhất, nhưng chỉ chấp nhận nếu điểm đó nằm trong một bán kính cho phép."""
    best = None
    for i, (x, y) in enumerate(path):
        d = math.hypot(x - hub[0], y - hub[1])
        if d <= radius and (best is None or d < best[0]):
            best = (d, i)
    return best[1] if best else None


def _split_path_at(path, idx, hub):
    """Split path at idx, projecting the cut to the hub."""
    """Tách đường tại idx, chiếu vết cắt đến trung tâm."""
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
    """Tách đường tại trung tâm giao cắt để mỗi nhánh dừng đúng tại điểm nối."""
    """
    (0,0) ---- (10,0) ---- (20,0) ---- (30,0) ---- (40,0)
                         *
                        hub
    Sau hàm này, nó sẽ thành:
    [
    [(0, 0), (10, 0), (20, 0), hub],
    [hub, (30, 0), (40, 0)]
]                   
    """
    
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
                idx = _nearest_index(path, hub, cut_radius) # tìm điểm trên path gần hub nhất, nhưng chỉ chấp nhận nếu điểm đó nằm trong một bán kính cho phép.
                if idx is None:
                    continue
                if 0 < idx < len(path) - 1: # Chỉ cắt nếu hub nằm ở giữa path
                    if split_idx is None or abs(idx - len(path) // 2) < abs(split_idx - len(path) // 2): # chọn hub gần trung tâm path nhất để cắt
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
    """Nếu đầu path đang gần hub nhưng chưa chạm hub,
        thì thêm một path ngắn nối từ endpoint đến hub."""
    """
    paths   = danh sách centerline path hiện tại
    hubs    = các điểm junction đã tìm được
    max_gap = khoảng cách tối đa được phép nối
    min_gap = khoảng cách tối thiểu để nối
    Trả về
    danh sách path cũ cộng thêm các đoạn nối mới.
    """
    if not hubs:
        return paths
    connectors = [] # các đoạn mới nối từ endpoint đến hub
    seen = set() # để tránh nối trùng, lưu lại các cặp endpoint-hub đã nối
    endpoints = [] # danh sách đầu/cuối của tất cả path
    for path in paths:
        if len(path) < 2:
            continue
        endpoints.append(path[0])
        endpoints.append(path[-1])

    for hub in hubs:
        hx, hy = hub
        for (ex, ey) in endpoints:
            gap = math.hypot(ex - hx, ey - hy) # tính khoảng cách từ endpoint đến hub
            if gap < min_gap or gap > max_gap:
                continue
            key = (round(ex, 1), round(ey, 1), round(hx, 1), round(hy, 1)) # tạo key ( endpoint-hub) để lưu vào seen, tránh nối trùng
            rev = (round(hx, 1), round(hy, 1), round(ex, 1), round(ey, 1)) # tạo key ngược ( hub-endpoint) để kiểm tra xem đã nối từ hub đến endpoint chưa
            if key in seen or rev in seen:
                continue
            seen.add(key)
            connectors.append([(ex, ey), (hx, hy)])
    return paths + connectors


def bridge_endpoint_gaps(paths, mask, w, h, hubs=(), min_gap=_BRIDGE_MIN_GAP, max_gap=_BRIDGE_MAX_GAP):
    """Connect nearby free endpoints (loop openings, stem gaps)."""
    """Nối điểm đầu tự do gần nhau (khe hở vòng, khoảng trống thân)."""
    """
    paths = danh sách centerline path
    mask  = ảnh nhị phân, vùng đen là 1
    w, h  = kích thước ảnh
    hubs  = danh sách junction/hub đã biết
    min_gap, max_gap = khoảng cách nhỏ/lớn nhất được phép nối
    """
    endpoints = []
    for pi, path in enumerate(paths):
        if len(path) < 2:
            continue
        endpoints.append((pi, 0, path[0])) # đầu path
        endpoints.append((pi, 1, path[-1])) # cuối path

    def _near_hub(pt, radius=14.0): # kiểm tra xem điểm pt có nằm gần hub không, nếu gần thì bỏ qua không nối endpoint này
        for hx, hy in hubs:
            if math.hypot(pt[0] - hx, pt[1] - hy) <= radius:
                return True
        return False

    used = set() # lưu lại các endpoint đã được nối, tránh nối trùng
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
            pj, sj, (x2, y2) = endpoints[j] # lấy endpoint thứ j
            if pi == pj or _near_hub((x2, y2)): # nếu endpoint j cùng path với endpoint i hoặc nằm gần hub thì bỏ qua
                continue
            gap = math.hypot(x2 - x1, y2 - y1)
            if gap < min_gap or gap > max_gap: # nếu khoảng cách từ endpoint i đến endpoint j không nằm trong khoảng [min_gap, max_gap] thì bỏ qua
                continue
            if not _segment_inside_mask(mask, w, h, x1, y1, x2, y2): # nếu đoạn nối từ endpoint i đến endpoint j không nằm hoàn toàn trong vùng đen của mask thì bỏ qua
                continue
            d1 = _end_dir(paths[pi], si, span=min(12.0, gap)) # lấy hướng tiếp tuyến tại endpoint i
            d2 = _end_dir(paths[pj], sj, span=min(12.0, gap)) # lấy hướng tiếp tuyến tại endpoint j
            ux, uy = (x2 - x1) / gap, (y2 - y1) / gap # vector đơn vị từ endpoint i đến endpoint j
            score = min(d1[0] * ux + d1[1] * uy, -(d2[0] * ux + d2[1] * uy)) # tính điểm số mức độ 2 endpoint có hướng về nhau và ngược chiều nhau
            if score < 0.35:
                continue
            if best is None or gap < best[0]: # nếu khoảng cách từ endpoint i đến endpoint j nhỏ hơn khoảng cách tốt nhất trước đó thì lưu lại cặp endpoint này
                best = (gap, j, x2, y2)
        if best is not None: # nếu tìm được endpoint j tốt nhất để nối với endpoint i thì thêm đoạn nối từ endpoint i đến endpoint j vào danh sách connectors và đánh dấu cả 2 endpoint đã được nối
            _, j, x2, y2 = best
            connectors.append([(x1, y1), (x2, y2)])
            used.add(i)
            used.add(j)
    return paths + connectors


def finalize_topology(paths, mask, dt, w, h, stroke_width=45):
    """Split arms at crossings and add junction connectors like the reference."""
    """Tách nhánh tại điểm giao và thêm đoạn nối điểm nối như bản tham chiếu."""
    """
    Input:
    nhiều path centerline đã được nối tương đối
    Output:
    path sạch hơn, có hub rõ ràng, có connector ngắn ở junction/gap
    """
    hubs = _path_crossing_hubs(paths) # tìm các hub từ giao điểm, endpoint gần nhau, endpoint gần path khác (T-junction)
    split = split_at_crossings(paths, hubs) if hubs else paths # nếu có hub thì tách path tại hub, nếu không thì giữ nguyên paths
    connected = add_junction_connectors(split, hubs, max_gap=stroke_width * 0.85) # nếu có hub thì thêm đoạn nối từ endpoint đến hub, nếu không thì giữ nguyên split
    connected = bridge_endpoint_gaps(connected, mask, w, h, hubs,
                                     max_gap=stroke_width * 0.38) # nối các endpoint gần nhau nhưng không nằm gần hub, khoảng cách tối đa được phép nối là 0.38 * stroke_width
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
    """Hợp nhất các đoạn thẳng hàng bị tách bởi xử lý điểm nối.

    Sau khi split_at_crossings cắt thân tại trung tâm điểm nối, hai nửa
    của mỗi thân thẳng hàng và điểm đầu của chúng gặp nhau gần (hoặc tại)
    trung tâm. add_junction_connectors có thể thêm đoạn cụt ngắn. Pha này
    hợp nhất các đoạn thẳng hàng có hướng tiếp tuyến điểm đầu gần song song
    (|dot| > cos_threshold), bỏ qua nhánh vuông góc (vd: thân dọc vào
    thanh ngang ở 90 độ) và đoạn cụt nối ngắn không thẳng hàng với hướng
    thân chính.
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
