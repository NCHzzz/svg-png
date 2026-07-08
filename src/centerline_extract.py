"""
Core centerline extraction via contour-based chord sampling
(not morphological thinning — the name is a misnomer kept for history).

binary mask → contour trace → perpendicular rays across the stroke
→ midpoints → chain midpoints by contour order.

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

# Tính toán biến đổi khoảng cách chamfer 3-4 hai lần ( thẳng cost 3, chéo cost 4). 
# Trả về một mảng 2D dt, trong đó dt[y][x] là khoảng cách chamfer 
# từ pixel (x, y) đến pixel foreground ( nền trắng ) gần nhất.
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
        for x in range(w): # từ trên xuống dưới, từ trái sang phải
            if mask[y][x] == 0:
                continue
            best = dt[y][x]
            if y > 0:
                if x > 0 and dt[y - 1][x - 1] + 4 < best: # trên trái ( đi chéo )
                    best = dt[y - 1][x - 1] + 4
                if dt[y - 1][x] + 3 < best: # trên ( đi thẳng )
                    best = dt[y - 1][x] + 3
                if x + 1 < w and dt[y - 1][x + 1] + 4 < best: # trên phải ( đi chéo )
                    best = dt[y - 1][x + 1] + 4
            if x > 0 and dt[y][x - 1] + 3 < best: # trái ( đi thẳng )
                best = dt[y][x - 1] + 3
            dt[y][x] = best
    for y in range(h - 1, -1, -1): # từ dưới lên trên, từ phải sang trái
        for x in range(w - 1, -1, -1):
            if mask[y][x] == 0:
                continue
            best = dt[y][x]
            if y + 1 < h:
                if dt[y + 1][x] + 3 < best: # dưới ( đi thẳng )
                    best = dt[y + 1][x] + 3
                if x + 1 < w and dt[y + 1][x + 1] + 4 < best: # dưới phải ( đi chéo )
                    best = dt[y + 1][x + 1] + 4
                if x > 0 and dt[y + 1][x - 1] + 4 < best: # dưới trái ( đi chéo )
                    best = dt[y + 1][x - 1] + 4
            if x + 1 < w and dt[y][x + 1] + 3 < best: # phải ( đi thẳng ) 
                best = dt[y][x + 1] + 3
            dt[y][x] = best
    return dt

# Kiểm tra xem pixel (x, y) có phải là biên của mask hay không. 
# Một pixel được coi là biên nếu nó là pixel foreground (mask[y][x] != 0) 
# và có ít nhất một pixel láng giềng là background (mask[ny][nx] == 0) 
# hoặc nằm ngoài biên của hình ảnh.
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

# Moore-style boundary following algorithm to trace all contours in a binary mask.
def trace_all_contours(mask, width, height):
    """Moore-style boundary following. Returns all contours (outer + holes)
    as ordered pixel lists."""
    w, h = width, height
    visited = bytearray(w * h)
    dirs = [
    (1, 0),    # 0: phải
    (1, 1),    # 1: phải-dưới
    (0, 1),    # 2: dưới
    (-1, 1),   # 3: trái-dưới
    (-1, 0),   # 4: trái
    (-1, -1),  # 5: trái-trên
    (0, -1),   # 6: trên
    (1, -1)    # 7: phải-trên
    ]
    contours = []
    
    for y0 in range(h):
        for x0 in range(w): # tìm kiếm điểm biên đầu tiên chưa được thăm
            if not _is_boundary(mask, x0, y0, w, h) or visited[y0 * w + x0]:
                continue
            start = (x0, y0)
            contour = [start]
            visited[start[1] * w + start[0]] = 1
            cx, cy = start # điểm hiện tại
            px, py = -1, -1 # điểm trước đó
            ss = 0 # start search direction ( bắt đầu hướng 3h )
            closed = False
            for _ in range(w * h): # vòng lặp ( số lần nhỏ hơn pixel count, tìm điểm tiếp theo )
                found = False
                for d in range(8):
                    di = (ss + d) % 8
                    dx, dy = dirs[di]
                    nx, ny = cx + dx, cy + dy
                    if 0 <= nx < w and 0 <= ny < h and _is_boundary(mask, nx, ny, w, h):
                        if (nx, ny) == (px, py):
                            continue
                        if (nx, ny) == start and len(contour) > 2:
                            closed = True
                            found = True
                            break
                        contour.append((nx, ny))
                        visited[ny * w + nx] = 1
                        px, py = cx, cy
                        cx, cy = nx, ny
                        ss = (di + 5) % 8 # băt đầu tìm kiếm từ hướng trước đó + 5 (ngược chiều kim đồng hồ)
                        found = True
                        break
                if closed or not found: # Nếu đã đóng vòng hoặc không tìm thấy điểm tiếp theo, thoát khỏi vòng lặp
                    break
            if len(contour) >= 10: # Chỉ giữ những đường biên đủ lớn, bỏ các đường biên quá nhỏ.
                contours.append(contour)
    return contours

# Tính toán vector tiếp tuyến tại điểm contour[i] 
# bằng cách sử dụng một cửa sổ các điểm lân cận.
def _tangent(contour, i, window):
    n = len(contour)
    im = (i - window) % n #  index của điểm phía sau i một đoạn window
    ip = (i + window) % n # index của điểm phía trước i một đoạn window
    dx = contour[ip][0] - contour[im][0]
    dy = contour[ip][1] - contour[im][1]
    mag = math.hypot(dx, dy)
    return (dx / mag, dy / mag) if mag > 1e-6 else (0.0, 0.0) # chuẩn hóa 0-1

# chọn hướng pháp tuyến đi vào trong shape 
def _inward_normal(mask, x, y, tx, ty, w, h):
    """Of the two perpendiculars to the tangent, pick the one pointing
    into the shape (deeper probe stays inside longer)."""
    # Từ tangent (tx, ty), xoay 90° để ra pháp tuyến. Công thức:
    # n1: xoay +90°:   n1x = -ty    n1y = tx
    # n2: xoay -90°:   n2x = ty     n2y = -tx
    n1x, n1y = -ty, tx
    n2x, n2y = ty, -tx
    # đếm độ sâu, trả về sô bước đi được trước khi gặp biên hoặc vượt ra ngoài hình ảnh.
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
    dt = chamfer_dt(mask, w, h) # Distance transform — mỗi pixel đen lưu khoảng cách đến pixel trắng gần nhất
    contours = trace_all_contours(mask, w, h) # Danh sách các đường biên, mỗi contour là 1 vòng pixel khép kín
    print(f"  Contours: {len(contours)} ({sum(len(c) for c in contours)} pts)")

    max_igap = sample_step * 8 # max khoảng cách theo index trên contour:
    max_sgap = 25.0 # ngưỡng để tránh nối sai 2 điểm centerline ở xa nhau.
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

            ox, oy = fx - nx, fy - ny  # last point still inside / điểm mép đối diện
            mx = (x + ox) / 2.0 # trung điểm
            my = (y + oy) / 2.0 # trung điểm
            chord = math.hypot(ox - x, oy - y) # đường xuyên qua stroke
            if chord < 4.0:
                continue
            mxi, myi = int(mx), int(my)
            if not (0 <= mxi < w and 0 <= myi < h and mask[myi][mxi]): # kiểm tra kĩ lại xem midpoint còn nằm trong shape không
                continue

            # Centeredness filter: clearance(midpoint) must match half chord.
            # Strict tolerance rejects diagonal chords at cap corners and
            # chords reaching into junction blobs; each stroke is sampled
            # from both sides, so losing one side's samples is harmless.
            # chứng minh điểm thẳng hàng
            half = chord / 2.0
            clearance = dt[myi][mxi] / 3.0
            if abs(half - clearance) > max(2.5, 0.15 * half): # cho phép sai số
                continue

            # Refine: the chord crosses the medial axis where clearance
            # peaks. The plain midpoint is biased on curved strokes (the
            # chord is perpendicular to one side only); slide to the max.
            # trượt dòng trên chord và chọn điểm có distance transform lớn nhất.
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
            samples.append((idx, bx, by, best_c)) #idx     = index của điểm contour gốc
                                                # bx, by  = điểm centerline tìm được
                                                # best_c  = clearance tốt nhất tại điểm đó
        if len(samples) < 3:
            continue

        # Chain consecutive samples; break on contour-index or spatial gaps.
        chains = [[samples[0]]]
        for k in range(1, len(samples)):
            pidx, px_, py_, _ = chains[-1][-1] # điểm trước node chain hiện tại
            cidx, cx_, cy_, _ = samples[k] # node chain hiện tại
            if (cidx - pidx) > max_igap or math.hypot(cx_ - px_, cy_ - py_) > max_sgap:
                chains.append([samples[k]]) # Tạo chain mới, chứa điểm hiện tại làm phần tử đầu tiên.
            else:
                chains[-1].append(samples[k]) # Không vượt ngưỡng nào → điểm hiện tại thuộc cùng chain với điểm trước → thêm vào chain cuối cùng.

        # The contour is a closed loop: the last chain may continue into the
        # first one across the index wrap-around.
        def _wrap_connects(tail, head):
            igap = head[0][0] + n - tail[-1][0] # Khoảng cách theo index từ điểm cuối của tail đi vòng qua cuối contour rồi về điểm đầu của head.
            sgap = math.hypot(head[0][1] - tail[-1][1], head[0][2] - tail[-1][2]) # sample dạng (idx, mx, my, best_c) => head[0][1] là mx, head[0][2] là my. Khoảng cách Euclidean giữa điểm cuối tail và điểm đầu head.
            return igap <= max_igap and sgap <= max_sgap # Nếu chain cuối và chain đầu cách nhau ít trên vòng contour
                                                            # và midpoint cuối/đầu cũng nằm gần nhau trong ảnh
                                                            # => cho phép nối.
        closed = False
        if len(chains) > 1 and _wrap_connects(chains[-1], chains[0]):
            chains[0] = chains.pop() + chains[0] # bốc chain cuối, góc vào đầu. Ví dụ A->B->C thì sau đó là C + A -> B
        elif len(chains) == 1 and _wrap_connects(chains[0], chains[0]): # Nếu chỉ có một chain, code kiểm tra xem cuối chain đó có nối lại được với đầu chain đó không.
            closed = True

        for ch in chains:
            if len(ch) < 3:
                continue
            pts = [(mx, my) for _, mx, my, _ in ch] # chỉ lấy tọa độ midpoint, bỏ index và clearance
            if closed:
                pts.append(pts[0])
            paths.append(pts)

    print(f"  Raw midpoint tracks: {len(paths)}")
    return paths, dt
