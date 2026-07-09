"""
Path fitting utilities: RDP simplify, straight-line detection, resample,
smooth, axis snap, spike removal. Pure Python, no third-party libs.
File này dùng để làm sạch, làm mượt, đơn giản hóa và chuẩn hóa path để SVG output đẹp hơn.
"""

import math


def resample_polyline(poly, spacing=4.0):
    """Resample polyline at regular spacing. Preserve endpoints."""
    if len(poly) < 2: # nếu polyline có ít hơn 2 điểm thì trả về polyline gốc
        return poly
    dists = [0.0]
    for i in range(1, len(poly)): # tính khoảng từ điểm đầu tới từng điểm trong polyline
        dx = poly[i][0] - poly[i - 1][0]
        dy = poly[i][1] - poly[i - 1][1]
        dists.append(dists[-1] + math.hypot(dx, dy)) # lưu khoảng cách từ điểm đầu tới điểm i vào dists
    total = dists[-1] # tổng chiều dài polyline
    if total < spacing: # đường ngắn hơn khoảng cách cần lấy mẫu, thì không cần thêm điểm giữa. 
        return [poly[0], poly[-1]]
    result = [poly[0]]
    pos = spacing # bắt đầu lấy mẫu từ khoảng cách spacing
    si = 0 # index của đoạn hiện tại trong polyline gốc. Ví dụ: si = 0: poly[0] -> poly[1] , si = 1: poly[1] -> poly[2]
    while pos < total:
        while si < len(dists) - 1 and dists[si + 1] < pos: # tìm đoạn hiện tại trong polyline gốc mà khoảng cách từ điểm đầu tới điểm si+1 lớn hơn pos
            si += 1
        if si >= len(poly) - 1:
            break
        # tính t, x, y của điểm mới trên đoạn polyline[si] -> polyline[si+1] sao cho khoảng cách từ điểm đầu tới điểm mới là pos
        t = (pos - dists[si]) / (dists[si + 1] - dists[si]) if dists[si + 1] > dists[si] else 0
        x = poly[si][0] + t * (poly[si + 1][0] - poly[si][0])
        y = poly[si][1] + t * (poly[si + 1][1] - poly[si][1])
        result.append((x, y))
        pos += spacing
    result.append(poly[-1])
    return result


def smooth_polyline(poly, window=3):
    """Moving-average smooth preserving endpoints."""
    """
    Hàm này dùng để làm mượt polyline bằng cách lấy trung bình các điểm xung quanh mỗi điểm.
    Nó giữ nguyên điểm đầu và điểm cuối, chỉ sửa các điểm ở giữa.
    """
    # window: số lượng điểm xung quanh mỗi điểm để lấy trung bình. Ví dụ window=3 thì lấy 1 điểm trước, 1 điểm sau và điểm hiện tại.
    if window <= 1 or len(poly) <= 2 or len(poly) <= window:
        return poly
    result = [poly[0]]
    half = window // 2
    for i in range(1, len(poly) - 1): # duyệt các điểm ở giữa
        s = max(0, i - half) # vị trí bắt đầu lấy mẫu, không vượt quá đầu polyline
        e = min(len(poly), i + half + 1) # vị trí kết thúc lấy mẫu, không vượt quá cuối polyline
        sx = sum(p[0] for p in poly[s:e])  # tính tổng x của các điểm trong cửa sổ
        sy = sum(p[1] for p in poly[s:e]) # tính tổng y của các điểm trong cửa sổ
        n = e - s # số lượng điểm trong cửa sổ
        result.append((sx / n, sy / n))
    result.append(poly[-1])
    return result

# tính khoảng cách từ điểm P(px, py) đến đoạn thẳng AB(ax, ay) -> (bx, by)
def _perp_distance(px, py, ax, ay, bx, by):
    """Perpendicular distance from point P to line AB."""
    vx, vy = bx - ax, by - ay
    wx, wy = px - ax, py - ay
    vv = vx * vx + vy * vy # bình phương độ dài đoạn AB
    if vv < 1e-12:
        return math.hypot(px - ax, py - ay) # nếu gần bằng 0, A và B gần như trùng nhau, trả về khoảng cách từ P đến A
    t = (wx * vx + wy * vy) / vv # t cho biết hình chiếu của P nằm ở đâu trên đoạn AB.
    """
    t = 0    tại A
    t = 0.5  giữa đoạn AB
    t = 1    tại B
    """
    if t <= 0:
        return math.hypot(px - ax, py - ay)
    if t >= 1:
        return math.hypot(px - bx, py - by)
    # Q = (qx, qy) là điểm trên AB gần P nhất
    qx = ax + t * vx
    qy = ay + t * vy
    return math.hypot(px - qx, py - qy) # trả về khoảng cách từ P đến Q

# thuật toán Ramer-Douglas-Peucker ( giảm số điểm của polyline nhưng vẫn giữ hình dạng chính )
def _rdp(poly, epsilon, start, end, result):
    """Recursive RDP."""
    dmax = 0 # khoảng cách lệnh lớn nhất tìm được
    imax = start # index của điểm có khoảng cách lớn nhất
    for i in range(start + 1, end):
        d = _perp_distance(poly[i][0], poly[i][1],
                           poly[start][0], poly[start][1],
                           poly[end][0], poly[end][1])
        if d > dmax:
            dmax = d
            imax = i
    if dmax > epsilon: # Nếu điểm lệch xa nhất vượt quá sai số cho phépthì đoạn này chưa đủ thẳng
        _rdp(poly, epsilon, start, imax, result) # xử lý đoạn từ start đến điểm lệch xa nhất
        _rdp(poly, epsilon, imax, end, result) # xử lý đoạn từ điểm lệch xa nhất đến end
    else:
        result.append(poly[end])


def rdp_simplify(poly, epsilon=2.0): # epsilon: ngưỡng sai số cho phép. epsilon càng nhỏ → giữ nhiều điểm hơn.
    """Ramer-Douglas-Peucker simplification."""
    if len(poly) < 3:
        return poly
    result = [poly[0]]
    _rdp(poly, epsilon, 0, len(poly) - 1, result)# xử lý điểm đầu đến điểm cuối toàn bộ polyline
    return result # trả về danh sách đã rút gọn


def is_almost_straight(poly, max_error=3.0):
    """Check if polyline approximates a straight line."""
    if len(poly) < 3:
        return True
    ax, ay = poly[0]
    bx, by = poly[-1]
    for i in range(1, len(poly) - 1):
        if _perp_distance(poly[i][0], poly[i][1], ax, ay, bx, by) > max_error:
            return False
    return True


def fit_line_if_straight(poly, max_error=3.0):
    """If polyline is nearly straight, replace with start-end line."""
    if len(poly) < 3:
        return poly
    if is_almost_straight(poly, max_error):
        return [poly[0], poly[-1]]
    return poly

# 
def _axis_span(poly):
    """Return (horizontal_span, vertical_span) of the polyline bbox.
        tính độ trải rộng của path theo trục x và trục y
    """
    xs = [p[0] for p in poly]
    ys = [p[1] for p in poly]
    return max(xs) - min(xs), max(ys) - min(ys)


def snap_axis_aligned(poly, ratio=4.0, max_error=6.0):
    """Snap nearly horizontal or vertical polylines to a constant axis."""
    """
    Hàm này kiểm tra path có gần ngang hoặc gần dọc không. Nếu có, nó ép thẳng lại.
    ratio=4.0: một chiều phải lớn hơn chiều kia ít nhất 4 lần.
    max_error=6.0: điểm được lệch tối đa 6px so với trục trung bình.
    """
    if len(poly) < 2:
        return poly
    hx, vy = _axis_span(poly) # tính độ trải rộng của path theo trục x và trục y
    if hx < 1e-6 and vy < 1e-6: # gần bằng 0, trả về
        return poly
    if hx >= ratio * vy: # gần ngang
        avg_y = sum(p[1] for p in poly) / len(poly) # tính trung bình y của tất cả các điểm
        if all(abs(p[1] - avg_y) <= max_error for p in poly): # Kiểm tra tất cả điểm có lệch khỏi avg_y quá nhiều không
            return [(p[0], avg_y) for p in poly] # Giữ nguyên x, nhưng ép toàn bộ y về avg_y.
    if vy >= ratio * hx:
        avg_x = sum(p[0] for p in poly) / len(poly)
        if all(abs(p[0] - avg_x) <= max_error for p in poly):
            return [(avg_x, p[1]) for p in poly]
    return poly


def remove_spikes(poly, min_leg=10.0, angle_thresh=0.45):
    """Drop short back-and-forth wiggles (common near junction blobs)."""
    if len(poly) < 4:
        return poly
    out = [poly[0]]
    for i in range(1, len(poly) - 1):
        ax, ay = out[-1]
        bx, by = poly[i]
        cx, cy = poly[i + 1]
        ab = math.hypot(bx - ax, by - ay)
        bc = math.hypot(cx - bx, cy - by)
        if ab < min_leg and bc < min_leg:
            v1x, v1y = bx - ax, by - ay # vector từ A đến B
            v2x, v2y = cx - bx, cy - by # vector từ B đến C
            m1 = math.hypot(v1x, v1y) # độ dài vector AB
            m2 = math.hypot(v2x, v2y) # độ dài vector BC
            if m1 > 1e-6 and m2 > 1e-6:
                cos_a = (v1x * v2x + v1y * v2y) / (m1 * m2) # tính cos góc giữa 2 vector AB và BC
                if cos_a < -angle_thresh:
                    continue
        out.append((bx, by))
    out.append(poly[-1])
    return out
