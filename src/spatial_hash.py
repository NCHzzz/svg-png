"""
Spatial hash grid for fast nearest-point and range queries on 2D coordinates.

Used throughout the pipeline: dedupe (two-track removal), fragment pruning,
T-connection detection, and junction endpoint registration.
"""

import math
from collections import defaultdict

_CELL = 8.0


class PointGrid:
    """Spatial hash for nearest-point and predicate queries."""
    """
    Chia nhỏ hình thành các hình nhỏ hơn
    Key là tọa độ của ô, value là list các điểm trong ô đó
    
    """
    def __init__(self):
        self.grid = defaultdict(list)

    def add(self, x, y, tag=None):
        self.grid[(int(x // _CELL), int(y // _CELL))].append((x, y, tag))

    # Trả về điểm gần nhất trong bán kính radius, hoặc None nếu không có điểm nào
    def near(self, x, y, radius, exclude_tag=None):
        """Return (dist, x, y, tag) of nearest point within radius, or None. 
        dist = khoảng cách từ (x, y) tới điểm tìm được
        ox, oy = tọa độ điểm tìm được
        tag = thông tin đi kèm điểm đó
        ví dụ: Tìm điểm đã lưu nào cách (100, 200) không quá 10px.
        Nếu tìm thấy: (3.5, 102, 201, tag)
        Nghĩa là điểm (102, 201) cách (100, 200) khoảng 3.5 px.
        """
        gx, gy = int(x // _CELL), int(y // _CELL)
        reach = int(radius // _CELL) + 1 # bán kính tìm kiếm tính theo số ô ( ví dụ radius = 20, _CELL = 8 thì reach = 3 ô )
        best = None
        for dgy in range(-reach, reach + 1):
            for dgx in range(-reach, reach + 1):
                for (ox, oy, tag) in self.grid.get((gx + dgx, gy + dgy), ()):
                    if exclude_tag is not None and tag == exclude_tag:
                        continue
                    d = math.hypot(ox - x, oy - y)
                    if d <= radius and (best is None or d < best[0]):
                        best = (d, ox, oy, tag)
        return best

    def any_match(self, x, y, radius, pred):
        """True if any point within radius satisfies pred(tag).
        Giải quyết câu hỏi: Có điểm nào gần (x, y) và tag của điểm đó thỏa điều kiện pred không?
        """
        gx, gy = int(x // _CELL), int(y // _CELL)
        reach = int(radius // _CELL) + 1
        rr = radius * radius
        for dgy in range(-reach, reach + 1):
            for dgx in range(-reach, reach + 1):
                for (ox, oy, tag) in self.grid.get((gx + dgx, gy + dgy), ()):
                    if (ox - x) ** 2 + (oy - y) ** 2 <= rr and pred(tag):
                        return True
        return False
