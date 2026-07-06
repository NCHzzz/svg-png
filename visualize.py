"""
Visual diagnostics (pure Python, writes PNGs):

  python visualize.py overlay [shape...]          mask gray + ref green + ours red
  python visualize.py zoom <shape> x0 y0 x1 y1    same, zoomed to a region
  python visualize.py cover <shape> [out|ref]     shape-reconstruction diff:
                                                  gray covered, red missed, blue overshoot
"""

import sys
import math
import zlib
import struct
from png_decode import decode_png
from centerline_extract import binarize
from test_compare import (parse_svg_paths, path_to_subpaths,
                          _mask_grid, _svg_full_stroke_grid)

SHAPES = ["letter_H", "letter_K", "arrow-turn-down-left",
          "arrow-pointer", "number_3", "number_6", "ampersand"]


def write_png_rgb(filepath, width, height, rows):
    def chunk(tag, data):
        c = struct.pack('>I', len(data)) + tag + data
        return c + struct.pack('>I', zlib.crc32(tag + data) & 0xFFFFFFFF)
    raw = b''.join(b'\x00' + bytes(row) for row in rows)
    png = (b'\x89PNG\r\n\x1a\n'
           + chunk(b'IHDR', struct.pack('>IIBBBBB', width, height, 8, 2, 0, 0, 0))
           + chunk(b'IDAT', zlib.compress(raw))
           + chunk(b'IEND', b''))
    with open(filepath, 'wb') as f:
        f.write(png)


def _render_region(name, x0, y0, x1, y1, out_file, size=512):
    sx = size / (x1 - x0)
    sy = size / (y1 - y0)
    w, h, px = decode_png(f"challenge_sample/{name}.png")
    mask = binarize(px, w, h)
    rows = [bytearray([255] * (size * 3)) for _ in range(size)]
    for y in range(size):
        my = int(y0 + y / sy)
        for x in range(size):
            mx = int(x0 + x / sx)
            if 0 <= mx < w and 0 <= my < h and mask[my][mx]:
                off = x * 3
                rows[y][off:off + 3] = b'\xd8\xd8\xd8'

    def draw(svg, color):
        for d in parse_svg_paths(svg):
            for poly in path_to_subpaths(d, step=2.0):
                for i in range(1, len(poly)):
                    ax, ay = (poly[i-1][0]-x0)*sx, (poly[i-1][1]-y0)*sy
                    bx, by = (poly[i][0]-x0)*sx, (poly[i][1]-y0)*sy
                    steps = max(1, int(math.hypot(bx-ax, by-ay)))
                    for s in range(steps + 1):
                        t = s / steps
                        gx, gy = int(ax + t*(bx-ax)), int(ay + t*(by-ay))
                        for dy in (-1, 0, 1):
                            for dx in (-1, 0, 1):
                                if 0 <= gx+dx < size and 0 <= gy+dy < size:
                                    off = (gx+dx) * 3
                                    rows[gy+dy][off:off+3] = bytes(color)

    draw(f"challenge_sample_results/{name}.svg", (0, 160, 0))
    draw(f"out/{name}.svg", (220, 0, 0))
    write_png_rgb(out_file, size, size, rows)
    print("wrote", out_file)


def overlay(names):
    for name in names:
        _render_region(name, 0, 0, 1024, 1024, f"viz_overlay_{name}.png")


def zoom(name, x0, y0, x1, y1):
    _render_region(name, x0, y0, x1, y1, f"viz_zoom_{name}.png")


def cover(name, which="out"):
    svg = f"out/{name}.svg" if which == "out" else f"challenge_sample_results/{name}.svg"
    size = 256
    m = _mask_grid(name, size)
    s = _svg_full_stroke_grid(svg, size)
    rows = []
    for y in range(size):
        row = bytearray()
        for x in range(size):
            mi, si = m[y * size + x], s[y * size + x]
            if mi and si:
                row += b'\xc8\xc8\xc8'
            elif mi:
                row += b'\xff\x00\x00'
            elif si:
                row += b'\x00\x00\xff'
            else:
                row += b'\xff\xff\xff'
        rows.append(row)
    out_file = f"viz_cover_{name}_{which}.png"
    write_png_rgb(out_file, size, size, rows)
    print("wrote", out_file)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else "overlay"
    if cmd == "overlay":
        overlay(sys.argv[2:] or SHAPES)
    elif cmd == "zoom":
        zoom(sys.argv[2], *map(int, sys.argv[3:7]))
    elif cmd == "cover":
        cover(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "out")
    else:
        print(__doc__)
