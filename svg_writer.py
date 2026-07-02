"""
Simplified SVG writer. No more connector hacks — topology fixed in graph phase.
"""

import math


def _path_d(edge):
    """Build SVG path d string from polyline points."""
    if len(edge) < 2:
        return ""
    parts = [f'M{edge[0][0]:.3f} {edge[0][1]:.3f}']
    for (x, y) in edge[1:]:
        parts.append(f'L{x:.3f} {y:.3f}')
    return ''.join(parts)


def write_svg(filepath, edges, width=1024, height=1024, stroke_width=45,
              stroke_color="#000000"):
    """Write SVG file with one <path> per centerline polyline.

    Graph phase already handles junctions and endpoint snapping.
    SVG writer just renders with round caps — topology is clean.
    """
    edges = [e for e in edges if len(e) >= 2]

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" preserveAspectRatio="xMidYMid meet">',
        ' <g id="centerline-shapes">',
    ]

    for idx, edge in enumerate(edges):
        lines.append(
            f'  <path id="path-{idx + 1}" d="{_path_d(edge)}" '
            f'stroke="{stroke_color}" stroke-width="{stroke_width}" '
            f'fill="none" stroke-linecap="round" stroke-linejoin="round" />'
        )

    lines.append(' </g>')
    lines.append('</svg>')

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
