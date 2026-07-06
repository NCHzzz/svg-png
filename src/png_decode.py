"""
Pure-Python PNG decoder (stdlib only — no PIL/pillow).

Supports common non-interlaced 8-bit PNG formats:
- color type 0: grayscale
- color type 2: RGB
- color type 4: grayscale + alpha
- color type 6: RGBA

Returns (width, height, pixels), where pixels is a list of grayscale rows.
Transparent pixels are composited over white background.
"""

import struct
import zlib


class PNGDecodeError(Exception):
    pass


def _paeth_predictor(a, b, c):
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def _unfilter_scanline(filter_type, scanline, prev_scanline, bpp):
    line = bytearray(scanline)
    if filter_type == 0:  # None
        return line
    if filter_type == 1:  # Sub
        for i in range(len(line)):
            left = line[i - bpp] if i >= bpp else 0
            line[i] = (line[i] + left) & 0xFF
        return line
    if filter_type == 2:  # Up
        for i in range(len(line)):
            up = prev_scanline[i] if prev_scanline is not None else 0
            line[i] = (line[i] + up) & 0xFF
        return line
    if filter_type == 3:  # Average
        for i in range(len(line)):
            left = line[i - bpp] if i >= bpp else 0
            up = prev_scanline[i] if prev_scanline is not None else 0
            line[i] = (line[i] + (left + up) // 2) & 0xFF
        return line
    if filter_type == 4:  # Paeth
        for i in range(len(line)):
            left = line[i - bpp] if i >= bpp else 0
            up = prev_scanline[i] if prev_scanline is not None else 0
            upleft = prev_scanline[i - bpp] if (prev_scanline is not None and i >= bpp) else 0
            line[i] = (line[i] + _paeth_predictor(left, up, upleft)) & 0xFF
        return line
    raise PNGDecodeError(f"Unknown filter type: {filter_type}")


def _bytes_per_pixel(color_type):
    if color_type == 0:
        return 1
    if color_type == 2:
        return 3
    if color_type == 4:
        return 2
    if color_type == 6:
        return 4
    raise PNGDecodeError(f"Unsupported color type: {color_type}")


def _to_gray(row, x, color_type):
    if color_type == 0:
        return row[x]
    if color_type == 2:
        off = x * 3
        return (row[off] + row[off + 1] + row[off + 2]) // 3
    if color_type == 4:
        off = x * 2
        g, a = row[off], row[off + 1]
        # Composite over white background.
        return (g * a + 255 * (255 - a)) // 255
    if color_type == 6:
        off = x * 4
        r, g, b, a = row[off], row[off + 1], row[off + 2], row[off + 3]
        gray = (r + g + b) // 3
        # Composite over white background.
        return (gray * a + 255 * (255 - a)) // 255
    raise PNGDecodeError(f"Unsupported color type: {color_type}")


def decode_png(filepath):
    with open(filepath, 'rb') as f:
        data = f.read()

    if data[:8] != b'\x89PNG\r\n\x1a\n':
        raise PNGDecodeError("Not a valid PNG file")

    pos = 8
    idat_chunks = []
    width = height = bit_depth = color_type = interlace = None

    while pos < len(data):
        if pos + 8 > len(data):
            raise PNGDecodeError("Truncated PNG chunk header")
        length = struct.unpack('>I', data[pos:pos + 4])[0]
        pos += 4
        chunk_type = data[pos:pos + 4].decode('ascii', errors='replace')
        pos += 4
        chunk_data = data[pos:pos + length]
        pos += length
        pos += 4  # CRC

        if chunk_type == 'IHDR':
            if length != 13:
                raise PNGDecodeError(f"IHDR length != 13: {length}")
            width = struct.unpack('>I', chunk_data[0:4])[0]
            height = struct.unpack('>I', chunk_data[4:8])[0]
            bit_depth = chunk_data[8]
            color_type = chunk_data[9]
            interlace = chunk_data[12]
            if bit_depth != 8:
                raise PNGDecodeError(f"Unsupported bit depth: {bit_depth} (expected 8)")
            if color_type not in (0, 2, 4, 6):
                raise PNGDecodeError(f"Unsupported color type: {color_type}")
            if interlace != 0:
                raise PNGDecodeError(f"Unsupported interlace: {interlace} (expected 0)")
        elif chunk_type == 'IDAT':
            idat_chunks.append(chunk_data)
        elif chunk_type == 'IEND':
            break

    if width is None or height is None:
        raise PNGDecodeError("No IHDR chunk found")

    try:
        raw_data = zlib.decompress(b''.join(idat_chunks))
    except zlib.error as e:
        raise PNGDecodeError(f"zlib decompression failed: {e}")

    bpp = _bytes_per_pixel(color_type)
    stride = 1 + width * bpp
    expected_size = height * stride
    if len(raw_data) < expected_size:
        raise PNGDecodeError(f"Decompressed data too short: {len(raw_data)} < {expected_size}")

    pixels = []
    prev_scanline = None
    for y in range(height):
        offset = y * stride
        filter_type = raw_data[offset]
        raw_scanline = raw_data[offset + 1:offset + stride]
        unfiltered = _unfilter_scanline(filter_type, raw_scanline, prev_scanline, bpp)
        gray_row = bytearray(width)
        for x in range(width):
            gray_row[x] = _to_gray(unfiltered, x, color_type)
        pixels.append(gray_row)
        prev_scanline = unfiltered

    return width, height, pixels
