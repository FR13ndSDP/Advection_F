import math
from ascii import colorize

pixels = ".,:;i1tfLCG08@"
for char in pixels:
    pixels = pixels.replace(char, chr(ord(char) + 0xFEE0))


def asciify(r, g, b, a):
    rawChar = getRawChar(r, g, b, a)

    return colorize.colorize(rawChar, (r, g, b))


def getRawChar(r, g, b, a):
    value = intensity(r, g, b, a)
    precision = 255 * 3 / (len(pixels) - 1)
    rawChar = pixels[int(round(value / precision))]

    return rawChar


def intensity(r, g, b, a):
    return (r + g + b) * a
