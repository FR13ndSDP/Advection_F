#!/usr/bin/env python
from ascii import asciify
from ascii import colorize


def rgb_to_hex(rgb):
    return "%02x%02x%02x" % rgb


# import image stuff
from PIL import Image


def loadFromFile(filename, columns=60):
    im = Image.open(filename)
    size = im.size
    rows = columns * size[1] / size[0]
    rows = int(round(rows))
    im = im.resize((columns, rows))
    px = im.load()
    size = im.size
    output = ""
    for y in range(0, size[1]):
        for x in range(0, size[0]):
            _px = px[x, y]
            _a = asciify.asciify(_px[0], _px[1], _px[2], 1)
            output = output + _a
        output = output + "\n"
    return output


def onePixel(r, g, b):
    return asciify.asciify(r, g, b, 1)


def color(code, rawChar):
    ansi = colorize.wrapAnsi256(code, 0)
    reset = "\x1b[0m"
    return ansi + rawChar + reset
