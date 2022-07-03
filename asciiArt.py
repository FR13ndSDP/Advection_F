#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import ascii
import sys

if __name__ == "__main__":
    if (len(sys.argv) == 3):
        output = ascii.loadFromFile(sys.argv[1], int(sys.argv[2])) # load a funny picture :)
    elif (len(sys.argv) == 2):
        output = ascii.loadFromFile(sys.argv[1]) # use default image size
    else:
        print("arguments needed")
        sys.exit()

    print(output)
