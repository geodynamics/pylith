#!/usr/bin/env nemesis

import sys
import pymupdf

if len(sys.argv) != 2:
    print("usage: pdf_to_svg.py FILENAME")
    sys.exit(1)

filename_pdf = sys.argv[1]
doc = pymupdf.open(filename_pdf) # open a document
page = doc[0]

filename_svg = filename_pdf.replace(".pdf", ".svg")
with open(filename_svg, "w", encoding="utf-8") as svg_out:
    svg_out.write(page.get_svg_image())
