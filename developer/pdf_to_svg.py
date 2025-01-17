#!/usr/bin/env nemesis

import sys
import pymupdf

if len(sys.argv) != 2:
    print("usage: pdf_to_svg.py FILENAME")
    sys.exit(1)

filename_pdf = sys.argv[1]
doc = pymupdf.open(filename_pdf) # open a document
for i_page, page in enumerate(doc):

    if len(doc) > 1:
        filename_single = filename_pdf.replace(".pdf", f"-{i_page}.pdf")
        single_pdf = pymupdf.open()
        single_pdf.insert_pdf(doc, from_page=i_page, to_page=i_page)
        print(f"Writing {filename_single}...")
        single_pdf.save(filename_single)

    filename_svg = filename_pdf.replace(".pdf", ".svg") if len(doc) == 1 else filename_pdf.replace(".pdf", f"-{i_page}.svg")
    with open(filename_svg, "w", encoding="utf-8") as svg_out:
        print(f"Writing {filename_svg}...")
        svg_out.write(page.get_svg_image())
