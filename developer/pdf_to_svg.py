#!/usr/bin/env nemesis

import re
import sys
import pymupdf

# PyMuPDF renders each rounded box as a masked <image> (colored fill + text)
# drawn on top of a <path> that has NO fill attribute. Such a path defaults to
# fill="black" in SVG and, clipped to the rounded-rect mask, shows up as the
# black border around the box. The drop shadow is a separate sibling path with
# fill="#666666", so leaving fills alone everywhere except those border paths
# removes the border while preserving the shadow.
#
# We target only the offending paths: a <path> with no fill attribute that is a
# direct child of a <g mask="url(#mask_...)"> group. Each such path gets
# fill="none".
BORDER_PATH_RE = re.compile(
    r'(<g mask="url\(#mask_\d+\)">\s*)(<path(?![^>]*\bfill=)[^>]*?)(/>)'
)


def remove_box_borders(svg):
    """Make box borders invisible while keeping drop shadows."""
    return BORDER_PATH_RE.sub(r'\1\2 fill="none"\3', svg)


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
        svg_out.write(remove_box_borders(page.get_svg_image()))
