#!/usr/bin/env nemesis

import subprocess
import pymupdf

subprocess.run(["pdflatex", "diagrams.tex"], check=True)

FILENAMES = (
    "geometry",
    "geometry-gmsh",
    "geometry-cubit",
    "step01-diagram",
)

doc = pymupdf.open("diagrams.pdf") # open a document
for i_page, page in enumerate(doc):

    name = FILENAMES[i_page]
    diagram = pymupdf.open()
    diagram.insert_pdf(doc, from_page=i_page, to_page=i_page)
    diagram.save(f"{name}.pdf")

    with open(f"{name}.svg", "w", encoding="utf-8") as svg_out:
      svg_out.write(page.get_svg_image())
