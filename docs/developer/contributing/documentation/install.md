# Building the documentation

You can build a local copy of the documentation using Sphinx.

## Prerequisites

You do not have to install PyLith to generate a local copy of the documentation.
Nevertheless, you must have Python and the following Python packages installed:

- **sphinx** (v3.5 or later)
- **myst-parser** (0.14.0 or later)
- **pydata-sphinx-theme** (0.6.2 or later)
- **sphinx-copybutton**
- **sphinxcontrib.bibtex**
- **sphinx_design**

These packages are all included in the PyLith development Docker container.

## Generating the documentation

Use the `build.sh` script in the `docs` directory to generate the documentation.
The default format is html, which will place the files in the `_build/html` directory.
PDF and epub formats can also be generated.

```{code-block} bash
cd docs

# Generate the documentation as html.
./build.sh

# Generate the documentation as a PDF.
./build.sh latex
cd _build/latex && make

# Generate the documentation as an epub.
./build.sh epub
```
