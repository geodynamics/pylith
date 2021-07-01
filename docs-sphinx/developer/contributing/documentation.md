# Documentation

:::{note}
These instructions apply to the new Sphinx-based documentation.
:::

## Building the documentation

### Prerequisites

You must have the following Python packages installed:

* **sphinx** (v3.5 or later)
* **myst-parser** (0.14.0 or later)
* **pydata-sphinx-theme** (0.6.2 or later)

### Generating the documentation

Use the `build.sh` script in the `docs-sphinx` directory to generate the documentation.
The default format is html, which will place the files in the `_build/html` directory.
PDF and epub formats can also be generated.

```{code-block} console
cd docs-sphinx

# Generate the documentation as html.
$ ./build.sh

# Generate the documentation as a PDF.
$ ./build.sh latex
$ cd _build/latex && make

# Generate the documentation as an epub.
$ ./build.sh epub
```

## Contributing to the documentation

### Style guide

1. Use Markedly Structure Text (MyST), not reStructure Text (rST).
2. Place each sentence on its own single line.
