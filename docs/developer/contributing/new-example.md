# Creating a new example

The workflow for creating a new example is similar to the workflow for creating a new feature.
A new example is added using a feature branch with additions to the `examples` directory and documentation.

## Checklist

1. Example files
    1. Example files: `README.md`, `Makefile.am`, parameter files, spatial databases, mesh generation scripts.
    2. Add directory in `examples/Makefile.am`
    3. Update `configure.ac` to include example Makefile.
2. User Guide
    1. Include sections: overview, meshing, common, steps, exercises
    2. Update `update-synopsis.sh` (generate list of features from example parameter files) to include example directory.
    3. Use include directive in example Markdown file to add list of features, etc in synopsis file.
    4. Add example to table and toc in `examples/index.md`.
    5. Add files to `docs/Makefile.am`
    6. Check that all Markdown files have each sentence on its own single line.
