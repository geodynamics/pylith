# Visual Studio Code

Visual Studio Code (VS Code) provides a extensive integrated development environment for writing code, source control, building, and running tests.
The PyLith repository is setup to include local settings for VS Code for building, running tests, and debugging.

## Setting up VS Code

The PyLith settings for VS Code make use of the following environment variables, which are automatically set in the PyLith development Docker container.
If you are working in VS Code with local source code, then we recommend setting these environment variables as part of activating the PyLith virtual environment.

| Environment variable | Decription                                              |
| :------------------- | :------------------------------------------------------ |
| `PETSC_DIR`          | Directory for PETSc                                     |
| `PETSC_ARCH`         | Build label for PETSc debugging configuration           |
| `PYLITH_BUILDDIR`    | Top-level directory where we build PyLith [^vscode]     |
| `PYLITH_DIR`         | Directory containing CIG-related dependencies [^vscode] |
| `PYLITHDEPS_DIR`     | Directory containing external dependencies [^vscode]    |
| `PYTHON_INCDIR`      | Directory containing Python header files [^vscode]      |
| `MPI_INCDIR`         | Directory containing MPI header files [^vscode]         |
| `PROJ_INCDIR`        | Directory containing Proj header files [^vscode]        |
| `CATCH_INCDIR`       | Directory containing Catch2 header files [^vscode]      |

## Starting VS Code

If you have PyLith source code on your local machine, then it is best to start VS Code from a terminal in which you have activated the PyLith virtual environment.
This helps VS Code know the location of dependencies.

```{code-block} console
# Change to top-level PyLith source directory
cd src/pylith

# Start VS Code and open the top-level directory
code .
```

If you are connecting to the PyLith development Docker container, then you can start VS Code using a terminal or the application icon.
Select `File`->`Open folder...` and navigate to the top-level PyLith source directory; this ensures VS Code loads the settings from `.vscode` in the top-level PyLith source directory.

## Extensions

We strongly recommend the following extensions to leverage all of the default settings included in the PyLith repository.

- C/C++
- Python
- Uncrustify
- C++ TestMate
- Test Explorer UI
- markdownlint
- MyST-Markdown

Other extensions that you may find useful include:

- Dev Containers
- Live Share
- autoconf
- Code Spell Checker
- Markdown all in One
- Material Icon Theme
- Remote-SSH
- Git Graph
- Git Lens
