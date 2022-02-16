#!/usr/bin/env nemesis

from importlib import import_module
import inspect
import types
import contextlib
import io
import os
import argparse
import pathlib

from pylith.apps.PetscApplication import PetscApplication
from pythia.pyre.components.Component import Component

class App():
    """Application to generate Markdown documentation for each PyLith component.
    """

    def __init__(self):
        """Constructor."""
        # PyLith app information must be generated manually.
        self.toc = {
            "apps": ["PyLithApp"]
        }

    def main(self):
        """Application entry point.
        """
        args = self._parse_command_line()
        output_path = pathlib.Path(args.output_path)
        src_path = pathlib.Path(args.src_path)
        for file_path in src_path.glob("**/*.py"):
            module_path = self._get_module_path(src_path, file_path)
            self.mod = import_module(module_path)
            cls_name = file_path.stem
            if App._has_class(self.mod, cls_name):
                filename = App._get_filename(output_path, module_path)
                print(f"Working on {filename}...")
                self._process_class(cls_name, filename)
        self._write_toc(output_path)

    def _write_toc(self, output_path):
        """Write table of contents."""
        for key,items in self.toc.items():
            filename = pathlib.Path(key) / "index.md"
            with open(filename, "w") as mfile:
                self._write_subtoc(mfile, key, items)
        filename = output_path / "index.md"
        with open(filename, "w") as mfile:
            self._write_maintoc(mfile, self.toc)
        
    @staticmethod
    def _write_subtoc(mfile, key, items):
        items = sorted(items)
        mfile.write(f"# {key}\n\n")
        mfile.write(":::{toctree}\n")
        mfile.write("---\nmaxdepth: 1\n---\n")
        for item in items:
            mfile.write(f"{item}.md\n")
        mfile.write(":::\n")

    @staticmethod
    def _write_maintoc(mfile, toc):
        mfile.write("(sec-user-components)=\n# PyLith Components\n\n")
        mfile.write(":::{toctree}\n")
        mfile.write("---\nmaxdepth: 2\n---\n")
        keys = sorted(toc.keys())
        for key in keys:
            mfile.write(f"{key}/index.md\n")
        mfile.write(":::\n")

    def _process_class(self, cls_name, filename):
        try:
            cls = getattr(self.mod, cls_name)
            if hasattr(cls, "components"):
                cls_obj = getattr(self.mod, cls_name)()
            else:
                print(f"Skipping {cls_name} object. Not a Pyre component.")
                return
        except:
            print(f"Skipping {cls_name} object. Could not construct object.")
            return
        if isinstance(cls_obj, Component) and not isinstance(cls_obj, PetscApplication):
            file_dir = pathlib.Path(filename).parents[0]
            file_dir.mkdir(parents=True, exist_ok=True)
            if not file_dir.name in self.toc:
                self.toc[file_dir.name] = []
            self.toc[file_dir.name].append(cls_name)
            with open(filename, "w") as mfile:
                self._write_markdown(mfile, cls_obj, cls_name, self.mod.__name__)
                
    @staticmethod
    def _write_markdown(mfile, cls_obj, cls_name, full_name):
        """Show documentation for class with name cls_name in module."""
        mfile.write(f"# {cls_name}\n\n")
        mfile.write(f"Full name: `{full_name}`\n\n")

        if cls_obj.__doc__:
            mfile.write(cls_obj.__doc__)
            mfile.write("\n\n")
        else:
            print(f"Component {cls_name} missing doc string.")

        if App._has_function(cls_obj, "showComponents"):
            with contextlib.redirect_stdout(io.StringIO()) as sout:
                cls_obj.showComponents()
            mfile.write("## Pyre Facilities\n\n")
            mfile.write("```{code-block} bash\n")
            mfile.write(sout.getvalue())
            mfile.write("```\n\n")
        if App._has_function(cls_obj, "showProperties"):
            with contextlib.redirect_stdout(io.StringIO()) as sout:
                cls_obj.showProperties()
            mfile.write("## Pyre Properties\n\n")
            mfile.write("```{code-block} bash\n")
            mfile.write(sout.getvalue())
            mfile.write("```\n")

    def _parse_command_line(self):
        """Parse command line arguments.
        """
        description = "Generate Markdown documentation for Pyre components."
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument("--src-path", action="store", dest="src_path", metavar="PATH", type=str, required=True,
            help="Path to search for Pyre components.")
        parser.add_argument("--output-path", action="store", dest="output_path", metavar="PATH", type=str, default=".",
            help="Path for Markdown output.")
        return parser.parse_args()

    @staticmethod        
    def _get_module_path(src_path, file_path):
        rel_path = file_path.relative_to(src_path.parents[0])
        module_path = ".".join(rel_path.parts)
        return module_path.replace(".py", "")

    @staticmethod
    def _get_filename(output_path, module_path):
        relative_path = module_path.split(".")[1:]
        file_path = output_path / pathlib.Path(*relative_path)
        filename = file_path.as_posix() + ".md"
        return filename

    @staticmethod
    def _has_function(obj, name):
        """Return True if obj has function name."""
        return hasattr(obj, name) and type(inspect.getattr_static(obj, name)) == types.FunctionType

    @staticmethod
    def _has_class(obj, name):
        """Return True if obj has class name."""
        return hasattr(obj, name) and inspect.isclass(inspect.getattr_static(obj, name))


if __name__ == "__main__":
    App().main()