#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

EXAMPLE_FILES = [
    "test",
]

# ----------------------------------------------------------------------
from collections import OrderedDict
COLUMNS = OrderedDict({
    "General": [
        "Dimension",
        "Coordinate system",
        "Mesh generator",
        "Cells",
        "Problem type",
        "Time dependence",
    ],
    "Boundary Conditions": [
        "Dirichlet",
        "Neumann",
        "Absorbing",
        "Point force",
    ],
    "Fault": [
        "Prescribed slip",
        "Slip time function",
        "Constitutive model",
        "Static friction",
        "Slip-weakening friction",
        "Time-weakening friction",
        "Rate-state friction w/ageing",
        "Traction perturbation",
    ],
    "Bulk Rheology": [
        "Linear elastic",
        "Linear Maxwell viscoelastic",
        "Generalixed Maxwell viscoelastic",
        "Powerlaw viscoelastic",
        "Drucker-Prager elastoplastic",
        "Incompressible linear elastic",
        "Porous linear elastic",
        "Stress/strain formulation",
        "Inertia",
        "Gravity",
        "Reference state",
    ],
    "Solver": [
        "Solver",
        "Preconditioner",
        "Time stepping",
    ],
    "Output": [
        "Format",
        "Domain output",
        "Surface output",
        "Point output",
        "State variable output",
        "ParaView",
        "Matplotlib",
    ],
    "Spatial Database": [
        "Uniform",
        "Simple",
        "Simple grid",
        "Composite",
        "Time history",
    ],
})

SCHEMA = {
    "$schema": "http://json-schema.org/schema#",
    "id": "PyLith examples schema",
    "title": "Example",
    "type": "object",
    "properties": {
        "General": {
            "type": "object",
            "properties": {
                "Dimension": { "enum": [2, 3], },
                "Coordinate system": { "enum": ["Cart", "Proj"], },
                "Mesh generator": { "enum": ["ASCII", "CUBIT", "LaGriT"], },
                "Cells": { "enum": ["Tri", "Quad", "Tet", "Hex"], },
                "Refinement": { "enum": ["2x", "4x", "8x"], },
                "Reordering": { "type": "boolean", },
                "Problem type": { "enum": ["TD", "GF"], },
                "Time dependence": { "enum": ["S", "QS", "D"], },
            },
            "required": [
                "Dimension", 
                "Coordinate system",
                "Mesh generator",
                "Cells",
                "Problem type",
                "Time dependence",
            ],
        },
        "Boundary Conditions": {
            "type": "object",
            "properties": {
                "Dirichlet": { "type": "integer", "minimum": 0, },
                "Neumann": { "type": "integer", "minimum": 0, },
                "Absorbing": { "type": "integer", "minimum": 0, },
                "Point force": { "type": "integer", "minimum": 0, },
            },
        },
        "Fault": {
            "type": "object",
            "properties": {
                "Prescribed slip": { "type": "integer", "minimum": 0, },
                "Slip time function": { "enum": ["Step", "Rate", "Liu", "Brune", "User"], },
                "Constitutive model": { "type": "integer", "minimum": 0, },
                "Static friction": { "type": "boolean" },
                "Slip-weakening friction": { "type": "boolean" },
                "Time-weakening friction": { "type": "boolean" },
                "Rate-state friction w/ageing": { "type": "boolean" },
                "Traction perturbation": { "type": "boolean" },
            },
        },
        "Bulk Rheology": {
            "type": "object",
            "properties": {
                "Linear elastic": { "type": "integer", "minimum": 0, },
                "Linear Maxwell viscoelastic": { "type": "integer", "minimum": 0, },
                "Generalized Maxwell viscoelastic": { "type": "integer", "minimum": 0, },
                "Powerlaw viscoelastic": { "type": "integer", "minimum": 0, },
                "Drucker-Prager elastoplastic": { "type": "integer", "minimum": 0, },
                "Incompressible linear elastic": { "type": "integer", "minimum": 0, },
                "Porous linear elastic": { "type": "integer", "minimum": 0, },
                "Stress/stress formulation": { "enum": ["Inf","Fin"], },
                "Inertia": { "type": "boolean" },
                "Reference state": { "type": "boolean" },
                "Gravity": { "type": "boolean" },
            },
        },
        "Solver": {
            "type": "object",
            "properties": {
                "Solver": { "enum": ["L","NL"], },
                "Preconditioner": { "enum": ["ILU", "ASM", "SCHUR", "CUST"], },
                "Time stepping": { "enum": ["BE", "FE", "RK"], },
            },
            "required": ["Solver"],
        },
        "Output": {
            "type": "object",
            "properties": {
                "Format": { "enum": ["VTK", "H5", "H5Ext"], },
                "Domain output": { "type": "integer", "minimum": 0, },
                "Surface output": { "type": "integer", "minimum": 0, },
                "Point output": { "type": "integer", "minimum": 0, },
                "State variable output": { "type": "integer", "minimum": 0, },
                "ParaView": { "type": "boolean", },
                "Matplotlib": { "type": "boolean", },
            },
            "required": ["Format"],
        },
        "Spatial Database": {
            "type": "object",
            "properties": {
                "Uniform": { "type": "integer", "minimum": 0, },
                "Simple": { "type": "integer", "minimum": 0, },
                "Simple grid": { "type": "integer", "minimum": 0, },
                "Time history": { "type": "integer", "minimum": 0, },
            },
        },
    },
}

# ----------------------------------------------------------------------
class Table(object):
    """Abstract base class for a table of example features.
    
    \usepackage[table]{xcolor}
    
    \rowcolors{2}{gray!25}{white}
    \begin{tabular}{cc}
    \rowcolor{gray!50}
    Table head & Table head\\
    Some values & Some values\\
    Some values & Some values\\
    Some values & Some values\\
    Some values & Some values\\
    \end{tabular}

    """
    
    def __init__(self, fout, columns):
        """Constructor.

        :param fout: File object for LaTeX output.
        """
        self.fout = fout
        self.columns = columns
        return

    def writeHeader(self):
        """Write table header.
        """
        f = self.fout

        f.write("\\rowcolors{2}{yellow!30}{white}\n")
        ctags = ["|l|%% Example"]
        for category in self.columns:
            cols = COLUMNS[category]
            ctags.append("*{%d}c|%% %s" % (len(cols), category))
        f.write("\\begin{tabular}{%s\n}\n" % "\n    ".join(ctags))
        f.write("\\hline\n")

        f.write("\\rowcolor{blue!10}\n")
        f.write("Example\n")
        for category in self.columns:
            cols = COLUMNS[category]
            f.write("& \multicolumn{%d}{c|}{%s}\n" % (len(cols), category))
        f.write("\\\\ \n")
        f.write("%%%%\n")
        f.write("\\hline\n")
        f.write("\\rowcolor{blue!10}\n")
        f.write("\n")
        for category in self.columns:
            f.write("%% %s\n" % category)
            cols = COLUMNS[category]
            for m in cols:
                f.write("& \\rlabel{%s}\n" % m)
        f.write("\\\\\n")
        f.write("\\hline\n")
        return

    def writeRow(self, label, example):
        """Write table row corresponding to features for a given example.
        
        :param label: Label for example.
        :param example: Features for a given example.
        """
        f = self.fout
        fromBoolean = lambda v: "\yes" if v else ""
        fromInt = lambda v: "%dx" % v

        f.write("%s\n" % label)
        for category in self.columns:
            cols = COLUMNS[category]
            if not category in example:
                for col in cols:
                    f.write("& ")
                continue
            for col in cols:
                if not col in example[category]:
                    f.write("& ")
                    continue
                value = example[category][col]
                colSchema = SCHEMA["properties"][category]["properties"][col]
                s = str(value)
                if "type" in colSchema:
                    if colSchema["type"] == "boolean":
                        s = fromBoolean(value)
                    elif colSchema["type"] == "integer":
                        s = fromInt(value)
                f.write("& %s " % s)
            f.write("\n")
        f.write("\\\\ \\hline\n")
        return

    def writeFooter(self):
        """Write table footer.
        """
        self.fout.write("\end{tabular}\n")
        return

    @staticmethod
    def writeDocHeader(fout):
        """Write document header.
        """
        fout.write("\\documentclass[10pt]{standalone}\n")
        fout.write("\\usepackage{pifont}\n")
        fout.write("\\usepackage{graphicx}\n")
        fout.write("\\usepackage[table]{xcolor}\n")
        fout.write("\\newcommand{\\rlabel}[1]{\\rotatebox[origin=l]{90}{#1}}\n")
        fout.write("\\newcommand{\\yes}{\\ding{52}}")
        fout.write("\\begin{document}\n")
        return

    @staticmethod
    def writeDocFooter(fout):
        """Write document footer.
        """
        fout.write("\\end{document}\n")
        return

# ----------------------------------------------------------------------
class App(object):
    """Application for generating information about examples.
    """

    def __init__(self):
        """Constructor.
        """
        self.examples = None
        return

    def initialize(self):
        """Initialize application.
        """
        import importlib
        self.examples = OrderedDict()
        for filename in EXAMPLE_FILES:
            m = importlib.import_module(filename)
            for label, example in m.examples.items():
                self.examples[label] = example
        return

    def validate(self):
        """Validate using JSON schema.
        """
        import jsonschema
        ok = True
        for label, example in self.examples.items():
            try:
                jsonschema.validate(example, SCHEMA)
                print("Example '%s': OK" % label)
            except jsonschema.exceptions.ValidationError as err:
                print("Example '%s': ERROR" % label)
                print(str(err))
                ok = False
        return ok

    def tabulate(self, filename, standalone=False):
        """Generate LaTeX table with list of features covered by each example.

        :param filename: Name of file for output with LaTeX code.
        """
        tableColumns = [
            [
                "General",
                "Solver",
                "Boundary Conditions",
                "Fault"
            ],
            [
                "Bulk Rheology",
                "Output",
                "Spatial Database",
            ],
        ]        
        with open(filename, "w") as fout:
            if standalone:
                Table.writeDocHeader(fout)
            for columns in tableColumns:
                table = Table(fout, columns)
                table.writeHeader()
                for label, example in self.examples.items():
                    table.writeRow(label, example)
                table.writeFooter()
            if standalone:
                Table.writeDocFooter(fout)
        return

    def listFeatures(self, exname):
        """Generate LaTeX itemized list of features covered by a given example.
        
        LaTeX code for itemized list is written to stdout.

        :param exame: Name of example in examples list.
        """
        features = self.examples[exname]

        print("EXAMPLE %s" % exname)
        print("\\begin{itemize}")
        for key, value in features.items():
            print("\\item %s %s" % (value, key))
        print("\\end{itemize}")
        return

# ======================================================================
if __name__ == "__main__":
    import argparse

    DESCRIPTION = "Application for generating LaTeX code for (1) a table of features covered in each example and (2) an itemized list of features covered in any given example."

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--validate-only", action="store_true", dest="validateOnly")
    parser.add_argument("--table", action="store_true", dest="table")
    parser.add_argument("--table-filename", action="store", dest="tableFilename", default="example_table.tex")
    parser.add_argument("--table-standalone", action="store_true", dest="tableStandalone")
    parser.add_argument("--list-features", action="store", dest="features", default=None)
    args = parser.parse_args()

    app = App()
    app.initialize()
    ok = app.validate()
    if not args.validateOnly and ok:
        if args.table:
            app.tabulate(args.tableFilename, args.tableStandalone)
        if args.features:
            fexs = args.features.split(",")
            for ex in fexs:
                app.listFeatures(ex)

# End of file
