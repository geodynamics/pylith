#!/usr/bin/env python3
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

from collections import OrderedDict

EXAMPLE_FILES = [
    "subduction3d",
]

# ----------------------------------------------------------------------


class Features(object):
    CS_CART = "Cart"
    CS_PROJ = "Proj"

    CELL_TRI = "Tri"
    CELL_QUAD = "Quad"
    CELL_TET = "Tet"
    CELL_HEX = "Hex"

    MESH_ASCII = "ASCII"
    MESH_CUBIT = "CUBIT"
    MESH_LAGRIT = "LaGriT"

    REFINE_2 = "2x"
    REFINE_4 = "4x"
    REFINE_8 = "8x"

    PROB_TIMEDEPENDENT = "TD"
    PROB_GREENSFNS = "GF"

    TIMEDEP_STATIC = "S"
    TIMEDEP_QUASISTATIC = "QS"
    TIMEDEP_DYNAMIC = "D"

    SOLVER_LINEAR = "L"
    SOLVER_NONLINEAR = "NL"

    TS_BWDEULER = "BE"
    TS_FWDEULER = "FE"
    TS_RUNGEKUTTA = "RK"

    PRECOND_ILU = "ILU"
    PRECOND_ASM = "ASM"
    PRECOND_SCHUR = "SCHUR"
    PRECOND_CUSTOM = "Cust"
    PRECOND_ML = "ML"
    PRECOND_GAMG = "GAMG"
    PRECOND_ML_CUSTOM = "ML+Cust"

    SLIPFN_STEP = "Step"
    SLIPFN_RATE = "Rate"
    SLIPFN_LIU = "Liu"
    SLIPFN_BRUNE = "Brune"
    SLIPFN_TIMEHISTORY = "User"

    STRAINFORM_INFINITESIMAL = "Inf"
    STRAINFORM_FINITE = "Fin"

    OUTPUT_VTK = "VTK"
    OUTPUT_HDF5 = "H5"
    OUTPUT_HDF5EXT = "H5Ext"

    from collections import OrderedDict
    CATEGORIES = [
        ("General",
         # Properties
         OrderedDict((
             ("Dimension", {"enum": [2, 3], }),
             ("Coordinate system", {
                 "enum": [CS_CART, CS_PROJ],
                 "description": "%s: Cartesian, %s: geographic projection" % (CS_CART, CS_PROJ),
             }),
             ("Mesh generator", {
                 "enum": [MESH_ASCII, MESH_CUBIT, MESH_LAGRIT],
                 "description": "%s: ASCII, %s: CUBIT/Trelis, %s: LaGriT" % (MESH_ASCII, MESH_CUBIT, MESH_LAGRIT),
             }),
             ("Cells", {
                 "enum": [CELL_TRI, CELL_QUAD, CELL_TET, CELL_HEX],
             }),
             ("Refinement", {"enum": [REFINE_2, REFINE_4, REFINE_8], }),
             ("Reordering", {"type": "boolean", }),
             ("Problem type", {
                 "enum": [PROB_TIMEDEPENDENT, PROB_GREENSFNS],
                 "description": "%s: time dependent, %s: Green's functions" % (PROB_TIMEDEPENDENT, PROB_GREENSFNS),
             }),
             ("Time dependence", {
                 "enum": [TIMEDEP_STATIC, TIMEDEP_QUASISTATIC, TIMEDEP_DYNAMIC],
                 "description": "%s: static, %s: quasi-static, %s: dynamic" % (TIMEDEP_STATIC, TIMEDEP_QUASISTATIC, TIMEDEP_DYNAMIC),
             }),
         )),
            # Required
            ["Dimension", "Coordinate system",
                "Mesh generator", "Cells", "Problem type"],
         ),
        ("Solver",
         # Properties
         OrderedDict((
             ("Solver", {
                     "enum": [SOLVER_LINEAR, SOLVER_NONLINEAR],
                     "description": "%s: linear, %s: nonlinear" % (SOLVER_LINEAR, SOLVER_NONLINEAR),
                     }),
             ("Preconditioner", {
                 "enum": [PRECOND_ILU, PRECOND_ASM, PRECOND_SCHUR, PRECOND_CUSTOM, PRECOND_ML, PRECOND_GAMG, PRECOND_ML_CUSTOM],
                 "description": "%s: ILU, %s: Additive Schwarz, %s: Schur complement, %s: custom, %s: ML algebraic multigrid, %s: geometric algebraic multigrid" % (PRECOND_ILU, PRECOND_ASM, PRECOND_SCHUR, PRECOND_CUSTOM, PRECOND_ML, PRECOND_GAMG),
             }),

             ("Time stepping", {
                 # "enum": [TS_BWDEULER, TS_FWDEULER, TS_RUNGEKUTTA],
                 # "description": "%s: Backward Euler, %s: Forward Euler, %s: Runge-Kutta" % (TS_BWDEULER, TS_FWDEULER, TS_RUNGEKUTTA),
                 "enum": [TS_BWDEULER, TS_FWDEULER],
                 "description": "%s: Backward Euler, %s: Forward Euler" % (TS_BWDEULER, TS_FWDEULER),
             }),
         )),
            # Required
            ["Solver"],
         ),
        ("Boundary Condition",
         # Properties
         OrderedDict((
             ("Dirichlet", {"type": "integer", "minimum": 0, }),
             ("Neumann", {"type": "integer", "minimum": 0, }),
             ("Absorbing", {"type": "integer", "minimum": 0, }),
             ("Point force", {"type": "integer", "minimum": 0, }),
         )),
            # Required
            None,
         ),
        ("Fault",
         # Properties
         OrderedDict((
             ("Prescribed slip", {"type": "integer", "minimum": 0, }),
             ("Slip time function", {"enum": [
                 SLIPFN_STEP, SLIPFN_RATE, SLIPFN_LIU, SLIPFN_BRUNE, SLIPFN_TIMEHISTORY], }),
             ("Constitutive model", {"type": "integer", "minimum": 0, }),
             ("Static friction", {"type": "boolean"}),
             ("Slip-weakening friction", {"type": "boolean"}),
             ("Time-weakening friction", {"type": "boolean"}),
             ("Rate-state friction w/ageing", {"type": "boolean"}),
             ("Traction perturbation", {"type": "boolean"}),
         )),
            # Required
            None,
         ),
        ("Bulk Rheology",
         # Properties
         OrderedDict((
             ("Linear elastic", {"type": "integer", "minimum": 0, }),
             ("Linear Maxwell viscoelastic", {
                 "type": "integer", "minimum": 0, }),
             ("Generalized Maxwell viscoelastic", {
                 "type": "integer", "minimum": 0, }),
             ("Powerlaw viscoelastic", {"type": "integer", "minimum": 0, }),
             ("Drucker-Prager elastoplastic",
              {"type": "integer", "minimum": 0, }),
             #("Incompressible linear elastic", { "type": "integer", "minimum": 0, }),
             #("Porous linear elastic", { "type": "integer", "minimum": 0, }),
             ("Stress/strain formulation", {
                 "enum": [STRAINFORM_INFINITESIMAL, STRAINFORM_FINITE],
                 "description": "%s: infinitesimal, %s: small, finite strain" % (STRAINFORM_INFINITESIMAL, STRAINFORM_FINITE),
             }),
             ("Inertia", {"type": "boolean"}),
             ("Reference state", {"type": "boolean"}),
             ("Gravity", {"type": "boolean"}),
         )),
            # Required
            ["Stress/strain formulation"],
         ),
        ("Output",
         # Properties
         OrderedDict((
             ("Format", {
                     "enum": [OUTPUT_VTK, OUTPUT_HDF5, OUTPUT_HDF5EXT],
                     "description": "%s: VTK, %s: HDF5, %s: HDF5 w/external datasets" % (OUTPUT_VTK, OUTPUT_HDF5, OUTPUT_HDF5EXT),
                     }),
             ("Domain output", {"type": "integer", "minimum": 0, }),
             ("Surface output", {"type": "integer", "minimum": 0, }),
             ("Point output", {"type": "integer", "minimum": 0, }),
             ("State variable output", {
                 "type": "integer", "minimum": 0, }),
             ("ParaView", {"type": "boolean", }),
             ("Matplotlib", {"type": "boolean", }),
         )),
            # Required
            ["Format"],
         ),
        ("Spatial Database",
         # Properties
         OrderedDict((
             ("Uniform", {"type": "integer", "minimum": 0, }),
             ("Simple", {"type": "integer", "minimum": 0, }),
             ("Simple grid", {"type": "integer", "minimum": 0, }),
             ("Composite", {"type": "integer", "minimum": 0, }),
             ("Time history", {"type": "integer", "minimum": 0, }),
         )),
            # Required
            None,
         )
    ]

    SCHEMA = {
        "$schema": "http://json-schema.org/schema#",
        "id": "PyLith examples schema",
        "title": "Example",
        "type": "object",
        "properties": {},
        "required": ["General", "Solver", "Boundary Conditions", "Bulk Rheology", "Output", "Spatial Database"],
    }
    for category, properties, required in CATEGORIES:
        SCHEMA["properties"][category] = {
            "type": "object",
            "properties": properties,
            "additionalProperties": False}
        if required:
            SCHEMA["properties"][category]["required"] = required


# ----------------------------------------------------------------------
class Table(object):
    """Object tabulating features in examples.
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

        # f.write("\\begin{table}[htbp]\n")
        f.write("\\rowcolors{2}{yellow!30}{white}\n")
        ctags = ["|l|%% Example"]
        for category in self.columns:
            cols = Features.SCHEMA["properties"][category]["properties"]
            ctags.append("*{%d}c|%% %s" % (len(cols), category))
        f.write("\\resizebox{\\textwidth}{!}{%\n")
        f.write("\\begin{tabular}{%s\n}\n" % "\n    ".join(ctags))
        f.write("\\hline\n")

        f.write("\\rowcolor{blue!10}\n")
        f.write("Example\n")
        for category in self.columns:
            cols = Features.SCHEMA["properties"][category]["properties"]
            f.write("& \multicolumn{%d}{c|}{%s}\n" % (len(cols), category))
        f.write("\\\\ \n")
        f.write("%%%%\n")
        f.write("\\hline\n")
        f.write("\\rowcolor{blue!10}\n")
        f.write("\n")
        for category in self.columns:
            f.write("%% %s\n" % category)
            cols = Features.SCHEMA["properties"][category]["properties"]
            for m in cols.keys():
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
        def fromBoolean(v): return "\yes" if v else ""
        def fromInt(v): return "x%d" % v

        f.write("%s\n" % label)
        for category in self.columns:
            cols = Features.SCHEMA["properties"][category]["properties"]
            if not category in example:
                for col in cols:
                    f.write("& ")
                continue
            for col, property in cols.items():
                if not col in example[category]:
                    f.write("& ")
                    continue
                value = example[category][col]
                s = str(value)
                if "type" in property:
                    if property["type"] == "boolean":
                        s = fromBoolean(value)
                    elif property["type"] == "integer":
                        s = fromInt(value)
                f.write("& %s " % s)
            f.write("\n")
        f.write("\\\\ \\hline\n")
        return

    def writeFooter(self):
        """Write table footer.
        """
        f = self.fout

        f.write("\\end{tabular}}\n")
        f.write("\\par\n")

        for category in self.columns:
            cols = Features.SCHEMA["properties"][category]["properties"]
            for label, col in cols.items():
                if "description" in col:
                    f.write("{\\bf %s} -- %s. " % (label, col["description"]))
        f.write("\\\\ \n")
        # f.write("\\end{table}")
        return

    @staticmethod
    def writeDocHeader(fout):
        """Write document header.
        """
        fout.write("\\documentclass[10pt]{article}\n")
        fout.write("\\usepackage{pifont}\n")
        fout.write("\\usepackage{graphicx}\n")
        fout.write("\\usepackage[table]{xcolor}\n")
        fout.write(
            "\\newcommand{\\rlabel}[1]{\\rotatebox[origin=l]{90}{#1}}\n")
        fout.write("\\newcommand{\\yes}{\\ding{52}}\n")
        fout.write("\n\\begin{document}\n")
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
                jsonschema.validate(example, Features.SCHEMA)
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
                "Spatial Database",
            ],
            [
                "Boundary Condition",
                "Fault",
                "Bulk Rheology",
                "Output",
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
    parser.add_argument("--validate-only",
                        action="store_true", dest="validateOnly")
    parser.add_argument("--table", action="store_true", dest="table")
    parser.add_argument("--table-filename", action="store",
                        dest="tableFilename", default="example_table.tex")
    parser.add_argument("--table-standalone",
                        action="store_true", dest="tableStandalone")
    parser.add_argument("--list-features", action="store",
                        dest="features", default=None)
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
