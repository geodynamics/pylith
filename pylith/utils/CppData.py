#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ----------------------------------------------------------------------
#

## @file pylith/utils/CppData.py
##
## @brief Python object to create C++ object holding data values.
##
## Useful in unit testing of C++ objects where data is generate with
## Python code.
##
## If parent property is set, we assume object is providing only the
## data, so data is private and object needs a constructor and
## destructor. Otherwise, object just has public data and no methods.
##
## Factory: cpp_data

from pyre.components.Component import Component

import string
import numpy

# CppData class
class CppData(Component):
  """
  Python objec to create C++ object holding data values.
  """

  class Inventory(Component.Inventory):
    """
    Python object for managing CppData facilities and properties.

    Useful in unit testing of C++ objects where data is generate with
    Python code.

    If parent property is set, we assume object is providing only the
    data, so data is private and object needs a constructor and
    destructor. Otherwise, object just has public data and no methods.

    Factory: cpp_data
    """

    ## @class Inventory
    ## Python object for managing CppData facilities and properties.
    ##
    ## \b Properties
    ## @li \b header Filename for header for C++ files
    ## @li \b objname Name of C++ object
    ## @li \b namespace Tuple of strings forming namespace for object
    ## @li \b parent Name of parent object
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    header = pyre.inventory.str("header", default="header.hh")
    header.meta['tip'] = "Filename for header for C++ files."

    objname = pyre.inventory.str("object", default="object")
    objname.meta['tip'] = "Name of C++ object"

    namespace = pyre.inventory.list("namespace", default=["pylith"])
    namespace.meta['tip'] = "Tuple of strings forming namspace for object."

    parent = pyre.inventory.str("parent", default="")
    parent.meta['tip'] = "Name of parent object."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="cppdata"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="cpp_data")
    self.scalars = []
    self.arrays = []
    self.app = ""
    return


  def write(self, app):
    """
    Write header and implementation file.
    """
    self.app = app
    self._writeHeader()
    self._writeSource()
    return


  def addScalar(self, vtype, name, value, format):
    """
    Add scalar to object's members.
    """
    data = {'type': vtype,
            'name': name,
            'value': value,
            'format': format}
    self.scalars.append(data)
    return


  def addArray(self, vtype, name, values, format, ncols):
    """
    Add array to object's members.
    """
    data = {'type': vtype,
            'name': name,
            'values': values,
            'format': format,
            'ncols': ncols}
    self.arrays.append(data)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.header = self.inventory.header
    self.objname = self.inventory.objname
    self.namespace = self.inventory.namespace
    self.parent = self.inventory.parent
    return


  def _writeHeader(self):
    """
    Write C++ header file.
    """

    filename = "%s.hh" % self.objname

    fileOut = open(filename, "w")
    self._insertHeader(fileOut)

    # Write define information
    hhname = "%s_%s_hh" % (string.join(self.namespace, "_"),
                           self.objname.lower())
    fileOut.write("#if !defined(%s)\n" % hhname)
    fileOut.write("#define %s\n" % hhname)
    fileOut.write("\n")

    if self.parent != "":
      fileOut.write("#include \"%s.hh\"\n\n" % self.parent)

    # Write namespace information
    level = 0
    for name in self.namespace:
      fileOut.write("%s" % string.join(["  "]*level))
      fileOut.write("namespace %s {\n" % name)
      level += 1
    fileOut.write("%s" % string.join(["  "]*level))
    fileOut.write("class %s;\n" % self.objname)
    level -= 1
    for name in self.namespace:
      fileOut.write("%s" % string.join(["  "]*level))
      fileOut.write("} // %s\n" % name)
      level -= 1
    fileOut.write("\n")

    # Write class opening
    fileOut.write("class %s::%s" % \
                  (string.join(self.namespace, "::"),
                   self.objname))
    if self.parent != "":
      fileOut.write(" : public %s\n" % self.parent)
    else:
      fileOut.write("\n")
    fileOut.write("{\n\n")

    # Write data
    if self.parent != "":
      fileOut.write("public: \n\n")
      fileOut.write("  /// Constructor\n")
      fileOut.write("  %s(void);\n\n" % self.objname)
      fileOut.write("  /// Destructor\n")
      fileOut.write("  ~%s(void);\n\n" % self.objname)

      fileOut.write("private:\n\n")
    else:
      fileOut.write("public:\n\n")

    # Write scalar information
    for scalar in self.scalars:
      fileOut.write("  static const %s %s;\n\n" % \
                    (scalar['type'], scalar['name']))

    # Write array information
    for array in self.arrays:
      fileOut.write("  static const %s %s[];\n\n" % \
                    (array['type'], array['name']))


    # Write class closing
    fileOut.write("};\n\n")
    fileOut.write("#endif // %s\n" % hhname)

    self._insertFooter(fileOut)
    fileOut.close()
    return


  def _writeSource(self):
    """
    Write C++ source file.
    """

    filename = "%s.cc" % self.objname

    fileOut = open(filename, "w")
    self._insertHeader(fileOut)

    fileOut.write("#include \"%s.hh\"\n" % self.objname)
    fileOut.write("\n")

    # Write scalar information
    for scalar in self.scalars:
      cppformat = "const %s %s::%s::%s = " + scalar['format'] + ";\n\n"
      fileOut.write(cppformat % \
                    (scalar['type'],
                     string.join(self.namespace, "::"), self.objname,
                     scalar['name'], scalar['value']))

    # Write array information
    for array in self.arrays:
      cppformat = "const %s %s::%s::%s[] = {\n"
      fileOut.write(cppformat % \
                    (array['type'],
                     string.join(self.namespace, "::"), self.objname,
                     array['name']))
      icol = 0
      for value in numpy.ravel(array['values']):
        cppformat = "%s," % array['format']
        fileOut.write(cppformat % value)
        icol += 1
        if icol == array['ncols']:
          fileOut.write("\n")
          icol = 0
      fileOut.write("};\n\n")

    if self.parent != "":
      self._writeLifecycle(fileOut)

    self._insertFooter(fileOut)
    fileOut.close()
    return


  def _writeLifecycle(self, fileOut):
    """
    Write default constructor and destructor.
    """
    # Constructor
    fileOut.write("%s::%s::%s(void)\n" % \
                  (string.join(self.namespace, "::"),
                   self.objname, self.objname))

    fileOut.write("{ // constructor\n")

    for scalar in self.scalars:
      n = scalar['name']
      fileOut.write("  %s = %s;\n" % (n[1:], n))
    for array in self.arrays:
      n = array['name']
      fileOut.write("  %s = const_cast<%s*>(%s);\n" % \
                    (n[1:], array['type'], n))
    fileOut.write("} // constructor\n\n")

    # Destructor
    fileOut.write("%s::%s::~%s(void)\n" % \
                  (string.join(self.namespace, "::"),
                   self.objname, self.objname))
    fileOut.write("{}\n\n")
    return


  def _insertHeader(self, fileOut):
    """
    Insert header file into output file.
    """
    fileIn = open(self.header, "r")
    for line in fileIn:
      fileOut.write(line)
    fileIn.close()

    fileOut.write("// DO NOT EDIT THIS FILE\n")
    fileOut.write("// This file was generated from python application "
                  "%s.\n\n" % self.app)
    return


  def _insertFooter(self, fileOut):
    """
    Insert footer into output file.
    """
    fileOut.write("\n// End of file\n")
    return


# FACTORIES ////////////////////////////////////////////////////////////

def cpp_data():
  """
  Factory associated with CppData.
  """
  return CppData()


# End of file 
