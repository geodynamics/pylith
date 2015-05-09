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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/utils/CppData.py
##
## @brief Python object to create C++ object holding test data.
##
## Useful in unit testing of C++ objects where test data is generated
## with Python code.
##
## If parent property is set, we assume the object is providing only the
## data, so data is private and object needs a constructor and
## destructor. Otherwise, the object just has public data and no methods.

import numpy

# CppData class
class CppData(object):
  """
  Python object to create C++ object holding test data.
  """
  scalars = []
  arrays = []


  def __init__(self, namespace, objname, parent=None, header="header.hh"):
    """
    Constructor.
    """
    self.namespace = namespace.split(",")
    self.objname = objname
    self.parent = parent
    self.creator = "Unknown"
    self.filenameHeader = header

    return


  def run(self):
    self._compute()
    self._collect()
    self._write()
    return


  # PROTECTED METHODS //////////////////////////////////////////////////

  def _compute(self):
    return


  def _collect(self):
    raise NotImplementedError("Implement _collect() in child class.")
    return


  def _write(self):
    """
    Write header and implementation file.
    """
    self._writeHeader()
    self._writeSource()
    return


  def _addScalar(self, vtype, name, value, format):
    """
    Add scalar to object's members.
    """
    data = {'type': vtype,
            'name': name,
            'value': value,
            'format': format}
    self.scalars.append(data)
    return


  def _addArray(self, vtype, name, values, format, ncols):
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

  def _writeHeader(self):
    """
    Write C++ header file.
    """

    filename = "%s.hh" % self.objname

    fileOut = open(filename, "w")
    self._insertHeader(fileOut)

    # Write define information
    hhname = "%s_%s_hh" % ("_".join(self.namespace), self.objname.lower())
    fileOut.write("#if !defined(%s)\n" % hhname)
    fileOut.write("#define %s\n" % hhname)
    fileOut.write("\n")

    if self.parent:
      fileOut.write("#include \"%s.hh\"\n\n" % self.parent)

    # Write namespace information
    level = 0
    for name in self.namespace:
      fileOut.write("  "*level)
      fileOut.write("namespace %s {\n" % name)
      level += 1
    fileOut.write("  "*level)
    fileOut.write("class %s;\n" % self.objname)
    level -= 1
    for name in self.namespace:
      fileOut.write("  "*level)
      fileOut.write("} // %s\n" % name)
      level -= 1
    fileOut.write("\n")

    # Write class opening
    fileOut.write("class %s::%s" % ("::".join(self.namespace), self.objname))
    if self.parent:
      fileOut.write(" : public %s\n" % self.parent)
    else:
      fileOut.write("\n")
    fileOut.write("{\n\n")

    # Write data
    if self.parent:
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
      fileOut.write("  static const %s %s;\n\n" % (scalar['type'], scalar['name']))

    # Write array information
    for array in self.arrays:
      if array['values'] is None:
        fileOut.write("  static const %s* %s;\n\n" % (array['type'], array['name']))
      else:
        fileOut.write("  static const %s %s[];\n\n" % (array['type'], array['name']))

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
      if not isinstance(scalar['value'], bool):
        value = scalar['value']
      else:
        value = str(scalar['value']).lower()
      fileOut.write(cppformat % (scalar['type'], "::".join(self.namespace), self.objname, scalar['name'], value))

    # Write array information
    for array in self.arrays:
      if array['values'] is None:
        cppformat = "const %s* %s::%s::%s = NULL;\n\n"
        fileOut.write(cppformat % (array['type'], "::".join(self.namespace), self.objname, array['name']))
      else:
        cppformat = "const %s %s::%s::%s[] = {\n"
        fileOut.write(cppformat % (array['type'], "::".join(self.namespace), self.objname, array['name']))
        icol = 0
        for value in numpy.ravel(array['values']):
          cppformat = "%s," % array['format']
          fileOut.write(cppformat % value)
          icol += 1
          if icol == array['ncols']:
            fileOut.write("\n")
            icol = 0
        fileOut.write("};\n\n")

    if self.parent:
      self._writeLifecycle(fileOut)

    self._insertFooter(fileOut)
    fileOut.close()
    return


  def _writeLifecycle(self, fileOut):
    """
    Write default constructor and destructor.
    """
    # Constructor
    fileOut.write("%s::%s::%s(void)\n" % ("::".join(self.namespace), self.objname, self.objname))

    fileOut.write("{ // constructor\n")

    for scalar in self.scalars:
      n = scalar['name']
      if "char*" != scalar['type']:
        fileOut.write("  %s = %s;\n" % (n[1:], n))
      else:
        fileOut.write("  %s = const_cast<char*>(%s);\n" % (n[1:], n))
    for array in self.arrays:
      n = array['name']
      fileOut.write("  %s = const_cast<%s*>(%s);\n" % (n[1:], array['type'], n))
    fileOut.write("} // constructor\n\n")

    # Destructor
    fileOut.write("%s::%s::~%s(void)\n" % ("::".join(self.namespace), self.objname, self.objname))
    fileOut.write("{}\n\n")
    return


  def _insertHeader(self, fileOut):
    """
    Insert header file into output file.
    """
    fileIn = open(self.filenameHeader, "r")
    for line in fileIn:
      fileOut.write(line)
    fileIn.close()

    fileOut.write("// DO NOT EDIT THIS FILE\n")
    fileOut.write("// This file was generated from python application "
                  "%s.\n\n" % self.creator)
    return


  def _insertFooter(self, fileOut):
    """
    Insert footer into output file.
    """
    fileOut.write("\n// End of file\n")
    return


# End of file 
