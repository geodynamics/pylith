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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/PointsList.py
##
## @brief Python object for reading list of points from a file.
##
## Factory: output_manager

from pyre.components.Component import Component

# Validator for filename
def validateFilename(value):
  """
  Validate filename with list of points.
  """
  if 0 == len(value):
    raise ValueError("Filename for list of points not specified.")
  return value


# PointsList class
class PointsList(Component):
  """
  Python object for reading a list of points from a file.

  @class Inventory
  Python object for managing PointsList facilities and properties.
  
  \b Properties
  @li \b filename Filename for list of points.
  @li \b comment_delimiter Delimiter for comments.
  @li \b value_delimiter Delimiter used to separate values.
  
  \b Facilities
  @li None

  Factory: points_list
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  filename = pyre.inventory.str("filename", default="", 
                                validator=validateFilename)
  filename.meta['tip'] = "Filename for list of points."

  commentDelimiter = pyre.inventory.str("comment_delimiter", default="#")
  commentDelimiter.meta['tip'] = "Delimiter for comments."

  valueDelimiter = pyre.inventory.str("value_delimiter", default="")
  valueDelimiter.meta['tip'] = "Delimiter used to separate values."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="pointslist"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return


  def read(self):
    """
    Read points from file.
    """
    import numpy
    if len(self.valueDelimiter) == 0:
      points = numpy.loadtxt(self.filename, comments=self.commentDelimiter)
    else:
      points = numpy.loadtxt(self.filename,
                             comments=self.commentDelimiter, 
                             delimiter=self.valueDelimiter)
    ndims = len(points.shape)
    if ndims == 1:
      spaceDim = points.shape[0]
      points = points.reshape((1,spaceDim))
    return points
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    try:
      Component._configure(self)
      self.filename = self.inventory.filename
      self.commentDelimiter = self.inventory.commentDelimiter
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring points list "
                       "(%s):\n%s" % (aliases, err.message))

    return


# FACTORIES ////////////////////////////////////////////////////////////

def points_list():
  """
  Factory associated with PointsList.
  """
  return PointsList()


# End of file 
