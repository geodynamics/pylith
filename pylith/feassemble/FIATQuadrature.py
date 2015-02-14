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

## @file pylith/feassemble/FIATQuadrature.py
##
## @brief Python object for special FIAT quadrature schemes.

from FIAT.quadrature import QuadratureRule
class CollocatedQuadratureRule(QuadratureRule):
  """
  Quadrature points colocated with vertices.
  """
  def __init__(self, ref_el, m):
    from FIAT.lagrange import Lagrange
    vertices = Lagrange(ref_el, m).dual.get_nodes()
    pts = [v.get_point_dict().keys()[0] for v in vertices]
    npts = len(pts)
    wts = (ref_el.volume()/npts,)*npts
    
    QuadratureRule.__init__(self, ref_el, pts, wts)
    return


# End of file 
