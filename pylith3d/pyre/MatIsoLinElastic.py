#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from Material import Material

# Class for isotropic linear elastic material
class MatIsoLinElastic(Material):
  """Python manager for an isotropic linear elastic material."""

  # PUBLIC METHODS /////////////////////////////////////////////////////////////

  def __init__(self, name="matisoelastic"):
    """Constructor."""
    Material.__init__(self, name)

    import pylith3d as pl3d

    # Define dimensions and function pointers.
    # Note that density is also included as a property, even though all material
    # models will have this.
    self.materialModel = {'materialType': "IsotropicLinearElastic",
                          'numProps': 3,
                          'propNames': ["Density", "YoungsModulus", "PoissonsRatio"],
                          'numStateVars': 12,
                          'stateVarNames': ["StressXX", "StressYY", "StressZZ", \
                                            "StressXY", "StressYZ", "StressXZ", \
                                            "StrainXX", "StrainYY", "StrainZZ", \
                                            "StrainXY", "StrainYZ", "StrainXZ"],
                          'numState0Vars': 6,
                          'state0VarNames': ["Stress0XX", "Stress0YY", "Stress0ZZ", \
                                            "Stress0XY", "Stress0YZ", "Stress0XZ"],
                          'fptrMatPrt': pl3d.mat_prt_1,
                          'fptrElasMat': pl3d.elas_mat_1,
                          'fptrElasStrs': pl3d.elas_strs_1,
                          'fptrTdMatinit': pl3d.td_matinit_1,
                          'fptrTdStrs': pl3d.td_strs_1,
                          'fptrTdStrsMat': pl3d.td_strs_mat_1,
                          'fptrTdPrestrMat': pl3d.prestr_mat_1,
                          'fptrGetState': pl3d.get_state_1,
                          'fptrUpdateState': pl3d.update_state_1}

# version
__id__ = "$Id$"

# End of file 
