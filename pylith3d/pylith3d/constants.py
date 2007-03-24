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


# Parameters that are invariant for this geometry type
numberSpaceDimensions = 3
numberDegreesFreedom = 3
stateVariableDimension = 6
materialMatrixDimension = 21
numberSkewDimensions = 2
numberSlipDimensions = 5
numberSlipNeighbors = 4

# self.listIddmat = [
#     1, 2, 3, 4, 5, 6,
#     2, 7, 8, 9,10,11,
#     3, 8,12,13,14,15,
#     4, 9,13,16,17,18,
#     5,10,14,17,19,20,
#     6,11,15,18,20,21]
# Changed this to correspond to BLAS packed symmetric matrix format.
listIddmat = [
     1, 2, 4, 7,11,16,
     2, 3, 5, 8,12,17,
     4, 5, 6, 9,13,18,
     7, 8, 9,10,14,19,
    11,12,13,14,15,20,
    16,17,18,19,20,21]


# Invariant parameters related to element type
maxElementNodes = 20
maxElementNodes2d = 4


# Invariant parameters related to material model
maxMaterialModels = 20
maxStateVariables = 30


numberAllowedVolumeElementTypes = 1


# end of file 
