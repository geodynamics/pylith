c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine matmod_def(infmatmod)
c
c...  quick and dirty method of defining global parameters for all
c     material models.
c
c     At present, the maximum number of material models is 20, but this
c     can easily be changed by editing materials.inc.  This routine
c     defines 6 parameters for each material model:
c
c         infmatmod(1,matmod) = definition flag
c                               0 = material is not defined
c                               1 = material is defined
c         infmatmod(2,matmod) = number of state variables.  If a
c                               material is defined, it will have at
c                               least 12 state variables (stress and
c                               strain).  At present the maximum number
c                               of state variables is 24, if both
c                               viscous and plastic strain are required.
c         infmatmod(3,matmod) = number of material properties needed to
c                               define the material.  The minimum value
c                               should be 3 for an isotropic linear
c                               elastic material (density is always the
c                               first property).
c         infmatmod(4,matmod) = material property variation flag.
c                               This flag indicates whether the material
c                               matrix varies as a function of the state
c                               variables.  If not, the same material
c                               matrix may be used for all elements in a
c                               material group.
c                               0 = material matrix independent of state
c                                   variables.
c                               1 = material matrix depends on state
c                                   variables.
c         infmatmod(5,matmod) = material property behavior flag.
c                               0 = time response of material is linear.
c                               1 = time response of material is
c                                   nonlinear.
c         infmatmod(6,matmod) = number of initial state variables.  At
c                               present, this should always be 6, since
c                               initial state variables are always
c                               assumed to be stresses, but this may
c                               change in the future.
c                 
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer infmatmod(6,nmatmodmax)
c
c...  local variables
c
      integer i,j
c
cdebug      write(6,*) "Hello from matmod_def_f!"
c
      call ifill(infmatmod,izero,6*nmatmodmax)
c
c...  Definitions for isotropic elastic material
c
      infmatmod(1,1) =ione
      infmatmod(2,1) =12
      infmatmod(3,1) =ithree
      infmatmod(6,1) =6
c
c...  Dummy definitions for additional elastic material models
c
      infmatmod(2,2)=12
      infmatmod(2,3)=12
      infmatmod(2,4)=12
      infmatmod(6,2) =6
      infmatmod(6,3) =6
      infmatmod(6,4) =6
c
c...  Definition for linear Maxwell viscoelastic material
c
      infmatmod(1,5) = ione
      infmatmod(2,5) = 18
      infmatmod(3,5) = ifour
      infmatmod(4,5) = ione
      infmatmod(5,5) = izero
      infmatmod(6,5) = 6
c
c...  Definition for power-law Maxwell viscoelastic material
c
      infmatmod(1,6) = ione
      infmatmod(2,6) = 18
      infmatmod(3,6) = ifive
      infmatmod(4,6) = ione
      infmatmod(5,6) = ione
      infmatmod(6,6) = 6
c
c...  Definition for generalized linear Maxwell viscoelastic material
c
      infmatmod(1,7) = ione
      infmatmod(2,7) = 30
      infmatmod(3,7) = 9
      infmatmod(4,7) = ione
      infmatmod(5,7) = izero
      infmatmod(6,7) = 6
c
c...  Dummy definitions for remaining materials
c
      do i=8,nmatmodmax
        do j=2,6
          infmatmod(j,i)=infmatmod(j,5)
        end do
      end do
c
      return
      end
c
c version
c $Id: matmod_def.f,v 1.5 2005/03/28 23:52:18 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
