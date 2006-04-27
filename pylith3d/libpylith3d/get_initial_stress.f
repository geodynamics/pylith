c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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
      subroutine get_initial_stress(
     & state,state0,ivfamily,nvfamilies,numelv,nstatesz,nstatesz0,      ! elemnt
     & infmatmod,                                                       ! materl
     & nen,ngauss)                                                      ! eltype
c
c...  routine to copy computed stresses into initial stress array
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nvfamilies,numelv,nstatesz,nstatesz0,nen,ngauss
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      double precision state(nstatesz),state0(nstatesz0)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer ifam,nelfamily,matmodel,indstate,indstate0,nstate,nstate0
      integer ielf,l
cdebug      integer idb,jdb
c
c...  included variable definitions
c
c
cdebug      write(6,*) "Hello from get_initial_stress_f!"
c
c
c...  loop over number of element families
c
      do ifam=1,nvfamilies
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indstate=ivfamily(3,ifam)
        indstate0=ivfamily(4,ifam)
        nstate=infmatmod(2,matmodel)
        nstate0=infmatmod(6,matmodel)
c
c...  loop over elements in a family and then gauss points for each
c     element
c
        do ielf=1,nelfamily
          do l=1,ngauss
            indstate=indstate+(l-1)*nstate
            indstate0=indstate0+(l-1)*nstate0
            call dcopy(nstr,state(indstate),ione,state0(indstate0),ione)
          end do
        end do
      end do
      return
      end
c
c version
c $Id: get_initial_stress.f,v 1.2 2005/03/21 19:48:57 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
