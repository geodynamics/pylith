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
      subroutine esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,
     & efsts,efstsi,gam,dlam,val,dval,deltp,alfap,iopd,plas)
c
c...program to compute the effective stress function and its
c   derivative for a given effective strain and initial stress state
c
c**** Note that this routine is currently set up for the previous
c     method of determining material behavior (only one material
c     type with the option of including viscous or plastic behavior).
c     In the next version of the code, there should be a different
c     routine like this for every material type (except for purely
c     elastic, which won't need it).  That way, it won't be necessary
c     to perform switches based on whether viscous or plastic behavior
c     is being used.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopd
      double precision ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,efsts,efstsi
      double precision gam,dlam,val,dval,deltp,alfap
      logical plas
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  local variables
c
      double precision strtau,a,dgam,ddlam,ddl1,ddl2
c
c*      write(6,*) "Hello from esf_f!"
      dlam=dl1*efsts+dl2
      if(dlam.le.zero.or.(.not.plas)) dlam=zero
c*      if(dlam.le.zero) dlam=zero
      strtau=(one-alfap)*efstsi+alfap*efsts
      gam=half*(strtau/emhu)**(anpwr-one)/emhu
      a=ae+t2*gam
      if(efsts.ne.zero) a=a+dlam/(t1*efsts)
      val=a*a*efsts*efsts-bs*gam*gam+c*gam-ds
      if(iopd.eq.1) return
      dgam=alfap*(anpwr-one)*gam/strtau
      ddlam=dl1
      if(dlam.eq.zero) ddlam=zero
      ddl1=zero
      ddl2=zero
      if(efsts.ne.zero) then
        ddl1=ddlam/(t1*efsts)
        ddl2=dlam/(t1*efsts*efsts)
      end if
      dval=two*a*a*efsts-two*bs*gam*dgam+c*dgam+
     & two*a*efsts*efsts*(ddl1-ddl2+deltp*dgam*alfap)
      return
      end
c
c version
c $Id: esf.f,v 1.1 2004/06/21 20:59:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
