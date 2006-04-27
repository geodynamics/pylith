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
      subroutine write_connect(iens,ivfamily,indxiel,nen,ngauss,numelv,
     & ietypev,nvfamilies,kw,kp,idout,idsk,ofile,pfile,ierr,errstrng)
c
c      this subroutine writes the element node array, as well as
c      material model and element type.
c
c      Error codes:
c          0:  No error
c          2:  Error opening output file
c          4:  Write error
c
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nen,ngauss,numelv,ietypev,nvfamilies,kw,kp,idout,idsk,ierr
      integer iens(nen,numelv),ivfamily(5,nvfamilies),indxiel(numelv)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  local variables
c
      character*4 head(20)
      data head/20*'node'/
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer ielg,ifam,nelfamily,matmod,ielf,ielgprev,j,npage
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from write_connect_f!"
c
c...  output plot info, if desired
c
      if(idsk.eq.ione) then
        call open_append(kp,pfile,"old","",
     &                   ierr,errstrng,"write_connect")
        if(ierr.ne.izero) return
        write(kp,1000,err=50) numelv,ietypev,nen,ngauss
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmod=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+ione
            ielgprev=indxiel(ielg)
            write(kp,1000) ielg,ielgprev,matmod,(iens(j,ielg),j=1,nen)
          end do
        end do
        close(kp)
      else if(idsk.eq.itwo) then
        call open_append(kp,pfile,"old","unformatted",
     &                   ierr,errstrng,"write_connect")
        if(ierr.ne.izero) return
        write(kp,err=50) numelv,ietypev,nen,ngauss
        write(kp,err=50) ivfamily
        write(kp,err=50) iens
        close(kp)
      end if
c
c...  output ascii info, if desired
c
      if(idout.gt.izero) then
        call open_append(kw,ofile,"old","",
     &                   ierr,errstrng,"write_connect")
        if(ierr.ne.izero) return
        npage=50
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmod=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+ione
            ielgprev=indxiel(ielg)
            if(ielg.eq.ione.or.mod(ielg,npage).eq.izero) then
              write(kw,2000)(head(j),j,j=1,nen)
              write(kw,2500)
            end if
            write(kw,3000) ielg,ielgprev,matmod,ietypev,ngauss,
     &       (iens(j,ielg),j=1,nen)
          end do
        end do
        write(kw,4000)
        close(kw)
      end if
c
c...  normal return
c
      return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_connect"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
1000  format(30i7)
2000  format(1x,///,
     x          ' e l e m e n t  d a t a ',//,5x,
     x          ' element    element    material  ',20(a4,i2,4x))
2500  format(5x,' reordered  initial     model    ',/)
3000  format(6x,i7,25(5x,i7))
4000  format(//)
c
      end
c
c version
c $Id: write_connect.F,v 1.1 2005/08/05 19:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
