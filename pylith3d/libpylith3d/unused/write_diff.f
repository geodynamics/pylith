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
      subroutine write_diff(diforc,nslip,idhist,numslp,numdif,numnp,
     & kw,idout,ofile,ierr,errstrng)
c
c...  prints differential forces applied to slippery nodes
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file (if numdif.ne.zero)
c         4:  Write error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numslp,numdif,numnp,kw,idout,ierr
      integer nslip(nsdim,numslp),idhist(numnp)
      double precision diforc(ndof,numnp)
      character ofile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer nlines,npage,i,j
      logical nonzed
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open input file
c
      ierr=izero
      if(numslp.eq.izero.or.numdif.eq.izero.or.idout.eq.izero) return
c
c...  loop over differential forces and output results, if desired
c
      nlines=izero
      npage=50
      open(kw,file=ofile,err=40,status="old",access="append")
      do i=1,numnp
        nonzed=.false.
        do j=1,ndof
          if(diforc(j,i).ne.zero) nonzed=.true.
        end do
        if(nonzed) then
          nlines=nlines+1
          if(nlines.eq.ione.or.mod(nlines,npage).eq.izero) then
            write(kw,6000) (labeld(j),j=1,ndof)
            write(kw,*) ' '
          end if
          write(kw,7000,err=50) i,idhist(i),(diforc(j,i),j=1,ndof)
        end if
      end do
      close(kw)
c
c...  normal return
c
      return
c
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="write_diff"
        close(kw)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_diff"
        close(kw)
        return
c
 6000 format(//' differential forces on slippery nodes'//
     & 1x,'  node #   hfac  ',7x,6(a4,11x))
 7000 format(1x,i7,3x,i7,3x,6(3x,1pe12.5))
      end
c
c version
c $Id: write_diff.f,v 1.1 2005/08/05 19:58:07 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
