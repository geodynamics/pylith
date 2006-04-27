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
      subroutine write_bc(bond,ibond,numnp,kw,idout,ofile,ierr,errstrng)
c
c...  subroutine to write boundary conditions to ascii output file.
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file
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
      integer numnp,kw,idout,ierr
      integer ibond(ndof,numnp)
      double precision bond(ndof,numnp)
      character ofile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  local constants
c
      character ltype(4)*4
      data ltype/'free','disp','velo','forc'/
c
c...  intrinsic functions
c
      intrinsic mod,index
c
c...  local variables
c
      integer ihist(3),itype(3)
      integer i,n,nlines,npage,imode
      logical nonzed
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
      if(idout.eq.izero) return
c
c...  output BC to ascii file, if requested
c
      call open_append(kw,ofile,"old","",ierr,errstrng,"write_bc")
      if(ierr.ne.izero) return
      write(kw,1000) (labeld(i),i=1,ndof)
      write(kw,2000)
      nlines=0
      npage=50
      do n=1,numnp
        nonzed=.false.
        do i=1,ndof
          if(ibond(i,n).ne.0) nonzed=.true.
        end do
        if(nonzed) then
          nlines=nlines+1
          if(mod(nlines,npage).eq.0) then
            write(kw,1000) (labeld(i),i=1,ndof)
            write(kw,2000)
          end if
          do i=1,ndof
            imode=ibond(i,n)
            ihist(i)=imode/10
            itype(i)=imode-10*ihist(i)
          end do
          write(kw,3000,err=20) n,(bond(i,n),ltype(itype(i)+1),ihist(i),
     &     i=1,ndof)
        end if
      end do
      close(kw)
c
c...  normal return
c
      return
c
c...  write error
c
20    continue
        ierr=4
        errstrng="write_bc"
        close(kw)
        return
c
1000  format(1x,///,'  n o d a l   f o r c e s   a n d   d i s p l a',
     1 ' c e m e n t s',//,
     2 '      key to boundary condition codes:',//,
     3 '          free = unconstrained degree of freedom; the',/,
     4 '                 numerical value is meaningless.',/,
     5 '          disp = fixed displacement',/,
     6 '          velo = constant velocity',/,
     7 '          forc = applied (constant) force',//,
     8 '      hfac = load history factor applied (0 if none)',///,
     9 ' node # ',3(4x,a4,6x,'type  hfac  '))
2000  format(' ')
3000  format(1x,i7,2x,3(1pe12.5,2x,a4,1x,i3,4x))
      end
c
c version
c $Id: write_bc.F,v 1.1 2005/08/05 19:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
