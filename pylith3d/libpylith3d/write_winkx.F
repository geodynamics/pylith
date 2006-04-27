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
      subroutine write_winkx(winkxdef,iwinkxdef,iwinkxid,nwinkxe,kw,
     & idout,ofile,ierr,errstrng)
c
c..  program for printing data on winkler restoring forces
c
c          winkdef(ndof,nwinke) = values of winkler restoring spring
c                             constant, force(i,j)=-wink(i,j)*disp(i,j)
c
c          iwinkdef(ndof,nwinke) = application mode:
c                               iwink = 0, no winkler forces,
c                               iwink = 1, applied throuthout computation
c                               iwink = -n, uses load history factor n
c          iwinkid(nwinke)     = node number for application of winkler
c                                force
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
c
c...  subroutine arguments
c
      integer nwinkxe,kw,idout,ierr
      integer iwinkxdef(ndof,nwinkxe),iwinkxid(nwinkxe)
      double precision winkxdef(ndof,nwinkxe)
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
      integer i,j,n,nlines,npage
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open output file
c
      ierr=izero
      if(idout.eq.izero.or.nwinkxe.eq.izero) return
      call open_append(kw,ofile,"old","",
     &                 ierr,errstrng,"write_winkx")
      if(ierr.ne.izero) return
c
c...  loop over winkler forces and output results
c
      nlines=izero
      npage=50
      do i=1,nwinkxe
        n=iwinkxid(i)
        if(mod(nlines,npage).eq.izero) then
          write(kw,1000,err=50)
          write(kw,3000,err=50) (labeld(j),j=1,ndof)
          write(kw,2005,err=50)
        end if
        nlines=nlines+1
        write(kw,2020,err=50) n,(winkxdef(j,i),iwinkxdef(j,i),j=1,ndof)
      end do
      close(kw)
c
c...  normal return
c
      return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_winkx"
        close(kw)
        return
c
 2005 format(/)
 2020 format(1x,i7,2x,3(1pe12.5,2x,i7,3x))
c 1000 format(///,'  w i n k l e r   r e s t o r i n g   f o r c e s',//)
 1000 format(///,'  d i f f e r e n t i a l   w i n k l e r   r e s',
     & ' t o r i n g   f o r c e s',//)
 3000 format(
     x '      key to application mode codes:',//
     x '          mode = 0, no winkler restoring forces; the',/,
     x '                 numerical value to the left is meaningless.',/,
     x '          mode = 1, forces applied at all times',/,
     x '          mode = -n, uses n^th load history factor',///,
     x 1x,'node # ',a4,' winklr     mode',5(2x,a4,' winklr     mode'),/)
      end
c
c version
c $Id: write_winkx.F,v 1.1 2005/08/05 19:58:07 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
