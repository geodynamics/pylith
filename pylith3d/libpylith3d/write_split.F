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
      subroutine write_split(fault,nfault,numfn,kw,kp,
     & idout,idsk,ofile,pfile,ierr,errstrng)
c
c...  prints data on split nodes
c
c      fault(ndof,numfn)  = amount of splitting for each dof
c      nfault(3,numfn)    = element, node, and history for each entry
c      numfn              = number of split node entries
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
      integer numfn,kw,kp,idout,idsk,ierr
      integer nfault(3,numfn)
      double precision fault(ndof,numfn)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer i,j,npage
c
c...  open output files
c
      if(idout.gt.izero) then
        call open_append(kw,ofile,"old","",ierr,errstrng,"write_split")
        if(ierr.ne.izero) return
      end if
      if(idsk.eq.ione) then
        call open_append(kp,pfile,"old","",ierr,errstrng,"write_split")
        if(ierr.ne.izero) return
      else if(idsk.eq.itwo) then
        call open_append(kp,pfile,"old","unformatted",
     &                   ierr,errstrng,"write_split")
        if(ierr.ne.izero) return
      end if
c
c...  output results if requested
c
      if(idsk.eq.ione) write(kp,"(i7)") numfn
      if(idsk.eq.itwo) write(kp) numfn
      if(idout.gt.izero.or.idsk.eq.ione) then
        npage=50
        do i=1,numfn
          if((i.eq.ione.or.mod(i,npage).eq.izero).and.idout.gt.izero)
     &     write(kw,3000,err=60)
          if(idout.gt.izero) write(kw,4000,err=60) (nfault(j,i),j=1,3),
     &     (fault(j,i),j=1,ndof)
          if(idsk.eq.ione) write(kp,5000,err=60) (nfault(j,i),j=1,2)
        end do
      end if
      if(idsk.eq.itwo.and.numfn.ne.izero) write(kp,err=60) nfault
      if(idsk.eq.itwo.and.numfn.ne.izero) write(kp,err=60) fault
      if(idout.gt.izero) close(kw)
      if(idsk.gt.izero) close(kp)
c
c...  normal return
c
      return
c
c...  error writing to output file
c
 60   continue
        ierr=4
        errstrng="write_split"
        if(idsk.gt.izero) close(kp)
        if(idout.gt.izero) close(kw)
        return
c
 3000 format(1x///1x,'s p l i t  n o d e   d a t a'//
     &'     ihist =  0, slip in elastic computation'/
     &'           = -1, slip at constant velocity '/
     &'           =  n, uses load history factor n'//
     &' elem',6x,'node',7x,'ihist',7x,'split',12x,'split',12x,'split'/
     &'  num       num',20x,'dof1',13x,'dof2',13x,'dof3'//)
 4000 format(i7,5x,i7,5x,i7,3(5x,1pe12.5))
 5000 format(2i7)
      end
c
c version
c $Id: write_split.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
