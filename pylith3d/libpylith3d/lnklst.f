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
      subroutine lnklst(neq,lm,lmx,numelv,nen,nee,indx,link,nbrs,iwork,
     & nsizea,nnz,numsn,ierr,errstrng)
c
c      program to create a linked list of nonzero row and column entries
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer neq,numelv,nen,nee,iwork,nsizea,nnz,numsn,ierr
      integer lm(nee,numelv),lmx(nee,numelv)
      integer indx(neq),link(iwork),nbrs(iwork)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer inext,ielg,neeq,i,irow,loc,j,icol
      integer loctmp,ival,locsav
c
cdebug      write(6,*) "Hello from lnklst_f!"
c
      call ifill(indx,izero,neq)
      call ifill(link,izero,iwork)
      call ifill(nbrs,izero,iwork)
c
      neeq=nee
      if(numsn.ne.izero) neeq=2*nee
      inext=ione
      do ielg=1,numelv
c
c      check that available storage is not exceeded
c
        if((inext+nee*nee).gt.iwork) then
          ierr=300
          write(errstrng,700) inext+nee*nee-iwork
          return
        end if
c
        do i=1,neeq
          if(i.le.nee) then
            irow=lm(i,ielg)
          else
            irow=abs(lmx(i-nee,ielg))
          end if
          if(irow.ne.izero) then
            loc=indx(irow)
            do j=1,neeq
              if(j.le.nee) then
                icol=lm(j,ielg)
              else
                icol=abs(lmx(j-nee,ielg))
              end if
              if(icol.ne.izero.and.icol.ne.irow) then
                if(loc.eq.izero) then
                  loc=inext
                  indx(irow)=inext
                  nbrs(inext)=icol
                  link(inext)=-irow
                  inext=inext+1
                else
                  loctmp=loc
 40               continue
                    ival=nbrs(loctmp)
                    if(icol.eq.ival) goto 30
                    locsav=loctmp
                    loctmp=link(loctmp)
                  if(loctmp.gt.izero) goto 40
                  link(locsav)=inext
                  nbrs(inext)=icol
                  link(inext)=-irow
                  inext=inext+1
                end if
              end if
 30           continue
            end do
          end if
        end do
      end do
      nsizea=inext-ione
      nnz=nsizea+neq+ione
 700  format("lnklst:  Working storage exceeded by ",i7)
      return
      end
c
c version
c $Id: lnklst.f,v 1.6 2005/03/25 23:04:26 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
