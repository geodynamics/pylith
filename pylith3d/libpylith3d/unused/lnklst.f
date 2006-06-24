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
      subroutine lnklst(neq,lm,lmx,infiel,nconsz,numelt,infetype,indx,
     & link,nbrs,iwork,nsizea,nnz,numsn,ierr,errstrng)
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
      integer neq,nconsz,numelt,iwork,nsizea,nnz,numsn,ierr
      integer lm(ndof*nconsz),lmx(ndof*nconsz),infiel(7,numelt)
      integer infetype(4,netypes),indx(neq),link(iwork),nbrs(iwork)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer inext,iel,indlm,ietype,nee,neeq,i,ii,irow,loc,j,jj,icol
      integer loctmp,ival,locsav
cdebug      integer idb,jdb,idbind
c
cdebug      write(6,*) "Hello from lnklst_f!"
c
cdebug      write(6,*) "neq,nconsz,numelt,iwork,nsizea,nnz,numsn:",neq,nconsz,
cdebug     & numelt,iwork,nsizea,nnz,numsn
cdebug      open(15,file="lnklst.info")
cdebug      write(15,*) (lm(idb),idb=1,ndof*nconsz)
cdebug      write(15,*) (lmx(idb),idb=1,ndof*nconsz)
cdebug      do idb=1,numelt
cdebug        write(15,*) (infiel(jdb,idb),jdb=1,6)
cdebug      end do
cdebug      do idb=1,netypes
cdebug        write(15,*) (infetype(jdb,idb),jdb=1,4)
cdebug      end do
cdebug      close(15)
      call ifill(indx,izero,neq)
      call ifill(link,izero,iwork)
      call ifill(nbrs,izero,iwork)
cdebug      do idb=1,netypes
cdebug        write(6,*) "infetype:",(infetype(jdb,idb),jdb=1,4)
cdebug      end do
cdebug      write(6,*) "lm:",(lm(idb),idb=1,ndof*nconsz)
cdebug      write(6,*) "lmx:",(lm(idb),idb=1,ndof*nconsz)
c
      inext=ione
      do iel=1,numelt
c
c      check that available storage is not exceeded
c
        if((inext+neemax*neemax).gt.iwork) then
          ierr=300
          write(errstrng,700) inext+neemax*neemax-iwork
cdebug          write(6,*) "Uncaught exception?"
          return
        end if
        indlm=ndof*(infiel(1,iel)-ione)+ione
        ietype=infiel(3,iel)
        nee=infetype(4,ietype)
        neeq=nee
        if(numsn.ne.izero) neeq=itwo*nee
cdebug        write(6,*) "iel,indlm,ietype,nee,neeq:",
cdebug     &   iel,indlm,ietype,nee,neeq
c
        do i=1,neeq
          ii=indlm+i-1
          irow=lm(ii)
          if(i.gt.nee) irow=abs(lmx(ii-nee))
          if(irow.ne.izero) then
            loc=indx(irow)
            do j=1,neeq
              jj=indlm+j-1
              icol=lm(jj)
              if(j.gt.nee) icol=abs(lmx(jj-nee))
              if(icol.ne.izero.and.icol.ne.irow) then
cdebug                write(6,*) "iel,i,j,ii,jj,irow,icol,loc,inext:",
cdebug     &           iel,i,j,ii,jj,irow,icol,loc,inext
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
cdebug      write(6,*) "After loop in lnklst!"
cdebug      call flush(6)
      nsizea=inext-ione
      nnz=nsizea+neq+ione
cdebug      write(6,*) "nsizea,nnz,neq:",nsizea,nnz,neq
cdebug      call flush(6)
cdebug      idbind=max(inext-100,1)
cdebug      write(6,*) "end of nbrs:",(nbrs(idb),idb=idbind,inext)
cdebug      write(6,*) "end of link:",(link(idb),idb=idbind,inext)
cdebug      idbind=max(neq-100,1)
cdebug      write(6,*) "end of indx:",(indx(idb),idb=idbind,neq)
cdebug      call flush(6)
cdebug      open(15,file="makemsr.info")
cdebug      write(15,*) neq,nnz,iwork
cdebug      write(15,*) (indx(idb),idb=1,neq)
cdebug      write(15,*) (link(idb),idb=1,iwork)
cdebug      write(15,*) (nbrs(idb),idb=1,iwork)
cdebug      close(15)
 700  format("lnklst:  Working storage exceeded by ",i7)
      return
      end
c
c version
c $Id: lnklst.f,v 1.1 2005/03/25 22:52:39 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
