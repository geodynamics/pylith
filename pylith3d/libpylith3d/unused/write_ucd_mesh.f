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
      subroutine write_ucd_mesh(
     & x,numnp,                                                         ! global
     & ien,ivfamily,numelv,nvfamilies,                                  ! elemnt
     & sh,nen,ngauss,ietype,                                            ! eltype
     & istatout,nstatout,                                               ! ioopts
     & kucd,iucd,ucdroot)                                               ! ioinfo
c
c...  Specialized routine to output mesh info for SCEC benchmarks.
c     This routine writes out the headers, coordinates, and
c     connectivities for the given mesh, and then writes out the
c     header and coordinates of the gauss points to a separate file.
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
      integer numnp,numelv,nvfamilies,nen,ngauss,ietype,kucd,iucd
      integer ien(nen,numelv),ivfamily(5,nvfamilies)
      integer istatout(nstatesmax,3),nstatout(3)
      character ucdroot*(*)
      double precision x(nsd,numnp),sh(nsd+1,nen,ngauss)
c
c...  local constants
c
      integer inducd(nenmax,10)
      data inducd/
     & 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 5, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 4, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,
     & 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18, 0, 0,
     & 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,14,13,15, 0, 0, 0, 0, 0,
     & 5, 1, 2, 3, 4,10,11,12,13, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, 0,
     & 4, 1, 2, 3, 8, 9,10, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      character eltype(10)*5
      data eltype/"hex","wrick","prism","pyr","tet",
     & "hex","wrick","prism","pyr","tet"/
c
c...  intrinsic functions
c
      intrinsic char,real
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer nnattr,neattr,nngattr,negattr,nmattr,nnlist,nglist,i,j
      integer iprestress,nstep,iopt
      integer ielg,ifam,ielf,nelfamily,matmodel
      integer numelg,l,ngpts,indtype,ictype,igpt
      integer ibyte,intlen,floatlen,istride
      integer itmp(20)
      double precision xl(nsd,nenmax),xg(nsd)
      character magnum*1
c
cdebug      write(6,*) "Hello from write_ucdmesh_f!"
c
      magnum=char(7)
      ibyte=ione
      intlen=ifour
      floatlen=ifour
      nnattr=itwo*ndof
      neattr=izero
      nngattr=nstatout(1)+nstatout(2)+nstatout(3)
      negattr=izero
      nmattr=izero
      numelg=izero
      ngpts=ngauss*numelv
      nnlist=numelv*nen
      nglist=izero
      iprestress=izero
      nstep=izero
      iopt=ione
c
c...  determine element type info
c
      indtype=ione
      if(ietype.eq.28) indtype=itwo
      if(ietype.eq.29) indtype=ithree
      if(ietype.eq.30) indtype=ifour
      if(ietype.eq.31) indtype=ifive
      if(ietype.gt.31) indtype=isix
      if(ietype.eq.59) indtype=iseven
      if(ietype.eq.60) indtype=ieight
      if(ietype.eq.61) indtype=inine
      if(ietype.eq.62) indtype=10
      if(indtype.eq.ione.or.indtype.eq.isix) ictype=iseven
      if(indtype.eq.itwo.or.indtype.eq.iseven) ictype=ieight
      if(indtype.eq.ithree.or.indtype.eq.ieight) ictype=isix
      if(indtype.eq.ifour.or.indtype.eq.inine) ictype=ifive
      if(indtype.eq.ifive.or.indtype.eq.10) ictype=ifour
c
c...  open UCD output file
c
      call open_ucd(kucd,iprestress,nstep,ucdroot,iopt,iucd)
c
c...  ascii ucd output
c
      if(iucd.eq.ione) then
c
c...  header info
c
        write(kucd,"(5i7)") numnp,numelv,nnattr,neattr,nmattr
c
c...  write nodal coordinates
c
        do i=1,numnp
          write(kucd,"(i7,3(2x,1pe15.8))") i,(x(j,i),j=1,nsd)
        end do
c
c...  write element connectivities
c
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmodel=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+1
            do j=1,nen
              itmp(j)=ien(j,ielg)
            end do
            write(kucd,"(2i7,2x,a4,20i7)") ielg,matmodel,
     &       eltype(indtype),(itmp(inducd(j,indtype)),j=1,nen)
          end do
        end do
        close(kucd)
c
c...  write Gauss point info
c
        iopt=itwo
        call open_ucd(kucd,iprestress,nstep,ucdroot,iopt,iucd)
        write(kucd,"(5i7)") ngpts,numelg,nngattr,negattr,nmattr
c
c...  write Gauss coordinates
c
        igpt=izero
        do ielg=1,numelv
          call lcoord(x,xl,ien(1,ielg),nen,numnp)
          do l=1,ngauss
            igpt=igpt+ione
            xg(1)=zero
            xg(2)=zero
            xg(3)=zero
            do j=1,nen
              xg(1)=xg(1)+xl(1,j)*sh(4,j,l)
              xg(2)=xg(2)+xl(2,j)*sh(4,j,l)
              xg(3)=xg(3)+xl(3,j)*sh(4,j,l)
            end do
            write(kucd,"(i7,3(2x,1pe15.8))") igpt,(xg(j),j=1,nsd)
          end do
        end do
        close(kucd)
c
c...  binary UCD output
c
      else if(iucd.eq.itwo) then
c
c...  header info
c
        write(kucd,rec=ibyte) magnum,numnp,numelv,nnattr,neattr,nmattr,
     &   nnlist
        ibyte=ibyte+ione+6*intlen
c
c...  element information
c
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmodel=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+1
            write(kucd,rec=ibyte) ielg,matmodel,nen,ictype
            ibyte=ibyte+4*intlen
          end do
        end do
c
c...  element connectivity
c
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmodel=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+1
            do j=1,nen
              itmp(j)=ien(j,ielg)
            end do
            write(kucd,rec=ibyte) (itmp(inducd(j,indtype)),j=1,nen)
            ibyte=ibyte+nen*intlen
          end do
        end do
c
c...  nodal coordinates
c
        do i=1,nsd
          write(kucd,rec=ibyte) (real(x(i,j)),j=1,numnp)
          ibyte=ibyte+numnp*floatlen
        end do
        close(kucd)
c
c...  write Gauss point info
c
        ibyte=ione
        iopt=itwo
        call open_ucd(kucd,iprestress,nstep,ucdroot,iopt,iucd)
        write(kucd,rec=ibyte) magnum,ngpts,numelg,nngattr,negattr,
     &   nmattr,nglist
        ibyte=ibyte+ione+6*intlen
c
c...  write Gauss coordinates
c
        istride=floatlen*ngpts
        do ielg=1,numelv
          call lcoord(x,xl,ien(1,ielg),nen,numnp)
          do l=1,ngauss
            xg(1)=zero
            xg(2)=zero
            xg(3)=zero
            do j=1,nen
              xg(1)=xg(1)+xl(1,j)*sh(4,j,l)
              xg(2)=xg(2)+xl(2,j)*sh(4,j,l)
              xg(3)=xg(3)+xl(3,j)*sh(4,j,l)
            end do
            write(kucd,rec=ibyte) real(xg(1))
            write(kucd,rec=ibyte+istride) real(xg(2))
            write(kucd,rec=ibyte+2*istride) real(xg(3))
            ibyte=ibyte+floatlen
          end do
        end do
        close(kucd)
      end if
c        
      return
      end
c
c version
c $Id: write_ucd_mesh.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
