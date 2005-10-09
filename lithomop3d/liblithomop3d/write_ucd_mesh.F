c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
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
c...  The current setup uses autoconf to determine whether the fortran
c     compiler supports old-style f77 unformatted direct access,
c     stream I/O, or neither.  If neither is supported, binary output
c     is turned off.
c
#if defined FORTRAN_F77_IO
        write(kucd,rec=ibyte) magnum,numnp,numelv,nnattr,neattr,nmattr,
     &   nnlist
#elif defined FORTRAN_STREAM_IO
        write(kucd,pos=ibyte) magnum,numnp,numelv,nnattr,neattr,nmattr,
     &   nnlist
#else
#warning "Binary UCD output not supported for this Fortran compiler!"
#endif
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
#if defined FORTRAN_F77_IO
            write(kucd,rec=ibyte) ielg,matmodel,nen,ictype
#elif defined FORTRAN_STREAM_IO
            write(kucd,pos=ibyte) ielg,matmodel,nen,ictype
#else
#warning "Binary UCD output not supported for this Fortran compiler!"
#endif
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
#if defined FORTRAN_F77_IO
            write(kucd,rec=ibyte) (itmp(inducd(j,indtype)),j=1,nen)
#elif defined FORTRAN_STREAM_IO
            write(kucd,pos=ibyte) (itmp(inducd(j,indtype)),j=1,nen)
#else
#warning "Binary UCD output not supported for this Fortran compiler!"
#endif
            ibyte=ibyte+nen*intlen
          end do
        end do
c
c...  nodal coordinates
c
        do i=1,nsd
#if defined FORTRAN_F77_IO
          write(kucd,rec=ibyte) (real(x(i,j)),j=1,numnp)
#elif defined FORTRAN_STREAM_IO
          write(kucd,pos=ibyte) (real(x(i,j)),j=1,numnp)
#else
#warning "Binary UCD output not supported for this Fortran compiler!"
#endif
          ibyte=ibyte+numnp*floatlen
        end do
        close(kucd)
c
c...  write Gauss point info
c
        ibyte=ione
        iopt=itwo
        call open_ucd(kucd,iprestress,nstep,ucdroot,iopt,iucd)
#if defined FORTRAN_F77_IO
        write(kucd,rec=ibyte) magnum,ngpts,numelg,nngattr,negattr,
     &   nmattr,nglist
#elif defined FORTRAN_STREAM_IO
        write(kucd,pos=ibyte) magnum,ngpts,numelg,nngattr,negattr,
     &   nmattr,nglist
#else
#warning "Binary UCD output not supported for this Fortran compiler!"
#endif
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
#if defined FORTRAN_F77_IO
            write(kucd,rec=ibyte) real(xg(1))
            write(kucd,rec=ibyte+istride) real(xg(2))
            write(kucd,rec=ibyte+2*istride) real(xg(3))
#elif defined FORTRAN_STREAM_IO
            write(kucd,pos=ibyte) real(xg(1))
            write(kucd,pos=ibyte+istride) real(xg(2))
            write(kucd,pos=ibyte+2*istride) real(xg(3))
#else
#warning "Binary UCD output not supported for this Fortran compiler!"
#endif
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
c $Id: write_ucd_mesh.F,v 1.1 2005/08/05 19:58:07 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
