c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
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
     & ien,infiel,numelt,nconsz,                                        ! elemnt
     & sh,infetype,                                                     ! eltype
     & istatout,                                                        ! ioopts
     & kucd,ucdroot)                                                    ! ioinfo
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
      integer numnp,numelt,nconsz,kucd
      integer ien(nconsz),infiel(6,numelt),infetype(4,netypes)
      integer istatout(2,nstatesmax)
      character ucdroot*(*)
      double precision x(nsd,numnp),sh(nsd+1,nenmax,ngaussmax,netypes)
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
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer nnattr,neattr,nngattr,negattr,nmattr,i,j,indien,imat
      integer ietype,ngauss,nen,numelg,l,ngpts,indtype,i1,i2,igpt
      integer itmp(20)
      double precision xl(nsd,nenmax),xg(nsd)
      character filenm*200
c
cdebug      write(6,*) "Hello from write_ucdmesh_f!"
c
      nnattr=ndof
      neattr=izero
      nngattr=izero
      negattr=izero
      nmattr=izero
      numelg=izero
      ngpts=izero
      do i=1,nstatesmax
        nngattr=nngattr+nstr*istatout(1,i)
        nngattr=nngattr+nstr*istatout(2,i)
      end do
      i1=nnblnk(ucdroot)
      i2=nchar(ucdroot)
c
c...  write mesh info
c
      filenm=ucdroot(i1:i2)//".mesh.inp"
      open(kucd,file=filenm,status="new")
      write(kucd,"(5i7)") numnp,numelt,nnattr,neattr,nmattr
c
c...  write nodal coordinates
c
      do i=1,numnp
        write(kucd,"(i7,3(2x,1pe15.8))") i,(x(j,i),j=1,nsd)
      end do
c
c...  write element connectivities
c
      do i=1,numelt
        indien=infiel(1,i)
        imat=infiel(2,i)
        ietype=infiel(3,i)
        nen=infetype(2,ietype)
        ngauss=infetype(1,ietype)
        ngpts=ngpts+ngauss
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
        do j=1,nen
          itmp(j)=ien(indien+j-1)
        end do
        write(kucd,"(2i7,2x,a4,20i7)") i,imat,eltype(indtype),
     &   (itmp(inducd(j,indtype)),j=1,nen)
      end do
      close(kucd)
c
c...  write Gauss point info
c
      filenm=ucdroot(i1:i2)//".gmesh.inp"
      open(kucd,file=filenm,status="new")
      write(kucd,"(5i7)") ngpts,numelg,nngattr,negattr,nmattr
c
c...  write Gauss coordinates
c
      igpt=izero
      do i=1,numelt
        indien=infiel(1,i)
        ietype=infiel(3,i)
        nen=infetype(2,ietype)
        ngauss=infetype(1,ietype)
        call lcoord(x,xl,ien(indien),nen,numnp)
        do l=1,ngauss
          igpt=igpt+ione
          xg(1)=zero
          xg(2)=zero
          xg(3)=zero
          do j=1,nen
            xg(1)=xg(1)+xl(1,j)*sh(4,j,l,ietype)
            xg(2)=xg(2)+xl(2,j)*sh(4,j,l,ietype)
            xg(3)=xg(3)+xl(3,j)*sh(4,j,l,ietype)
          end do
          write(kucd,"(i7,3(2x,1pe15.8))") igpt,(xg(j),j=1,nsd)
        end do
      end do
      close(kucd)
      return
      end
c
c version
c $Id: write_ucd_mesh.f,v 1.4 2004/08/31 19:14:54 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
