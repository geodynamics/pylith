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
      subroutine write_ucdmesh(
     & x,numnp,                                                         ! global
     & ien,infiel,numelt,nconsz,                                        ! elemnt
     & gauss,sh,infetype,                                               ! eltype
     & istatout,                                                        ! ioopts
     & idout,idsk,kw,kp,fileroot)                                       ! ioinfo
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
      integer numnp,numelt,nconsz,idout,idsk,kw,kp
      integer ien(nconsz),infiel(6,numelt),infetype(4,netypes)
      integer istatout(2,nstatesmax)
      character fileroot*(*)
      double precision x(nsd,numnp),gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
c
c...  local constants
c
      character eltype(10)*4
      data eltype/"hex","wrk","wdg","pyr","tet",
     & "qhex","qwrk","qwdg","qpyr","qtet"/
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
      i1=nnblnk(fileroot)
      i2=nchar(fileroot)
c
c...  write mesh info
c
      filenm=fileroot(i1:i2)//"mesh.inp"
      open(kw,file=filenm,status="new")
      write(kw,"(5i7)") numnp,numelt,nnattr,neattr,nmattr
c
c...  write nodal coordinates
c
      do i=1,numnp
        write(kw,"(i7,3(2x,1pe15.8))") i,(x(j,i),j=1,nsd)
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
        write(kw,"(2i7,2x,a4,20i7)") i,imat,eltype(indtype),
     &   (ien(j),j=indien,indien+nen-1)
      end do
      close(kw)
c
c...  write Gauss point info
c
      filenm=fileroot(i1:i2)//"gmesh.inp"
      open(kw,file=filenm,status="new")
      write(kw,"(5i7)") ngpts,numelg,nngattr,negattr,nmattr
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
          write(kw,"(i7,3(2x,1pe15.8))") igpt,(xg(j),j=1,nsd)
        end do
      end do
      close(kw)
      return
      end
c
c version
c $Id: write_ucdmesh.f,v 1.1 2004/08/23 18:48:15 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
