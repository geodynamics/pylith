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
      subroutine read_connect(neni,infetype,infmat,infmatmod,ien,infiel,
     & indmat,imgrp,nconsz,numelt,numnp,numat,nstatesz,ndmatsz,kr,kw,kp,
     & idout,idsk,ifile,ofile,pfile,ierr,errstrng)
c
c      this subroutine reads element types and connectivities, material
c      types, and infinite element info.
c      The initial element type and the infinite element info are used
c      together to form the final element types, which includes all
c      infinite element permutations.
c
c      Error codes:
c          0:  No error
c          1:  Error opening input file
c          3:  Read error
c        106:  Illegal element type
c        101:  Attempt to use undefined material model
c        107:  Illegal material type
c        108:  Undefined node used for connectivity
c
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nconsz,numelt,numnp,numat,nstatesz,ndmatsz,kr,kw,kp
      integer idout,idsk,ierr
      integer neni(netypesi),infetype(4,netypes),infmat(3,numat)
      integer infmatmod(5,nmatmodmax),ien(nconsz),infiel(6,numelt)
      integer indmat(numat),imgrp(numat)
      character ifile*(*),ofile*(*),pfile*(*),errstrng*(*)
c
c...  local variables
c
      character*4 head(20)
      data head/20*'node'/
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer imat,npage,matmod,i,n,j,inf,ietypei,ngauss,nstate,imatvar
      integer ietype,nen,i1,i2
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from read_connect_f!"
cdebug      write(6,*) nconsz,numelt,numnp,numat,nstatesz,ndmatsz,kr,kw,kp
cdebug      write(6,*) (neni(idb),idb=1,netypesi)
cdebug      write(6,*) "infetype:"
cdebug      do idb=1,netypes
cdebug        write(6,*) (infetype(jdb,idb),jdb=1,4)
cdebug      end do
cdebug      write(6,*) "infmat:"
cdebug      do idb=1,numat
cdebug        write(6,*) (infmat(jdb,idb),jdb=1,3)
cdebug      end do
cdebug      write(6,*) "infmatmod:"
cdebug      do idb=1,nmatmodmax
cdebug        write(6,*) (infmatmod(jdb,idb),jdb=1,5)
cdebug      end do
      call ifill(ien,izero,nconsz)
      call ifill(indmat,izero,numat)
      call ifill(imgrp,izero,numat)
      call ifill(infiel,izero,isix*numelt)
      nstatesz=izero
      ndmatsz=izero
      ierr=izero
c
c...  set initial pointers for different material types
c
      indmat(1)=ione
      do imat=2,numat
        indmat(imat)=indmat(imat-1)+infmat(2,imat-1)
      end do
cdebug      write(6,*) "indmat:",(indmat(idb),idb=1,numat)
c
c...  read connectivity, material number, and infinite element info
c
      infiel(1,1)=ione
      open(kr,file=ifile,status="old",err=20)
      call pskip(kr)
      do i=1,numelt
        if(i.ne.ione) infiel(1,i)=infiel(1,i-1)+neni(ietypei)
        read(kr,*,end=30,err=30) n,ietypei,imat,inf,
     &   (ien(j),j=infiel(1,i),infiel(1,i)+neni(ietypei)-1)
        infiel(2,i)=imat
c
c...  set global element type, after first checking for illegal
c     element type
c
        if(ietypei.le.izero.or.ietypei.gt.netypesi) then
          ierr=106
          errstrng="read_connect"
          return
        end if
        call infcmp(ietypei,infiel(3,i),inf)
c
c...  index element into its material type group, after first checking
c     for illegal material type or material model
c
        matmod=infmat(1,imat)
        if(imat.le.izero.or.imat.gt.numat) then
          ierr=107
          errstrng="read_connect"
          return
        end if
        if(matmod.gt.nmatmodmax.or.infmatmod(1,matmod).eq.izero) then
          ierr=101
          errstrng="read_connect"
          return
        end if
cdebug        write(6,*) "i,imat,indmat(imat):",i,imat,indmat(imat)
        infiel(4,indmat(imat))=i
        indmat(imat)=indmat(imat)+ione
c
c...  set state variable and material matrix indices, and increment the
c     total size of the state variable and material matrices.
c
        infiel(5,i)=nstatesz+ione
        ngauss=infetype(1,infiel(3,i))
        nstate=infmatmod(2,matmod)
        imatvar=infmatmod(4,matmod)
        nstatesz=nstatesz+nstate*ngauss
        if(imatvar.eq.izero) then
          if(imgrp(imat).eq.izero) then
            infiel(6,i)=ndmatsz+ione
            ndmatsz=ndmatsz+ngaussmax
            imgrp(imat)=infiel(6,i)
          else
            infiel(6,i)=imgrp(imat)
cdebug            infiel(6,i)=infiel(6,infiel(4,i-1))
          end if
        else
          infiel(6,i)=ndmatsz+ione
          ndmatsz=ndmatsz+ngauss
        end if
cdebug        write(6,*) i,(infiel(jdb,i),jdb=1,6)
      end do
      close(kr)
cdebug      write(6,*) "infiel:"
cdebug      do idb=1,numelt
cdebug        write(6,*) idb,(infiel(jdb,idb),jdb=1,6)
cdebug      end do
c
c...  output plot info, if desired
c
      if(idsk.eq.izero) then
        open(kp,file=pfile,err=40,status="old",access="append")
        write(kp,1000,err=50) numelt
        do i=1,numelt
          imat=infiel(2,i)
          ietype=infiel(3,i)
          ngauss=infetype(1,ietype)
          nen=infetype(2,ietype)
          i1=infiel(1,i)
          i2=infiel(1,i)+nen-1
          write(kp,1000) i,imat,ietype,ngauss,(ien(j),j=i1,i2)
        end do
        close(kp)
      else if(idsk.eq.ione) then
        open(kp,file=pfile,err=40,status="old",access="append",
     &   form="unformatted")
        write(kp,err=50) numelt,nconsz
        write(kp,err=50) infetype
        write(kp,err=50) infiel
        write(kp,err=50) ien
        close(kp)
      end if
c
c...  output ascii info, if desired
c
      if(idout.gt.izero) then
        open(kw,file=ofile,status="old",access="append")
        npage=50
        do i=1,numelt
          imat=infiel(2,i)
          ietype=infiel(3,i)
          ngauss=infetype(1,ietype)
          nen=infetype(2,ietype)
          i1=infiel(1,i)
          i2=infiel(1,i)+nen-1
          if(i.eq.ione.or.mod(n,npage).eq.izero) then
            write(kw,2000)(head(j),j,j=1,nen)
            write(kw,2500)
          end if
          write(kw,3000) i,imat,ietype,ngauss,(ien(j),j=i1,i2)
        end do
        write(kw,4000)
        close(kw)
      end if
c
c......test for zero or out-of-bound entries in ien array
c
      do i=1,nconsz
        if(ien(i).le.izero.or.ien(i).gt.numnp) then
          ierr=108
          errstrng="read_connect"
          return
        end if
      end do
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_connect"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_connect"
        close(kr)
        return
c
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="read_connect"
        close(kw)
        close(kp)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="read_connect"
        close(kw)
        close(kp)
        return
c
1000  format(30i7)
2000  format(1x,///,
     x          ' e l e m e n t  d a t a ',//,5x,
     x          ' element    material  element  number  ',20(a4,i2,4x))
2500  format(5x,'             number     type   gauss  ',/,
     x       5x,'                               points  ',/)
3000  format(6x,i7,25(5x,i7))
4000  format(//)
c
      end
c
c version
c $Id: read_connect.f,v 1.4 2004/08/12 02:24:21 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
