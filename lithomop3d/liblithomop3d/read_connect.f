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
      subroutine read_connect(infmatmod,ivflist,nen,ngauss,ien,ivfamily,
     & ivftmp,maxvfamilies,nvfamilies,nprestrflag,numelv,numnp,
     & nstatesz,nstatesz0,npropsz,kr,kw,kp,idout,idsk,ifile,ofile,pfile,
     & ierr,errstrng)
c
c      this subroutine reads element types and connectivities, material
c      types, and infinite element info, and sorts things into element
c      families.
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
      integer nen,ngauss,maxvfamilies,nvfamilies,nprestrflag,numelv
      integer numnp,nstatesz,nstatesz0,npropsz,kr,kw,kp,idout,idsk,ierr
      integer infmatmod(6,nmatmodmax),ivflist(3,maxvfamilies)
      integer ivfamily(5,nvfamilies),ivftmp(nvfamilies),ien(nen,numelv)
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
      integer ifam,incstate0,nstate0max,nelfamily,j,matmod,nstate
      integer nstate0,nprop,i,n,ietypev,imat,inf,iloc,npage,ielf,ielg
      integer itmp(nenmax)
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from read_connect_f!"
cdebug      write(6,*) "nen,ngauss,maxvfamilies,nvfamilies,nprestrflag:"
cdebug      write(6,*) nen,ngauss,maxvfamilies,nvfamilies,nprestrflag
cdebug      write(6,*) "numelv,numnp,nstatesz,nstatesz0,npropsz,kr,kw,kp:"
cdebug      write(6,*) numelv,numnp,nstatesz,nstatesz0,npropsz,kr,kw,kp
cdebug      write(6,*) "idout,idsk,ierr:"
cdebug      write(6,*) idout,idsk,ierr
      call ifill(ien,izero,nen*numelv)
      call ifill(ivfamily,izero,5*nvfamilies)
      ifam=izero
      nstatesz=izero
      nstatesz0=izero
      npropsz=izero
      incstate0=izero
      nstate0max=izero
      nelfamily=izero
      if(nprestrflag.ne.izero) incstate0=ione
      ierr=izero
c
c...  set up element family information array (ivfamily)
c
      do j=1,maxvfamilies
cdebug        write(6,*) "j,ivflist:",j,ivflist(1,j),ivflist(2,j)
        if(ivflist(1,j).ne.izero) then
          ifam=ifam+ione
          if(ifam.eq.ione) then
            ivftmp(ifam)=ione
          else
            ivftmp(ifam)=ivftmp(ifam-ione)+nelfamily
          end if
          nelfamily=ivflist(1,j)
          matmod=ivflist(3,j)
          ivfamily(1,ifam)=nelfamily
          ivfamily(2,ifam)=matmod
          nstate=infmatmod(2,matmod)
          nstate0=infmatmod(6,matmod)
          nstate0max=max(nstate0max,nstate0)
          nprop=infmatmod(3,matmod)
          ivfamily(3,ifam)=nstatesz+ione
          ivfamily(4,ifam)=nstatesz0+ione
          ivfamily(5,ifam)=npropsz+ione
          nstatesz=nstatesz+nstate*ngauss*nelfamily
          nstatesz0=nstatesz0+incstate0*nstate0*ngauss*nelfamily
          npropsz=npropsz+nprop
cdebug          write(6,*) "ifam,ivftmp(ifam),nelfamily,matmod:"
cdebug          write(6,*) ifam,ivftmp(ifam),nelfamily,matmod
        end if
      end do
c
c...  set size of state0 array if there are no prestresses
c
      if(nprestrflag.eq.izero) nstatesz0=ngauss*nstate0max
c
c...  read connectivity, material number, and infinite element info
c
      open(kr,file=ifile,status="old",err=20)
      call pskip(kr)
      do i=1,numelv
        read(kr,*,end=30,err=30) n,ietypev,imat,inf,(itmp(j),j=1,nen)
        ifam=imat
        matmod=ivfamily(2,ifam)
cdebug        write(6,*) "i,n,ietypev,imat,inf,ifam,matmod:"
cdebug        write(6,*) i,n,ietypev,imat,inf,ifam,matmod
c
c...  check for illegal element type
c
        if(ietypev.le.izero.or.ietypev.gt.netypesi) then
          ierr=106
          errstrng="read_connect"
          return
        end if
clater        call infcmp(ietypev,infiel(3,i),inf)
c
c...  check for illegal material type or material model
c
        if(imat.le.izero.or.imat.gt.nvfamilies) then
          ierr=107
          errstrng="read_connect"
          return
        end if
        if(matmod.gt.nmatmodmax.or.infmatmod(1,matmod).eq.izero) then
          ierr=101
          errstrng="read_connect"
          return
        end if
c
c...  store element nodes in the proper family and update the pointer
c     array for that family
c
        iloc=ivftmp(imat)
        call iquate(ien(1,iloc),itmp,nen)
        ivftmp(imat)=ivftmp(imat)+ione
      end do
      close(kr)
c
c...  output plot info, if desired
c
      if(idsk.eq.ione) then
        open(kp,file=pfile,err=40,status="old",access="append")
        write(kp,1000,err=50) numelv
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmod=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+ione
            write(kp,1000) ielg,matmod,ietypev,ngauss,
     &       (ien(j,ielg),j=1,nen)
          end do
        end do
        close(kp)
      else if(idsk.eq.itwo) then
        open(kp,file=pfile,err=40,status="old",access="append",
     &   form="unformatted")
        write(kp,err=50) numelv,nen,ngauss
        write(kp,err=50) ivfamily
        write(kp,err=50) ien
        close(kp)
      end if
c
c...  output ascii info, if desired
c
      if(idout.gt.izero) then
        open(kw,file=ofile,status="old",access="append")
        npage=50
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmod=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+ione
            if(ielg.eq.ione.or.mod(ielg,npage).eq.izero) then
              write(kw,2000)(head(j),j,j=1,nen)
              write(kw,2500)
            end if
            write(kw,3000) ielg,matmod,ietypev,ngauss,
     &       (ien(j,ielg),j=1,nen)
          end do
        end do
        write(kw,4000)
        close(kw)
      end if
c
c......test for zero or out-of-bound entries in ien array
c
      do n=1,numelv
        do j=1,nen
          if(ien(j,n).le.izero.or.ien(j,n).gt.numnp) then
            ierr=108
            errstrng="read_connect"
            return
          end if
        end do
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
c $Id: read_connect.f,v 1.8 2005/04/01 23:24:41 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
