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
      subroutine read_connect(ien,mat,ienind,nen,neni,ietype,numel,
     & ngauss,netypes,netypesi,nconsize,numelt,numnp,numat,kr,kw,kp,
     & idout,idsk,ierr,ifile,ofile,pfile)
c
c      this subroutine reads element connectivities, material types, and
c      infinite element info
c      The initial element type and the infinite element info are used
c      together to form the final element types, which includes all
c      infinite element permutations.
c
c      Error codes:
c          0:  No error
c          1:  Error opening input file
c          2:  Units not specified (not applicable for this routine)
c          3:  Read error
c          4:  Bad material number
c          5:  Bad connectivity
c
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer netypes,netypesi,nconsize,numelt,numnp,numat,kr,kw,kp
      integer idout,idsk,ierr
      integer ien(nconsize),mat(numelt),ienind(numelt),nen(netypes)
      integer neni(netypesi),ietype(numelt),numel(netypes)
      integer ngauss(netypes)
      character ifile*(*),ofile*(*),pfile*(*)
c
c...  defined constants
c
      include "nconsts.inc"
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
      integer npage,n,i,j,inf,ietypei,i1,i2
c
c...  read connectivity, material number, and infinite element info
c
      call ifill(numel,izero,netypes)
      ierr=izero
      ienind(1)=ione
      open(kr,file=ifile,status="old",err=20)
      call pskip(kr)
      do i=1,numelt
        if(i.ne.1) ienind(i)=ienind(i-1)+neni(ietypei)-1
        read(kr,*,end=30,err=30) n,ietypei,mat(i),inf,
     &   (ien(j),j=ienind(i),ienind(i)+neni(ietypei)-1)
        call infcmp(ietypei,ietype(i),inf)
        numel(ietype(i))=numel(ietype(i))+1
      end do
      close(kr)
c
c...  output plot info, if desired
c
      if(idsk.eq.0) then
        open(kp,file=pfile,status="old",access="append")
        write(kp,1000) netypes,numelt
        write(kp,1000) (numel(j),j=1,netypes)
        write(kp,1000) (ngauss(j),j=1,netypes)
        do i=1,numelt
          i1=ienind(i)
          i2=ienind(i)+nen(ietype(i))-1
          write(kp,1000) i,mat(i),ietype(i),i1,(ien(j),j=i1,i2)
        end do
        close(kp)
      else if(idsk.eq.1) then
        open(kp,file=pfile,status="old",access="append",
     &   form="unformatted")
        write(kp),netypes,numelt
        write(kp) numel
        write(kp) ngauss
        write(kp) mat
        write(kp) ietype
        write(kp) ienind
        write(kp) ien
        close(kp)
      end if
c
c...  output ascii info, if desired
c
      if(idout.gt.0) then
        open(kw,file=ofile,status="old",access="append")
        npage=50
        do n=1,numel
          if(n.eq.1.or.mod(n,npage).eq.0) then
            write(kw,2000)(head(i),i,i=1,nen(ietype(n)))
            write(kw,2500)
          end if
          i1=ienind(n)
          i2=ienind(n)+nen(ietype(n))-1
          write(kw,3000) n,mat(n),ietype(n),(ien(i),i=i1,i2)
        end do
        write(kw,4000)
        close(kw)
      end if
c
c......test for zero or out-of-bound entries in ien array and
c      zero or out-of bounds material numbers
c
      do n=1,numel
        if((mat(n).le.0).or.(mat(n).gt.numat)) then
          ierr=4
          return
        end if
        do i=ienind(n),ienind(n)+nen(ietype(n))-1
          if((ien(i).le.0).or.(ien(i).gt.numnp)) then
            ierr=5
            return
          end if
        end do
      end do
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        close(kr)
        return
c
1000  format(30i7)
2000  format(1x,///,
     x' e l e m e n t  d a t a ',//,5x,
     x' element    material  infinite',20(a4,i2,4x))
2500  format(5x,'  number    element',/,
     x       5x,'              code',/)
3000  format(6x,i5,10(5x,i7))
4000  format(//)
6000  format('ien array generation failure in element # ',i7,
     & ' at node # ',i7)
7000  format('material number of element # ',i7,
     & ' is ',i5, ', which is not defined')
c
      end
c
c version
c $Id: read_connect.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
