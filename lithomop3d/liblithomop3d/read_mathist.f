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
      subroutine read_mathist(mhist,infmat,numat,npropsz,nhist,
     & kr,kw,kp,idout,idsk,mhfile,ofile,pfile,ierr,errstrng)
c
c...     reads material history information.
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file
c         3:  Read error
c         4:  Write error
c       100:  Attempt to use undefined history
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numat,npropsz,nhist,kr,kw,kp,idout,idsk,ierr
      integer mhist(npropsz),infmat(3,numat)
      character mhfile*(*),ofile*(*),pfile*(*),errstrng*(*)
c
c...  local variables
c
      integer imat,iprop,imhist,indprop,indpropg,i
c
      ierr=izero
      call ifill(mhist,izero,npropsz)
c
c...  open files
c
      open(kr,file=mhfile,status="old",err=10)
      if(idout.gt.izero) open(kw,file=ofile,err=40,status="old",
     & access="append")
      if(idsk.eq.ione) open(kp,file=pfile,err=40,status="old",
     & access="append")
      if(idsk.eq.itwo) open(kp,file=pfile,err=40,status="old",
     & access="append",form="unformatted")
      if(idout.gt.izero) write(kw,800)
c
c...  read material histories, if present
c
 70   continue
        call pskip(kr)
        read(kr,*,end=60,err=30) imat,iprop,imhist
        if(imhist.le.izero.or.imhist.gt.nhist) then
          ierr=100
          errstrng="read_mathist"
          return
        end if
        indprop=infmat(3,imat)
        indpropg=indprop+iprop-ione
        mhist(indpropg)=imhist
        if(idout.gt.izero) write(kw,810,err=50) imat,iprop,imhist
        go to 70
 60   continue
      close(kr)
      if(idout.gt.izero) close(kw)
      if(idsk.eq.ione) write(kp,820,err=50) (mhist(i),i=1,npropsz)
      if(idsk.eq.itwo) write(kp,err=50) mhist
      close(kp)
c
c...  normal return
c
 10   continue
      return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_mathist"
        close(kr)
        return
c
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="read_mathist"
        if(idout.gt.izero) close(kw)
        close(kp)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="read_mathist"
        if(idout.gt.izero) close(kw)
        close(kp)
        return
c
800   format(//,
     & ' material property histories',/,
     & ' The specified properties for the following materials follow',/,
     & ' the given time history:',//,5x,
     & '      Material #     Property #     History #',/)
 810  format(6x,i5,8x,i5,8x,i5)
820   format(16i5)
c
      end
c
c version
c $Id: read_mathist.f,v 1.3 2004/08/25 01:12:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
