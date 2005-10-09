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
      subroutine write_connect(iens,ivfamily,indxiel,nen,ngauss,numelv,
     & ietypev,nvfamilies,kw,kp,idout,idsk,ofile,pfile,ierr,errstrng)
c
c      this subroutine writes the element node array, as well as
c      material model and element type.
c
c      Error codes:
c          0:  No error
c          2:  Error opening output file
c          4:  Write error
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
      integer nen,ngauss,numelv,ietypev,nvfamilies,kw,kp,idout,idsk,ierr
      integer iens(nen,numelv),ivfamily(5,nvfamilies),indxiel(numelv)
      character ofile*(*),pfile*(*),errstrng*(*)
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
      integer ielg,ifam,nelfamily,matmod,ielf,ielgprev,j,npage
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from write_connect_f!"
c
c...  output plot info, if desired
c
      if(idsk.eq.ione) then
c
c...  The current setup uses autoconf to determine whether the fortran
c     compiler supports access="append", position="append", or neither.
c     If neither is supported an error is returned.
c
#if defined FORTRAN_ACCESS_APPEND
        open(kp,file=pfile,err=40,status="old",access="append")
#elif defined FORTRAN_POSITION_APPEND
        open(kp,file=pfile,err=40,status="old",position="append")
#else
#error "Cannot determine method for appending to Fortran files!"
#endif
        write(kp,1000,err=50) numelv,ietypev,nen,ngauss
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmod=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+ione
            ielgprev=indxiel(ielg)
            write(kp,1000) ielg,ielgprev,matmod,(iens(j,ielg),j=1,nen)
          end do
        end do
        close(kp)
      else if(idsk.eq.itwo) then
#if defined FORTRAN_ACCESS_APPEND
        open(kp,file=pfile,err=40,status="old",access="append",
     &   form="unformatted")
#elif defined FORTRAN_POSITION_APPEND
        open(kp,file=pfile,err=40,status="old",position="append",
     &   form="unformatted")
#else
#error "Cannot determine method for appending to Fortran files!"
#endif
        write(kp,err=50) numelv,ietypev,nen,ngauss
        write(kp,err=50) ivfamily
        write(kp,err=50) iens
        close(kp)
      end if
c
c...  output ascii info, if desired
c
      if(idout.gt.izero) then
#if defined FORTRAN_ACCESS_APPEND
        open(kw,file=ofile,status="old",access="append")
#elif defined FORTRAN_POSITION_APPEND
        open(kw,file=ofile,status="old",position="append")
#else
#error "Cannot determine method for appending to Fortran files!"
#endif
        npage=50
        ielg=izero
        do ifam=1,nvfamilies
          nelfamily=ivfamily(1,ifam)
          matmod=ivfamily(2,ifam)
          do ielf=1,nelfamily
            ielg=ielg+ione
            ielgprev=indxiel(ielg)
            if(ielg.eq.ione.or.mod(ielg,npage).eq.izero) then
              write(kw,2000)(head(j),j,j=1,nen)
              write(kw,2500)
            end if
            write(kw,3000) ielg,ielgprev,matmod,ietypev,ngauss,
     &       (iens(j,ielg),j=1,nen)
          end do
        end do
        write(kw,4000)
        close(kw)
      end if
c
c...  normal return
c
      return
c
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="write_connect"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_connect"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
1000  format(30i7)
2000  format(1x,///,
     x          ' e l e m e n t  d a t a ',//,5x,
     x          ' element    element    material  ',20(a4,i2,4x))
2500  format(5x,' reordered  initial     model    ',/)
3000  format(6x,i7,25(5x,i7))
4000  format(//)
c
      end
c
c version
c $Id: write_connect.F,v 1.1 2005/08/05 19:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
