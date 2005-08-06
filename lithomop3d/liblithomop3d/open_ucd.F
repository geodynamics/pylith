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
      subroutine open_ucd(kucd,iprestress,nstep,ucdroot,iopt,iucd)
c
c...  Routine to determine filename for UCD output and open the file.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer kucd,iprestress,nstep,iopt,iucd
      character ucdroot*(*)
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer i1,i2
      character filenm*200,cstep*5
c
c...  At present, it is assumed that all f77 compilers support
c     unformatted direct access with a record length of 1 (byte).
c     If exceptions are found we will have to use FORTRAN_COMPILER_NAME
c     to differentiate.
c
#if defined FORTRAN_TYPE_F77
      integer len
      len=ione
#endif
c
cdebug      write(6,*) "Hello from open_ucd_f!"
c
      i1=nnblnk(ucdroot)
      i2=nchar(ucdroot)
      if(iopt.eq.1) then
        filenm=ucdroot(i1:i2)//".mesh.inp"
      else if(iopt.eq.2) then
        filenm=ucdroot(i1:i2)//".gmesh.inp"
      else if(iopt.eq.3) then
        if(iprestress.eq.izero) then
          write(cstep,"(i5.5)") nstep
        else
          cstep="prest"
        end if
        filenm=ucdroot(i1:i2)//".mesh.time."//cstep//".inp"
      else if(iopt.eq.4) then
        if(iprestress.eq.izero) then
          write(cstep,"(i5.5)") nstep
        else
          cstep="prest"
        end if
        filenm=ucdroot(i1:i2)//".gmesh.time."//cstep//".inp"
      end if
      if(iucd.eq.ione) then
        open(kucd,file=filenm,status="new")
      else if(iucd.eq.itwo) then
c
c...  The g95 compiler supports stream i/o, while the gfortran compiler
c     does not.  At present, it is assumed that all f95 compilers other
c     than gfortran support this feature.  Undoubtedly, exceptions will
c     be found and the list below will need to be expanded.
c
#if defined FORTRAN_TYPE_F77
        open(kucd,file=filenm,status="new",access="direct",recl=len,
     &   form="unformatted")
#elif defined FORTRAN_NAME_GFORTRAN
#warning "Binary UCD output not supported with gfortran!"
#else
        open(kucd,file=filenm,status="new",access="stream",
     &   form="unformatted")
#endif
      end if
      return
      end
c
c version
c $Id: open_ucd.F,v 1.1 2005/08/05 19:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
