c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  Copyright 2005 Rensselaer Polytechnic Institute.
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
      subroutine write_props(prop,grav,ivfamily,infmatmod,nvfamilies,
     & npropsz,idout,idsk,kw,kp,ofile,pfile,ierr,errstrng)
c
c...  program to print properties for each material type.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nvfamilies,npropsz,idout,idsk,kw,kp,ierr
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      double precision prop(npropsz),grav(ndof)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  external routines
c
      include "mat_prt_ext.inc"
c
c...  local variables
c
      integer numat,ifam,matmodel,indprop,nprop,i
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
cdebug      write(6,*) "Hello from write_props_f!"
c
      ierr=izero
      numat=nvfamilies
c
c...  open output files and output number of materials
c
      if(idout.gt.izero) then
        open(kw,file=ofile,err=10,status="old",access="append")
        write(kw,1000,err=20) numat
      end if
      if(idsk.eq.ione) then
        open(kp,file=pfile,err=10,status="old",access="append")
        write(kp,"(i5)",err=20) numat
      end if
      if(idsk.eq.itwo) then
        open(kp,file=pfile,err=10,status="old",access="append",
     &   form="unformatted")
        write(kp,err=20) numat
      end if
c
c...  loop over number of element families, which presently correspond
c     to material types
c
      do ifam=1,nvfamilies
        matmodel=ivfamily(2,ifam)
        indprop=ivfamily(5,ifam)
        nprop=infmatmod(3,matmodel)
        if(matmodel.eq.1) then
          call mat_prt_1(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.2) then
          call mat_prt_2(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.3) then
          call mat_prt_3(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.4) then
          call mat_prt_4(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.5) then
          call mat_prt_5(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.6) then
          call mat_prt_6(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.7) then
          call mat_prt_7(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.8) then
          call mat_prt_8(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.9) then
          call mat_prt_9(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.10) then
          call mat_prt_10(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.11) then
          call mat_prt_11(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.12) then
          call mat_prt_12(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.13) then
          call mat_prt_13(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.14) then
          call mat_prt_14(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.15) then
          call mat_prt_15(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.16) then
          call mat_prt_16(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.17) then
          call mat_prt_17(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.18) then
          call mat_prt_18(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.19) then
          call mat_prt_19(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        else if(matmodel.eq.20) then
          call mat_prt_20(prop(indprop),nprop,ifam,idout,idsk,kw,kp,
     &     ierr,errstrng)
        end if
        if(ierr.ne.izero) return
      end do
c
c...  output gravitational acceleration
c
      if(idout.gt.izero) then
        write(kw,1100,err=20)
        do i=1,ndof
          write(kw,1200,err=20) labeld(i),grav(i)
        end do
      end if
c
c...  normal return
c
      if(idout.gt.izero) close(kw)
      close(kp)
      return
c
c...  error opening output file
c
 10   continue
        ierr=2
        errstrng="write_props"
        if(idout.gt.izero) close(kw)
        close(kp)
        return
c
c...  error writing to output file
c
 20   continue
        ierr=4
        errstrng="write_props"
        if(idout.gt.izero) close(kw)
        close(kp)
        return
c
 1000 format(1x,///,
     & " m a t e r i a l   s e t   d a t a                    ",//,5x,
     & " number of material sets . . . . . . . . . . (numat) =",i5)
 1100 format(//,
     & ' a c c e l e r a t i o n   o f   g r a v i t y         ',//)
 1200 format(5x,
     & a4,'-direction  . . . . . . . . . . . . . . . . =',1pe17.5/)
c
      end
c
c version
c $Id: write_props.f,v 1.5 2005/03/22 19:59:57 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
