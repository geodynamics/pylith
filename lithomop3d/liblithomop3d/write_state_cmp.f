c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c 
c 	  All worldwide rights reserved.  A license to use, copy, modify and
c         distribute this software for non-commercial research purposes only
c         is hereby granted, provided that this copyright notice and
c         accompanying disclaimer is not modified or removed from the software.
c     
c         DISCLAIMER:  The software is distributed "AS IS" without any express
c         or implied warranty, including but not limited to, any implied
c         warranties of merchantability or fitness for a particular purpose
c         or any warranty of non-infringement of any current or pending patent
c         rights.  The authors of the software make no representations about
c         the suitability of this software for any particular purpose.  The
c         entire risk as to the quality and performance of the software is with
c         the user.  Should the software prove defective, the user assumes the
c         cost of all necessary servicing, repair or correction.  In
c         particular, neither Rensselaer Polytechnic Institute, nor the authors
c         of the software are liable for any indirect, special, consequential,
c         or incidental damages related to the software, to the maximum extent
c         the law permits.
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine write_state_cmp(
     & state,dstate,nelfamily,nstate,ielg,                              ! elemfamily
     & get_state,                                                       ! materl
     & ngauss,                                                          ! eltype
     & delti,nstep,                                                     ! timdat
     & itmp1,itmp2,nstatestot,nstatesinc,nstatesrate,nstatesout,        ! ioopts
     & idout,idsk,iucd,kw,kp,kucd)                                      ! ioinfo
c
c...  computation routine to update state variables within an element family
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nelfamily,nstate,ielg,ngauss,nstep,nstatestot,nstatesinc
      integer nstatesrate,nstatesout,idout,idsk,iucd,kw,kp,kucd
      integer itmp1(nstatesmax,3),itmp2(nstatesout)
      double precision state(nstate,ngauss,nelfamily)
      double precision dstate(nstate,ngauss,nelfamily)
      double precision delti
c
c...  external routines
c
      external get_state
c
c...  included dimension and type statements
c
      include "labels_dim.inc"
c
c...  local constants
c
      integer npage
      data npage/50/
c
c...  local variables
c
      integer ielf,l,nout
      double precision sout(3*nstatesmax)
c
c...  included variable definitions
c
      include "labels_dim.inc"
c
      nout=izero
c
c...  loop over elements in a family and write requested output
c
      do ielf=1,nelfamily
	do l=1,ngauss
          nout=nout+ione
	  call get_state(state,dstate,sout,nstate)
          call daxpy(nstatesmax,delti,sout(nstatesmax+ione),ione,
     &     sout(2*nstatesmax+ione),ione)
	end do
        ielg=ielg+ione
      end do
c
      return
      end
c
c version
c $Id: write_state_cmp.f,v 1.1 2005/03/23 02:44:50 willic3 Exp $
c
c End of file 
