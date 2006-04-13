c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  Program segment to define material model 5:
c
c       Model number:                      5
c       Model name:                        IsotropicLinearMaxwellViscoelastic
c       Number material properties:        4
c       Number state variables:            3
c       Tangent matrix varies with state:  True
c       Material properties:               Density
c                                          Young's modulus
c                                          Poisson's ratio
c                                          Viscosity
c
      subroutine mat_prt_5(prop,nprop,matnum,idout,idsk,kw,kp,
     & ierr,errstrng)
c
c...  subroutine to output material properties for material model 5.
c
c     Error codes:
c         4:  Write error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nprop,matnum,idout,idsk,kw,kp,ierr
      double precision prop(nprop)
      character errstrng*(*)
c
c...  local constants
c
      character labelp(4)*15,modelname*37
      data labelp/"Density",
     &            "Young's modulus",
     &            "Poisson's ratio",
     &            "Viscosity"/
      parameter(modelname="Isotropic Linear Maxwell Viscoelastic")
      integer mattype
      parameter(mattype=5)
c
c...  local variables
c
      integer i
c
cdebug      write(6,*) "Hello from mat_prt_5_f!"
c
c...  output plot results
c
      if(idsk.eq.ione) then
	write(kp,"(3i7)",err=10) matnum,mattype,nprop
	write(kp,"(1pe15.8,20(2x,1pe15.8))",err=10) (prop(i),i=1,nprop)
      else if(idsk.eq.itwo) then
	write(kp,err=10) matnum,mattype,nprop
	write(kp,err=10) prop
      end if
c
c...  output ascii results, if desired
c
      if(idout.gt.izero) then
	write(kw,700,err=10) matnum,modelname,nprop
	do i=1,nprop
	  write(kw,710,err=10) labelp(i),prop(i)
        end do
      end if
c
      return
c
c...  error writing to output file
c
 10   continue
        ierr=4
        errstrng="mat_prt_5"
        return
c
 700  format(/,5x,"Material number:       ",i7,/,5x,
     &            "Material type:         ",a37,/,5x,
     &            "Number of properties:  ",i7,/)
 710  format(15x,a15,3x,1pe15.8)
c
      end
c
c
      subroutine elas_mat_5(dmat,prop,iddmat,nprop,ierr,errstrng)
c
c...  subroutine to form the material matrix for an integration point
c     for the elastic solution.  The material matrix is assumed to be
c     independent of the state variables in this case.
c     Note also that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision dmat(nddmat),prop(nprop)
c
c...  local variables
c
      integer i,j
      double precision e,pr,pr1,pr2,pr3,fac,dd,od,ss
c
cdebug      write(6,*) "Hello from elas_mat_5_f!"
c
      call fill(dmat,zero,nddmat)
      e=prop(2)
      pr=prop(3)
      pr1=one-pr
      pr2=one+pr
      pr3=one-two*pr
      fac=e/(pr2*pr3)
      dd=pr1*fac
      od=pr*fac
      ss=half*pr3*fac
      do i=1,3
        dmat(iddmat(i,i))=dd
        dmat(iddmat(i+3,i+3))=ss
        do j=i+1,3
          dmat(iddmat(i,j))=od
        end do
      end do
      return
      end
c
c
      subroutine elas_strs_5(state,state0,ee,scur,dmat,nstate,nstate0,
     & ierr,errstrng)
c
c...  subroutine to compute stresses for the elastic solution.  For this
c     material, there are 3 state variables:  total stress, total
c     strain, and viscous strain.  The current total strain is contained
c     in ee.
c
c     state(1:6)   = Cauchy stress
c     state(7:12)  = linear strain
c     state(13:18) = viscous strain
c
c     The state0 array contains initial stresses.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,nstate0,ierr
      double precision state(nstate),state0(nstate0),ee(nstr),scur(nstr)
      double precision dmat(nddmat)
      character errstrng*(*)
c
cdebug      write(6,*) "Hello from elas_strs_5_f!"
c
      call dcopy(nstr,ee,ione,state(7),ione)
      call dcopy(nstr,state0,ione,state,ione)
      call dspmv("u",nstr,one,dmat,state(7),ione,one,state,ione)
      call dcopy(nstr,state,ione,scur,ione)
      return
      end
c
c
      subroutine td_matinit_5(state,dstate,state0,dmat,prop,rtimdat,
     & rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,ierr,
     & errstrng)
c
c...  subroutine to form the material matrix for an integration point
c     for the time-dependent solution.  This routine is meant to be
c     called at the beginning of a time step, before strains have been
c     computed.  Thus, for some time-dependent materials, the material
c     matrix will only be an approximation of the material matrix for
c     the current iteration.
c     As this is a linear viscoelastic material, the tangent matrix is
c     actually independent of the current stresses and strains, so they
c     are not used.
c     Note that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,nstate0,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision state(nstate),dstate(nstate)
      double precision state0(nstate0),dmat(nddmat),prop(nprop),tmax
      logical matchg
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  local variables
c
      double precision e,pr,vis,f1,f2
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_matinit_5_f!"
c
      call fill(dmat,zero,nddmat)
      tmax=big
      e=prop(2)
      pr=prop(3)
      vis=prop(4)
      f1=third*e/(one-two*pr)
      f2=third/((one+pr)/e+half*alfap*deltp/vis)
      dmat(iddmat(1,1))=f1+two*f2
      dmat(iddmat(2,2))=dmat(iddmat(1,1))
      dmat(iddmat(3,3))=dmat(iddmat(1,1))
      dmat(iddmat(1,2))=f1-f2
      dmat(iddmat(1,3))=dmat(iddmat(1,2))
      dmat(iddmat(2,3))=dmat(iddmat(1,2))
      dmat(iddmat(4,4))=half*three*f2
      dmat(iddmat(5,5))=dmat(iddmat(4,4))
      dmat(iddmat(6,6))=dmat(iddmat(4,4))
      return
      end
c
c
      subroutine td_strs_5(state,dstate,state0,ee,scur,dmat,prop,
     & rtimdat,rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,
     & ierr,errstrng)
c
c...  subroutine to compute the current values for stress, total strain,
c     and viscous strain increment.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,nstate0,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      logical matchg
      double precision state(nstate),dstate(nstate),state0(nstate0)
      double precision ee(nstr),scur(nstr),dmat(nddmat),prop(nprop),tmax
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  local constants
c
      double precision diag(6)
      data diag/one,one,one,zero,zero,zero/
c
c...  local variables
c
      integer i
      double precision e,pr,vis,rmu,f1,f2,emean,smean,eet,stau,smeantp
      double precision sdev,sdevtp,smean0,sdev0,fac1,fac2,rvis
      double precision eett(6)
c
cdebug      integer idb
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_strs_5!"
c
      call dcopy(nstr,ee,ione,eett,ione)
      eett(4)=half*eett(4)
      eett(5)=half*eett(5)
      eett(6)=half*eett(6)
      e=prop(2)
      pr=prop(3)
      vis=prop(4)
      rmu=half*e/(one+pr)
      rvis=one/vis
      tmax=two*vis/rmu
      fac1=half*deltp*rvis
      fac2=(one+pr)/e
      f1=fac1*(one-alfap)
      f2=one/(fac2+fac1*alfap)
      emean=third*(ee(1)+ee(2)+ee(3))
      smean=e*emean/(one-two*pr)
      smean0=third*(state0(1)+state0(2)+state0(3))
      smeantp=third*(state(1)+state(2)+state(3))
cdebug      write(6,*) "state-begin:",(state(idb),idb=1,nstate)
cdebug      write(6,*) "dstate-begin:",(dstate(idb),idb=1,nstate)
cdebug      write(6,*) "ee:",(ee(idb),idb=1,nstr)
cdebug      write(6,*) "eett:",(eett(idb),idb=1,nstr)
cdebug      write(6,*)
cdebug     & "e,pr,rmu,rvis,tmax,fac1,fac2,f1,f2,emean,smean,smean0,smeantp:"
cdebug      write(6,*)
cdebug     & e,pr,rmu,rvis,tmax,fac1,fac2,f1,f2,emean,smean,smean0,smeantp
      do i=1,nstr
        eet=eett(i)-diag(i)*emean-state(i+12)
        sdevtp=state(i)-diag(i)*smeantp
        sdev0=state0(i)-diag(i)*smean0
        sdev=f2*(eet-f1*sdevtp+fac2*sdev0)
        dstate(i)=sdev+diag(i)*(smean+smean0)
        stau=(one-alfap)*sdevtp+alfap*sdev
        dstate(i+12)=fac1*stau
        dstate(i+6)=ee(i)
        scur(i)=dstate(i)
cdebug        write(6,*) "i,eet,eett(i),state(i+12),sdevtp,sdev0,sdev:"
cdebug        write(6,*) i,eet,eett(i),state(i+12),sdevtp,sdev0,sdev
cdebug        write(6,*) "dstate(i),stau,dstate(i+6),dstate(i+12),scur(i):"
cdebug        write(6,*) dstate(i),stau,dstate(i+6),dstate(i+12),scur(i)
      end do
cdebug      write(6,*) "state-end:",(state(idb),idb=1,nstate)
cdebug      write(6,*) "dstate-end:",(dstate(idb),idb=1,nstate)
c
      return
      end
c
c
      subroutine td_strs_mat_5(state,dstate,state0,ee,scur,dmat,prop,
     & rtimdat,rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,
     & ierr,errstrng)
c
c...  subroutine to compute the current stress and updated material
c     matrix for the time-dependent solution.
c     Note that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,nstate0,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      logical matchg
      double precision state(nstate),dstate(nstate),state0(nstate0)
      double precision ee(nstr),scur(nstr),dmat(nddmat),prop(nprop),tmax
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_strs_mat_5_f!"
c
      call td_matinit_5(state,dstate,state0,dmat,prop,rtimdat,rgiter,
     & ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,ierr,errstrng)
      call td_strs_5(state,dstate,state0,ee,scur,dmat,prop,rtimdat,
     & rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,ierr,
     & errstrng)
c
      return
      end
c
c
      subroutine prestr_mat_5(dmat,prop,tpois,tyoungs,iddmat,ipauto,
     & nprop,ierr,errstrng)
c
c...  subroutine to form the material matrix for an integration point
c     for prestress computation.  The material matrix is assumed to be
c     independent of the state variables in this case.
c     Note also that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer ipauto,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision tpois,tyoungs,dmat(nddmat),prop(nprop)
c
c...  local variables
c
      double precision ptmp(10)
c
cdebug      write(6,*) "Hello from prestr_mat_5_f!"
c
      call dcopy(nprop,prop,ione,ptmp,ione)
      if(ipauto.eq.ione) then
        ptmp(2)=tyoungs
        ptmp(3)=tpois
      end if
      call elas_mat_5(dmat,ptmp,iddmat,nprop,ierr,errstrng)
      return
      end
c       
c       
      subroutine get_state_5(state,dstate,sout,nstate)
c
c...  routine to transfer state variables into sout
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nstate
      double precision state(nstate),dstate(nstate),sout(3*nstatesmax)
c
cdebug      write(6,*) "Hello from get_state_5_f!"
c
      call dcopy(nstate,state,ione,sout,ione)
      call dcopy(nstate,dstate,ione,sout(nstatesmax+ione),ione)
      return
      end
c
c
      subroutine update_state_5(state,dstate,nstate)
c
c...  routine to update state variables at the end of a time step.
c     After updating, state should contain the current total values
c     and dstate should contain the incremental changes since the
c     previous time step.
c     On input, dstate contains the current stress and strain values and
c     state contains the values from the previous time step.
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
      integer nstate
      double precision state(nstate),dstate(nstate)
c
c...  local data
c
      double precision sub
      data sub/-1.0d0/
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from update_state_5_f!"
c
cdebug      write(6,*) "state-begin:",(state(idb),idb=1,nstate)
cdebug      write(6,*) "dstate-begin:",(dstate(idb),idb=1,nstate)
      call daxpy(nstate-6,sub,state,ione,dstate,ione)
      call daxpy(nstate,one,dstate,ione,state,ione)
cdebug      write(6,*) "state-end:",(state(idb),idb=1,nstate)
cdebug      write(6,*) "dstate-end:",(dstate(idb),idb=1,nstate)
c
      return
      end
c
c
c version
c $Id: mat_5.f,v 1.8 2005/04/01 23:07:43 willic3 Exp $

c Generated automatically by Fortran77Mill on Tue May 18 14:18:50 2004

c End of file 
