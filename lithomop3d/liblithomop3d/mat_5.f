c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c...  output plot results
c
      if(idsk.eq.izero) then
	write(kp,"(3i7)",err=10) matnum,mattype,nprop
	write(kp,"(1pe15.8,20(2x,1pe15.8))",err=10) (prop(i),i=1,nprop)
      else if(idsk.eq.ione) then
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
cdebug      integer idb
c
cdebug      write(6,*) "Hello from elas_mat_5_f!"
c
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
cdebug      write(6,*) "dmat:",(dmat(idb),idb=1,nddmat)
      return
      end
c
c
      subroutine elas_strs_5(state,ee,dmat,nstate,ierr,errstrng)
c
c...  subroutine to compute stresses for the elastic solution.  For this
c     material, there are 3 state variables:  total stress, total
c     strain, and viscous strain.  The current total strain is contained
c     in ee.
c
c     state(nstr,1) = Cauchy stress
c     state(nstr,2) = linear strain
c     state(nstr,3) = viscous strain
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
      integer nstate,ierr
      double precision state(nstr,nstate),ee(nstr),dmat(nddmat)
      character errstrng*(*)
c
      call dcopy(nstr,ee,ione,state(1,2),ione)
      call dspmv("u",nstr,one,dmat,state(1,2),ione,zero,state(1,1),ione)
      return
      end
c
c
      subroutine td_matinit_5(state,dstate,dmat,prop,rtimdat,rgiter,
     & ntimdat,iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
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
      integer nstate,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision state(nstr,nstate),dstate(nstr,nstate)
      double precision dmat(nddmat),prop(nprop),tmax
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
cdebug      integer idb
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_matinit_5_f!"
c
      tmax=big
      e=prop(2)
      pr=prop(3)
      vis=prop(4)
      f1=third*e/(one-two*pr)
      f2=one/((one+pr)/e+half*alfap*deltp/vis)
      dmat(iddmat(1,1))=f1+two*third*f2
      dmat(iddmat(2,2))=dmat(iddmat(1,1))
      dmat(iddmat(3,3))=dmat(iddmat(1,1))
      dmat(iddmat(1,2))=f1-third*f2
      dmat(iddmat(1,3))=dmat(iddmat(1,2))
      dmat(iddmat(2,3))=dmat(iddmat(1,2))
      dmat(iddmat(4,4))=half*f2
      dmat(iddmat(5,5))=dmat(iddmat(4,4))
      dmat(iddmat(6,6))=dmat(iddmat(4,4))
cdebug      write(6,*) "dmat:",(dmat(idb),idb=1,nddmat)
      return
      end
c
c
      subroutine td_strs_5(state,dstate,ee,dmat,prop,rtimdat,rgiter,
     & ntimdat,iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
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
      integer nstate,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      logical matchg
      double precision state(nstr,nstate),dstate(nstr,nstate),ee(nstr)
      double precision dmat(nddmat),prop(nprop),tmax
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
      double precision sdev,sdevtp
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_strs_5!"
c
      e=prop(2)
      pr=prop(3)
      vis=prop(4)
      rmu=half*e/(one+pr)
      tmax=two*vis/rmu
      f1=half*deltp*(one-alfap)/vis
      f2=one/((one+pr)/e+half*alfap*deltp/vis)
      emean=third*(ee(1)+ee(2)+ee(3))
      smean=e*emean/(one-two*pr)
      smeantp=third*(state(1,1)+state(2,1)+state(3,1))
cdebug      write(6,*) "e,pr,vis,rmu,tmax,f1,f2,emean,smean,smeantp:",
cdebug     & e,pr,vis,rmu,tmax,f1,f2,emean,smean,smeantp
      do i=1,nstr
        eet=ee(i)-diag(i)*emean-state(i,3)
        sdevtp=state(i,1)-diag(i)*smeantp
        sdev=f2*(eet-f1*sdevtp)
        dstate(i,1)=sdev+diag(i)*smean
        stau=(one-alfap)*sdevtp+alfap*sdev
        dstate(i,3)=half*deltp*stau/vis
        dstate(i,2)=ee(i)
cdebug        write(6,*) "i,eet,sdevtp,sdev,ds1,ds2,ds3,s1,s2,s3,stau:",
cdebug     &   i,eet,sdevtp,sdev,dstate(i,1),dstate(i,2),dstate(i,3),
cdebug     &   state(i,1),state(i,2),state(i,3),stau
      end do
c
      return
      end
c
c
      subroutine td_strs_mat_5(state,dstate,ee,dmat,prop,rtimdat,rgiter,
     & ntimdat,iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
c
c...  subroutine to compute the current stress and updated material
c     matrix for the time-dependent solution.  Since this is a purely
c     elastic material, the material matrix should not change unless the
c     material properties have changed for a time step.
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
      integer nstate,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      logical matchg
      double precision state(nstr,nstate),dstate(nstr,nstate),ee(nstr)
      double precision dmat(nddmat),prop(nprop),tmax
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
      if(matchg) call td_matinit_5(state,dstate,dmat,prop,rtimdat,
     & rgiter,ntimdat,iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
      call td_strs_5(state,dstate,ee,dmat,prop,rtimdat,rgiter,ntimdat,
     & iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
c
      return
      end
c       

c version
c $Id: mat_5.f,v 1.3 2004/08/12 02:01:55 willic3 Exp $

c Generated automatically by Fortran77Mill on Tue May 18 14:18:50 2004

c End of file 
