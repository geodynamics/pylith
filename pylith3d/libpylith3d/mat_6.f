c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  Program segment to define material model 6:
c
c       Model number:                      6
c       Model name:                        IsotropicPowerLawMaxwellViscoelastic
c       Number material properties:        5
c       Number state variables:            18
c       Tangent matrix varies with state:  True
c       Material properties:               Density
c                                          Young's modulus
c                                          Poisson's ratio
c                                          Power-law exponent
c                                          Viscosity coefficient
c
      subroutine mat_prt_6(prop,nprop,matnum,idout,idsk,kw,kp,
     & ierr,errstrng)
c
c...  subroutine to output material properties for material model 6.
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
      character labelp(5)*21,modelname*40
      data labelp/"Density",
     &            "Young's modulus",
     &            "Poisson's ratio",
     &            "Power-law exponent",
     &            "Viscosity coefficient"/
      parameter(modelname="Isotropic Power-Law Maxwell Viscoelastic")
      integer mattype
      parameter(mattype=6)
c
c...  local variables
c
      integer i
c
cdebug      write(6,*) "Hello from mat_prt_6_f!"
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
        errstrng="mat_prt_6"
        return
c
 700  format(/,5x,"Material number:       ",i7,/,5x,
     &            "Material type:         ",a40,/,5x,
     &            "Number of properties:  ",i7,/)
 710  format(15x,a21,3x,1pe15.8)
c
      end
c
c
      subroutine elas_mat_6(dmat,prop,iddmat,nprop,ierr,errstrng)
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
cdebug      write(6,*) "Hello from elas_mat_6_f!"
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
      subroutine elas_strs_6(prop,nprop,state,state0,ee,scur,dmat,tmax,
     & nstate,nstate0,ierr,errstrng)
c
c...  subroutine to compute stresses for the elastic solution. For this
c     material, there are 3 sets of state variables:  total stress,
c     total strain, and viscous strain. The Maxwell time is computed,
c     even though this is the elastic solution, as an aid in determining
c     the proper time step size for the next step.
c     The current total strain is contained in ee and the computed
c     total stress should be copied to scur.
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
      integer nprop,nstate,nstate0,ierr
      double precision prop(nprop),state(nstate),state0(nstate0)
      double precision ee(nstr),scur(nstr)
      double precision dmat(nddmat),tmax
      character errstrng*(*)
c
c...  local variables
c
      double precision e,pr,anpwr,emhu,rmu
      double precision sdev(nstr),sinv1,steff
c
cdebug      write(6,*) "Hello from elas_strs_6_f!"
c
      call dcopy(nstr,ee,ione,state(7),ione)
      call dcopy(nstr,state0,ione,state,ione)
      call dspmv("u",nstr,one,dmat,state(7),ione,one,state,ione)
      call dcopy(nstr,state,ione,scur,ione)
c
c...  compute Maxwell time for current stress state
c
      e=prop(2)
      pr=prop(3)
      anpwr=prop(4)
      emhu=prop(5)
      rmu=half*e/(one+pr)
      call invar(sdev,sinv1,steff,scur)
      tmax=big
      if(steff.ne.zero) tmax=((emhu/steff)**(anpwr-one))*emhu/rmu
      return
      end
c
c
      subroutine td_matinit_6(state,dstate,state0,dmat,prop,rtimdat,
     & rgiter,iopt,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,ierr,
     & errstrng)
c
c...  subroutine to form the material matrix for an integration point
c     for the time-dependent solution.  This routine is meant to be
c     called at the beginning of a time step, before strains have been
c     computed.  Thus, for some time-dependent materials, the material
c     matrix will only be an approximation of the material matrix for
c     the current iteration.
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
      integer iopt,nstate,nstate0,nprop,ierr
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
      double precision e,pr,anpwr,emhu,rmu
      double precision sdev(nstr),sinv1,steff
      double precision sigma,sigmap,sxx,syy,szz,sxy,syz,sxz
      double precision cmat(nddmat)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_matinit_6_f!"
c
c
c...  get stress invariants and compute elastic material matrix
c
      if(iopt.eq.1) then
        call invar(sdev,sinv1,steff,state)
      else
        call invar(sdev,sinv1,steff,dstate)
      end if
      call elas_mat_6(dmat,prop,iddmat,nprop,ierr,errstrng)
      tmax=big
c
c...  if second deviatoric invariant is zero, use only elastic solution
c
      if(steff.eq.zero) then
        return
c
c...  otherwise, invert elastic matrix and augment it by the viscous
c     Jacobian.
c
      else
        call dpptrf('u',nstr,dmat,ierr)
        if(ierr.ne.izero) then
          errstrng="td_matinit_6(1)"
          if(ierr.lt.izero) then
            ierr=117
          else
            ierr=118
          end if
        end if
        call dpptri('u',nstr,dmat,ierr)
        if(ierr.ne.izero) then
          errstrng="td_matinit_6(2)"
          if(ierr.lt.izero) then
            ierr=117
          else
            ierr=119
          end if
        end if
        call fill(cmat,zero,nddmat)
c
c...  define material properties
c
        e=prop(2)
        pr=prop(3)
        anpwr=prop(4)
        emhu=prop(5)
        rmu=half*e/(one+pr)
c
c...  get stress invariants and compute factors for Jacobian
c
        sigma=(steff/emhu)**(anpwr-one)
        sigmap=(anpwr-one)*sigma
        sigma=half*alfap*deltp*sigma/emhu
        sigmap=fourth*alfap*deltp*sigmap/emhu
        sxx=sdev(1)/steff
        syy=sdev(2)/steff
        szz=sdev(3)/steff
        sxy=two*sdev(4)/steff
        syz=two*sdev(5)/steff
        sxz=two*sdev(6)/steff
c
c...  compute viscous Jacobian
c
        cmat(iddmat(1,1))=sigma*two*third+sigmap*sxx*sxx
        cmat(iddmat(1,2))=-sigma*third   +sigmap*sxx*syy
        cmat(iddmat(1,3))=-sigma*third   +sigmap*sxx*szz
        cmat(iddmat(1,4))=                sigmap*sxx*sxy
        cmat(iddmat(1,5))=                sigmap*sxx*syz
        cmat(iddmat(1,6))=                sigmap*sxx*sxz
        cmat(iddmat(2,2))=sigma*two*third+sigmap*syy*syy
        cmat(iddmat(2,3))=-sigma*third   +sigmap*syy*szz
        cmat(iddmat(2,4))=                sigmap*syy*sxy
        cmat(iddmat(2,5))=                sigmap*syy*syz
        cmat(iddmat(2,6))=                sigmap*syy*sxz
        cmat(iddmat(3,3))=sigma*two*third+sigmap*szz*szz
        cmat(iddmat(3,4))=                sigmap*szz*sxy
        cmat(iddmat(3,5))=                sigmap*szz*syz
        cmat(iddmat(3,6))=                sigmap*szz*sxz
        cmat(iddmat(4,4))=sigma*two      +sigmap*sxy*sxy
        cmat(iddmat(4,5))=                sigmap*sxy*syz
        cmat(iddmat(4,6))=                sigmap*sxy*sxz
        cmat(iddmat(5,5))=sigma*two      +sigmap*syz*syz
        cmat(iddmat(5,6))=                sigmap*syz*sxz
        cmat(iddmat(6,6))=sigma*two      +sigmap*sxz*sxz
c
c...  add this to inverted elastic matrix and invert again
c
        call daxpy(nddmat,one,cmat,ione,dmat,ione)
        call dpptrf('u',nstr,dmat,ierr)
        if(ierr.ne.izero) then
          errstrng="td_matinit_6(3)"
          if(ierr.lt.izero) then
            ierr=117
          else
            ierr=118
          end if
        end if
        call dpptri('u',nstr,dmat,ierr)
        if(ierr.ne.izero) then
          errstrng="td_matinit_6(4)"
          if(ierr.lt.izero) then
            ierr=117
          else
            ierr=119
          end if
        end if
        tmax=((emhu/steff)**(anpwr-one))*emhu/rmu
      end if
      return
      end
c
c
      subroutine td_strs_6(state,dstate,state0,ee,scur,dmat,prop,
     & rtimdat,rgiter,iter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,
     & matchg,ierr,errstrng)
c
c...  subroutine to compute the current values for stress, total strain,
c     and viscous strain increment for the time-dependent solution.
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
      integer iter,nstate,nstate0,nprop,ierr
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
c... external functions
c
      double precision sprod
      external sprod
c
c...  local constants
c
c
c...  local variables
c
      integer i
      double precision e,pr,anpwr,emhu,rmu
      double precision sinv1t,sefft,sinv1tdt,sefftdt,beta,betat,betatdt
      double precision rmt,rmtdt
      double precision sdevt(nstr),sdevtdt(nstr),eetdt(nstr)
      double precision dtmp(nddmat)

c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_strs_6!"
c
      ierr=0
c
c...  define material properties
c
      e=prop(2)
      pr=prop(3)
      anpwr=prop(4)
      emhu=prop(5)
      rmu=half*e/(one+pr)
      tmax=big
c
c...  for first iteration, beta is computed using only the stresses at
c     time t.  Otherwise, the stress estimates at time t + dt are also
c     used.
c
      call invar(sdevt,sinv1t,sefft,state)
      if(iter.eq.ione) then
        do i=1,nstr
          beta=half*(sefft/emhu)**(anpwr-one)*sdevt(i)/emhu
          dstate(i+12)=state(i+12)+deltp*beta
          eetdt(i)=ee(i)-dstate(i+12)
        end do
      else
        rmt=one-alfap
        rmtdt=alfap
        call invar(sdevtdt,sinv1tdt,sefftdt,dstate)
        do i=1,nstr
          betat=half*(sefft/emhu)**(anpwr-one)*sdevt(i)/emhu
          betatdt=half*(sefftdt/emhu)**(anpwr-one)*sdevtdt(i)/emhu
          dstate(i+12)=state(i+12)+deltp*(rmt*betat+rmtdt*betatdt)
          eetdt(i)=ee(i)-dstate(i+12)
        end do
      end if
c
c...  compute current stress estimates and copy state variables into the
c     appropriate spot
c
      call elas_mat_6(dtmp,prop,iddmat,nprop,ierr,errstrng)
      call dcopy(nstr,ee,ione,dstate(7),ione)
      call dcopy(nstr,state0,ione,dstate,ione)
      call dspmv("u",nstr,one,dtmp,eetdt,ione,one,dstate,ione)
      call dcopy(nstr,dstate,ione,scur,ione)
c
c...  compute current Maxwell time
c
      call invar(sdevtdt,sinv1tdt,sefftdt,dstate)
      if(sefftdt.ne.zero) tmax=((emhu/sefftdt)**(anpwr-one))*emhu/rmu
c
      return
      end
c
c
      subroutine td_strs_mat_6(state,dstate,state0,ee,scur,dmat,prop,
     & rtimdat,rgiter,iter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,
     & matchg,ierr,errstrng)
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
      integer iter,nstate,nstate0,nprop,ierr
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
c...  external routines
c
c
c...  local constants
c
c
c...  local variables
c
      integer iopt
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
      ierr=0
      iopt=2
c
c...  compute current stress estimates for time-dependent solution
c
      call td_strs_6(state,dstate,state0,ee,scur,dmat,prop,
     & rtimdat,rgiter,iter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,
     & matchg,ierr,errstrng)
c
c...  compute tangent stiffness corresponding to new stress estimates
c
      call td_matinit_6(state,dstate,state0,dmat,prop,rtimdat,
     & rgiter,iopt,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,ierr,
     & errstrng)
c
      return
      end
c
c
      subroutine prestr_mat_6(dmat,prop,tpois,tyoungs,iddmat,ipauto,
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
      call dcopy(nprop,prop,ione,ptmp,ione)
      if(ipauto.eq.ione) then
        ptmp(2)=tyoungs
        ptmp(3)=tpois
      end if
      call elas_mat_6(dmat,ptmp,iddmat,nprop,ierr,errstrng)
      return
      end
c
c
      subroutine get_state_6(state,dstate,sout,nstate)
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
cdebug      write(6,*) "Hello from get_state_6_f!"
c
      call dcopy(nstate,state,ione,sout,ione)
      call dcopy(nstate,dstate,ione,sout(nstatesmax+ione),ione)
c
      return
      end
c
c
      subroutine update_state_6(state,dstate,nstate)
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
cdebug      write(6,*) "Hello from update_state_6_f!"
c
      call daxpy(nstate,sub,state,ione,dstate,ione)
      call daxpy(nstate,one,dstate,ione,state,ione)
c
      return
      end
c
c       
c version
c $Id: mat_6.f,v 1.5 2005/04/01 23:10:17 willic3 Exp $

c Generated automatically by Fortran77Mill on Tue May 18 14:18:50 2004

c End of file 
