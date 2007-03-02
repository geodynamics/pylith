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
c
c...  invert elastic matrix
c
        call invsymp(nddmat,nstr,dmat,ierr,errstrng)
        if(ierr.ne.izero) return
c
c...  define material properties
c
        e=prop(2)
        pr=prop(3)
        anpwr=prop(4)
        emhu=prop(5)
        rmu=half*e/(one+pr)
c
c...  compute viscous Jacobian and add it to inverted elastic matrix
c
        call dbds_6(cmat,nddmat,iddmat,sdev,nstr,steff,anpwr,emhu,
     &   alfap,deltp)
        call daxpy(nddmat,one,cmat,ione,dmat,ione)
c
c...  invert augmented matrix
c
        call invsymp(nddmat,nstr,dmat,ierr,errstrng)
        if(ierr.ne.izero) return
        tmax=((emhu/steff)**(anpwr-one))*emhu/rmu
      end if
      return
      end
c
c
      subroutine jaccmp_6(scur,n,fvec,NP,fjac,rpar,nrpar,ipar,nipar,
     & ierr,errstrng)
c
c...  subroutine to form the Jacobian used in finding the zero of the
c     stress function.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
      include "parmat_6.inc"
c
c...  subroutine arguments
c
      integer n,NP,nrpar,nipar,ierr
c***  note that the dimension below is a bit kludgy and does not allow
c***  anything else to be put in the ipar array
      integer ipar(nstr,nstr)
      double precision scur(n),fvec(n),fjac(NP,NP),rpar(nrpar)
      character errstrng*(*)
c
c...  local variables
c
      integer i,j
ctest      integer iddmat(nstr,nstr)
ctest      equivalence(ipar(indiddmat),iddmat)
      double precision dtmp(nddmat),cmat(nddmat)
      double precision sdev(nstr),sinv1,steff
c
cdebug      write(6,*) "Hello from jaccmp_6_f!"
c
c
c...  get stress invariants and a copy of inverted elasticity matrix
c
      call invar(sdev,sinv1,steff,scur)
      call dcopy(nddmat,rpar(inddmati),ione,dtmp,ione)
c
c...  if second deviatoric invariant is zero, use only elastic solution
c
      if(steff.eq.zero) then
        call invsymp(nddmat,nstr,dtmp,ierr,errstrng)
        if(ierr.ne.izero) return
        do i=1,nstr
          do j=1,nstr
            fjac(i,j)=dtmp(ipar(i,j))
          end do
        end do
        return
c
c...  otherwise, invert elastic matrix and augment it by the viscous
c     Jacobian.
c
      else
c
c...  compute viscous Jacobian and add it to inverted elastic matrix
c
        call dbds_6(cmat,nddmat,ipar,sdev,nstr,steff,
     &   rpar(indanpwr),rpar(indemhu),rpar(indalfap),rpar(inddeltp))
        call daxpy(nddmat,one,cmat,ione,dtmp,ione)
c
c...  invert augmented matrix and store results in fjac
c
        call invsymp(nddmat,nstr,dtmp,ierr,errstrng)
        do i=1,nstr
          do j=1,nstr
            fjac(i,j)=dtmp(ipar(i,j))
          end do
        end do
      end if
      return
      end
c
c
      subroutine dbds_6(cmat,nddmat,iddmat,sdev,nstr,steff,anpwr,emhu,
     & alfap,deltp)
c
c...  subroutine to compute derivatives of viscous strain rate with
c     respect to stress.
c
      implicit none
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nddmat,nstr
      integer iddmat(nstr,nstr)
      double precision cmat(nddmat),sdev(nstr),steff,anpwr,emhu
      double precision alfap,deltp
c
c...  local variables
c
      double precision sigma,sigmap,sxx,syy,szz,sxy,syz,sxz
c
cdebug      write(6,*) "Hello from dbds_6!"
c
      call fill(cmat,zero,nddmat)
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
      return
      end
c
c
      subroutine betacmp_6(strs,beta,sdev,seff,sinv1,emhu,anpwr,nstr)
c
c...  subroutine to compute viscous strain rate vector beta
c     Routine also returns stress invariants
c
      implicit none
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstr
      double precision strs(nstr),beta(nstr),sdev(nstr),seff,sinv1
      double precision emhu,anpwr
c
c...  local variables
c
      double precision rm
c
c... compute stress invariants and constants for loop
c
      call fill(beta,zero,nstr)
      call invar(sdev,sinv1,seff,strs)
      rm=half*(seff/emhu)**(anpwr-one)/emhu
c
c...  compute strain rate components
c
      call daxpy(nstr,rm,sdev,ione,beta,ione)
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
      include "parmat_6.inc"
      integer nrpar,nipar
      parameter(nrpar=31,nipar=36)
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
      double precision sprod,dasum
      external sprod,funcv_6,jaccmp_6,dasum
c
c...  local constants
c
      double precision sub
      data sub/-1.0d0/
c
c...  local variables
c
      integer i
      double precision e,pr,anpwr,emhu,rmu
      double precision test0,rm,rmt,rmtdt
      double precision sdev(nstr),betat(nstr),betatdt(nstr),sinv1,seff
      double precision rpar(nrpar)
      logical check
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
      test0=dasum(nstate0,state0,ione)
c
c...  define material properties
c
      e=prop(2)
      pr=prop(3)
      anpwr=prop(4)
      emhu=prop(5)
      rmu=half*e/(one+pr)
      tmax=big
      rmt=deltp*(one-alfap)
      rm=sub*rmt
      rmtdt=deltp*alfap
c
c...  for first iteration, use stresses from previous step as initial
c     guess, otherwise use current estimate.
c
      if(iter.eq.ione) call dcopy(nstr,state,ione,dstate,ione)
c
c...  compute constant part of iterative solution equation.
c
c...  current total strain
      call dcopy(nstr,ee,ione,rpar(indconst),ione)
c...  inverted elasticity matrix times initial stresses
      call elas_mat_6(rpar(inddmati),prop,iddmat,nprop,ierr,errstrng)
      call invsymp(nddmat,nstr,rpar(inddmati),ierr,errstrng)
      if(ierr.ne.izero) return
      if(test0.ne.zero) call dspmv("u",nstr,one,rpar(inddmati),state0,
     & ione,one,rpar(indconst),ione)
c...  viscous strain from previous time step
      call daxpy(nstr,sub,state(13),ione,rpar(indconst),ione)
c...  contribution to viscous strain from stresses at beginning of step
      call betacmp_6(state,betat,sdev,seff,sinv1,emhu,anpwr,nstr)
      call daxpy(nstr,rm,betat,ione,rpar(indconst),ione)
c
c...  finish defining parameters for stress computation then call Newton
c     routine to compute stresses.
c
      rpar(indalfap)=alfap
      rpar(inddeltp)=deltp
      rpar(indemhu)=emhu
      rpar(indanpwr)=anpwr
      call newt(dstate,nstr,rpar,nrpar,iddmat,nipar,funcv_6,jaccmp_6,
     & check,ierr,errstrng)
      if(ierr.ne.izero) return
c
c...  update state variables for current stress estimates.
c
      call betacmp_6(dstate,betatdt,sdev,seff,sinv1,emhu,anpwr,nstr)
      call dcopy(nstr,ee,ione,dstate(7),ione)
      call dcopy(nstr,state(13),ione,dstate(13),ione)
      call daxpy(nstr,rmt,betat,ione,dstate(13),ione)
      call daxpy(nstr,rmtdt,betatdt,ione,dstate(13),ione)
      call dcopy(nstr,dstate,ione,scur,ione)
c
c...  compute current Maxwell time
c
      call invar(sdev,sinv1,seff,dstate)
      if(seff.ne.zero) tmax=((emhu/seff)**(anpwr-one))*emhu/rmu
c
      return
      end
c
c
      subroutine funcv_6(n,scur,r,rpar,nrpar,ipar,nipar,
     & ierr,errstrng)
c
c...  routine to compute the vector of functions(r) evaluated at the
c     given stress state (contained in scur).
c
c...  The rpar array contains:
c     1-6:      Constant part of vector, which includes current total
c               strain, viscous strain from previous step, stress-related
c               quantities from prefious time step, and prestress term.
c     7-27:     Inverted elastic material matrix.
c     28:       Integration parameter alfap.
c     29:       Time step size deltp.
c     30:       Viscosity coefficient emhu.
c     31:       Power-law exponent anpwr.
c
c...  The ipar array contains:
c     1-36:     The iddmat index array
c
      implicit none
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
      include "parmat_6.inc"
c
c...  subroutine arguments
c
      integer n,nrpar,nipar,ierr
      integer ipar(nipar)
      double precision scur(n),r(n),rpar(nrpar)
      character errstrng*(*)
c
c...  local data
c
      double precision sub
      data sub/-1.0d0/
c
c...  local variables
c
      double precision betatdt(nstr),rmtdt
      double precision sdev(nstr),seff,sinv1
c
c...  viscous strain rate multiplier
c
      rmtdt=-rpar(indalfap)*rpar(inddeltp)
c
c...  compute current viscous strain rate
c
      call betacmp_6(scur,betatdt,sdev,seff,sinv1,
     & rpar(indemhu),rpar(indanpwr),nstr)
c
c...  compute contributions to r and accumulate
c
      call dcopy(nstr,rpar(indconst),ione,r,ione)
      call daxpy(nstr,rmtdt,betatdt,ione,r,ione)
      call dspmv("u",nstr,sub,rpar(inddmati),scur,ione,one,r,ione)
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
