c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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
c...  Program segment to define material model 1:
c
c       Model number:                      1
c       Model name:                        IsotropicLinearElastic
c       Number material properties:        3
c       Number state variables:            2
c       Tangent matrix varies with state:  False
c       Material properties:               Density
c                                          Young's modulus
c                                          Poisson's ratio
c
      subroutine mat_prt_1(prop,nprop,matnum,idout,idsk,kw,kp,
     & ierr,errstrng)
c
c...  subroutine to output material properties for material model 1.
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
      character labelp(3)*15,modelname*24
      data labelp/"Density",
     &            "Young's modulus",
     &            "Poisson's ratio"/
      parameter(modelname="Isotropic Linear Elastic")
      integer mattype
      parameter(mattype=1)
c
c...  local variables
c
      integer i
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
        errstrng="mat_prt_1"
        return
c
 700  format(/,5x,"Material number:       ",i7,/,5x,
     &            "Material type:         ",a24,/,5x,
     &            "Number of properties:  ",i7,/)
 710  format(15x,a15,3x,1pe15.8)
c
      end
c
c
      subroutine elas_mat_1(dmat,prop,iddmat,nprop,ierr,errstrng)
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
      subroutine elas_strs_1(state,state0,ee,scur,dmat,nstate,nstate0,
     & ierr,errstrng)
c
c...  subroutine to compute stresses for the elastic solution.  For this
c     material, there are just 2 sets of state variables:  total stress
c     and total strain.  The current total strain is contained in ee and
c     the computed total stress should be copied to scur.
c
c     state(1:6)  = Cauchy stress
c     state(7:12) = linear strain
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
      call dcopy(nstr,ee,ione,state(7),ione)
      call dcopy(nstr,state0,ione,state,ione)
      call dspmv("u",nstr,one,dmat,state(7),ione,one,state,ione)
      call dcopy(nstr,state,ione,scur,ione)
      return
      end
c
c
      subroutine td_matinit_1(state,dstate,state0,dmat,prop,rtimdat,
     & rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,ierr,
     & errstrng)
c
c...  subroutine to form the material matrix for an integration point
c     for the time-dependent solution.  This routine is meant to be
c     called at the beginning of a time step, before strains have been
c     computed.  Thus, for some time-dependent materials, the material
c     matrix will only be an approximation of the material matrix for
c     the current iteration.
c     Since this is an elastic material, the material matrix remains
c     unchanged unless a material property has changed (indicated by the
c     matchg flag).
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
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
      tmax=big
      if(.not.matchg) return
c
      call elas_mat_1(dmat,prop,iddmat,nprop,ierr,errstrng)
      return
      end
c
c
      subroutine td_strs_1(state,dstate,state0,ee,scur,dmat,prop,
     & rtimdat,rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,
     & ierr,errstrng)
c
c...  subroutine to compute the current stress for the time-dependent
c     solution.
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
      tmax=big
      call dcopy(nstr,ee,ione,dstate(7),ione)
      call dcopy(nstr,state0,ione,dstate,ione)
      call dspmv("u",nstr,one,dmat,dstate(7),ione,one,dstate,ione)
      call dcopy(nstr,dstate,ione,scur,ione)
c
      return
      end
c
c
      subroutine td_strs_mat_1(state,dstate,state0,ee,scur,dmat,prop,
     & rtimdat,rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,
     & ierr,errstrng)
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
      if(matchg) call elas_mat_1(dmat,prop,iddmat,nprop,ierr,errstrng)
      call td_strs_1(state,dstate,state0,ee,scur,dmat,prop,rtimdat,
     & rgiter,ntimdat,iddmat,tmax,nstate,nstate0,nprop,matchg,
     & ierr,errstrng)
c
      return
      end
c
c
      subroutine prestr_mat_1(dmat,prop,tpois,tyoungs,iddmat,ipauto,
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
      call elas_mat_1(dmat,ptmp,iddmat,nprop,ierr,errstrng)
      return
      end
c
c
      subroutine get_state_1(state,dstate,sout,nstate)
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
      call dcopy(nstate,state,ione,sout,ione)
      call dcopy(nstate,dstate,ione,sout(nstatesmax+ione),ione)
      return
      end
c
c
      subroutine update_state_1(state,dstate,nstate)
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
c...  local variables
c
      call daxpy(nstate,sub,state,ione,dstate,ione)
      call daxpy(nstate,one,dstate,ione,state,ione)
c
      return
      end
c
c       

c version
c $Id: mat_1.f,v 1.19 2005/04/01 23:10:16 willic3 Exp $

c Generated automatically by Fortran77Mill on Tue May 18 14:18:50 2004

c End of file 
