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
      subroutine matprt1(prop,nprop,matnum,idout,idsk,kw,kp,ierr,ofile,
     & pfile)
c
c...  subroutine to output material properties for material model 1.
c
c     Error codes:
c         0:  No error
c         1:  Error opening output file
c         2:  Not used
c         3:  Write error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nprop,matnum,idout,idsk,kw,kp,ierr
      double precision prop(nprop)
      character ofile*(*),pfile*(*)
c
c...  parameters
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
c...  open output files
c
      ierr=0
      if(idout.gt.0) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.0) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.1) open(kp,file=pfile,status="old",access="append",
     & form="unformatted")
c
c...  output plot results
c
      if(idsk.eq.0) then
	write(kp,"(3i7)") matnum,mattype,nprop
	write(kp,"(1pe15.8,20(2x,1pe15.8))") (prop(i),i=1,nprop)
      else if(idsk.eq.1) then
	write(kp) matnum,mattype,nprop
	write(kp) prop
      end if
c
c...  output ascii results, if desired
c
      if(idout.gt.0) then
	write(kw,700) matnum,modelname,nprop
	do i=1,nprop
	  write(kw,710) labelp(i),prop(i)
        end do
      end if
c
 700  format("Material number:   ",i7,/,10x,a80,/,
     &       "Number of properties:  ",i7,/)
 710  format(15x,a15,3x,1pe15.8)
c
      return
      end
c
c
      subroutine prtstrs1(state,nstr,nstate,iel,igauss,idout,idsk,kw,kp)
c
c...  subroutine to print out stress components for material type 1.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nstate,iel,igauss,idout,idsk,kw,kp
      double precision state(nstr,nstate)
c
c...  parameters
c
      integer ione
      parameter(ione=1)
c
c...  local variables
c
      double precision str(6)
      integer i
c
      if(idout.gt.0) write(kw,700) iel,igauss,(state(i,1),i=1,nstr)
      if(idsk.eq.0) write(kp,710) (state(i,1),i=1,nstr)
      if(idsk.eq.1) then
        call dcopy(nstr,state(1,1),ione,str,ione)
        write(kp) str
      end if
c
 700  format(2i7,6(2x,1pe15.8))
 710  format(1pe15.8,5(2x,1pe15.8))
c
      return
      end
c
c
      subroutine prtstrn1(state,nstr,nstate,iel,igauss,idout,idsk,kw,kp)
c
c...  subroutine to print out strain components for material type 1.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nstate,iel,igauss,idout,idsk,kw,kp
      double precision state(nstr,nstate)
c
c...  parameters
c
      integer ione
      parameter(ione=1)
c
c...  local variables
c
      double precision str(6)
      integer i
c
      if(idout.gt.0) write(kw,700) iel,igauss,(state(i,2),i=1,nstr)
      if(idsk.eq.0) write(kp,710) (state(i,2),i=1,nstr)
      if(idsk.eq.1) then
        call dcopy(nstr,state(1,2),ione,str,ione)
        write(kp) str
      end if
c
 700  format(2i7,6(2x,1pe15.8))
 710  format(1pe15.8,5(2x,1pe15.8))
c
      return
      end
c
c
      subroutine prtbeta1(state,nstr,nstate,iel,igauss,idout,idsk,kw,kp)
c
c...  subroutine to print out viscous strain components for material
c     type 1.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nstate,iel,igauss,idout,idsk,kw,kp
      double precision state(nstr,nstate)
c
c...  parameters
c
      double precision zero
      parameter(zero=0.0d0)
c
c...  local variables
c
      double precision str(6)
      integer i
c
      call fill(str,zero,nstr)
      if(idout.gt.0) write(kw,700) iel,igauss,(str(i),i=1,nstr)
      if(idsk.eq.0) write(kp,710) (str(i),i=1,nstr)
      if(idsk.eq.1) write(kp) str
c
 700  format(2i7,6(2x,1pe15.8))
 710  format(1pe15.8,5(2x,1pe15.8))
c
      return
      end
c
c
      subroutine prtbetb1(state,nstr,nstate,iel,igauss,idout,idsk,kw,kp)
c
c...  subroutine to print out plastic strain components for material
c     type 1.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nstate,iel,igauss,idout,idsk,kw,kp
      double precision state(nstr,nstate)
c
c...  parameters
c
      double precision zero
      parameter(zero=0.0d0)
c
c...  local variables
c
      double precision str(6)
      integer i
c
      call fill(str,zero,nstr)
      if(idout.gt.0) write(kw,700) iel,igauss,(str(i),i=1,nstr)
      if(idsk.eq.0) write(kp,710) (str(i),i=1,nstr)
      if(idsk.eq.1) write(kp) str
c
 700  format(2i7,6(2x,1pe15.8))
 710  format(1pe15.8,5(2x,1pe15.8))
c
      return
      end
c
c
      subroutine elasmat1(dmat,prop,iddmat,nstr,nddmat,nprop)
c
c...  subroutine to form the material matrix for an integration point
c     for the elastic solution.  The material matrix is assumed to be
c     independent of the state variables in this case.
c     Note also that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nddmat,nprop,nstr
      integer iddmat(nstr,nstr)
      double precision dmat(nddmat),prop(nprop)
c
c...  parameters
c
      integer inc
      double precision zero
      parameter(inc=1,zero=0.0d0)
c
c...  local variables
c
      integer i,j
      double precision e,pr,pr1,pr2,pr3,fac,dd,od,ss
c
      e=prop(2)
      pr=prop(3)
      pr1=1.0d0-pr
      pr2=1.0d0+pr
      pr3=1.0d0-2.0d0*pr
      fac=e/(pr2*pr3)
      dd=pr1*fac
      od=pr*fac
      ss=0.5d0*pr3*fac
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
      subroutine elasstrs1(state,ee,dmat,nstr,nstate,nddmat)
c
c...  subroutine to compute stresses for the elastic solution.  For this
c     material, there are just 2 state variables:  total stress and
c     total strain.  The current total strain is contained in ee.
c
c     state(nstr,1) = Cauchy stress
c     state(nstr,2) = linear strain
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nstate,nddmat
      double precision state(nstr,nstate),ee(nstr),dmat(nddmat)
c
c...  parameters
c
      double precision zero,one
      integer ione
      parameter(zero=0.0d0,one=1.0d0,ione=1)
c
      call dcopy(nstr,ee,ione,state(1,2),ione)
      call dspmv("u",nstr,one,dmat,state(1,2),ione,zero,state(1,1),ione)
      return
      end
c
c
      subroutine tdmat1(state,dmat,prop,iddmat,nstr,nstate,nddmat,nprop,
     & matchg)
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
c...  subroutine arguments
c
      integer nstr,nstate,nddmat,nprop
      integer iddmat(nstr,nstr)
      double precision state(nstr,nstate),dmat(nddmat),prop(nprop)
      logical matchg
c
c...  parameters
c
      integer inc
      double precision zero
      parameter(inc=1,zero=0.0d0)
c
c...  local variables
c
      integer i,j
      double precision e,pr,pr1,pr2,pr3,fac,dd,od,ss
c
      if(.not.matchg) return
c
      call fill(dmat,zero,nddmat)
      e=prop(2)
      pr=prop(3)
      pr1=1.0d0-pr
      pr2=1.0d0+pr
      pr3=1.0d0-2.0d0*pr
      fac=e/(pr2*pr3)
      dd=pr1*fac
      od=pr*fac
      ss=0.5d0*pr3*fac
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
      subroutine tdmatstrs1(state,dmat,prop,iddmat,nstr,nstate,nddmat,
     & nprop)
c
c...  subroutine to compute the current stress and updated material
c     matrix for the time-dependent solution.  Since this is a purely
c     elastic material, the material matrix should not change.  If the
c     material properties have changed for a time step, this should
c     already be reflected in the current dmat.
c     Note that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nstate,nddmat,nprop
      integer iddmat(nstr,nstr)
      double precision state(nstr,nstate),dmat(nddmat),prop(nprop)
c
c...  parameters
c
      double precision zero,one
      integer ione
      parameter(zero=0.0d0,one=1.0d0,ione=1)
c
      call dspmv("u",nstr,one,dmat,state(1,2),ione,zero,state(1,1),ione)
c
      return
      end
c       

c version
c $Id: mat1.f,v 1.3 2004/05/24 21:03:25 willic3 Exp $

c Generated automatically by Fortran77Mill on Tue May 18 14:18:50 2004

c End of file 
