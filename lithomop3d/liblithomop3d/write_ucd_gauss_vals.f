c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine write_ucd_gauss_vals(
     & state,dstate,infiel,nstatesz,numelt,                             ! elemnt
     & infmat,infmatmod,ismatmod,numat,                                 ! materl
     & infetype,                                                        ! eltype
     & delt,nstep,                                                      ! timdat
     & istatout,                                                        ! ioopts
     & kucd,ucdroot,iprestress)                                         ! ioinfo
c
c...  Specialized routine to output element info for SCEC benchmarks.
c     This routine creates the nodal value portion of the UCD file for
c     each time step, including the header information.  Note that a
c     complete UCD file is formed by concatenating the Gauss output from
c     write_ucdmesh with the file created by this routine.
c     The state variables output are determined from the istatout array.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstatesz,numelt,numat,nstep,kucd,iprestress
      integer infiel(6,numelt),infmat(3,numat),infmatmod(5,nmatmodmax)
      integer ismatmod(nstatesmax,nmatmodmax),infetype(4,netypes)
      integer istatout(2,nstatesmax)
      character ucdroot*(*)
      double precision delt
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
c
c...  local constants
c
      integer ival(48)
      data ival/48*1/
      character sout(6,12)*19
      data sout/"Stress-xx,Pa","Stress-yy,Pa","Stress-zz,Pa",
     &          "Stress-xy,Pa","Stress-yz,Pa","Stress-xz,Pa",
     &          "Strain-xx,ND","Strain-yy,ND","Strain-zz,ND",
     &          "Strain-xy,ND","Strain-yz,ND","Strain-xz,ND",
     &          "Vstrain-xx,ND","Vstrain-yy,ND","Vstrain-zz,ND",
     &          "Vstrain-xy,ND","Vstrain-yz,ND","Vstrain-xz,ND",
     &          "Pstrain-xx,ND","Pstrain-yy,ND","Pstrain-zz,ND",
     &          "Pstrain-xy,ND","Pstrain-yz,ND","Pstrain-xz,ND",
     &          "Stress-rate-xx,Pa-s","Stress-rate-yy,Pa-s",
     &          "Stress-rate-zz,Pa-s","Stress-rate-xy,Pa-s",
     &          "Stress-rate-yz,Pa-s","Stress-rate-xz,Pa-s",
     &          "Strain-rate-xx,1/s","Strain-rate-yy,1/s",
     &          "Strain-rate-zz,1/s","Strain-rate-xy,1/s",
     &          "Strain-rate-yz,1/s","Strain-rate-xz,1/s",
     &          "Vstrain-rate-xx,1/s","Vstrain-rate-yy,1/s",
     &          "Vstrain-rate-zz,1/s","Vstrain-rate-xy,1/s",
     &          "Vstrain-rate-yz,1/s","Vstrain-rate-xz,1/s",
     &          "Pstrain-rate-xx,1/s","Pstrain-rate-yy,1/s",
     &          "Pstrain-rate-zz,1/s","Pstrain-rate-xy,1/s",
     &          "Pstrain-rate-yz,1/s","Pstrain-rate-xz,1/s",
     &          "Stress-incr-xx,Pa","Stress-incr-yy,Pa",
     &          "Stress-incr-zz,Pa","Stress-incr-xy,Pa",
     &          "Stress-incr-yz,Pa","Stress-incr-xz,Pa",
     &          "Strain-incr-xx,ND","Strain-incr-yy,ND",
     &          "Strain-incr-zz,ND","Strain-incr-xy,ND",
     &          "Strain-incr-yz,ND","Strain-incr-xz,ND",
     &          "Vstrain-incr-xx,ND","Vstrain-incr-yy,ND",
     &          "Vstrain-incr-zz,ND","Vstrain-incr-xy,ND",
     &          "Vstrain-incr-yz,ND","Vstrain-incr-xz,ND",
     &          "Pstrain-incr-xx,ND","Pstrain-incr-yy,ND",
     &          "Pstrain-incr-zz,ND","Pstrain-incr-xy,ND",
     &          "Pstrain-incr-yz,ND","Pstrain-incr-xz,ND"/
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer matmodpt(nstatesmax,nmatmodmax)
      integer i1,i2,nnvals,i,j,ist,k,ngpts,iel,imat,ietype,ngauss
      integer matmodel,nstate,inds,indstate,l,indstateg
      character filenm*200,cstep*5
      double precision stmp(48,ngaussmax),tmult
c
cdebug      write(6,*) "Hello from write_ucd_gauss_vals!"
c
      i1=nnblnk(ucdroot)
      i2=nchar(ucdroot)
      if(iprestress.eq.izero) then
        write(cstep,"(i5.5)") nstep
      else
        cstep="prest"
      end if
c
c...  determine how many values will be written
c
      nnvals=izero
      do i=1,nstatesmax
        if(istatout(1,i).ne.izero) nnvals=nnvals+nstr
        if(istatout(2,i).ne.izero) nnvals=nnvals+nstr
      end do
c
c...  define pointer array for location of state variables for each
c     material model.
c
      do i=1,nmatmodmax
        matmodpt(1,i)=izero
        do j=2,nstatesmax
          matmodpt(j,i)=matmodpt(j-1,i)+ismatmod(j-1,i)
        end do
      end do
c
c...  write mesh info
c
      filenm=ucdroot(i1:i2)//".gmesh.time."//cstep//".inp"
      open(kucd,file=filenm,status="new")
      write(kucd,"(50i5)") nnvals,(ival(i),i=1,nnvals)
      do j=1,2
        do i=1,nstatesmax
          ist=izero
          if(istatout(j,i).ne.izero) ist=i
          if(j.eq.itwo.and.istatout(j,i).eq.ione) ist=i+ifour
          if(j.eq.itwo.and.istatout(j,i).eq.itwo) ist=i+ieight
          if(ist.ne.izero) then
            do k=1,nstr
              write(kucd,"(a19)") sout(k,ist)
            end do
          end if
        end do
      end do
c
c...  write state variables at Gauss points
c
      ngpts=izero
      do iel=1,numelt
        imat=infiel(2,iel)
        ietype=infiel(3,iel)
        ngauss=infetype(1,ietype)
        matmodel=infmat(1,imat)
        nstate=infmatmod(2,matmodel)
        inds=1
        do i=1,nstatesmax
          if(istatout(1,i).ne.izero) then
            indstate=infiel(5,iel)+matmodpt(i,matmodel)
            do l=1,ngauss
              indstateg=indstate+(l-1)*nstate
              if(ismatmod(i,matmodel).eq.izero) then
                call fill(stmp(inds,l),zero,nstr)
              else
                call dcopy(nstr,state(1,indstateg),ione,stmp(inds,l),
     &           ione)
              end if
            end do
            inds=inds+nstr
          end if
        end do
        do i=1,nstatesmax
          if(istatout(2,i).ne.izero) then
            tmult=one
            if(istatout(2,i).eq.ione.and.delt.gt.zero) tmult=one/delt
            indstate=infiel(5,iel)+matmodpt(i,matmodel)
            do l=1,ngauss
              indstateg=indstate+(l-1)*nstate
              call fill(stmp(inds,l),zero,nstr)
              if(ismatmod(i,matmodel).ne.izero) call daxpy(nstr,tmult,
     &         dstate(1,indstateg),ione,stmp(inds,l),ione)
            end do
            inds=inds+nstr
          end if
        end do
        do l=1,ngauss
          ngpts=ngpts+ione
          write(kucd,"(i7,48(2x,1pe15.8))") ngpts,(stmp(j,l),
     &     j=1,inds-ione)
        end do
      end do
      close(kucd)
      return
      end
c
c version
c $Id: write_ucd_gauss_vals.f,v 1.3 2005/01/19 20:34:25 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
