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
      subroutine write_ucd_node_vals(d,deld,deltp,nstep,numnp,kucd,
     & ucdroot,iprestress)
c
c...  Specialized routine to output displacement info for SCEC
c     benchmarks.
c     This routine creates the nodal value portion of the UCD file for
c     each time step, including the header information.  Note that a
c     complete UCD file is formed by concatenating the mesh output from
c     write_ucdmesh with the file created by this routine.
c     At present, only the total displacements are written.  In the near
c     future, the variable idispout should be used to determine whether
c     to output displacements, displacement increments, and/or
c     velocities.
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
      integer nstep,numnp,kucd,iprestress
      character ucdroot*(*)
      double precision d(ndof,numnp),deld(ndof,numnp),deltp
c
c...  local constants
c
      integer nnvals
      parameter(nnvals=3)
      integer ival(nnvals)
      data ival/1,1,1/
      character dout(3)*9
      data dout/"X-Displ,m","Y-Displ,m","Z-Displ,m"/
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer i,j,i1,i2
      character filenm*200,cstep*5
c
cdebug      write(6,*) "Hello from write_ucd_node_vals!"
c
      i1=nnblnk(ucdroot)
      i2=nchar(ucdroot)
      if(iprestress.eq.izero) then
        write(cstep,"(i5.5)") nstep
      else
        cstep="prest"
      end if
c
c...  write mesh info
c
      filenm=ucdroot(i1:i2)//".mesh.time."//cstep//".inp"
      open(kucd,file=filenm,status="new")
      write(kucd,"(7i7)") nnvals,(ival(i),i=1,nnvals)
      do i=1,nnvals
        write(kucd,"(a9)") dout(i)
      end do
c
c...  write nodal displacements
c
      do i=1,numnp
        write(kucd,"(i7,3(2x,1pe15.8))") i,(d(j,i),j=1,ndof)
      end do
      close(kucd)
      return
      end
c
c version
c $Id: write_ucd_node_vals.f,v 1.3 2005/01/19 20:28:04 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
