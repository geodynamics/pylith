c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine write_ucd_node_vals(d,deld,tfault,dfault,nfault,numfn,
     & dx,deldx,idslp,numsn,deltp,nstep,numnp,kucd,iucd,ucdroot,
     & iprestress)
c
c...  Specialized routine to output displacement info for SCEC
c     benchmarks.
c     This routine creates the nodal value portion of the UCD file for
c     each time step, including the header information.  Note that a
c     complete UCD file is formed by concatenating the mesh output from
c     write_ucdmesh with the file created by this routine.
c     At present, displacements and velocities are written.  In the near
c     future, the variable idispout should be used to determine whether
c     to output displacements, displacement increments, and/or
c     velocities.
c     Note also that the split node output is not really a UCD file,
c     since it contains both node and element information.
c     Either binary or ascii UCD files are written, depending on the
c     value of the iucd parameter.
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
      integer numfn,numsn,nstep,numnp,kucd,iucd,iprestress
      integer nfault(3,numfn),idslp(numsn)
      double precision deltp
      double precision d(ndof,numnp),deld(ndof,numnp)
      double precision tfault(ndof,numfn),dfault(ndof,numfn)
      double precision dx(ndof,numnp),deldx(ndof,numnp)
      character ucdroot*(*)
c
c...  local constants
c
      integer nnvals
      parameter(nnvals=6)
      integer ival(nnvals)
      data ival/1,1,1,1,1,1/
      character dout(nnvals)*9
      data dout/"X-Displ,m","Y-Displ,m","Z-Displ,m",
     &          "X-Vel,m/s","Y-Vel,m/s","Z-Vel,m/s"/
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer i,j,i1,i2,inode
      integer ibyte,intlen,floatlen,len
      double precision rmult
      double precision vmin(6),vmax(6)
      character filenm*200,cstep*5
      character nlabels*1024,nunits*1024
c
cdebug      write(6,*) "Hello from write_ucd_node_vals!"
c
      if(iucd.eq.izero) return
c
c...  set up labels and byte-counting info for binary output
c
      nlabels="X-Displacement.Y-Displacement.Z-Displacement."//
     &        "X-Velocity.Y-Velocity.Z-Velocity."
      nlabels(1024:1024)="0"
      nunits="meters.meters.meters."//
     &       "meters/second.meters/second.meters/second."
      nunits(1024:1024)="0"
      ibyte=ione
      len=ione
      intlen=ifour
      floatlen=ifour
c
c...  determine filename
c
      if(deltp.ne.zero) then
        rmult=one/deltp
      else
        rmult=one
      end if
      i1=nnblnk(ucdroot)
      i2=nchar(ucdroot)
      if(iprestress.eq.izero) then
        write(cstep,"(i5.5)") nstep
      else
        cstep="prest"
      end if
      filenm=ucdroot(i1:i2)//".mesh.time."//cstep//".inp"
c
c...  ascii UCD output
c
      if(iucd.eq.ione) then
c
c...  write mesh info
c
        open(kucd,file=filenm,status="new")
        write(kucd,"(7i7)") nnvals,(ival(i),i=1,nnvals)
        do i=1,nnvals
          write(kucd,"(a9)") dout(i)
        end do
c
c...  write nodal displacements
c
        do i=1,numnp
          write(kucd,"(i7,6(2x,1pe15.8))") i,(d(j,i),j=1,ndof),
     &     (rmult*deld(j,i),j=1,ndof)
        end do
        close(kucd)
c
c...  write split node displacements, if there are any
c
        if(numfn.ne.0) then
          filenm=ucdroot(i1:i2)//".mesh.split.time."//cstep//".inp"
          open(kucd,file=filenm,status="new")
          write(kucd,"(7i7)") nnvals,(ival(i),i=1,nnvals)
          do i=1,nnvals
            write(kucd,"(a9)") dout(i)
          end do
          do i=1,numfn
            write(kucd,"(2i7,6(2x,1pe15.8))") nfault(1,i),nfault(2,i),
     &       (tfault(j,i),j=1,ndof),(rmult*dfault(j,i),j=1,ndof)
          end do
          close(kucd)
        end if
c
c...  write slippery node displacements, if there are any
c
        if(numsn.ne.0) then
          filenm=ucdroot(i1:i2)//".mesh.slip.time."//cstep//".inp"
          open(kucd,file=filenm,status="new")
          write(kucd,"(7i7)") nnvals,(ival(i),i=1,nnvals)
          do i=1,nnvals
            write(kucd,"(a9)") dout(i)
          end do
          do i=1,numsn
            inode=idslp(i)
            write(kucd,"(i7,6(2x,1pe15.8))") idslp(inode),
     &       (dx(j,inode),j=1,ndof),(rmult*deldx(j,inode),j=1,ndof)
          end do
          close(kucd)
        end if
c
c...  binary UCD output
c
      else if(iucd.eq.itwo) then
        open(kucd,file=filenm,status="new",access="direct",recl=len,
     &   form="unformatted")
c
c...  determin min/max of data values
c
        do i=1,6
          vmin(i)=1.0d30
          vmax(i)=-1.0d30
        end do
        do i=1,numnp
          vmin(1)=min(vmin(1),d(1,i))
          vmin(2)=min(vmin(2),d(2,i))
          vmin(3)=min(vmin(3),d(3,i))
          vmin(4)=min(vmin(4),rmult*deld(1,i))
          vmin(5)=min(vmin(5),rmult*deld(2,i))
          vmin(6)=min(vmin(6),rmult*deld(3,i))
          vmax(1)=max(vmax(1),d(1,i))
          vmax(2)=max(vmax(2),d(2,i))
          vmax(3)=max(vmax(3),d(3,i))
          vmax(4)=max(vmax(4),rmult*deld(1,i))
          vmax(5)=max(vmax(5),rmult*deld(2,i))
          vmax(6)=max(vmax(6),rmult*deld(3,i))
        end do
c
c...  write header info
c
        write(kucd,rec=ibyte) nlabels,nunits
        ibyte=ibyte+2048
        write(kucd,rec=ibyte) nnvals,(ival(i),i=1,nnvals)
        ibyte=ibyte+intlen*(1+nnvals)
        write(kucd,rec=ibyte) (real(vmin(i)),i=1,nnvals)
        ibyte=ibyte+floatlen*nnvals
        write(kucd,rec=ibyte) (real(vmax(i)),i=1,nnvals)
        ibyte=ibyte+floatlen*nnvals
c
c...  write node data
c
        do j=1,ndof
          write(kucd,rec=ibyte) (real(d(j,i)),i=1,numnp)
          ibyte=ibyte+numnp*floatlen
        end do
        do j=1,ndof
          write(kucd,rec=ibyte) (real(rmult*deld(j,i)),i=1,numnp)
          ibyte=ibyte+numnp*floatlen
        end do
c
c...  write terminating values and close file
c
        write(kucd,rec=ibyte) (real(vmax(i)),i=1,nnvals)
        close(kucd)
c
c...  NOTE:  for now, binary UCD files are not written for split or
c     slippery nodes, since there does not appear to be a use for them.
c
      end if
      return
      end
c
c version
c $Id: write_ucd_node_vals.f,v 1.1 2005/08/05 19:58:09 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
