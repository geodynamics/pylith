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
      subroutine presurql(pres,pdir,xl,p,ien,iside,infin,n,nen,
     & idout,kto,kw)
c
c...subroutine to compute local load vector for traction BC
c
c     note defintions of side numbers:
c
c      iside = 1-4 sides of brick containing edge i,i+1
c      iside = 5  front face (nodes 1,2,3,4)
c      iside = 6  back face  (nodes 5,6,7,8)
c
c     pres is the traction value, and is specified at each node on the
c     face.
c     pdir gives the direction cosines for the traction direction if
c     iside < 0.  If iside is positive, the traction is assumed to be a
c     pressure load, and the direction at each node is taken to be the
c     plane normal computed for nodes i-1, i, and i+1.  Positive values
c     of pres indicate a compressive stress when iside > 0.
c
c**  Note that this routine currently assumed 2x2 quadrature over an
c    element face, regardless of the integration order being used in
c    the rest of the code.  This may be altered in the future.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
c
c...  subroutine arguments
c
      integer iside,infin,n,nen,idout,kto,kw
      integer ien(nen)
      double precision pres(nen/2),pdir(npdir),xl(nsd,nen),p(ndof,nen)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
      integer iefc(4,6)
      data iefc/1,5,6,2,
     &          2,6,7,3,
     &          3,7,8,4,
     &          4,8,5,1,
     &          1,2,3,4,
     &          5,8,7,6/
c
      double precision rg(8),sg(8),tg(8),rf(6),sf(6),tf(6),ri(6),si(6)
      double precision ti(6)
      data rg/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0/
      data sg/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0/
      data tg/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0/
      data rf/ 0d0, 1d0, 0d0,-1d0, 0d0, 0d0/
      data sf/-1d0, 0d0, 1d0, 0d0, 0d0, 0d0/
      data tf/ 0d0, 0d0, 0d0, 0d0, 1d0,-1d0/
      data ri/ 1d0, 0d0, 1d0, 0d0, 1d0, 1d0/
      data si/ 0d0, 1d0, 0d0, 1d0, 1d0, 1d0/
      data ti/ 1d0, 1d0, 1d0, 1d0, 0d0, 0d0/
c
c...  intrinsic functions
c
      intrinsic sqrt,abs
c
c...  user-defined functions
c
      double precision dnrm2
      external dnrm2
c
c...  local variables
c
      integer iiside,if1,if2,i,im,ip,nm,ne,np,j,iopt,ipnt,ii,k,inode,l
      double precision g,vmag,ds,det
      double precision sh(4,8),ptmp(3,4),v1(3),v2(3),vc(3),rgs(3)
c
cdebug      write(6,*) "Hello from presurql_f!"
c
      g=root3i
c
c...  Compute direction cosines for pressure BC and corresponding
c     pressure loads, if requested.  Otherwise, compute the component
c     of traction in each direction, as specified by pdir.
c
      iiside=abs(iside)
      if1=1
      if2=3
      if(iiside.eq.2.or.iiside.eq.4) then
        if1=2
        if2=3
      end if
      if(iiside.eq.5.or.iiside.eq.6) then
        if1=1
        if2=2
      end if
      call fill(p,zero,ndof*nen)
      call fill(ptmp,zero,ndof*nen/2)
      if(iside.gt.0) then
        do i=1,nen/2
          im=i-1
          if(im.eq.0) im=4
          ip=i+1
          if(ip.eq.5) ip=1
          nm=iefc(im,iiside)
          ne=iefc(i,iiside)
          np=iefc(ip,iiside)
          do j=1,3
            v1(j)=xl(j,nm)-xl(j,ne)
            v2(j)=xl(j,np)-xl(j,ne)
          end do
          call cross(vc,v1,v2)
          vmag=pres(i)/dnrm2(nsd,vc,ione)
          call daxpy(nsd,vmag,vc,ione,ptmp(1,i),ione)
        end do
      else
        do i=1,nen/2
          call daxpy(nsd,pres(i),pdir,ione,ptmp(1,i),ione)
        end do
      end if
c
c...  perform numerical integration over element face to obtain
c     equivalent nodal loads
c
      iopt=1
      do l=1,nen/2
        ipnt=iefc(l,iiside)
        rgs(1)=rf(iiside)+ri(iiside)*g*rg(ipnt)
        rgs(2)=sf(iiside)+si(iiside)*g*sg(ipnt)
        rgs(3)=tf(iiside)+ti(iiside)*g*tg(ipnt)
clater        call shapql(rgs,xl,det,sh,ien,nen,infin,iopt,n,idout,kto,kw)
c
c...  compute surface element for use in integration
c
        do i=1,3
          v1(i)=zero
          v2(i)=zero
          do j=1,nen/2
            ii=iefc(j,iiside)
            v1(i)=v1(i)+sh(if1,ii)*xl(i,ii)
            v2(i)=v2(i)+sh(if2,ii)*xl(i,ii)
          end do
        end do
        call cross(vc,v1,v2)
        ds=dnrm2(nsd,vc,ione)
        do k=1,nen/2
          inode=iefc(k,iiside)
          p(1,inode)=p(1,inode)+ptmp(1,l)*ds*sh(4,inode)
          p(2,inode)=p(2,inode)+ptmp(2,l)*ds*sh(4,inode)
          p(3,inode)=p(3,inode)+ptmp(3,l)*ds*sh(4,inode)
        end do
      end do
      return
      end
c
c version
c $Id: presurql.f,v 1.3 2004/08/12 22:53:10 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
