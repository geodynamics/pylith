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
      subroutine td_strs_cmp_ss(
     & bintern,neq,                                                     ! force
     & x,d,numnp,iddmat,                                                ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,               ! elemfamily
     & nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,               ! elemfamily
     & ielg,                                                            ! elemnt
     & prop,nprop,elas_strs,td_strs,matchg,tminmax,                     ! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & rtimdat,ntimdat,rgiter,                                          ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...program to compute the total stress and strain for the current
c   iteration for a given material model
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,numnp,numslp,numfn,nelfamily,nstate,nstate0,n0states
c
c... the ielg variable below contains the global element number, while
c    the local variable ielf contains the element number within a given
c    family.
c
      integer ielg
      integer nprestrflag,ipstrs,ipauto,nprop,nen,ngauss,nee,numrot,ierr
      integer iddmat(nstr,nstr)
      integer ien(nen,nelfamily),lm(ndof*nen,nelfamily)
      integer lmx(ndof*nen,nelfamily),lmf(nen,nelfamily)
      character errstrng*(*)
      logical matchg
      double precision bintern(neq),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp)
      double precision tfault(ndof,numfn)
      double precision state(nstate,ngauss,nelfamily)
      double precision dstate(nstate,ngauss,nelfamily)
      double precision state0(nstate0,ngauss,n0states)
      double precision dmat(nddmat,ngauss,nelfamily),prop(nprop),tminmax
      double precision gauss(nsd+1,ngauss),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss),skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  external routines
c
      external elas_strs,td_strs,getshape,bmatrix
c
c...  local variables
c
      integer ielf,incstate0,indstate0,l
      double precision tmax
      double precision dl(60),xl(60),scur(162),ee(162),p(60),det(27)
c
cdebug      integer idb,jdb
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_strs_cmp_ss_f!"
c
      incstate0=izero
      indstate0=ione
      if(ipstrs.ne.izero) then
        incstate0=ione
        indstate0=izero
      end if
c
c...  loop over elements in a family
c
      do ielf=1,nelfamily
        indstate0=indstate0+incstate0
c
c...  localize coordinates and displacements
c
        call lcoord(x,xl,ien(1,ielf),nen,numnp)
        call ldisp(dl,d,ien(1,ielf),nen,numnp)
        if(numfn.ne.0) call adfldp(dl,lmf(1,ielf),tfault,nen,numfn)
        if(numslp.ne.0) call addsn(dl,dx,ien(1,ielf),lmx(1,ielf),nen,
     &   numnp)
c
c...  compute strains
c
        call bdeld_ss(xl,dl,sh,shj,ee,det,gauss,ielg,nen,nee,ngauss,
     &   getshape,bmatrix,ierr,errstrng)
        if(ierr.ne.izero) return
c
c...  loop over gauss points, compute stresses, and transfer them into
c     scur
c
        do l=1,ngauss
          call td_strs(state(1,l,ielf),dstate(1,l,ielf),
     &     state0(1,l,indstate0),ee(nstr*(l-1)+1),scur(nstr*(l-1)+1),
     &     dmat(1,l,ielf),prop,rtimdat,rgiter,ntimdat,iddmat,tmax,
     &     nstate,nstate0,nprop,matchg,ierr,errstrng)
          if(ierr.ne.izero) return
          tminmax=min(tmax,tminmax)
        end do
c
c...  compute equivalent nodal loads
c
        call fill(p,zero,nee)
        call eforce(xl,sh,shj,det,gauss,scur,p,ielg,nen,ngauss,getshape,
     &   bmatrix,ierr,errstrng)
        if(ierr.ne.izero) return
        if(numrot.ne.izero) call rpforc(p,skew,ien(1,ielf),numnp,nen)
        call addfor(bintern,p,lm(1,ielf),lmx(1,ielf),neq,nee)
        ielg=ielg+ione
      end do
      return
      end
c
c version
c $Id: td_strs_cmp_ss.f,v 1.10 2005/04/08 00:37:25 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
