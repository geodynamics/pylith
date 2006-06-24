c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine stresn(x,b,d,dx,tfault,stn,eps,beta,betb,scur,st0,
     & dbeta,dbetb,skew,ien,lm,lmx,lmf,dmat,mat,prop,histry,infin,gauss,
     & rtimdat,stol,iddmat,nen,numel,ndof,nsd,numnp,neq,nee,
     & nstr,ngauss,nppts,ngem,nskdim,nhist,nprop,numat,numfn,numslp,
     & numrot,lastep,nstep,lgdefp,ibbarp,ivisc,iplas,imhist,ipstrs,
     & nprestr,nddmat,ndmat,idebug,idout,kto,kw,fulout)
c
c...program to compute the total stress and strain for the current
c   iteration
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nen,numel,ndof,nsd,numnp,neq,nee,nstr,ngauss,nppts,ngem
      integer nskdim,nhist,nprop,numat,numfn,numslp,numrot,lastep,nstep
      integer lgdefp,ibbarp,ivisc,iplas,imhist,ipstrs,nprestr,nddmat
      integer ndmat,idebug,idout,kto,kw
      integer ien(nen,numel),lm(ndof,nen,numel),lmx(ndof,nen,numel)
      integer lmf(nen,numel),mat(numel),infin(numel)
      double precision stol
      double precision x(nsd,numnp),b(neq),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn),stn(nstr,ngauss,numel)
      double precision eps(nstr,ngauss,numel),beta(nstr,ngauss,numel)
      double precision betb(nstr,ngauss,numel),scur(nstr,ngauss,numel)
      double precision st0(nstr,nppts,numel),dbeta(nstr,ngauss,numel)
      double precision dbetb(nstr,ngauss,numel),skew(nskdim,numnp)
      double precision dmat(nddmat,ngauss,ndmat),prop(nprop,numat)
      double precision histry(nhist,lastep+1),gauss(nsd+1,ngauss)
      logical fulout
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "rtimdat_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
      double precision big
      parameter(big=1.0d30)
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
cdebug      integer idb,jdb,ldb,ndb
      integer ldtmp,io1,npage,n,m,imat,l,ll,nstart,nend,k,nmin
      double precision tmin,emhu,anpwr,rmu,tau,sinv1,steff
      double precision dl(24),xl(24),ee(48),p(24),det(8),wr(24),ptmp(30)
      double precision eet(6),sdev(6)
      logical debug
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
c
cdebug      write(6,*) "Hello from stresn_f!"
cdebug2      write(6,*) nen,numel,ndof,nsd,numnp,neq,nee,nstr,ngauss,nppts
cdebug2      write(6,*) ngem,nskdim,nhist,nprop,numat,numfn,numslp,numrot
cdebug2      write(6,*) lastep,nstep,lgdefp,ibbarp,ivisc,iplas,imhist,ipstrs
cdebug2      write(6,*) nprestr,nddmat,ndmat,idebug,idout,kto,kw,fulout,stol
cdebug2      write(6,*) "x:",((x(idb,ndb),idb=1,nsd),ndb=1,numnp)
cdebug      write(6,*) "b:",(b(idb),idb=1,neq)
cdebug2      write(6,*) "d:",((d(idb,ndb),idb=1,ndof),ndb=1,numnp)
cdebug2      write(6,*) "dx:",((dx(idb,ndb),idb=1,ndof),ndb=1,numnp)
cdebug2      write(6,*) "tfault:",((tfault(idb,ndb),idb=1,ndof),ndb=1,numfn)
cdebug2      write(6,*) "stn:",(((stn(idb,ldb,ndb),idb=1,nstr),ldb=1,ngauss),
cdebug2     & ndb=1,numel)
cdebug2      write(6,*) "eps:",(((eps(idb,ldb,ndb),idb=1,nstr),ldb=1,ngauss),
cdebug2     & ndb=1,numel)
cdebug2      write(6,*) "beta:",(((beta(idb,ldb,ndb),idb=1,nstr),ldb=1,ngauss),
cdebug2     & ndb=1,numel)
cdebug      write(6,*) "betb:",(((betb(idb,ldb,ndb),idb=1,nstr),ldb=1,ngauss),
cdebug     & ndb=1,numel)
cdebug2      write(6,*) "scur:",(((scur(idb,ldb,ndb),idb=1,nstr),ldb=1,ngauss),
cdebug2     & ndb=1,numel)
cdebug      write(6,*) "st0:",(((st0(idb,ldb,ndb),idb=1,nstr),ldb=1,nppts),
cdebug     & ndb=1,numel)
cdebug2      write(6,*) "dbeta:",(((dbeta(idb,ldb,ndb),idb=1,nstr),
cdebug2     & ldb=1,ngauss),ndb=1,numel)
cdebug      write(6,*) "dbetb:",(((dbetb(idb,ldb,ndb),idb=1,nstr),
cdebug     & ldb=1,ngauss),ndb=1,numel)
cdebug      write(6,*) "skew:",((skew(idb,ndb),idb=1,nskdim),ndb=1,numnp)
cdebug2      write(6,*) "dmat:",(((dmat(idb,ldb,ndb),idb=1,nddmat),
cdebug2     & ldb=1,ngauss),ndb=1,ndmat)
cdebug2      write(6,*) "prop:",((prop(idb,ndb),idb=1,nprop),ndb=1,numat)
cdebug      write(6,*) "From stresn, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & ((gauss(idb,ndb),idb=1,nsd+1),ndb=1,ngauss)
cdebug2      write(6,*) "iddmat:",((iddmat(idb,ldb),idb=1,6),ldb=1,6)
c
      ldtmp=lgdefp
      tmin=big
      if(ipstrs.eq.1.and.nstep.eq.0) ldtmp=0
      io1=0
      if(ldtmp.gt.0) io1=1
      debug=(idebug.eq.1).and.(idout.gt.1)
      call fill(b,zero,neq)
      npage=50
      do n=1,numel
        m=mat(n)
        imat=n
        if((ivisc.eq.0).and.(iplas.eq.0)) imat=m
cdebug        write(6,*) "From point 1 in stresn, n, imat:",n,imat
        call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &   idout,kto,kw,imhist)
        emhu=ptmp(4)
        anpwr=ptmp(5)
        rmu=half*ptmp(1)/(one+ptmp(2))
        call lcoord(x,xl,ien(1,n),nen,nsd,numnp)
        call ldupdat(d,dx,tfault,dl,xl,ien(1,n),lmx(1,1,n),lmf(1,n),
     &   ndof,nsd,nen,numnp,numfn,numslp,io1,ldtmp)
cdebug        if(n.lt.20) then
cdebug          write(6,*) "xl:",(xl(idb),idb=1,24)
cdebug          write(6,*) "dl:",(dl(idb),idb=1,24)
cdebug        end if
cdebug        write(6,*) "From point 2 in stresn, n, imat:",n,imat
        call bdeldql(xl,dl,ee,wr,det,gauss,ien(1,n),infin(n),n,nen,nee,
     &   nsd,ndof,nstr,ngauss,ibbarp,idout,kto,kw)
cdebug        write(6,*) "From point 3 in stresn, n, imat:",n,imat
cdebug2        write(6,*) "From stresn after bdeldql:"
cdebug2        write(6,*) "xl:",(xl(idb),idb=1,24)
cdebug2        write(6,*) "dl:",(dl(idb),idb=1,24)
cdebug2        write(6,*) "ee:",(ee(idb),idb=1,48)
cdebug2        write(6,*) "wr:",(wr(idb),idb=1,24)
cdebug2        write(6,*) "det:",(det(idb),idb=1,8)
cdebug        write(6,*) "gauss:",((gauss(jdb,idb),jdb=1,nsd+1),idb=1,ngauss)
cdebug2        write(6,*) "ien:",(ien(jdb,n),jdb=1,nen)
cdebug2        write(6,*) "infin,n,nen,nee,nsd,ndof,nstr,ngauss,ibbarp,idout,"
cdebug2        write(6,*) "kto,kw:"
cdebug2        write(6,*) infin(n),n,nen,nee,nsd,ndof,nstr,ngauss,ibbarp,idout
cdebug2        write(6,*) kto,kw
cdebug        write(6,*) (ee(idb),idb=1,48)
cdebug        if(n.lt.20) write(6,*) "ee:",(ee(idb),idb=1,48)
        do l=1,ngauss
          ll=l
          if(nppts.eq.1) ll=1
cdebug          write(6,*) "From point 4 in stresn, n, l, ll:",n,l,ll
          if(ldtmp.ge.1) then
            call stresld(stn(1,l,n),scur(1,l,n),st0(1,ll,n),eps(1,l,n),
     &       beta(1,l,n),dbeta(1,l,n),betb(1,l,n),dbetb(1,l,n),
     &       dmat(1,l,imat),ee(nstr*(l-1)+1),wr(3*(l-1)+1),ptmp,
     &       rtimdat,stol,iddmat,n,nstr,ndof,nddmat,nprop,
     &       ipstrs,nprestr,nstep,lgdefp,idout,kto,kw,ivisc,iplas)
cdebug2            write(6,*) "Test 1"
          else
            call dcopy(nstr,ee(nstr*(l-1)+1),ione,eps(1,l,n),ione)
cdebug2            write(6,*) "From stresn, n,l,ll:",n,l,ll
            call dcopy(nstr,eps(1,l,n),ione,eet,ione)
cdebug            write(6,*) "From point 5 in stresn, n, l, ll:",n,l,ll
cdebug            write(6,*) "Point 2",n,l,ll
            if(nstep.gt.0.and.((ivisc.eq.1).or.(iplas.eq.1))) then
              nstart=ndof+1
              nend=nstr
              if(ngem.eq.0.or.ngem.eq.1) nstart=ndof+2
              do k=nstart,nend
                eet(k)=half*eet(k)
              end do
cdebug              if(n.lt.20) then
cdebug2                write(6,*) "Before esfcomp, n,l,ll",n,l,ll
cdebug2                write(6,*) "stn:",(stn(idb,l,n),idb=1,nstr)
cdebug2                write(6,*) "scur:",(scur(idb,l,n),idb=1,nstr)
cdebug2                write(6,*) "eet:",(eet(idb),idb=1,nstr)
cdebug2                write(6,*) "beta:",(beta(idb,l,n),idb=1,nstr)
cdebug2                write(6,*) "dbeta:",(dbeta(idb,l,n),idb=1,nstr)
cdebug                write(6,*) "betb:",(betb(idb,l,n),idb=1,nstr)
cdebug2                write(6,*) "ptmp:",(ptmp(idb),idb=1,10)
cdebug                write(6,*) "rtimdat:",(rtimdat(idb),idb=1,3)
cdebug2                write(6,*) stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp
cdebug2                write(6,*) idout,kto,kw,ivisc,iplas
cdebug                write(6,*) "dbetb:",(dbetb(idb,l,n),idb=1,nstr)
cdebug              end if
cdebug              write(6,*) "From point 6 in stresn, n, l, ll:",n,l,ll
              call esfcomp(stn(1,l,n),scur(1,l,n),eet,beta(1,l,n),
     &         dbeta(1,l,n),betb(1,l,n),dbetb(1,l,n),ptmp,
     &         rtimdat,stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,
     &         kto,kw,ivisc,iplas)
cdebug              write(6,*) "From point 7 in stresn, n, l, ll:",n,l,ll
cdebug2              write(6,*) "Right after esfcomp in stresn!"
cdebug2              write(6,*) "stn:",(stn(idb,l,n),idb=1,nstr)
cdebug2              write(6,*) "scur:",(scur(idb,l,n),idb=1,nstr)
cdebug2              write(6,*) "dbeta:",(dbeta(idb,l,n),idb=1,nstr)
cdebug              write(6,*) "Point 4",n,l,ll
            else
              call dcopy(nstr,eps(1,l,n),ione,eet,ione)
              do k=1,ndof
                eet(k)=eps(k,l,n)-third*ptmp(nprop-1)
              end do
	      call dspmv("l",nstr,one,dmat(1,l,imat),eet,ione,zero,
     &         scur(1,l,n),ione)
            end if
            if(nprestr.ne.0) call daxpy(nstr,one,st0(1,ll,n),ione,
     &       scur(1,l,n),ione)
          end if
          call invar(sdev,sinv1,steff,scur(1,l,n))
cdebug          write(6,*) "From point 8 in stresn, n, l, ll:",n,l,ll
          tau=tmin
          if(steff.ne.zero) tau=root3*(emhu/steff)**(anpwr-one)*emhu/rmu
          if(tau.lt.tmin) then
            tmin=tau
            nmin=n
          end if
        end do
c
c...compute equivalent nodal forces corresponding to stresses
c
        call fill(p,zero,nee)
cdebug        write(6,*) "From point 9 in stresn, n, infin(n):",n,infin(n)
        call eforceql(scur(1,1,n),p,xl,ptmp,gauss,ien(1,n),infin(n),n,
     &   nen,nee,nsd,nstr,ngauss,ibbarp,nprop,idout,kto,kw)
cdebug        write(6,*) "From point 10 in stresn, n, infin(n):",n,infin(n)
        if(numrot.ne.0) call rpforc(p,skew,ien(1,n),ndof,numnp,nen,
     &   nskdim)
        if(debug) then
          if((n.eq.1.or.mod(n,npage).eq.0)) write(kw,1000)
          call prntforc(n,p,ien(1,n),nen,ndof,idout,kw)
        end if
        call addfor(b,p,lm(1,1,n),lmx(1,1,n),neq,nee)
cdebug        write(6,*) "From point 11 in stresn, n, infin(n):",n,infin(n)
      end do
      if(tmin.lt.deltp) then
        write(kto,800) nmin,tmin,deltp
        if(idout.gt.1) write(kw,800) nmin,tmin,deltp
      end if
      if(ivisc.eq.1) then
        write(kto,810) tmin,nmin
        if(idout.gt.1) write(kw,810) tmin,nmin
      end if
cdebug      write(6,*) "From end of stresn, b:",(b(idb),idb=1,neq)
 800  format(//,' WARNING!',
     & '  Computed Maxwell time is less than step size for element',i8,
     & '!',/,
     & '  Minimum Maxwell time:  ',1pe15.8,/,
     & '  Time step size:        ',1pe15.8)
 810  format(//,'  Smallest Maxwell time for current iteration:',
     & 1pe15.8,/,
     & '  in element #',i8,//)
1000  format(//' local forces computed from stress field'//)
      return
      end
c
c version
c $Id: stresn.f,v 1.1 2004/07/12 20:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
