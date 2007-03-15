      SUBROUTINE newt(x,n,rpar,nrpar,ipar,nipar,funcv,jaccmp,check,
     & ierr,errstrng)
c
c...  routine to find zero of a function using Newton's method with
c     line search.
c     Adapted from Numerical Recipes.
c**   Note that in this version, the Jacobian is assumed to always
c     be symmetric and stored in packed triangular form, with dimensions
c     of 6x6.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
      include "ndimens.inc"
      integer NP,MAXITS
      double precision TOLF,TOLMIN,TOLX,STPMX
      parameter(NP=6,MAXITS=100,TOLF=1.0d-8,TOLMIN=1.0d-10,TOLX=1.0d-12,
     & STPMX=100.0d0)
c
c...  subroutine arguments
c
      integer n,nrpar,nipar,ierr
      integer ipar(nipar)
      double precision x(n),rpar(nrpar)
      logical check
      character errstrng*(*)
c
c...  external routines
c
      double precision ddot,dnrm2
      external funcv,jaccmp,ddot,dnrm2
c
c...  local variables
c
      double precision fvec(NP)
CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp
      INTEGER i,its,j,indx(NP)
      double precision d,den,f,fold,stpmax,sum,temp,test,fjac(nddmat)
      double precision g(NP),p(NP),xold(NP),fmin
      external fmin
c
      f=fmin(x,n,fvec,rpar,nrpar,ipar,nipar,funcv,ierr,errstrng)
      if(ierr.ne.izero) return
      test=zero
      do i=1,n
        test=max(test,abs(fvec(i)))
      end do
      if(test.lt.0.01d0*TOLF)return
      sum=dnrm2(n,x,ione)
      stpmax=STPMX*max(sum,dble(n))
      do its=1,MAXITS
        call jaccmp(x,n,fvec,nddmat,fjac,rpar,nrpar,ipar,nipar,
     &   ierr,errstrng)
        call dspmv("u",n,one,fjac,fvec,ione,zero,g,ione)
        call dcopy(n,x,ione,xold,ione)
        fold=f
        do i=1,n
          p(i)=-fvec(i)
        end do
        call dppsv('u',n,1,fjac,p,n,ierr)
        call lnsrch(n,xold,fold,g,p,x,f,fvec,stpmax,check,fmin,
     &   rpar,nrpar,ipar,nipar,funcv,ierr,errstrng)
        if(ierr.ne.izero) return
        test=zero
        do i=1,n
          test=max(test,abs(fvec(i)))
        end do
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          test=zero
          den=max(f,half*dble(n))
          do i=1,n
            temp=abs(g(i))*max(abs(x(i)),one)/den
            test=max(test,temp)
          end do
          if(test.lt.TOLMIN)then
            check=.true.
          else
            check=.false.
          endif
          return
        endif
        test=zero
        do i=1,n
          temp=(abs(x(i)-xold(i)))/max(abs(x(i)),one)
          test=max(test,temp)
        end do
        if(test.lt.TOLX)return
      end do
      ierr=120
      errstrng="newt"
      return
      END
