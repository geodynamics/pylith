      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,fvec,stpmax,check,func,rpar,
     & nrpar,ipar,nipar,funcv,ierr,errstrng)
c
c...  routine to perform a line search along the direction given by p.
c     Adapted from Numerical Recipes.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
      integer NP
      double precision ALF,TOLX
      PARAMETER (ALF=1.d-8,TOLX=1.d-12,NP=6)
c
c...  subroutine arguments
c
      integer n,nrpar,nipar,ierr
      integer ipar(nipar)
      double precision xold(n),fold,g(n),p(n),x(n),f,fvec(n),stpmax,func
      double precision rpar(nrpar)
      external func,funcv
      logical check
      character errstrng*(*)
CU    USES func
c
c...  external routines
c
      double precision dnrm2,ddot
      external dnrm2,ddot
c
c...  local variables
c
      INTEGER i
      double precision a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2
      double precision slope,sum,temp,test,tmplam,da
c
      check=.false.
      sum=dnrm2(n,p,ione)
      if(sum.gt.stpmax)then
        da=stpmax/sum
        call dscal(n,da,p,ione)
      endif
      slope=ddot(n,g,ione,p,ione)
      test=zero
      do i=1,n
        temp=abs(p(i))/max(abs(xold(i)),one)
        test=max(temp,test)
      end do
      alamin=TOLX/test
      alam=one
1     continue
        do i=1,n
          x(i)=xold(i)+alam*p(i)
        end do
        f=func(x,n,fvec,rpar,nrpar,ipar,nipar,funcv,ierr,errstrng)
        if(alam.lt.alamin)then
          call dcopy(n,xold,ione,x,ione)
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.one)then
            tmplam=-slope/(two*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam*alam-rhs2/(alam2*alam2))/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.zero)then
              tmplam=-slope/(two*b)
            else
              disc=b*b-three*a*slope
              tmplam=(-b+sqrt(disc))/(three*a)
            endif
            if(tmplam.gt.half*alam)tmplam=half*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,0.1d0*alam)
      goto 1
      END
