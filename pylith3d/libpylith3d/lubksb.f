      SUBROUTINE lubksb(a,n,np,indx,b)
c
c...  routine to perform backsubstitution on an LU-decomposed matrix.
c     Adapted from Numerical Recipes.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n,np
      INTEGER indx(n)
      double precision a(np,np),b(n)
c
c...  local variables
c
      INTEGER i,ii,j,ll
      double precision sum
c
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.0d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
