      SUBROUTINE ludcmp(a,n,np,indx,d)
c
c...  subroutine to perform an LU decomposition on a matrix.
c     Adapted from Numerical Recipes.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      integer NMAX
      double precision TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
c
c...  subroutine arguments
c
      integer n,np
      INTEGER indx(n)
      double precision d,a(np,np)
c
c...  local variables
c
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d=1.0d0
      do 12 i=1,n
        aamax=0.0d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.0d0) pause 'singular matrix in ludcmp'
        vv(i)=1.0d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.0d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.0d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.0d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
