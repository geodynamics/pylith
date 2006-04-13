      program tetcmp
c
c...  program to compute determinant and equivalent tetrahedral volume
c     for a given 4x4 matrix.
c
      implicit none
      double precision a(4,4),cof(3,3)
      integer kti,kto
      data kti,kto/5,6/
      integer i,j,ii
      double precision det3,onem,sgn,det,vol
      external det3
      integer ind(3,4)
      data ind/2,3,4,1,3,4,1,2,4,1,2,3/
      write(kto,*) "Enter matrix row by row:"
      do i=1,4
        read(kti,*) (a(i,j),j=1,4)
      end do
      onem=-1.0d0
      sgn=onem
      det=0.0d0
      do i=1,4
        sgn=sgn*onem
        do ii=1,3
          do j=2,4
            cof(j-1,ii)=a(j,ind(ii,i))
          end do
        end do
        det=det+sgn*a(1,i)*det3(cof)
      end do
      vol=det/6.0d0
      write(kto,700) det,vol
700   format("Determinant:  ",1pe15.8,/,
     &       "Volume:       ",1pe15.8)
      stop
      end
c
c
      function det3(x)
c
c...  function to compute determinant of a 3x3 matrix
c
      double precision det3
      double precision x(3,3)
      det3=x(1,1)*x(2,2)*x(3,3)-x(1,1)*x(2,3)*x(3,2)+
     &     x(1,2)*x(2,3)*x(3,1)-x(1,2)*x(2,1)*x(3,3)+
     &     x(1,3)*x(2,1)*x(3,2)-x(1,3)*x(2,2)*x(3,1)
      return
      end
