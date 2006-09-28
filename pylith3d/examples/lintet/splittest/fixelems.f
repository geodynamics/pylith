      program fixelems
c
c...  quick and dirty program to replace element node numbers with
c     new numbers based on a given key.
c
      implicit none
      integer nen,nel,nnode
      parameter(nen=4,nel=41,nnode=25)
      integer n,ntype,mat,inf,ien(nen)
      integer key(nnode)
      data key/11, 10, 25, 19, 23, 7, 18, 1, 24, 14,
     &         20, 22, 13, 16,  5, 21, 6, 2,  4, 17,
     &         12,  3, 15,  8,  9/
      character ifile*20,ofile*21
      data ifile/"element-ordering.txt"/
      data ofile/"element-reordered.txt"/
      integer i,j
c
      open(10,file=ifile,status="old")
      open(11,file=ofile,status="new")
      do i=1,nel
        read(10,*) n,ntype,mat,inf,(ien(j),j=1,nen)
        do j=1,nen
          ien(j)=key(ien(j))
        end do
        write(11,"(8i5)") n,ntype,mat,inf,(ien(j),j=1,nen)
      end do
      close(10)
      close(11)
      stop
      end
