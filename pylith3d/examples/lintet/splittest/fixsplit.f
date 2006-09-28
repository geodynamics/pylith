      program fixelems
c
c...  quick and dirty program to replace element node numbers with
c     new numbers based on a given key.
c
      implicit none
      integer nen,nel,nnode,nsplit
      parameter(nen=4,nel=41,nnode=25,nsplit=76)
      integer elem,node,hist
      double precision split(3)
      integer nkey(nnode)
      data nkey/11, 10, 25, 19, 23, 7, 18, 1, 24, 14,
     &         20, 22, 13, 16,  5, 21, 6, 2,  4, 17,
     &         12,  3, 15,  8,  9/
      integer ekey(nel)
      data ekey/ 1,  3,  8, 10, 12, 15, 17, 11, 18, 14,
     &          20, 24, 19, 16,  2,  9, 28, 26,  5, 22,
     &           6, 33,  7, 36, 25, 32, 31, 38, 29,  4,
     &          30, 34, 35, 40, 27, 37, 39, 13, 41, 21,
     &          23/

      character ifile*20,ofile*21
      data ifile/"splittest.split"/
      data ofile/"splittest-reordered.split"/
      integer i,j
c
      open(10,file=ifile,status="old")
      open(11,file=ofile,status="new")
      do i=1,nsplit
        read(10,*) elem,node,hist,(split(j),j=1,3)
        elem=ekey(elem)
        node=nkey(node)
        write(11,"(3i5,3(2x,1pe15.8))") elem,node,hist,(split(j),j=1,3)
      end do
      close(10)
      close(11)
      stop
      end
