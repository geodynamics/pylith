      program fixqt
c
c...  Program to take nodal coordinates and connectivities for a
c     quadratic tet mesh and insure that all non-vertex nodes occur
c     at mid-side.
c
      implicit none
      integer maxnodes,maxelmts,nsd,nen
      parameter(maxnodes=500000,maxelmts=500000,nsd=3,nen=10)
      double precision eps
      parameter(eps=1.0d-8)
      integer ien(nen)
      double precision xi(nsd,maxnodes)
c
      intrinsic iargc
      integer kni,kno,kci,nargs,iargc,i,j,numnp,n,numel,ietype
      integer mat,inf,k,node
      integer ind(2,6)
      character nifile*100,cifile*100,nofile*100,filenm*100
      data ind/1,2,2,3,3,1,1,4,2,4,3,4/
c
c...  open i/o files
c
      kni=10
      kno=11
      kci=12
      nargs=iargc()
      do i=1,nargs
        call getarg(i,filenm)
        if(index(filenm,'ni=').ne.0) then
          j=lnblnk(filenm)
          nifile=filenm(4:j)
        else if(index(filenm,'ci=').ne.0) then
          j=lnblnk(filenm)
          cifile=filenm(4:j)
        else if(index(filenm,'no=').ne.0) then
          j=lnblnk(filenm)
          nofile=filenm(4:j)
        end if
      end do
      open(kni,file=nifile,status="old")
      open(kno,file=nofile,status="new")
      open(kci,file=cifile,status="old")
c
c...  read nodal coordinates
c
      numnp=0
      do i=1,maxnodes
        read(kni,*,end=10) n,(xi(j,i),j=1,nsd)
        numnp=numnp+1
      end do
10    continue
c
c...  read connectivities and make sure non-vertex nodes are at mid-side
c
      numel=0
      do i=1,maxelmts
        read(kci,*,end=20) n,ietype,mat,inf,(ien(j),j=1,nen)
        numel=numel+1
        do j=1,6
          node=ien(j+4)
          do k=1,nsd
            xi(k,node)=0.5d0*(xi(k,ien(ind(1,j)))+xi(k,ien(ind(2,j))))
          end do
        end do
      end do
20    continue
c
c...  output nodal coordinates
c
      do i=1,numnp
        write(kno,'(i7,3(2x,1pe15.8))') i,(xi(j,i),j=1,nsd)
      end do
c
      stop
      end
