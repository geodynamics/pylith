      program lt2qt
c
c...  Program to take nodal coordinates and connectivities for a linear
c     tet mesh and convert it to a quadratic (10-node) tet mesh, by
c     inserting new nodes between the old ones.
c
      implicit none
      integer maxnodes,maxelmts,nsd,nenl,nenq,ietypeq
      parameter(maxnodes=500000,maxelmts=500000,nsd=3,nenl=4,nenq=10,
     & ietypeq=10)
      double precision eps
      parameter(eps=1.0d-8)
      integer ien(nenq)
      double precision xi(nsd,maxnodes),xn(nsd,maxnodes)
c
      intrinsic iargc
      integer kni,kno,kci,kco,nargs,iargc,i,j,ndi,ndn,n,numel,ietype
      integer mat,inf,k,l,numnp,node
      integer ind(2,6)
      double precision xtmp(3),dx
      character nifile*100,cifile*100,nofile*100,cofile*100,filenm*100
      logical new
      data ind/1,2,2,3,3,1,4,1,4,2,4,3/
c
c...  open i/o files
c
      kni=10
      kno=11
      kci=12
      kco=13
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
        else if(index(filenm,'co=').ne.0) then
          j=lnblnk(filenm)
          cofile=filenm(4:j)
        end if
      end do
      open(kni,file=nifile,status="old")
      open(kno,file=nofile,status="new")
      open(kci,file=cifile,status="old")
      open(kco,file=cofile,status="new")
c
c...  read nodal coordinates
c
      ndi=0
      ndn=0
      do i=1,maxnodes
        read(kni,*,end=10) n,(xi(j,i),j=1,nsd)
        write(kno,'(i7,3(2x,1pe15.8))') n,(xi(j,i),j=1,nsd)
        ndi=ndi+1
      end do
10    continue
      numnp=ndi
c
c...  read connectivities and insert new nodes if they don't already
c     exist
c
      numel=0
      do i=1,maxelmts
        read(kci,*,end=20) n,ietype,mat,inf,(ien(j),j=1,nenl)
        numel=numel+1
        do j=1,6
          do k=1,nsd
            xtmp(k)=0.5d0*(xi(k,ien(ind(1,j)))+xi(k,ien(ind(2,j))))
          end do
c
c...  see if new node has already been produced
c
          new=.true.
          do l=1,ndn
            dx=0.0d0
            do k=1,nsd
              dx=dx+abs(xtmp(k)-xn(k,l))
            end do
            if(dx.lt.eps) then
              new=.false.
              node=l
              go to 30
            end if
          end do
30        continue
          if(new) then
            ndn=ndn+1
            numnp=numnp+1
            ien(j+4)=numnp
            do k=1,nsd
              xn(k,ndn)=xtmp(k)
            end do
          else
            ien(j+4)=ndi+l
          end if
        end do
        write(kco,'(24i7)') i,ietypeq,mat,inf,(ien(j),j=1,nenq)
      end do
20    continue
c
c...  output additional nodal coordinates
c
      do i=1,ndn
        write(kno,'(i7,3(2x,1pe15.8))') i+ndi,(xn(j,i),j=1,nsd)
      end do
c
      stop
      end
