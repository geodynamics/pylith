      program addflt
c
c...  Quick and dirty program to take lithomop input/output files and
c     create alternate output that has separate nodes on either side of
c     faults.
c
      implicit none
      integer maxnodes,maxelmts,nsd,nen
      parameter(maxnodes=600000,maxelmts=600000,nenmax=20,netypes=10,
     & nsd=3)
      double precision eps
      parameter(eps=1.0d-8)
      integer ien(nenmax,maxelmts),mat(maxelmts)
      double precision x(nsd,maxnodes)
c
      integer kti,kto,kr,kw
      integer nchar,nnblnk
      external nchar,nnblnk
      integer ifb,ife,ilb,ile
      integer numnp,numel,nnout,neout,ngout
      integer n,i,j
      integer nen(netypes)
      data nen/8,7,6,5,4,20,18,15,13,10/
      double precision dx,dy,dz
      character fileroot*100,filenm*100,eltype(5)*5
      data eltype/"hex","wrick","prism","pyr","tet"/
      character line*500
c
c...  get file root and fault offset values
c
      kti=5
      kto=6
      kr=10
      kw=11
      write(kto,*) "Enter root filename:"
      read(kti,"(a80)") fileroot
      write(kto,*) "Enter x, y, and z offsets for new fault nodes:"
      read(kti,*) dx,dy,dz
      ifb=nnblnk(fileroot)
      ife=nchar(fileroot)
c
c...  read initial nodal coordinates
c
      filenm=fileroot(ifb:ife)//".mesh.inp"
      open(kr,file=filenm,status="old")
      read(kr,*) numnp,numel,nnout,neout,ngout
      do i=1,numnp
        read(kr,*) n,(x(j,i),j=1,nsd)
      end do
c
c...  read initial connectivities
c
      read(kr,"(a500)") line
      do i=1,5
        ilb=index(line,eltype(i))
        if(ilb.ne.0) then
          ietype=i+5
          go to 10
        end if
      end do
10    continue
      ile=ilb+4
      read(line,*) n,mat(1)
      read(line(ile+1:500),*,err=20) (ien(j,1),j=1,nen(ietype))
      go to 30
20    continue
      ietype=ietype-5
30    continue
      do i=2,numel
        read(kr,"(a500)") line
        read(line,*) n,mat(i)
        read(line(ile+1:500),*) (ien(j,i),j=1,nen(ietype))
      end do
      close(kr)
c
c...  read any split node descriptions
c
      nf=0
      numfn=0
      nnp=numnp
      split=.false.
      filenm=fileroot(ifb:ife)//".split"
      open(kr,file=filenm,status="old",err=40)
      split=.true.
      call pskip(kr)
50    continue
        read(kr,*,end=40) (nfault(j,numfn+1),j=1,3),(fault(j),j=1,ndof)
        numfn=numfn+1
        used=.false.
        do i=1,nf
          if(nfault(2,numfn).eq.ifs(1,i)) used=.true.
        end do
        if(.not.used) then
          nf=nf+1
          nnp=nnp+1
          ifs(1,nf)=nfault(2,numfn)
          ifs(2,nf)=nnp
          x(1,nnp)=x(1,ifs(1,nf))+dx
          x(2,nnp)=x(2,ifs(1,nf))+dy
          x(3,nnp)=x(3,ifs(1,nf))+dz
        end if
        sgn=1.0d0
        if(fault(1).ne.0.0d0) then
          sgn=sign(1.0d0,fault(1))
        else if(fault(2).ne.0.0d0) then
          sgn=sign(1.0d0,fault(2))
        else if(fault(3).ne.0.0d0) then
          sgn=sign(1.0d0,fault(3))
        end if
        go to 50
      close(kr)
40    continue
c
c...  read any slippery node descriptions
c
      ns=0
      numslp=0
      slip=.false.
      filenm=fileroot(ifb:ife)//".slip"
      open(kr,file=filenm,status="old",err=60)
      slip=.true.
      call pskip(kr)
70    continue
        read(kr,*,end=60) (nslip(j,numslp+1),j=1,5)
        numslp=numslp+1
        used=.false.
        do i=1+nf,ns
          if(slip(2,numslp).eq.ifs(1,i)) used=.true.
        end do
        if(.not.used) then
          ns=ns+1
          nnp=nnp+1
          ifs(1,nf+ns)=nslip(2,numslp)
          ifs(2,nf+ns)=nnp
          x(1,nnp)=x(1,ifs(1,nf+ns))+dx
          x(2,nnp)=x(2,ifs(1,nf+ns))+dy
          x(3,nnp)=x(3,ifs(1,nf+ns))+dz
        end if
        go to 70
      close(kr)
60    continue
c
c...  loop over time steps looking for UCD files
c
      ctime="prest"
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
