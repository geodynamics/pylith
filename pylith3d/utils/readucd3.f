      program readucd3
c
c...  quick and dirty code to translate AVS UCD format to a format
c     suitable for lithomop or pylith.  All boundaries for which BC
c     are to be applied should be flagged with a boundary condition
c     code.
c
c     Note that this code no longer deals with faults.  Instead, they
c     are meant to be handled by an auxiliary code that reads
c     LaGriT-specific fault info.
c
c     At present, code is only set up for linear tetrahedra, although
c     this would be easy to change.
c
      implicit none
c
c...  parameters
c
      integer nsd,ndof,maxnodes,maxelmts,nen
      integer maxbnds,maxnattr,ietyp,inf
      parameter(nsd=3,ndof=3,maxnodes=1000000,maxelmts=5400000,nen=4,
     & maxbnds=20,maxnattr=20,ietyp=5,inf=0)
      integer kti,kto,kr,kw
      parameter(kti=5,kto=6,kr=10,kw=11)
      double precision eps
      parameter(eps=1.0d-7)
c
c...  local constants
c
      integer icflip(4)
      data icflip/1,3,2,4/
      integer kww
      data kww/12/
      integer izero
      data izero/0/
      double precision zero
      data zero/0.0d0/
c
c...  parameters read from parameter file
c
      integer nbc,iconopt
      integer ibfield(maxbnds),ibcode(maxbnds),ibc(nsd,maxbnds)
      integer iac(maxbnds),isn(nsd)
      double precision bc(nsd,maxbnds)
      double precision cscale
      character cunits*20,dunits*20,vunits*20,funits*20
c
c...  filenames
c
      character fileroot*200,pfile*200,ifile*200,nfile*200,cfile*200
      character bcfile*200,afile*200
c
c...  parameters and variables read from UCD file
c
      integer numnp,numel,nnattr,neattr,nmattr
      integer ien(nen,maxelmts),mat(maxelmts)
      double precision x(nsd,maxnodes),attrn(maxnattr)
c
c...  values read from auxiliary file
c
      double precision bca(ndof,maxnodes)
      double precision xa,ya,za
c
c...  external routines
c
      double precision tetcmp
      integer nnblnk,nchar
      external nnblnk,nchar,tetcmp
c
c...  local variables
c
      double precision det,sgn,xl(nsd,nen),xtmp(nsd),bct(3)
      integer ientmp(nen)
      integer itmp(maxnattr),idir(nsd)
      integer ibcnode(maxbnds+1,maxnodes),ibctmp(3)
      integer i,j,k,l,n,jj,i1,i2,j1,j2,nenl,nsdl,ndofl,nf
      integer iattr,kk
      integer elem,node,nflip
      integer numbnd,ibct
      double precision bctmp(3)
      character cstring*14,etype*3,cfnum*2,descr*10
      character dstring*21,vstring*17,fstring*14
      data cstring/"coord_units = "/
      data dstring/"displacement_units = "/
      data vstring/"velocity_units = "/
      data fstring/"force_units = "/
      character stout*50
      logical aux
c
      aux=.false.
      nenl=nen
      nsdl=nsd
      ndofl=ndof
      write(kto,*) "Enter root name for all files.  Both input and"
      write(kto,*) "output files will all have this prefix:"
      read(kti,"(a200)") fileroot
      i1=nnblnk(fileroot)
      i2=nchar(fileroot)
      call ifill(ibcnode,izero,(maxbnds+1)*maxnodes)
c
c...  read parameter file
c
      pfile=fileroot(i1:i2)//".par"
      open(file=pfile,unit=kr,status="old")
c
c...  coordinate scaling factor and units
c
      call pskip(kr)
      read(kr,*) cscale,(idir(i),i=1,nsd)
      if(cscale.eq.0.0d0) cscale=1.0d0
      call pskip(kr)
      read(kr,*) cunits
c
c...  bc info
c
      call pskip(kr)
      read(kr,*) nbc
      do i=1,nbc
        call pskip(kr)
        read(kr,*) ibfield(i),ibcode(i),
     &   (ibc(j,ibcode(i)),bc(j,ibcode(i)),j=1,nsd),iac(ibcode(i))
        if(iac(ibcode(i)).ne.0) aux=.true.
      end do
      call pskip(kr)
      read(kr,*) dunits
      call pskip(kr)
      read(kr,*) vunits
      call pskip(kr)
      read(kr,*) funits
c
c...  connectivity order option
c
      call pskip(kr)
      read(kr,*) iconopt
      close(kr)
c
c...  read and output nodal coordinates
c
      ifile=fileroot(i1:i2)//".inp"
      open(file=ifile,unit=kr,status="old")
      nfile=fileroot(i1:i2)//".coord"
      open(file=nfile,unit=kw,status="new")
      stout=cstring//cunits
      write(kw,"(a50)") stout
      call pskip(kr)
      read(kr,*) numnp,numel,nnattr,neattr,nmattr
      do i=1,numnp
        read(kr,*) n,(xtmp(j),j=1,nsd)
        do j=1,nsd
          x(j,i)=cscale*xtmp(idir(j))
        end do
        write(kw,"(i7,3(2x,1pe15.8))") i,(x(j,i),j=1,nsd)
      end do
      close(kw)
c
c...  read and output connectivity info
c
      cfile=fileroot(i1:i2)//".connect"
      open(file=cfile,unit=kw,status="new")
      nflip=0
      if(iconopt.eq.1) then
        do i=1,numel
          read(kr,*) n,mat(i),etype,(ien(j,i),j=1,nen)
          write(kw,"(i7,3i4,4i7)") i,ietyp,mat(i),inf,(ien(j,i),j=1,nen)
        end do
      else if(iconopt.eq.2) then
        do i=1,numel
          read(kr,*) n,mat(i),etype,(ientmp(j),j=1,nen)
          do j=1,nen
            ien(j,i)=ientmp(icflip(j))
          end do
          write(kw,"(i7,3i4,4i7)") i,ietyp,mat(i),inf,(ien(j,i),j=1,nen)
        end do
      else if(iconopt.eq.3) then
        do i=1,numel
          read(kr,*) n,mat(i),etype,(ientmp(j),j=1,nen)
          call lcoord(ientmp,x,xl,nsdl,nenl,numnp)
          det=tetcmp(xl)
          if(det.lt.0.0d0) then
            nflip=nflip+1
            do j=1,nen
              ien(j,i)=ientmp(icflip(j))
            end do
          else
            do j=1,nen
              ien(j,i)=ientmp(j)
            end do
          end if
          write(kw,"(i7,3i4,4i7)") i,ietyp,mat(i),inf,(ien(j,i),j=1,nen)
        end do
      end if
      close(kw)
c
c...  read nodal attributes to determine which nodes are associated
c     with each boundary condition code.
c
      read(kr,*) nf,(itmp(i),i=1,nf)
      do i=1,nf
        read(kr,*) descr
      end do
      do i=1,numnp
        read(kr,*) n,(attrn(j),j=1,nnattr)
        do j=1,nbc
          iattr=nint(attrn(ibfield(j)))
          if(iattr.eq.ibcode(j)) then
            ibcnode(1,i)=ibcnode(1,i)+1
            ibcnode(ibcnode(1,i)+1,i)=iattr
          end if
        end do
      end do
      close(kr)
c
c...  if auxiliary file is being used, read BC from it
c
      if(aux) then
        afile=fileroot(i1:i2)//".aux"
        open(file=afile,unit=kr,status="old")
        do i=1,numnp
          read(kr,*) xa,ya,za,(bca(j,i),j=1,ndofl)
        end do
        close(kr)
      end if
c
c...  output BC and BC coordinates
c
      bcfile=fileroot(i1:i2)//".bc"
      open(file=bcfile,unit=kw,status="new")
      stout=dstring//dunits
      write(kw,"(a50)") stout
      stout=vstring//vunits
      write(kw,"(a50)") stout
      stout=fstring//funits
      write(kw,"(a50)") stout
      do i=1,numnp
        numbnd=ibcnode(1,i)
        if(numbnd.ne.0) then
          ibctmp(1)=0
          ibctmp(2)=0
          ibctmp(3)=0
          bctmp(1)=0.0d0
          bctmp(2)=0.0d0
          bctmp(3)=0.0d0
          do k=1,numbnd
            ibct=ibcnode(k+1,i)
            if(iac(ibct).eq.0) then
              do j=1,3
                bct(j)=bc(j,ibct)
              end do
            else
              do j=1,3
                bct(j)=bca(j,i)
              end do
            end if
            do j=1,3
              if(ibctmp(j).eq.0) then
                ibctmp(j)=ibc(j,ibct)
                bctmp(j)=bct(j)
              end if
            end do
          end do
          write(kw,"(i7,3(2x,i3),3(2x,1pe15.8))") i,
     &     (ibctmp(k),k=1,nsd),(bctmp(k),k=1,nsd)
        end if
      end do
      close(kw)
      write(kto,700) nflip
700   format("Number of connectivities flipped:  ",i7)
      stop
      end
c
c
      subroutine fill(arr,val,nlen)
c
c...  subroutine to fill a double precision array with a given value
c
      implicit none
      integer nlen
      double precision val,arr(nlen)
      integer i
      do i=1,nlen
        arr(i)=val
      end do
      return
      end
c
c
      subroutine ifill(iarr,ival,nlen)
c
c...  subroutine to fill an integer array with a given value
c
      implicit none
      integer ival,nlen
      integer iarr(nlen)
      integer i
      do i=1,nlen
        iarr(i)=ival
      end do
      return
      end
c
c
      subroutine lcoord(ien,x,xl,nsd,nen,numnp)
c
c...  subroutine to localize element coordinates
c
      implicit none
      integer nsd,nen,numnp
      integer ien(nen)
      double precision x(nsd,numnp),xl(nsd,nen)
c
      integer i,j,ii
c
      do i=1,nen
        ii=ien(i)
        do j=1,nsd
          xl(j,i)=x(j,ii)
        end do
      end do
      return
      end
c
c
      function nchar(string)
c
c...  determines the minimum nonblank length of a string
c
      implicit none
c
c...  parameter definitions
c
      character blank*1
      parameter(blank=' ')
c
c...  function arguments
c
      integer nchar
      character*(*) string
c
c...  intrinsic functions
c
      intrinsic len
c
c...  local variables
c
      integer nmax,i,itest
c
      nmax=len(string)
      nchar=0
      do i=1,nmax
        itest=nmax-i+1
        if(string(itest:itest).ne.blank) then
          nchar=itest
          return
        end if
      end do
      return
      end
c
c
      function nnblnk(string)
c
c       determines the position of the first nonblank entry
c       of a string (returns 1 if the first character is
c       not blank)
c
      implicit none
c
c...  parameter definitions
c
      character blank*1
      parameter(blank=' ')
c
c...  function arguments
c
      integer nnblnk
      character*(*) string
c
c... intrinsic functions
c
      intrinsic len
c
c...  local variables
c
      integer nmax,i
c
      nmax=len(string)
      nnblnk=nmax
      do i=1,nmax
        if(string(i:i).ne.blank) then
          nnblnk=i
          return
        end if
      end do
      return
      end
c
c
      subroutine pskip(iunit)
c
c      routine to skip lines beginning with the string # and blank
c      lines.
c      this routine ignores leading blanks before the key string.
c
c
      implicit none
c
c...  subroutine arguments
c
      integer iunit
c
c...  local constants
c
      character leader*1
      data leader/'#'/
c
c...  intrinsic functions
c
      intrinsic index
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer inblnk
      character string*80
c
 10   continue
        read(iunit,"(a80)",end=20) string
        if(nchar(string).eq.0) goto 10
        inblnk=nnblnk(string)
        if(index(string,leader).eq.inblnk) goto 10
      backspace(iunit)
 20   continue
      return
      end
c
c
      function tetcmp(xl)
c
c...  program to compute determinant and equivalent tetrahedral volume
c     for a given 4x4 matrix.
c
      implicit none
      double precision tetcmp,xl(3,4)
      double precision a(4,4),cof(3,3)
      data a(1,1),a(1,2),a(1,3),a(1,4)/1.0d0,1.0d0,1.0d0,1.0d0/
      integer i,j,ii
      double precision det3,onem,sgn
      external det3
      integer ind(3,4)
      data ind/2,3,4,1,3,4,1,2,4,1,2,3/
      do i=2,4
        do j=1,4
          a(i,j)=xl(i-1,j)
        end do
      end do
      onem=-1.0d0
      sgn=onem
      tetcmp=0.0d0
      do i=1,4
        sgn=sgn*onem
        do ii=1,3
          do j=2,4
            cof(j-1,ii)=a(j,ind(ii,i))
          end do
        end do
        tetcmp=tetcmp+sgn*a(1,i)*det3(cof)
cdebug        write(6,*) "i,a,tetcmp:",i,(a(j,i),j=1,4)
      end do
      return
      end
c
c
      function det3(x)
c
c...  function to compute determinant of a 3x3 matrix
c
      implicit none
      double precision det3
      double precision x(3,3)
      det3=x(1,1)*x(2,2)*x(3,3)-x(1,1)*x(2,3)*x(3,2)+
     &     x(1,2)*x(2,3)*x(3,1)-x(1,2)*x(2,1)*x(3,3)+
     &     x(1,3)*x(2,1)*x(3,2)-x(1,3)*x(2,2)*x(3,1)
      return
      end
