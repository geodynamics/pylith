      program readucd2
c
c...  quick and dirty code to translate AVS UCD format to a format
c     suitable for lithomop or pylith.  All boundaries for which BC
c     are to be applied should be flagged with a boundary condition
c     code.  This also applies to any faults in the model.
c
c     At present, code is only set up for linear tetrahedra, although
c     this would be easy to change.
c
      implicit none
c
c...  parameters
c
      integer nsd,ndof,maxnodes,maxelmts,nen,maxflts,maxfnodes
      integer maxbnds,maxfelems,maxnattr,maxeattr,ietyp,inf
      parameter(nsd=3,ndof=3,maxnodes=1000000,maxelmts=5400000,nen=4,
     & maxflts=6,maxfnodes=1000000,maxfelems=1000000,maxbnds=30,
     & maxnattr=20,maxeattr=20,ietyp=5,inf=0)
      integer kti,kto,kr,kw
      parameter(kti=5,kto=6,kr=10,kw=11)
      double precision eps
      parameter(eps=1.0d-7)
c
c...  local constants
c
      integer icflip(4)
      data icflip/1,3,2,4/
      integer kww,kwfb(maxflts),kwfc(maxflts)
      data kww/12/
      data kwfb/12,13,14,15,16,17/
      data kwfc/18,19,20,21,22,23/
      integer izero
      data izero/0/
      double precision zero
      data zero/0.0d0/
c
c...  parameters read from parameter file
c
      integer nbc,iconopt,numflt,ibfield
      integer ibcode(maxbnds),ibc(nsd,maxbnds),iac(maxbnds),isn(nsd)
      integer iftype(maxflts),iffield(maxflts),ifcode(maxflts)
      integer ifefield(maxflts),ifnorm(3)
      double precision bc(nsd,maxbnds),fsplit(nsd,2,maxflts)
      double precision cscale
      character cunits*20,dunits*20,vunits*20,funits*20
c
c...  filenames
c
      character fileroot*200,pfile*200,ifile*200,nfile*200,cfile*200
      character bcfile*200,afile*200
      character fbcfile*200,fbccfile*200
c
c...  parameters and variables read from UCD file
c
      integer numnp,numel,nnattr,neattr,nmattr
      integer ien(nen,maxelmts),mat(maxelmts)
      double precision x(nsd,maxnodes),attrn(maxnattr),attre(maxeattr)
      double precision fnorm(nsd,maxnodes)
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
      integer ifnodes(maxflts,maxfnodes),ifelems(2,maxflts,maxfelems)
      integer nfnodes(maxflts),nfelems(2,maxflts)
      integer ientmp(nen)
      integer ifhist(maxflts),itmp(maxnattr),idir(nsd)
      integer ibcnode(maxnodes)
      integer i,j,k,l,n,jj,i1,i2,j1,j2,nenl,nsdl,ndofl,nfltnodes,nf
      integer iattr,kk
      integer elem,node,nflip
      integer ibct
      character cstring*14,etype*3,cfnum*2,descr*10
      character dstring*21,vstring*17,fstring*14
      data cstring/"coord_units = "/
      data dstring/"displacement_units = "/
      data vstring/"velocity_units = "/
      data fstring/"force_units = "/
      character stout*50
      logical aux,getnorm
c
      aux=.false.
      getnorm=.false.
      nenl=nen
      nsdl=nsd
      ndofl=ndof
      write(kto,*) "Enter root name for all files.  Both input and"
      write(kto,*) "output files will all have this prefix:"
      read(kti,"(a200)") fileroot
      i1=nnblnk(fileroot)
      i2=nchar(fileroot)
      call fill(fnorm,zero,3*maxnodes)
      call ifill(ibcnode,izero,maxnodes)
      call ifill(nfnodes,izero,maxflts)
      call ifill(nfelems,izero,2*maxflts)
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
      read(kr,*) nbc,ibfield
      do i=1,nbc
        call pskip(kr)
        read(kr,*) ibcode(i),(ibc(j,ibcode(i)),bc(j,ibcode(i)),j=1,nsd),
     &   iac(ibcode(i))
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
c
c...  fault definitions
c
      call pskip(kr)
      read(kr,*) numflt,(ifnorm(i),i=1,3)
      if(ifnorm(1).gt.0.and.ifnorm(2).gt.0.and.ifnorm(3).gt.0)
     & getnorm=.true.
      do i=1,numflt
        call pskip(kr)
        read(kr,*) iftype(i),ifhist(i),ifcode(i),iffield(i),ifefield(i)
        call pskip(kr)
        read(kr,*) (fsplit(j,1,i),j=1,nsd)
        call pskip(kr)
        read(kr,*) (fsplit(j,2,i),j=1,nsd)
      end do
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
c     with each fault and boundary condition code.  Get node normals
c     if they are available.
c
      read(kr,*) nf,(itmp(i),i=1,nf)
      do i=1,nf
        read(kr,*) descr
      end do
      nfltnodes=0
      do i=1,numnp
        read(kr,*) n,(attrn(j),j=1,nnattr)
        if(getnorm) then
          fnorm(1,i)=attrn(ifnorm(1))
          fnorm(2,i)=attrn(ifnorm(2))
          fnorm(3,i)=attrn(ifnorm(3))
        end if
        iattr=nint(attrn(ibfield))
        do j=1,nbc
          if(iattr.eq.ibcode(j)) ibcnode(i)=iattr
        end do
c
        do j=1,numflt
          iattr=nint(attrn(iffield(j)))
          if(iattr.eq.ifcode(j)) then
            nfnodes(j)=nfnodes(j)+1
            nfltnodes=nfltnodes+1
            ifnodes(j,nfnodes(j))=n
          end if
        end do
      end do
      write(kto,*) (nfnodes(j),j=1,numflt)
c
c...  read element attributes to determine which elements are adjacent
c     to each fault
c
      read(kr,*) nf,(itmp(i),i=1,nf)
      do i=1,nf
        read(kr,*) descr
      end do
      do i=1,numel
        read(kr,*) n,(attre(j),j=1,neattr)
        do j=1,numflt
          iattr=nint(attre(ifefield(j)))
          if(iattr.eq.1) then
            nfelems(1,j)=nfelems(1,j)+1
            ifelems(1,j,nfelems(1,j))=n
          else if(iattr.eq.-1) then
            nfelems(2,j)=nfelems(2,j)+1
            ifelems(2,j,nfelems(2,j))=n
          end if
        end do
      end do
      write(kto,*) (nfelems(1,j),j=1,numflt)
      write(kto,*) (nfelems(2,j),j=1,numflt)
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
        ibct=ibcnode(i)
        if(ibct.ne.0) then
          if(iac(ibct).eq.0) then
            do j=1,3
              bct(j)=bc(j,ibct)
            end do
          else
            do j=1,3
              bct(j)=bca(j,i)
            end do
          end if
          write(kw,"(i7,3(2x,i3),3(2x,1pe15.8))") i,
     &     (ibc(k,ibct),k=1,nsd),(bct(k),k=1,nsd)
        end if
      end do
      close(kw)
c
c...  output fault info after determining which fault each node lies on
c
      do i=1,numflt
        write(cfnum,"(i2)") i
        j1=nnblnk(cfnum)
        j2=nchar(cfnum)
        fbcfile=fileroot(i1:i2)//"."//cfnum(j1:j2)//".fbc"
        fbccfile=fileroot(i1:i2)//"."//cfnum(j1:j2)//".fcoord"
        open(file=fbcfile,unit=kwfb(i),status="new")
        open(file=fbccfile,unit=kwfc(i),status="new")
      end do
c
c...  Loop over positive sides of fault, then negative.
c
      do kk=1,2
        do i=1,numflt
          do j=1,nsd
            sgn=sign(1.0d0,fsplit(j,kk,i))
            isn(j)=nint(sgn)
            if(fsplit(j,kk,i).eq.0.0d0) isn(j)=0
          end do
c
c...  First loop over faulted elements, then nodes attached to each element.
c
          do j=1,nfelems(kk,i)
            elem=ifelems(kk,i,j)
            do l=1,nen
              do k=1,nfnodes(i)
                node=ifnodes(i,k)
                if(ien(l,elem).eq.node) then
                  if(iftype(i).eq.1) then
                    write(kwfb(i),"(2i7,i4,3(2x,1pe15.8))") elem,
     &               node,ifhist(i),(fsplit(jj,kk,i),jj=1,nsd)
                  else
                    write(kwfb(i),"(2i7,3i4)") elem,node,
     &               (isn(jj),jj=1,nsd)
                  end if
                  write(kwfc(i),"(3i7,6(2x,1pe15.8))") 
     &             elem,node,kk,(x(jj,node),jj=1,nsd),
     &             (fnorm(jj,node),jj=1,nsd)
                end if
              end do
            end do
          end do
        end do
      end do
      do i=1,numflt
        close(kwfb(i))
        close(kwfc(i))
      end do
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
