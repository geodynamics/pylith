      program readnetgen
c
c...  quick and dirty code to translate netgen neutral format to a
c     format suitable for lithomop.  All boundaries for which BC
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
      integer nsd,maxnodes,maxelmts,nen,maxnbrs,maxflts,maxfnodes
      integer maxbnds,maxmatf,ietyp,inf
      parameter(nsd=3,maxnodes=500000,maxelmts=500000,nen=4,
     & maxnbrs=100,maxflts=6,maxfnodes=100000,
     & maxbnds=30,maxmatf=5,ietyp=5,inf=0)
      integer kti,kto,kr,kw
      parameter(kti=5,kto=6,kr=10,kw=11)
      double precision eps
      parameter(eps=1.0d-7)
c
c...  local constants
c
      integer icflip(4)
      data icflip/1,3,2,4/
      integer kwb,kwc,kwfb(maxflts),kwfc(maxflts)
      data kwb/12/
      data kwc/13/
      data kwfb/12,13,14,15,16,17/
      data kwfc/18,19,20,21,22,23/
      integer izero
      data izero/0/
c
c...  parameters read from parameter file
c
      integer nbc,iconopt
      integer ibcode(maxbnds),ibc(nsd,maxbnds),isn(nsd)
      integer iftype(maxflts),nmatf(2,maxflts),ifmat(2,maxflts,maxmatf)
      double precision bc(nsd,maxbnds),fsplit(nsd,2,maxflts)
      double precision cscale
c
c...  filenames
c
      character fileroot*200,pfile*200,ifile*200,nfile*200,cfile*200
      character bcfile*200,bccfile*200
      character fbcfile*200,fbccfile*200
c
c...  parameters and variables read from netgen file
c
      integer numnp,numel,numflt,nbcfac,ibcflt,ibcfac
      integer ien(nen,maxelmts),mat(maxelmts),ibcvert(3)
      double precision x(nsd,maxnodes)
c
c...  external routines
c
      double precision tetcmp
      integer nnblnk,nchar
      external nnblnk,nchar,tetcmp
c
c...  local variables
c
      double precision det,sgn,xl(nsd,nen),xtmp(nsd)
      integer iadjf(maxflts*maxfnodes,maxnbrs),inodef(maxflts*maxfnodes)
      integer ientmp(nen),nelsf(maxflts*maxfnodes),ifltnode(maxnodes)
      integer ifhist(maxflts),idir(nsd),ibcnode(maxnodes,maxbnds)
      integer nbcnodes(maxbnds)
      integer i,j,k,l,jj,i1,i2,j1,j2,nenl,nsdl,nfltnodes,kk
      integer iel,nflip
      character cstring*15,cbound(maxbnds)*3,cfnum*2
      data cstring/"coord_units = m"/
      data cbound/
     & "b01","b02","b03","b04","b05","b06","b07","b08","b09","b10",
     & "b11","b12","b13","b14","b15","b16","b17","b18","b19","b20",
     & "b21","b22","b23","b24","b25","b26","b27","b28","b29","b30"/
c
      nenl=nen
      nsdl=nsd
      write(kto,*) "Enter root name for all files.  Both input and"
      write(kto,*) "output files will all have this prefix:"
      read(kti,"(a200)") fileroot
      i1=nnblnk(fileroot)
      i2=nchar(fileroot)
      call ifill(ibcnode,izero,maxnodes*maxbnds)
      call ifill(nbcnodes,izero,maxbnds)
c
c...  read parameter file
c
      pfile=fileroot(i1:i2)//".par"
      open(file=pfile,unit=kr,status="old")
c
c...  coordinate scaling factor
c
      call pskip(kr)
      read(kr,*) cscale,(idir(i),i=1,nsd)
      if(cscale.eq.0.0d0) cscale=1.0d0
c
c...  bc info
c
      call pskip(kr)
      read(kr,*) nbc
      do i=1,nbc
        call pskip(kr)
        read(kr,*) ibcode(i),(ibc(j,i),bc(j,i),j=1,nsd)
      end do
c...  connectivity order option
c
      call pskip(kr)
      read(kr,*) iconopt
c
c...  fault definitions
c
      call pskip(kr)
      read(kr,*) numflt,ibcflt
      do i=1,numflt
        call pskip(kr)
        read(kr,*) iftype(i),ifhist(i),nmatf(1,i),nmatf(2,i)
        call pskip(kr)
        read(kr,*) (fsplit(j,1,i),j=1,nsd),(ifmat(1,i,j),j=1,nmatf(1,i))
        call pskip(kr)
        read(kr,*) (fsplit(j,2,i),j=1,nsd),(ifmat(2,i,j),j=1,nmatf(2,i))
      end do
      close(kr)
c
c...  read and output nodal coordinates
c
      ifile=fileroot(i1:i2)//".netgen"
      open(file=ifile,unit=kr,status="old")
      nfile=fileroot(i1:i2)//".coord"
      open(file=nfile,unit=kw,status="new")
      write(kw,"(a15)") cstring
      read(kr,*) numnp
      do i=1,numnp
        ifltnode(i)=0
        read(kr,*) (xtmp(j),j=1,nsd)
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
      read(kr,*) numel
      if(iconopt.eq.1) then
        do i=1,numel
          read(kr,*) mat(i),(ien(j,i),j=1,nen)
          write(kw,"(i7,3i4,4i7)") i,ietyp,mat(i),inf,(ien(j,i),j=1,nen)
        end do
      else if(iconopt.eq.2) then
        do i=1,numel
          read(kr,*) mat(i),(ientmp(j),j=1,nen)
          do j=1,nen
            ien(j,i)=ientmp(icflip(j))
          end do
          write(kw,"(i7,3i4,4i7)") i,ietyp,mat(i),inf,(ien(j,i),j=1,nen)
        end do
      else if(iconopt.eq.3) then
        do i=1,numel
          read(kr,*) mat(i),(ientmp(j),j=1,nen)
          call lcoord(ientmp,x,xl,nsdl,nenl,numnp)
          det=tetcmp(xl)
cdebug          write(kto,*) "i,det:",i,det
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
c...  read face BC codes to determine which nodes are associated
c     with each boundary and fault.
c
      read(kr,*) nbcfac
      do i=1,nbcfac
        read(kr,*) ibcfac,(ibcvert(j),j=1,3)
        if(ibcfac.eq.ibcflt) then
          do j=1,3
            ifltnode(ibcvert(j))=1
          end do
        else
          do j=1,nbc
            if(ibcfac.eq.ibcode(j)) then
              do k=1,3
                ibcnode(ibcvert(k),j)=ibcfac
              end do
            end if
          end do
        end if
      end do
      nfltnodes=0
      do i=1,numnp
        if(ifltnode(i).eq.1) then
          nfltnodes=nfltnodes+1
          inodef(nfltnodes)=i
          nelsf(nfltnodes)=0
        end if
      end do
c
c...  output BC and BC coordinates
c
      do i=1,nbc
        bcfile=fileroot(i1:i2)//"."//cbound(i)//".bc"
        bccfile=fileroot(i1:i2)//"."//cbound(i)//".coord"
        open(file=bcfile,unit=kwb,status="new")
        open(file=bccfile,unit=kwc,status="new")
        do j=1,numnp
          if(ibcnode(j,i).ne.izero) then
            write(kwb,"(i7,3(2x,i3),3(2x,1pe15.8))") j,
     &       (ibc(k,i),k=1,nsd),(bc(k,i),k=1,nsd)
            write(kwc,"(i7,3(2x,1pe15.8))") j,(x(k,j),k=1,nsd)
          end if
        end do
        close(kwb)
        close(kwc)
      end do
c
c...  determine which elements contain each node on the faults
c
      do i=1,numel
        do j=1,nfltnodes
          do k=1,nen
            if(ien(k,i).eq.inodef(j)) then
              nelsf(j)=nelsf(j)+1
              iadjf(j,nelsf(j))=i
            end if
          end do
        end do
      end do
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
c...  First find all elements on one side of the fault, then the other
c
      do kk=1,2
        do i=1,numflt
          do j=1,nsd
            sgn=sign(1.0d0,fsplit(j,kk,i))
            isn(j)=nint(sgn)
          end do
          do j=1,nfltnodes
            do k=1,nelsf(j)
              do l=1,nmatf(kk,i)
                iel=iadjf(j,k)
                if(mat(iel).eq.ifmat(kk,i,l)) then
                  if(iftype(i).eq.1) then
                    write(kwfb(i),"(2i7,i4,3(2x,1pe15.8))") iel,
     &               inodef(j),ifhist(i),(fsplit(jj,kk,i),jj=1,nsd)
                  else
                    write(kwfb(i),"(2i7,3i4)") iel,inodef(j),
     &               (isn(jj),jj=1,nsd)
                  end if
                  write(kwfc(i),"(3i7,3(2x,1pe15.8))") 
     &             iel,inodef(j),kk,(x(jj,inodef(j)),jj=1,nsd)
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
