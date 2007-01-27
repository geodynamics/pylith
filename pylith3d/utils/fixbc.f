      program fixbc
c
c...  Code to modify a set of boundary conditions along a specified
c     set of lines.  The code assumes a boundary condition input file
c     suitable for use by PyLith or LithoMop and outputs a files in the
c     same format.
c
      implicit none
c
c...  parameters
c
      integer maxnodes,maxlines
      parameter(maxnodes=5000000,maxlines=1000)
c
c...  global arrays
c
      double precision x(3,maxnodes),xl(3,2,maxlines),bcc(3,maxlines)
      double precision dc(3,maxlines)
      integer ibcc(3,maxlines),iuse(maxlines)
c
c...  dimensions
c
      integer nnodes,nlines
c
c...  filenames and unit numbers
c
      integer kti,kto,kp,kc,kbi,kbo
      common/units/kti,kto,kp,kc,kbi,kbo
      character pfile*200,cfile*200,bifile*200,bofile*200
c
c...  local variables
c
      integer i,j,k,node
      integer ibc(3)
      double precision eps,bc(3),rmag
      character string*80
      logical nused
c
c...  get filenames
c
      call files(pfile,cfile,bifile,bofile)
c
c...  read parameters
c
      open(kp,file=pfile,status="old")
      call pskip(kp)
      read(kp,*) nlines,eps
      do i=1,nlines
        call pskip(kp)
        read(kp,*) (xl(j,1,i),j=1,3),(xl(j,2,i),j=1,3)
        call pskip(kp)
        read(kp,*) (ibcc(j,i),j=1,3),(bcc(j,i),j=1,3)
        rmag=0.0d0
        nused=.false.
        do j=1,3
          dc(j,i)=xl(j,2,i)-xl(j,1,i)
          rmag=rmag+dc(j,i)*dc(j,i)
          if(abs(dc(j,i)).gt.eps.and.(.not.nused)) then
            iuse(i)=j
            nused=.true.
          end if
        end do
        do j=1,3
          dc(j,i)=dc(j,i)/sqrt(rmag)
        end do
      end do
      close(kp)
c
c...  read nodal coordinates
c
      open(kc,file=cfile,status="old")
      nnodes=0
      call pskip(kc)
      read(kc,*) string
      call pskip(kc)
 10   continue
        read(kc,*,end=20) node,(x(j,node),j=1,3)
        nnodes=nnodes+1
        go to 10
 20   continue
      close(kc)
c
c...  loop over BC entries
c
      open(kbi,file=bifile,status="old")
      open(kbo,file=bofile,status="new")
      call pskip(kbi)
      read(kbi,"(a80)") string
      write(kbo,"(a80)") string
      read(kbi,"(a80)") string
      write(kbo,"(a80)") string
      read(kbi,"(a80)") string
      write(kbo,"(a80)") string
      call pskip(kbi)
 30   continue
        read(kbi,*,end=40) node,(ibc(j),j=1,3),(bc(j),j=1,3)
        call getvals(x(1,node),ibc,bc,ibcc,bcc,xl,dc,iuse,eps,nlines)
        write(kbo,"(4i7,3(2x,1pe15.8))") node,(ibc(j),j=1,3),
     &   (bc(j),j=1,3)
        go to 30
 40   continue
      close(kbi)
      close(kbo)
      stop
      end
c
c
      subroutine files(pfile,cfile,bifile,bofile)
c
c...  subroutine to get filenames from command-line
c
      implicit none
      character pfile*(*),cfile*(*),bifile*(*),bofile*(*)
c
c...  intrinsic functions
c
      intrinsic iargc,index
c
c...  external functions
c
      integer nchar
      external nchar
c
c...  unit numbers
c
      integer kti,kto,kp,kc,kbi,kbo
      common/units/kti,kto,kp,kc,kbi,kbo
c
c...  local variables
c
      integer nargs,i,j
      character string*2000
      logical fflag(4)
c
      kti=5
      kto=6
      kp=10
      kc=11
      kbi=12
      kbo=14
c
      do i=1,4
        fflag(i)=.false.
      end do
c
      nargs=iargc()
      do i=1,nargs
        call getarg(i,string)
        if(index(string,'p=').ne.0) then
          j=nchar(string)
          pfile=string(3:j)
          fflag(1)=.true.
        else if(index(string,'c=').ne.0) then
          j=nchar(string)
          cfile=string(3:j)
          fflag(2)=.true.
        else if(index(string,'i=').ne.0) then
          j=nchar(string)
          bifile=string(3:j)
          fflag(3)=.true.
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          bofile=string(3:j)
          fflag(4)=.true.
        end if
      end do
c
      do i=1,4
        if(.not.fflag(i)) then
          write(kto,800)
          stop
        end if
      end do
c
 800  format("Usage:",/,
     & "    fixbc p=<parameter_file> c=<coord_file>",/,
     & "    i=<boundary_condition_input_file>",/,
     & "    o=<boundary_condition_output_file>")
      return
      end
c
c
      subroutine getvals(x,ibc,bc,ibcc,bcc,xl,dc,iuse,eps,nlines)
c
c...  subroutine to modify BC if they fall on one of the specified lines
c
      implicit none
c
c...  subroutine arguments
c
      integer nlines
      integer ibc(3),ibcc(3,nlines),iuse(nlines)
      double precision x(3),bc(3),bcc(3,nlines),xl(3,2,nlines)
      double precision dc(3,nlines),eps
      logical inline
c
c...  intrinsic functions
c
      intrinsic sign,sqrt
c
c...  local variables
c
      integer i,j,ind
      double precision t,xp,xmin,xmax
c
c...  loop over line segments to see if point lies on any of them
c
      do i=1,nlines
        ind=iuse(i)
        t=(x(ind)-xl(ind,1,i))/dc(ind,i)
        inline=.true.
        do j=1,3
          xp=xl(j,1,i)+t*dc(j,i)
          if(abs(xp-x(j)).gt.eps) inline=.false.
          xmin=min(xl(j,1,i),xl(j,2,i))
          xmax=max(xl(j,1,i),xl(j,2,i))
          if(x(j).lt.xmin.or.x(j).gt.xmax) inline=.false.
        end do
        if(inline) then
          do j=1,3
            if(ibcc(j,i).ge.0) then
              ibc(j)=ibcc(j,i)
              bc(j)=bcc(j,i)
            end if
          end do
          return
        end if
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
