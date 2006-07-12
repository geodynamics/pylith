      program scaleucd
c
c...  Program to scale the coordinates produced by netgen and output
c     them to a separate file.
c
      implicit none
c
c...  filenames and unit numbers
c
      integer kti,kto,kr,kw
      common/units/kti,kto,kr,kw
      character ifile*200,ofile*200
c
c...  local variables
c
      integer numnp,numel,nnattr,neattr,nmattr,n,i,j
      double precision x(3),scale
c
c...  get filenames and scale factor
c
      call files(scale,ifile,ofile)
      open(kr,file=ifile,status="old")
      open(kw,file=ofile,status="new")
c
c...  get number of nodes and loop over them
c
      read(kr,*) numnp,numel,nnattr,neattr,nmattr
      do i=1,numnp
        read(kr,*) n,(x(j),j=1,3)
        write(kw,"(3(2x,1pe15.8))") (scale*x(j),j=1,3)
      end do
      close(kr)
      close(kw)
      stop
      end
c
c
      subroutine files(scale,ifile,ofile)
c
c...  subroutine to get filenames and scale factor from command-line
c
      implicit none
      double precision scale
      character ifile*(*),ofile*(*)
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
      integer kti,kto,kr,kw
      common/units/kti,kto,kr,kw
c
c...  local variables
c
      integer nargs,i,j
      character string*2000
      logical fflag(3)
c
      kti=5
      kto=6
      kr=10
      kw=11
c
      do i=1,3
        fflag(i)=.false.
      end do
c
      nargs=iargc()
      do i=1,nargs
        call getarg(i,string)
        if(index(string,'i=').ne.0) then
          j=nchar(string)
          ifile=string(3:j)
          fflag(1)=.true.
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          ofile=string(3:j)
          fflag(2)=.true.
        else if(index(string,'s=').ne.0) then
          j=nchar(string)
          read(string(3:j),*) scale
          fflag(3)=.true.
        end if
      end do
c
      do i=1,3
        if(.not.fflag(i)) then
          write(kto,800)
          stop
        end if
      end do
c
 800  format("Usage:",/,
     & "    scaleucd i=<input_file> o=<output_file>",/,
     & "    s=<scale_factor>")
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
