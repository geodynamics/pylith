      program prescomp
c
c...  quick and dirty program to compute pressure BC over a given
c     set of faces.  The required inputs are a set of nodal
c     coordinates (PyLith/LithoMop format) and an existing traction
c     input file in PyLith format, with the vertices for each face
c     ordered CCW when looking at the face.
c
      implicit none
c
c...  parameters
c
      integer maxnodes,maxtrac,nsd
      parameter(maxnodes=2000000,maxtrac=2000000,nsd=3)
c
c...  global arrays
c
      integer nvert
      double precision x(nsd,maxnodes)
      double precision pval
c
c...  unit numbers
c
      integer kti,kto,kn,kt,kw,ka
      common/units/kti,kto,kn,kt,kw,ka
c
c...  local variables
c
      integer i,j,n,numnp
      integer ivert(4)
      double precision pvec(nsd),fcent(nsd)
      character dummy*80
      logical aux
c
c...  get command-line arguments
c
      call files(nvert,pval,aux)
c
c...  read nodal coordinates
c
      call pskip(kn)
      read(kn,"(a80)") dummy
      call pskip(kn)
      numnp=0
      do i=1,maxnodes
        read(kn,*,end=10) n,(x(j,i),j=1,nsd)
        numnp=numnp+1
      end do
 10   continue
      close(kn)
c
c...  read tractions and output computed results
c
      call pskip(kt)
      read(kt,"(a80)") dummy
      write(kw,"(a80)") dummy
      call pskip(kt)
      do i=1,maxtrac
        read(kt,*,end=20) (ivert(j),j=1,nvert)
        call getvec(x,pvec,fcent,ivert,pval,nsd,numnp,nvert,aux)
        write(kw,"(3i8,3(2x,1pe15.8))") (ivert(j),j=1,nvert),
     &   (pvec(j),j=1,nsd)
        if(aux) write(ka,"(3i8,3(2x,1pe15.8))") (ivert(j),j=1,nvert),
     &   (fcent(j),j=1,nsd)
      end do
 20   continue
      close(kt)
      close(kw)
      if(aux) close(ka)
      stop
      end
c
c
      subroutine cross(x1,x2,x3,pvec)
c
c...  routine to compute vector cross product given 3 points
c
      implicit none
c
c...  subroutine arguments
c
      double precision x1(3),x2(3),x3(3),pvec(3)
c
c...  local variables
c
      integer i
      double precision v(3),u(3)
c
c...  compute vectors
c
      do i=1,3
        u(i)=x1(i)-x2(i)
        v(i)=x3(i)-x2(i)
      end do
c
c...  compute cross product
c
      pvec(1)=v(2)*u(3)-v(3)*u(2)
      pvec(2)=v(3)*u(1)-v(1)*u(3)
      pvec(3)=v(1)*u(2)-v(2)*u(1)
      return
      end
c
c
      subroutine getvec(x,pvec,fcent,ivert,pval,nsd,numnp,nvert,aux)
c
c...  routine to compute stress vector for a given face
c
      implicit none
c
c...  subroutine arguments
c
      integer nsd,numnp,nvert
      integer ivert(nvert)
      double precision x(nsd,numnp),pvec(nsd),fcent(nsd),pval
      logical aux
c
c...  local variables
c
      integer i,j
      double precision v(3),vmag,fac
c
c...  compute cross product to get normal.
c     If nvert == 4, the average of 2 vectors is used.
c
      call cross(x(1,ivert(3)),x(1,ivert(2)),x(1,ivert(1)),pvec)
      if(nvert.eq.4) then
        call cross(x(1,ivert(1)),x(1,ivert(4)),x(1,ivert(3)),v)
        do i=1,nsd
          pvec(i)=0.5d0*(pvec(i)+v(i))
        end do
      end if
c
c...  scale vector
c
      vmag=0.0d0
      do i=1,nsd
        vmag=vmag+pvec(i)*pvec(i)
        fcent(i)=0.0d0
      end do
      vmag=sqrt(vmag)
      do i=1,nsd
        pvec(i)=pval*pvec(i)/vmag
      end do
c
c...  compute centroid of face
c
      do i=1,nvert
        do j=1,nsd
          fcent(j)=fcent(j)+x(j,ivert(i))
        end do
      end do
      fac=1.0d0/dble(nvert)
      fcent(1)=fcent(1)*fac
      fcent(2)=fcent(2)*fac
      fcent(3)=fcent(3)*fac
      return
      end
c
c
      subroutine files(nvert,pval,aux)
c
c...  subroutine to get command-line arguments
c
      implicit none
c
c...  subroutine arguments
c
      integer nvert
      double precision pval
      logical aux
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
      integer kti,kto,kn,kt,kw,ka
      common/units/kti,kto,kn,kt,kw,ka
c
c...  local variables
c
      integer nargs,i,j
      character cfile*500,tfile*500,ofile*500,afile*500,num*20
      character string*2000
      logical fflag(5)
c
      kti=5
      kto=6
      kn=10
      kt=11
      kw=12
c
      do i=1,5
        fflag(i)=.false.
      end do
c
      nargs=iargc()
      do i=1,nargs
        call getarg(i,string)
        if(index(string,'c=').ne.0) then
          j=nchar(string)
          cfile=string(3:j)
          fflag(1)=.true.
        else if(index(string,'t=').ne.0) then
          j=nchar(string)
          tfile=string(3:j)
          fflag(2)=.true.
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          ofile=string(3:j)
          fflag(3)=.true.
        else if(index(string,'a=').ne.0) then
          j=nchar(string)
          afile=string(3:j)
          aux=.true.
        else if(index(string,'n=').ne.0) then
          j=nchar(string)
          num=string(3:j)
          read(num,*) nvert
          fflag(4)=.true.
        else if(index(string,'p=').ne.0) then
          j=nchar(string)
          num=string(3:j)
          read(num,*) pval
          fflag(5)=.true.
        end if
      end do
c
      do i=1,5
        if(.not.fflag(i)) then
          write(kto,800)
          stop
        end if
      end do
      open(kn,file=cfile,status="old")
      open(kt,file=tfile,status="old")
      open(kw,file=ofile,status="new")
      if(aux) open(ka,file=afile,status="new")
c
 800  format("Usage:",/,
     & "    prescomp c=<coordinate_file> t=<traction_input_file>",/,
     & "    o=<traction_output_file> [a=<auxiliary_output_file>]",/,
     & "    n=<num_vertices_per_face> p=<pressure_value>")
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
