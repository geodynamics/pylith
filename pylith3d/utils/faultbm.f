      program faultbm
c
c...  quick and dirty program to compute fault boundary conditions,
c     given a set of points.  The user specifies coordinate ranges
c     for a number of regions, along with a set of coefficients
c     for a general equation of second degree in x, y, and z.
c     This program is suitable for specifying split node BC for the
c     program LithoMop.
c
c     This version of the code assumes that output has already been
c     produced by a code such as readnetgen or readucd, such that
c     fault information files (*.fcoord, *.fbc) have already been
c     produced.
c
c     This code has been modified from the original faultcalc2 code
c     to deal with overlapping regions in a special way such that
c     the minimum computed slip is used.  This is a specialized feature
c     that is really only useful for the SCEC benchmarks, which have
c     been set up this way.
c
      implicit none
c
c...  parameters
c
      integer maxsurf
      parameter(maxsurf=200)
c
c...  global arrays
c
      double precision scoef(10,3,maxsurf)
      double precision smin(3,maxsurf),smax(3,maxsurf)
c
c...  dimensions
c
      integer nsurf,nfnodes
c
c...  filenames and unit numbers
c
      integer kti,kto,kp,kn,ksi,kso
      common/units/kti,kto,kp,kn,ksi,kso
      character pfile*200,fnfile*200,sifile*200,sofile*200
c
c...  local variables
c
      integer i,j,k,nodes,elems,hists,nodec,elemc,isidec
      double precision eps,ux,uy,uz
      double precision x(3)
      character string*80
      logical insurf,usedc
c
c...  get filenames
c
      call files(pfile,fnfile,sifile,sofile)
c
c...  read parameters
c
      open(kp,file=pfile,status="old")
      call pskip(kp)
      read(kp,*) nsurf,eps
      do i=1,nsurf
        do j=1,3
          call pskip(kp)
          read(kp,*) smin(j,i),smax(j,i)
        end do
        do k=1,10
          call pskip(kp)
          read(kp,*) (scoef(k,j,i),j=1,3)
        end do
      end do
      close(kp)
c
c...  loop over fault entries
c
      open(kn,file=fnfile,status="old")
      open(ksi,file=sifile,status="old")
      open(kso,file=sofile,status="new")
 30   continue
        read(ksi,*,end=40) elems,nodes,hists,ux,uy,uz
        read(kn,*,end=40) elemc,nodec,isidec,(x(j),j=1,3)
        call getvals(x,ux,uy,uz,scoef,smin,smax,eps,nsurf,insurf)
        if(insurf) write(kso,"(3i7,3(2x,1pe15.8))") elems,nodes,hists,
     &   ux,uy,uz
        go to 30
 40   continue
      close(ksi)
      close(kso)
      close(kn)
      stop
      end
c
c
      function compd2(x,scoef)
c
c...  function to compute the value of an equation of second degree
c
      implicit none
      double precision compd2
      double precision x(3),scoef(10)
c
c...  compute value for given point and set of coefficients
c
      compd2=      scoef(1)*x(1)*x(1)+2.0d0*scoef(2)*x(1)*x(2)+
     &             scoef(3)*x(2)*x(2)+2.0d0*scoef(4)*x(1)*x(3)+
     &       2.0d0*scoef(5)*x(2)*x(3)+      scoef(6)*x(3)*x(3)+
     &       2.0d0*scoef(7)*x(1)     +2.0d0*scoef(8)*x(2)     +
     &       2.0e0*scoef(9)*x(3)     +      scoef(10)
      return
      end
c
c
      subroutine files(pfile,fnfile,sifile,sofile)
c
c...  subroutine to get filenames from command-line
c
      implicit none
      character pfile*(*),fnfile*(*),sifile*(*),sofile*(*)
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
      integer kti,kto,kp,kn,ksi,kso
      common/units/kti,kto,kp,kn,ksi,kso
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
      kn=11
      ksi=12
      kso=14
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
        else if(index(string,'n=').ne.0) then
          j=nchar(string)
          fnfile=string(3:j)
          fflag(2)=.true.
        else if(index(string,'i=').ne.0) then
          j=nchar(string)
          sifile=string(3:j)
          fflag(3)=.true.
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          sofile=string(3:j)
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
     & "    faultbm p=<parameter_file> n=<split_node_coord_file>",/,
     & "    i=<split_node_input_file> o=<split_node_output_file>")
      return
      end
c
c
      subroutine getvals(x,ux,uy,uz,scoef,smin,smax,eps,nsurf,insurf)
c
c...  subroutine to get x, y, and z dislocation values and return them
c     in ux, uy, and uz.  The initial signs of ux, uy, and uz are used
c     to determine the sign of the returned values.
c
      implicit none
c
c...  subroutine arguments
c
      integer nsurf
      double precision x(3),ux,uy,uz,scoef(10,3,nsurf)
      double precision smin(3,nsurf),smax(3,nsurf),eps
      logical insurf
c
c...  intrinsic functions
c
      intrinsic sign,sqrt
c
c...  external functions
c
      double precision compd2
      external compd2
c
c...  local variables
c
      integer i,j
      double precision xsign,ysign,zsign,u,uxt,uyt,uzt,ut
c
c...  loop over surfaces to see if point falls in any of them
c
      insurf=.false.
      u=1.0d30
      do i=1,nsurf
        do j=1,3
          if(x(j).lt.smin(j,i).or.x(j).gt.smax(j,i)) go to 10
        end do
        insurf=.true.
        xsign=sign(1.0d0,ux)
        ysign=sign(1.0d0,uy)
        zsign=sign(1.0d0,uz)
        uxt=xsign*compd2(x,scoef(1,1,i))
        uyt=ysign*compd2(x,scoef(1,2,i))
        uzt=zsign*compd2(x,scoef(1,3,i))
        ut=sqrt(ux*ux+uy*uy+uz*uz)
        if(ut.lt.u) then
          ux=uxt
          uy=uyt
          uz=uzt
          u=ut
        end if
        if(u.lt.eps) insurf=.false.
        return
 10     continue
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
