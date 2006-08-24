      program blockrot
c
c...  quick and dirty program to compute fault boundary conditions,
c     given a set of points and a rotation pole that specifies the
c     movement of the 'positive' side of the fault with respect to
c     the 'negative' side.  Positive and negative sides are determined
c     by the side number (1=positive) in the split node coordinate file.
c     The user must specify the names of the initial split node file
c     as well as a file containing split node coordinates and possibly
c     fault normal vectors (as produced by readucd2, for example).
c     The user must also specify values for x-position, y-position,
c     and rotation (degrees CCW) for the rotation pole, and may
c     optionally specify a coordinate scaling value and a dip-slip
c     cutoff value (degrees from vertical).
c     This program is suitable for specifying split node BC for the
c     LithoMop or PyLith codes.
c
c     This version of the code assumes that output has already been
c     produced by a code such as readnetgen, readucd, or readucd2,
c     such that fault information files (*.fcoord, *.fbc) have already
c     been produced.
c
      implicit none
c
c...  parameters
c
c
c...  global arrays
c
      double precision dval,scale,pole1(3),pole2(3)
      logical dflag,hasnorm
c
c...  intrinsic functions
c
      intrinsic acos,sin,abs
c
c...  filenames and unit numbers
c
      integer kti,kto,ksi,ksic,kso
      common/units/kti,kto,ksi,ksic,kso
      character sifile*200,sicfile*200,sofile*200
c
c...  local variables
c
      integer j,nodes,elems,hists,nodec,elemc,isidec
      double precision ux,uy,uz,pi,d2r
      double precision x(3),fnorm(3)
c
      pi=acos(-1.0d0)
      d2r=pi/180.0d0
c
c...  get filenames and other runtime info
c
      call files(sifile,sicfile,sofile,dflag,dval,scale,pole1,pole2)
c
c...  compute cutoff info in terms of projection onto z-axis
c
      if(dflag) dval=abs(sin(d2r*dval))
c
c...  determine whether fault normals are included
c
      open(ksi,file=sifile,status="old")
      open(ksic,file=sicfile,status="old")
      open(kso,file=sofile,status="new")
      hasnorm=.false.
      read(ksic,*,err=10) elemc,nodec,isidec,(x(j),j=1,3),
     & (fnorm(j),j=1,3)
      hasnorm=.true.
 10   continue
      backspace(ksic)
c
c...  loop over fault entries
c
      if(hasnorm) then
 20     continue
          read(ksi,*,end=30) elems,nodes,hists,ux,uy,uz
          read(ksic,*,end=30) elemc,nodec,isidec,(x(j),j=1,3),
     &     (fnorm(j),j=1,3)
          call getvals(x,fnorm,isidec,scale,dval,dflag,pole1,pole2,
     &     ux,uy,uz)
          write(kso,"(3i7,3(2x,1pe15.8))") elems,nodes,hists,ux,uy,uz
          go to 20
 30     continue
      else if(.not.dflag) then
 40     continue
          read(ksi,*,end=50) elems,nodes,hists,ux,uy,uz
          read(ksic,*,end=50) elemc,nodec,isidec,(x(j),j=1,3)
          call getvals(x,fnorm,isidec,scale,dval,dflag,pole1,pole2,
     &     ux,uy,uz)
          write(kso,"(3i7,3(2x,1pe15.8))") elems,nodes,hists,ux,uy,uz
          go to 40
 50     continue
      else
        write(kto,800)
        stop
      end if
      close(ksi)
      close(ksic)
      close(kso)
 800  format("ERROR!",/,
     &       "Fault normals must be present to use dip-slip cutoff!")
      stop
      end
c
c
      subroutine files(sifile,sicfile,sofile,dflag,dval,scale,pole1,
     & pole2)
c
c...  subroutine to get filenames and other runtime info from
c     command-line
c
c
c...  routine arguments
c
      implicit none
      double precision dval,scale,pole1(3),pole2(3)
      character sifile*(*),sicfile*(*),sofile*(*)
      logical dflag
c
c...  intrinsic functions
c
      intrinsic iargc,index,acos
c
c...  external functions
c
      integer nchar
      external nchar
c
c...  unit numbers
c
      integer kti,kto,ksi,ksic,kso
      common/units/kti,kto,ksi,ksic,kso
c
c...  local variables
c
      integer nargs,i,j
      double precision pi,d2r
      character string*2000
      logical fflag(9)
c
      kti=5
      kto=6
      ksi=10
      ksic=11
      kso=12
c
      pi=acos(-1.0d0)
      d2r=pi/180.0d0
c
      dflag=.false.
      dval=0.0d0
      scale=1.0d0
c
      do i=1,9
        fflag(i)=.false.
      end do
c
      nargs=iargc()
      do i=1,nargs
        call getarg(i,string)
        if(index(string,'i=').ne.0) then
          j=nchar(string)
          sifile=string(3:j)
          fflag(1)=.true.
        else if(index(string,'c=').ne.0) then
          j=nchar(string)
          sicfile=string(3:j)
          fflag(2)=.true.
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          sofile=string(3:j)
          fflag(3)=.true.
        else if(index(string,'x1=').ne.0) then
          j=nchar(string)
          read(string(4:j),*) pole1(1)
          fflag(4)=.true.
        else if(index(string,'y1=').ne.0) then
          j=nchar(string)
          read(string(4:j),*) pole1(2)
          fflag(5)=.true.
        else if(index(string,'r1=').ne.0) then
          j=nchar(string)
          read(string(4:j),*) pole1(3)
          fflag(6)=.true.
          pole1(3)=d2r*pole1(3)
        else if(index(string,'x2=').ne.0) then
          j=nchar(string)
          read(string(4:j),*) pole2(1)
          fflag(7)=.true.
        else if(index(string,'y2=').ne.0) then
          j=nchar(string)
          read(string(4:j),*) pole2(2)
          fflag(8)=.true.
        else if(index(string,'r2=').ne.0) then
          j=nchar(string)
          read(string(4:j),*) pole2(3)
          fflag(9)=.true.
          pole2(3)=d2r*pole2(3)
        else if(index(string,'d=').ne.0) then
          j=nchar(string)
          read(string(3:j),*) dval
          dflag=.true.
        else if(index(string,'s=').ne.0) then
          j=nchar(string)
          read(string(3:j),*) scale
        end if
      end do
c
      do i=1,9
        if(.not.fflag(i)) then
          write(kto,800)
          stop
        end if
      end do
c
 800  format("Usage:",/,
     & "    blockrot i=<split_node_input_file>",/,
     & "    c=<split_node_coord_file> o=<split_node_output_file>",/,
     & "    x1=<rotation_pole1_xcoord> y1=<rotation_pole1_ycoord>",/,
     & "    r1=<rotation1_amount> x2=<rotation_pole2_xcoord>",/,
     & "    y2=<rotation_pole2_ycoord> r2=<rotation2_amount>",/,
     & "    [d=<dip_slip_cutoff> s=<scale_factor>]",//,
     & "    Block 1 is the reference block, and rotations are in",/,
     & "    degrees.")
      return
      end
c
c
      subroutine getvals(x,fnorm,isidec,scale,dval,dflag,pole1,pole2,
     & ux,uy,uz)
c
c...  subroutine to get x, y, and z dislocation values and return them
c     in ux, uy, and uz.  The isidec value is used to determine the sign
c     of the returned values.
c
      implicit none
c
c...  parameters
c
      double precision eps
      parameter(eps=1.0d-3)
c
c...  subroutine arguments
c
      integer isidec
      double precision x(3),fnorm(3),pole1(3),pole2(3),scale,dval
      double precision ux,uy,uz
      logical dflag
c
c...  intrinsic functions
c
      intrinsic sqrt,abs,sin,acos
c
c...  local variables
c
      double precision uxl,uyl,uzl,uxt,uyt,uzt,uxn,uyn,uzn
      double precision sgn,dx1,dy1,dx2,dy2,dr,dfault,ssp,dsp,fac
      double precision ss(3),ds(3)
c
      sgn=1.0d0
      if(isidec.eq.2) sgn=-1.0d0
c
c...  Compute rotation.
c
      dx1=x(1)-pole1(1)
      dy1=x(2)-pole1(2)
      dx2=x(1)-pole2(1)
      dy2=x(2)-pole2(2)
      uxl=-pole1(3)*dy1+pole2(3)*dy2
      uyl=pole1(3)*dx1-pole2(3)*dx2
      uzl=0.0d0
c
c...  if cutoff values are being used, compute strike-slip and dip-slip
c     components.
c
      if(dflag) then
c
c...  make sure normal is normalized
c
        dr=sqrt(fnorm(1)*fnorm(1)+fnorm(2)*fnorm(2)+fnorm(3)*fnorm(3))
        fnorm(1)=fnorm(1)/dr
        fnorm(2)=fnorm(2)/dr
        fnorm(3)=fnorm(3)/dr
        dfault=sqrt(1.0d0-fnorm(3)*fnorm(3))
        if(dfault.lt.dval) then
c
c...  strike-slip unit vector is along-strike (khat X fnorm) unless the
c     fault is horizontal, in which case it is E-W (along x-direction).
c     Dipslip is cross product of normal with strike-slip, which should
c     be updip.
c
          if(abs(fnorm(3)).lt.(1.0d0-eps)) then
            dr=sqrt(fnorm(1)*fnorm(1)+fnorm(2)*fnorm(2))
            ss(1)=-fnorm(2)/dr
            ss(2)=fnorm(1)/dr
            ss(3)=0.0d0
          else
            ss(1)=1.0d0
            ss(2)=0.0d0
            ss(3)=0.0d0
          end if
          ds(1)=fnorm(2)*ss(3)-fnorm(3)*ss(2)
          ds(2)=fnorm(3)*ss(1)-fnorm(1)*ss(3)
          ds(3)=fnorm(1)*ss(2)-fnorm(2)*ss(1)
c
c...  project slip onto strike-slip direction and add contribution to total slip.
c
          ssp=uxl*ss(1)+uyl*ss(2)+uzl*ss(3)
          uxt=ssp*ss(1)
          uyt=ssp*ss(2)
          uzt=ssp*ss(3)
c
c...  compute normal components of slip, and compute dot product with
c     dipslip.  Scale vector so that horizontal components are equal
c     to block-normal motion.
c
          uxn=uxl-uxt
          uyn=uyl-uyt
          uzn=uzl-uzt
          dr=sqrt(uxn*uxn+uyn*uyn)
          dsp=uxn*ds(1)+uyn*ds(2)+uzn*ds(3)
          uxn=dsp*ds(1)
          uyn=dsp*ds(2)
          uzn=dsp*ds(3)
          fac=dr/sqrt(uxn*uxn+uyn*uyn)
          uxl=uxt+uxn*fac
          uyl=uyt+uyn*fac
          uzl=uzt+uzn*fac
        end if
      end if
c
c...  final slip is multiplied by sign corresponding to fault side, and
c     divided by 2 since half will be applied to each side.
c
      ux=sgn*uxl*0.5d0
      uy=sgn*uyl*0.5d0
      uz=sgn*uzl*0.5d0
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
