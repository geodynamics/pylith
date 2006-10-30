      program blockrot3
c
c...  This code computes split node boundary conditions due to the
c     rotation of one or more blocks about local Cartesian Euler poles.
c     The required inputs are:
c     1. A parameter file describing block and fault information.
c     2. Fault files (one per fault, produced by LaGriT) that
c        describe element/node pairs lying on faults.
c     The resulting output is a split node specification file suitable
c     for use by PyLith or LithoMop.
c
c     Code is presently just set up for linear tets, but could easily be
c     made more general.
c
      implicit none
c
c...  parameters
c
      integer maxfaults,maxblocks,maxentries,maxattached,maxdefs
      integer nsd
      double precision eps
      parameter(maxfaults=10000,maxblocks=100,maxentries=10000000,
     & maxattached=2000,maxdefs=10,nsd=3,eps=1.0d-3)
c
c...  global arrays
c
      integer nfault(2,maxentries)
      double precision split(nsd)
c
c...  info read or deduced from parameter file
c
      integer nfaults,numdefs,faultdef(3,maxdefs),blocknum(maxblocks)
      double precision cscale,dipcut,pole(3,maxblocks)
      logical dflag
c
c...  intrinsic functions
c
      intrinsic acos,sin,abs
c
c...  filenames and unit numbers
c
      integer kti,kto,kp,kf,kso
      common/units/kti,kto,kp,kf,kso
c
c...  local variables
c
      integer i,j,k
      integer node
      integer nblocks,nentries,nattached,blk1,blk2,faultnum,ind
      integer elems(maxattached),colors(maxattached)
      double precision pi,d2r
      double precision x(nsd),xtmp(nsd),fnorm(nsd)
      logical pairused
      character ffile*500
c
      pi=acos(-1.0d0)
      d2r=pi/180.0d0
      dflag=.false.
c
c...  get filenames and other runtime info
c
      call files()
c
c...  read info from parameter file
c
      call pskip(kp)
c
c...  global info
c
      read(kp,*) nblocks,nfaults,cscale,dipcut
      if(dipcut.lt.90.0d0-eps) dflag=.true.
      if(cscale.eq.0.0d0) cscale=1.0d0
c
c...  compute cutoff info in terms of projection onto z-axis
c
      if(dflag) dipcut=abs(sin(d2r*dipcut))
c
c...  block rotation poles
c
      do i=1,nblocks
        call pskip(kp)
        read(kp,*) ind,(pole(j,i),j=1,3)
        blocknum(ind)=i
        pole(3,i)=d2r*pole(3,i)
      end do
c
c...  fault definition files
c
      nentries=0
      call ifill(nfault,0,2*maxentries)
      do i=1,nfaults
        call pskip(kp)
        read(kp,*) numdefs
        do j=1,numdefs
          read(kp,*) (faultdef(k,j),k=1,3)
        end do
        read(kp,*) ffile
        open(kf,file=ffile,status="old")
c
c...  loop over entries in file
c
        call pskip(kf)
 10     continue
          read(kf,*,end=20) node,(xtmp(j),j=1,3),(fnorm(j),j=1,3)
          read(kf,*,end=20) nattached
          read(kf,*,end=20) (elems(j),j=1,nattached)
          read(kf,*,end=20) (colors(j),j=1,nattached)
          do j=1,3
            x(j)=cscale*xtmp(j)
          end do
c
c...  see if pair has been used previously.  If not, create a new entry
c     if the fault has been defined.
c
          do j=1,nattached
            pairused=.false.
            do k=1,nentries
              if(elems(j).eq.nfault(1,k).and.node.eq.nfault(2,k)) then
                pairused=.true.
                go to 30
              end if
            end do
 30         continue
            if(.not.pairused) then
c
c...  for now, choose first fault in definition list
c
              blk1=0
              blk2=0
              faultnum=0
              do k=1,numdefs
                if(colors(j).eq.faultdef(1,k)) then
                  blk1=colors(j)
                  blk2=faultdef(2,k)
                  faultnum=k
                  go to 40
                else if(colors(j).eq.faultdef(2,k)) then
                  blk1=colors(j)
                  blk2=faultdef(1,k)
                  faultnum=k
                  go to 40
                end if
              end do
 40           continue
              if(blk1.eq.0) then
                write(kto,810) elems(j),node,colors(j)
              else
                nentries=nentries+1
                nfault(1,nentries)=elems(j)
                nfault(2,nentries)=node
                call getvals(x,fnorm,dipcut,dflag,
     &           pole(1,blocknum(blk1)),pole(1,blocknum(blk2)),split)
c
c...  output split node values
c
                write(kso,"(3i8,3(2x,1pe15.8))")
     &           (nfault(k,nentries),k=1,2),faultdef(3,faultnum),
     &           (split(k),k=1,3)
              end if
            end if
          end do
          go to 10
 20     continue
        close(kf)
      end do
      close(kp)
      close(kso)
 810  format("WARNING!  No appropriate fault definition found for:",/,
     &       "Element:  ",i8,/,
     &       "Node:     ",i8,/,
     &       "Color:    ",i8)
c
      stop
      end
c
c
      subroutine files()
c
c...  subroutine to get filenames and other runtime info from
c     command-line
c
c
c...  routine arguments
c
      implicit none
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
      integer kti,kto,kp,kf,kso
      common/units/kti,kto,kp,kf,kso
c
c...  local variables
c
      integer nargs,i,j
      character string*2000
      character pfile*500,ofile*500
      logical fflag(2)
c
      kti=5
      kto=6
      kp=10
      kf=12
      kso=13
c
      do i=1,2
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
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          ofile=string(3:j)
          fflag(2)=.true.
        end if
      end do
c
      do i=1,2
        if(.not.fflag(i)) then
          write(kto,800)
          stop
        end if
      end do
c
      open(kp,file=pfile,status="old")
      open(kso,file=ofile,status="new")
c
 800  format("Usage:",/,
     & "    blockrot2 p=<parameter_file>",/,
     & "    o=<split_node_output_file>")
      return
      end
c
c
      subroutine getvals(x,fnorm,dipcut,dflag,pole1,pole2,split)
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
      double precision x(3),fnorm(3),pole1(3),pole2(3),dipcut
      double precision split(3)
      logical dflag
c
c...  intrinsic functions
c
      intrinsic sqrt,abs,sin,acos
c
c...  local variables
c
      double precision uxl,uyl,uzl,uxt,uyt,uzt,uxn,uyn,uzn,dnm
      double precision dx1,dy1,dx2,dy2,dr,dfault,ssp,dsp,fac
      double precision ss(3),ds(3)
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
        if(dfault.lt.dipcut) then
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
          dnm=sqrt(uxn*uxn+uyn*uyn)
          if(dnm.gt.0.0d0) then
            fac=dr/dnm
          else
            fac=0.0d0
          end if
          uxl=uxt+uxn*fac
          uyl=uyt+uyn*fac
          uzl=uzt+uzn*fac
        end if
      end if
c
c...  final slip is divided by 2 since half will be applied to each side.
c     we no longer worry about fault 'sign', since this should be taken
c     taken care of by the normal direction.
c
      split(1)=uxl*0.5d0
      split(2)=uyl*0.5d0
      split(3)=uzl*0.5d0
      return
      end
c
c
      subroutine ifill(iarr,ival,ilen)
c
c...  subroutine to fill an integer array with a given value
c
      implicit none
c
c...  subroutine arguments
c
      integer ival,ilen
      integer iarr(ilen)
c
c...  local variables
c
      integer i
c
      do i=1,ilen
        iarr(i)=ival
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
