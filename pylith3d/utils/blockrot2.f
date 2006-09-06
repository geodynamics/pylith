      program blockrot2
c
c...  This code computes split node boundary conditions due to the
c     rotation of one or more blocks about local Cartesian Euler poles.
c     The required inputs are:
c     1. A parameter file describing block and fault information.
c     2. A UCD file (as produced by LaGriT) describing the coordinates
c        and connectivities for a finite element mesh.
c     3. An auxiliary file (also produced by LaGriT) that describes
c        elements and associated faces lying on faults, including
c        fault normal information.
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
      integer maxblocks,maxfaults,maxnodes,maxelems,maxentries
      integer nelemnodes,nsd
      double precision eps
      parameter(maxblocks=1000,maxfaults=10000,maxnodes=2000000,
     & maxelems=10000000,maxentries=1000000,nelemnodes=4,nsd=3,
     & eps=1.0d-3)
c
c...  global arrays
c
      integer ien(nelemnodes,maxelems),nfault(4,maxentries)
      double precision x(nsd,maxnodes),split(nsd,maxentries)
      integer iface(3,4)
      data iface/2,3,4,
     &           1,4,3,
     &           1,2,4,
     &           1,3,2/
c
c...  info read or deduced from parameter file
c
      integer nblocks,nfaults,faultdef(3,maxfaults)
      double precision cscale,dipcut,pole(3,maxblocks)
      logical dflag
c
c...  intrinsic functions
c
      intrinsic acos,sin,abs
c
c...  filenames and unit numbers
c
      integer kti,kto,kp,ku,ka,kso
      common/units/kti,kto,kp,ku,ka,kso
c
c...  local variables
c
      integer i,j,k,n,mat
      integer elem,face,node,elemc,oelemc
      integer nen,nnodes,nelems,nnattr,neattr,nmattr
      integer nentries,blk1,blk2,faultnum,entrynum
      double precision pi,d2r
      double precision xtmp(nsd),fnorm(nsd,3)
      character etype*3
      logical pairused,highrank
c
      pi=acos(-1.0d0)
      d2r=pi/180.0d0
      dflag=.false.
      nen=nelemnodes
c
c...  get filenames and other runtime info
c
      call files()
c
c...  read info from parameter file
c
      call pskip(kp)
c...  global info
      read(kp,*) nblocks,nfaults,cscale,dipcut
      if(dipcut.lt.90.0d0-eps) dflag=.true.
      if(cscale.eq.0.0d0) cscale=1.0d0
c...  block rotation poles
      do i=1,nblocks
        call pskip(kp)
        read(kp,*) (pole(j,i),j=1,3)
        pole(3,i)=d2r*pole(3,i)
      end do
c...  fault rankings/definitions
      do i=1,nfaults
        call pskip(kp)
        read(kp,*) (faultdef(j,i),j=1,3)
      end do
      close(kp)
c
c...  compute cutoff info in terms of projection onto z-axis
c
      if(dflag) dipcut=abs(sin(d2r*dipcut))
c
c...  read info from UCD file
c
      read(ku,*) nnodes,nelems,nnattr,neattr,nmattr
      do i=1,nnodes
        read(ku,*) n,(xtmp(j),j=1,nsd)
        do j=1,nsd
          x(j,i)=cscale*xtmp(j)
        end do
      end do
c
      do i=1,nelems
        read(ku,*) n,mat,etype,(ien(j,i),j=1,nen)
      end do
      close(ku)
c
c...  read fault definitions and compute split node values
c
      nentries=0
      call ifill(nfault,0,4*maxentries)
      do i=1,maxentries
        read(ka,*,end=30) elem,face,elemc,oelemc,
     &   ((fnorm(j,k),j=1,nsd),k=1,3)
        blk1=min(elemc,oelemc)
        blk2=max(elemc,oelemc)
        faultnum=0
        do j=1,nfaults
          if(blk1.eq.faultdef(1,j).and.blk2.eq.faultdef(2,j)) then
            faultnum=j
            go to 10
          end if
        end do
 10     continue
        if(faultnum.eq.0) then
          write(kto,*) "Fault not defined for blocks:",blk1,blk2
          stop
        end if
c
c...  loop over nodes per face and number of entries so far to determine
c     whether element-node pair has already been used and whether the
c     current one is of higher rank.
c
        do j=1,3
          node=ien(iface(j,face),elem)
          pairused=.false.
          highrank=.false.
          do k=1,nentries
            if(elem.eq.nfault(1,k).and.node.eq.nfault(2,k)) then
              pairused=.true.
              if(nfault(4,k).gt.faultnum) highrank=.true.
              entrynum=k
              go to 20
            end if
          end do
          nentries=nentries+1
          entrynum=nentries
 20       continue
          if(.not.pairused.or.highrank) then
            nfault(1,entrynum)=elem
            nfault(2,entrynum)=node
            nfault(3,entrynum)=faultdef(3,faultnum)
            nfault(4,entrynum)=faultnum
            call getvals(x(1,node),fnorm(1,j),dipcut,dflag,
     &       pole(1,elemc),pole(1,oelemc),split(1,entrynum))
          end if
        end do
      end do
      write(kto,*) "Maximum number of entries exceeded!"
      stop
 30   continue
      close(ka)
c
c...  output results to split node file
c
      do i=1,nentries
        write(kso,"(3i8,3(2x,1pe15.8))") (nfault(j,i),j=1,3),
     & (split(j,i),j=1,3)
      end do
      close(kso)
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
      integer kti,kto,kp,ku,ka,kso
      common/units/kti,kto,kp,ku,ka,kso
c
c...  local variables
c
      integer nargs,i,j
      character string*2000
      character pfile*500,ufile*500,afile*500,ofile*500
      logical fflag(4)
c
      kti=5
      kto=6
      kp=10
      ku=11
      ka=12
      kso=13
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
        else if(index(string,'u=').ne.0) then
          j=nchar(string)
          ufile=string(3:j)
          fflag(2)=.true.
        else if(index(string,'a=').ne.0) then
          j=nchar(string)
          afile=string(3:j)
          fflag(3)=.true.
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          ofile=string(3:j)
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
      open(kp,file=pfile,status="old")
      open(ku,file=ufile,status="old")
      open(ka,file=afile,status="old")
      open(kso,file=ofile,status="new")
c
 800  format("Usage:",/,
     & "    blockrot2 p=<parameter_file> u=<UCD_input_file>",/,
     & "    a=<auxiliary_input_file> o=<split_node_output_file>")
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
      double precision uxl,uyl,uzl,uxt,uyt,uzt,uxn,uyn,uzn
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
          fac=dr/sqrt(uxn*uxn+uyn*uyn)
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
