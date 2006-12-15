      program makeucd
c
c...  Simple code to take UCD pieces created by LithoMop/PyLith and
c     create a UCD file for a single timestep.  Command-line options
c     allow a user to:
c       1.  Interpolate Gauss point values to nodes.
c       2.  Average Gauss point values to element centroids.
c       3.  Include nodal values in the output file.
c       4.  Include centroid values in the output file.
c
c     The code is presently set up to deal with either linear tets or
c     linear hexes, and it is assumed that full integration order is
c     used for both (1-point for tets, 2x2x2 for hexes).
c     Present version only reads ASCII files.  This should be updated.
c
      implicit none
c
c...  parameters
c
      integer nsd,maxnodes,maxelmts,maxnen,maxgauss,maxattr,maxeattr
      integer maxoattr
      parameter(nsd=3,maxnodes=1000000,maxelmts=1000000,maxnen=8,
     & maxgauss=8,maxattr=6,maxeattr=24,maxoattr=30)
c
c...  input/output and logical flag parameters
c
      integer kti,kto,km,kp,kg,ko
      common/units/kti,kto,km,kp,kg,ko
      logical passigned,gassigned,ninterp,nout,cout
c
c...  parameters and variables read/deduced from mesh UCD file
c
      integer numnp,numel,nen,ngauss,ngpts
      integer ien(maxnen*maxelmts),mat(maxelmts)
      double precision x(nsd,maxnodes)
c
c...  parameters and variables read/deduced from nodal values file
c
      integer nnattr,isnattr(maxattr),ntnattr
      double precision vnodes(maxattr*maxnodes)
      character nvnames(maxattr)*30
c
c...  parameters and variables read/deduced from element values file
c
      integer neattr,iseattr(maxeattr),nteattr
      double precision velemt(maxeattr*maxelmts*maxgauss)
      character evnames(maxeattr)*30
c
c...  output values
c
      integer nnattrout,neattrout,ntnattrout,nteattrout
      integer isnattrout(maxoattr),iseattrout(maxoattr)
      double precision vnout(maxoattr*maxnodes),vavg(maxoattr)
      character nvnamesout(maxoattr)*30,evnamesout(maxoattr)*30
c
c...  external routines
c
      integer nnblnk,nchar
      external nnblnk,nchar
c
c...  local constants
c
      integer indu(maxnen,2)
      data indu/1, 2, 3, 4, 5, 6, 7, 8,
     &          4, 1, 2, 3, 0, 0, 0, 0/
      integer indp(maxnen,2)
      data indp/1, 2, 3, 4, 5, 6, 7, 8,
     &          2, 3, 4, 1, 0, 0, 0, 0/
c
c...  local variables
c
      integer nsdl
      integer idum,i,j,ind,n,ietype
      integer ientmp(maxnen)
      character tstring*500,etype*3
      double precision smop(maxnodes)
      double precision sh((nsd+1)*maxnen*maxgauss)
      double precision gauss((nsd+1)*maxgauss),xl(nsd*maxnen)
c
      nsdl=nsd
c
c...  get command-line arguments, open files, and assign logical
c     variables
c
      call files(passigned,gassigned,ninterp,nout,cout)
c
c...  read mesh UCD file
c
      read(km,*) numnp,numel,idum,idum,idum
      do i=1,numnp
        read(km,*) n,(x(j,n),j=1,nsd)
      end do
      read(km,"(a500)") tstring
      ietype=0
      if(index(tstring,"hex").ne.0) then
        ietype=1
        nen=8
        ngauss=8
      else if(index(tstring,"tet").ne.0) then
        ietype=2
        nen=4
        ngauss=1
      else
        write(kto,*) "Unknown element type!"
        stop
      end if
      ngpts=ngauss*numel
      backspace(km)
      ind=0
      do i=1,numel
        read(km,*) n,mat(n),etype,(ientmp(j),j=1,nen)
        do j=1,nen
          ien(j+ind)=ientmp(indp(j,ietype))
        end do
        ind=ind+nen
      end do
      close(km)
c
c...  read nodal values file
c
      nnattr=0
      ntnattr=0
      if(passigned) then
        read(kp,*) nnattr,(isnattr(i),i=1,nnattr)
        do i=1,nnattr
          read(kp,"(a30)") nvnames(i)
          ntnattr=ntnattr+isnattr(i)
        end do
        ind=1
        do i=1,numnp
          read(kp,*) n,(vnodes(j),j=ind,ind+ntnattr-1)
          ind=ind+ntnattr
        end do
        close(kp)
      end if
c
c...  read Gauss point values file
c
      neattr=0
      nteattr=0
      if(gassigned) then
        read(kg,*) neattr,(iseattr(i),i=1,neattr)
        do i=1,neattr
          read(kg,"(a30)") evnames(i)
          nteattr=nteattr+iseattr(i)
        end do
        ind=1
        do i=1,ngpts
          read(kg,*) n,(velemt(j),j=ind,ind+nteattr-1)
          ind=ind+nteattr
        end do
        close(kg)
      end if
c
c... determine number of nodal and element values to be output and write
c    to file
c
      nnattrout=0
      neattrout=0
      ntnattrout=0
      nteattrout=0
      if(nout) then
        nnattrout=nnattr
        ntnattrout=ntnattr
        do i=1,nnattr
          isnattrout(i)=isnattr(i)
          nvnamesout(i)=nvnames(i)
        end do
        if(ninterp) then
          do i=1,neattr
            isnattrout(i+nnattrout)=iseattr(i)
            nvnamesout(i+nnattrout)=evnames(i)
          end do
          nnattrout=nnattrout+neattr
          ntnattrout=ntnattrout+nteattr
        end if
      end if
      if(cout) then
        neattrout=neattr
        nteattrout=nteattr
        do i=1,neattr
          iseattrout(i)=iseattr(i)
          evnamesout(i)=evnames(i)
        end do
      end if
      idum=0
      write(ko,"(5i9)") numnp,numel,nnattrout,neattrout,idum
c
c...  output coordinates and connectivities
c
      do i=1,numnp
        write(ko,"(i8,3(2x,1pe15.8))") i,(x(j,i),j=1,nsd)
      end do
      ind=0
      do i=1,numel
        do j=1,nen
          ientmp(j)=ien(ind+j)
        end do
        write(ko,"(2i8,2x,a3,10i9)") i,mat(i),etype,
     &   (ientmp(indu(j,ietype)),j=1,nen)
        ind=ind+nen
      end do
c
c...  output nodal values if requested, interpolating Gauss point
c     values to nodes if desired.
c
      if(nout) call nodout(
     & x,ien,numnp,numel,nsdl,nen,ngauss,ninterp,
     & vnodes,isnattr,nnattr,ntnattr,
     & velemt,iseattr,neattr,nteattr,
     & vnout,isnattrout,nnattrout,ntnattrout,nvnamesout,
     & smop,sh,gauss,xl)
c
c...  output values at element centroids if requested
c
      if(cout) call centout(
     & numel,ngauss,
     & velemt,neattr,nteattr,
     & vavg,iseattrout,neattrout,nteattrout,evnamesout)
c
      close(ko)
c
      stop
      end
c
c
      subroutine abtrans(a,b,c,nra,nca,nrb,lda,ldb,ldc)
c
c...subroutine to multiply the matrix a by b-transpose.  Results
c   are stored in matrix c.
c
      implicit none
      integer nra,nca,nrb,lda,ldb,ldc
      double precision a(lda,nca),b(ldb,nca),c(ldc,nrb)
c
      integer i,j,k
      double precision sum
c
      do i=1,nra
        do j=1,nrb
          sum=0.0d0
          do k=1,nca
            sum=sum+a(i,k)*b(j,k)
          end do
          c(i,j)=sum
        end do
      end do
      return
      end
c
c
      subroutine centout(
     & numel,ngauss,
     & velemt,neattr,nteattr,
     & vavg,iseattrout,neattrout,nteattrout,evnamesout)
c
c...  routine to output values at element centroids
c
      implicit none
      integer numel,ngauss
      integer neattr,nteattr,neattrout,nteattrout
      integer iseattrout(neattrout)
      double precision velemt(nteattr,ngauss,numel)
      double precision vavg(nteattrout)
      character evnamesout(neattrout)*(*)
c
c... i/o units
c
      integer kti,kto,km,kp,kg,ko
      common/units/kti,kto,km,kp,kg,ko
c
c...  local variables
c
      integer i,j,k
      double precision wt
c
c... output cell values header
c
      write(ko,"(30i5)") neattrout,(iseattrout(i),i=1,neattrout)
      do i=1,neattrout
        write(ko,"(a30)") evnamesout(i)
      end do
c
c...  compute values averaged over an element and output them
c
      wt=1.0d0/dble(ngauss)
      do i=1,numel
        call fill(vavg,0.0d0,nteattrout)
        do j=1,ngauss
          do k=1,nteattr
            vavg(k)=vavg(k)+wt*velemt(k,j,i)
          end do
        end do
        write(ko,"(i7,30(2x,1pe15.8))") i,(vavg(j),j=1,nteattrout)
      end do
c
      return
      end
c
c
      subroutine files(passigned,gassigned,ninterp,nout,cout)
c
c...  subroutine to set up i/o and run options
c
      implicit none
      logical ninterp,nout,cout
c
      integer kti,kto,km,kp,kg,ko
      common/units/kti,kto,km,kp,kg,ko
c
      intrinsic index,iargc
c
      integer nnblnk,nchar
      external nnblnk,nchar
c
      character mfile*500,pfile*500,gfile*500,ofile*500,string*500
      logical massigned,passigned,gassigned,oassigned
      integer nargs,i,j,ival
c
c...  unit numbers
c
      kti=5
      kto=6
      km=10
      kp=11
      kg=12
      ko=13
c
c...  logical flags
c
      ninterp=.true.
      nout=.true.
      cout=.false.
      massigned=.false.
      passigned=.false.
      gassigned=.false.
      oassigned=.false.
c
c... loop over command-line arguments
c
      nargs=iargc()
      do i=1,nargs
        call getarg(i,string)
        if(index(string,'m=').ne.0) then
          massigned=.true.
          j=len_trim(string)
          mfile=string(3:j)
        else if(index(string,'p=').ne.0) then
          passigned=.true.
          j=len_trim(string)
          pfile=string(3:j)
        else if(index(string,'g=').ne.0) then
          gassigned=.true.
          j=len_trim(string)
          gfile=string(3:j)
        else if(index(string,'o=').ne.0) then
          oassigned=.true.
          j=len_trim(string)
          ofile=string(3:j)
        else if(index(string,'n=').ne.0) then
          j=len_trim(string)
          read(string(3:j),*) ival
          if(ival.eq.0) then
            ninterp=.false.
          else if(ival.eq.1) then
            ninterp=.true.
          else
            write(kto,800)
            stop
          end if
        else if(index(string,'nf=').ne.0) then
          j=len_trim(string)
          read(string(4:j),*) ival
          if(ival.eq.0) then
            nout=.false.
          else if(ival.eq.1) then
            nout=.true.
          else
            write(kto,800)
            stop
          end if
        else if(index(string,'cf=').ne.0) then
          j=len_trim(string)
          read(string(4:j),*) ival
          if(ival.eq.0) then
            cout=.false.
          else if(ival.eq.1) then
            cout=.true.
          else
            write(kto,800)
            stop
          end if
        end if
      end do
c
c...  adjust flags for inconsistencies
c
      if(.not.nout) ninterp=.false.
c
c...  open files
c
      if(massigned) then
        open(km,file=mfile,status="old")
      else
        write(kto,800)
        stop
      end if
      if(passigned) open(kp,file=pfile,status="old")
      if(gassigned) open(kg,file=gfile,status="old")
      if(oassigned) then
        open(ko,file=ofile,status="new")
      else
        write(kto,800)
        stop
      end if
c
 800  format("Usage:",/,
     & "makeucd m=<mesh_file> o=<output_file> [p=<nodal_value_file>]",/,
     & "[g=<gauss_value_file>] [n=<gauss_interpolation_option]",/,
     & "[nf=<nodal_output_flag>] [cf=<element_output_flag>]",//,
     & "All options have values of 0 (no) or 1 (yes).",/,
     & "Default option values are:",/,
     & "n = 1 (interpolate Gauss values to nodes)",/,
     & "nf = 1 (output values at nodes)",/,
     & "cf = 0 (do not output values at element centroids)")
c
      return
      end
c
c
      subroutine fill(arr,val,len)
c
c...  routine to fill a double-precision array with a given value
c
      implicit none
      integer len
      double precision arr(len),val
c
      integer i
c
      do i=1,len
        arr(i)=val
      end do
      return
      end
c
c
      subroutine getjac(x,xs,det,shj,nsd,nen,iel)
c
c...  subroutine to compute the jacobian determinant given the element
c     coordinates and the shape functions in natural coordinates.
c
c       shj(1,nen),sh(2,nen),sh(3,nen) = x,y,and z derivatives
c                                        of shape functions
c       xs(nsd,nsd)                   = jacobian matrix
c       det                           = determinant of jacobian matrix
c       x(nsd,nen)                    = local nodal coordinates
c
      implicit none
c
c...  subroutine arguments
c
      integer nsd,nen,iel
      double precision x(nsd,nen),xs(nsd,nsd),det,shj(nsd+1,nen)
c
c... i/o units
c
      integer kti,kto,km,kp,kg,ko
      common/units/kti,kto,km,kp,kg,ko
c
c...calculate jacobian matrix for (x,y,z) to (r,s,t) transformation
c
cblas      call dgemm("n","t",nsd,nsd,nen,1.0d0,x,nsd,shj,nsd+1,0.0d0,xs,nsd)
      call abtrans(x,shj,xs,nsd,nen,nsd,nsd,nsd+1,nsd)
c
c...form determinant of jacobian matrix and check for error condition
c
      det=xs(1,1)*xs(2,2)*xs(3,3)+xs(1,2)*xs(2,3)*xs(3,1)+xs(1,3)
     & *xs(2,1)*xs(3,2)-xs(1,3)*xs(2,2)*xs(3,1)-xs(1,2)*xs(2,1)
     & *xs(3,3)-xs(1,1)*xs(2,3)*xs(3,2)
c
      if(det.le.0.0d0) then
        write(kto,700) iel,det
        stop
      end if
c
 700  format("getjac:  element # ",i7,2x,1pe15.8)
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
      subroutine nodout(
     & x,ien,numnp,numel,nsd,nen,ngauss,ninterp,
     & vnodes,isnattr,nnattr,ntnattr,
     & velemt,iseattr,neattr,nteattr,
     & vnout,isnattrout,nnattrout,ntnattrout,nvnamesout,
     & smop,sh,gauss,xl)
c
c...  subroutine to interpolate Gauss point values to nodes
c     (if requested) and output values at nodes.
c
      implicit none
c
c...  subroutine arguments
c
      integer numnp,numel,nsd,nen,ngauss
      integer nnattr,ntnattr
      integer neattr,nteattr
      integer nnattrout,ntnattrout
      integer ien(nen,numel)
      integer isnattr(nnattr),iseattr(neattr),isnattrout(nnattrout)
      double precision x(nsd,numnp)
      double precision vnodes(ntnattr,numnp)
      double precision velemt(nteattr,ngauss,numel)
      double precision vnout(ntnattrout,numnp)
      double precision smop(numnp),sh(nsd+1,nen,ngauss)
      double precision gauss(nsd+1,ngauss),xl(nsd,nen)
      logical ninterp
      character nvnamesout(nnattrout)*(*)
c
c... i/o units
c
      integer kti,kto,km,kp,kg,ko
      common/units/kti,kto,km,kp,kg,ko
c
c...  local variables
c
      integer i,j,k,l,node
      double precision xs(3,3),det
c
      if(ninterp) then
        if(nen.eq.4) then
          call plintet(sh,gauss,nsd,nen,ngauss)
        else
          call plinhex(sh,gauss,nsd,nen,ngauss)
        end if
c
c...  loop over elements and compute smoothing operator and weighted
c     state variables at nodes.
c
        call fill(smop,0.0d0,numnp)
        do i=1,numel
c
c...  localize nodal coordinates for element
c
          call lcoord(ien(1,i),x,xl,nsd,nen,numnp)
c
c...  loop over Gauss points
c
          do j=1,ngauss
            call getjac(xl,xs,det,sh,nsd,nen,i)
            do k=1,nen
              node=ien(k,i)
              smop(node)=smop(node)+det
              do l=1,nteattr
                vnout(l,node)=vnout(l,node)+velemt(l,j,i)*det
              end do
            end do
          end do
        end do
c
c...  invert smoothing operator and compute smoothed state variables
c
        do i=1,numnp
          smop(i)=1.0d0/smop(i)
          do j=1,nteattr
            vnout(j,i)=vnout(j,i)*smop(i)
          end do
        end do
      end if
c
c... output nodal values header
c
      write(ko,"(31i5)") nnattrout,(isnattrout(i),i=1,nnattrout)
      do i=1,nnattrout
        write(ko,"(a30)") nvnamesout(i)
      end do
c
c...  write nodal values
c
      if(.not.ninterp) then
        do i=1,numnp
          write(ko,"(i7,6(2x,1pe15.8))") i,(vnodes(j,i),j=1,ntnattr)
        end do
      else
        do i=1,numnp
          write(ko,"(i7,30(2x,1pe15.8))") i,(vnodes(j,i),j=1,ntnattr),
     &     (vnout(j,i),j=1,nteattr)
        end do
      end if
      return
      end
c
c
      subroutine plinhex(sh,gauss,nsd,nen,ngauss)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a trilinear hexahedron.
c
      implicit none
c
c...  subroutine arguments
c
      integer nsd,nen,ngauss
      double precision sh(nsd+1,nen,ngauss)
      double precision gauss(nsd+1,ngauss)
c
c...  local constants
c
      double precision r(8),s(8),t(8)
      data r/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0/
      data s/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0/
      data t/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0/
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l
      double precision rr,ss,tt,root3i,eighth,one
c
c...  definitions
c
      one=1.0d0
      root3i=sqrt(1.0d0/3.0d0)
      eighth=1.0d0/8.0d0
c
c...  Linear hex definition
c
      do l=1,ngauss
        gauss(1,l)=r(l)*root3i
        gauss(2,l)=s(l)*root3i
        gauss(3,l)=t(l)*root3i
        gauss(4,l)=one
      end do
c
      do l=1,ngauss
        do i=1,nen
          rr=one+r(i)*gauss(1,l)
          ss=one+s(i)*gauss(2,l)
          tt=one+t(i)*gauss(3,l)
          sh(4,i,l)=eighth*rr*ss*tt
          sh(1,i,l)=eighth*r(i)*ss*tt
          sh(2,i,l)=eighth*s(i)*rr*tt
          sh(3,i,l)=eighth*t(i)*rr*ss
        end do
      end do
c
      return
      end
c
c
      subroutine plintet(sh,gauss,nsd,nen,ngauss)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a linear tetrahedron.
c
      implicit none
c
c...  subroutine arguments
c
      integer nsd,nen,ngauss
      double precision sh(nsd+1,nen,ngauss)
      double precision gauss(nsd+1,ngauss)
c
c...  local constants
c
      double precision r(4),s(4),t(4),u(4)
      data u/ 1d0, 0d0, 0d0, 0d0/
      data r/ 0d0, 1d0, 0d0, 0d0/
      data s/ 0d0, 0d0, 1d0, 0d0/
      data t/ 0d0, 0d0, 0d0, 1d0/
c
c...  intrinsic functions
c
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i
      double precision rr,ss,tt,uu
      double precision tetvol,fourth
c
c...  definitions
c
      tetvol=1.0d0/6.0d0
      fourth=0.25d0
c
c...  Linear hex definition
c     One-point integration is used in all cases.
c
      gauss(1,1)=fourth
      gauss(2,1)=fourth
      gauss(3,1)=fourth
      gauss(4,1)=tetvol
c
      do i=1,nen
        rr=r(i)*gauss(1,1)
        ss=s(i)*gauss(2,1)
        tt=t(i)*gauss(3,1)
        uu=u(i)*gauss(3,1)
        sh(4,i,1)=rr+ss+tt+uu
        sh(1,i,1)=r(i)-u(i)
        sh(2,i,1)=s(i)-u(i)
        sh(3,i,1)=t(i)-u(i)
      end do
c
      return
      end
c
c version
c $Id: plintet.f,v 1.4 2005/03/22 04:45:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
