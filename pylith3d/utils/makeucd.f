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
      parameter(nsd=3,maxnodes=500000,maxelmts=500000,maxnen=8,
     & maxgauss=8,maxattr=6,maxeattr=24,maxoattr=30)
c
c...  input/output and logical flag parameters
c
      integer kti,kto,km,kp,kg,ko
      common/kti,kto,km,kp,kg,ko/
      logical passigned,gassigned,ninterp,gavg,nout,cout
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
      double precision vnodes(maxattr*naxnodes)
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
      double precision vnout(maxoattr*maxnodes),veout(maxoattr*maxelmts)
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
      integer idum,i,nen,ngauss,ind
      integer ientmp(maxnen)
      character tstring*500,etype*3
      double precision sh((nsd+1)*maxnen*maxgauss)
      double precision gauss((nsd+1)*maxgauss)
c
c...  get command-line arguments, open files, and assign logical
c     variables
c
      call files(passigned,gassigned,ninterp,gavg,nout,cout)
c
c...  read mesh UCD file
c
      read(km,*) numnp,numel,idum,idum,idum
      do i=1,numnp
        read(km,*) n,(x(j,n),j=1,nsd)
      end do
      read(km,*) tstring
      if(index(tstring,"hex").ne.0) then
        nen=8
        ngauss=8
      else if(index(tstring,"tet").ne.0) then
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
          ien(j+ind)=ientmp(indp(j))
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
        if(ninterp) nnattrout=nnattrout+neattr
      end if
      if(gavg) neattrout=neattr
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
     &   (ientmp(indu(j)),j=1,nen)
        ind=ind+nen
      end do

      nenl=nen
      nsdl=nsd
      write(kto,*) "Enter root name for all files.  Both input and"
      write(kto,*) "output files will all have this prefix:"
      read(kti,"(a200)") fileroot
      i1=nnblnk(fileroot)
      i2=nchar(fileroot)
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
c...  box dimensions and bc
c
      call pskip(kr)
      read(kr,*) ((xlim(j,i),j=1,2),i=1,nsd)
      do i=1,nsides
        call pskip(kr)
        read(kr,*) (ibc(j,i),bc(j,i),j=1,nsd)
      end do
c...  connectivity order option
c
      call pskip(kr)
      read(kr,*) iconopt
c
c...  fault definitions
c
      call pskip(kr)
      read(kr,*) numflt,ifattrn,ifattrv
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
c...  read and output nodal coordinates and bc info
c
      ifile=fileroot(i1:i2)//".inp"
      open(file=ifile,unit=kr,status="old")
      nfile=fileroot(i1:i2)//".coord"
      open(file=nfile,unit=kw,status="new")
      do i=1,nsides
        j1=nnblnk(cside(i))
        j2=nchar(cside(i))
        bcfile=fileroot(i1:i2)//"."//cside(i)(j1:j2)//".bc"
        bccfile=fileroot(i1:i2)//"."//cside(i)(j1:j2)//".coord"
        open(file=bcfile,unit=kwb(i),status="new")
        open(file=bccfile,unit=kwc(i),status="new")
      end do
      write(kw,"(a15)") cstring
      call pskip(kr)
      read(kr,*) numnp,numel,nnattr,neattr,nmattr
      do i=1,numnp
        read(kr,*) n,(xtmp(j),j=1,nsd)
        do j=1,nsd
          x(j,i)=cscale*xtmp(idir(j))
        end do
        write(kw,"(i7,3(2x,1pe15.8))") i,(x(j,i),j=1,nsd)
        jj=0
        do j=1,nsd
          do k=1,2
            jj=jj+1
            dr=abs(x(j,i)-xlim(k,j))
            if(dr.lt.eps) then
              write(kwb(jj),"(i7,3i4,3(2x,1pe15.8))") i,
     &         (ibc(l,jj),l=1,nsd),(bc(l,jj),l=1,nsd)
              write(kwc(jj),"(i7,3(2x,1pe15.8))") i,
     &         (x(l,i),l=1,nsd)
            end if
          end do
        end do
      end do
      close(kw)
      do i=1,nsides
        close(kwb(i))
        close(kwc(i))
      end do
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
c...  read nodal attributes to determine which nodes are associated
c     with each fault
c
      read(kr,*) nf,(itmp(i),i=1,nf)
      do i=1,nf
        read(kr,*) descr
      end do
      nfltnodes=0
      do i=1,numnp
        read(kr,*) n,(fattr(j),j=1,nnattr)
        iattr=nint(fattr(ifattrn))
        if(iattr.eq.ifattrv) then
          nfltnodes=nfltnodes+1
          inodef(nfltnodes)=i
          nelsf(nfltnodes)=0
        end if
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
      do i=1,nsides
        close(kwfb(i))
        close(kwfc(i))
      end do
      write(kto,700) nflip
700   format("Number of connectivities flipped:  ",i7)
      stop
      end
c
c
      subroutine files(passigned,gassigned,ninterp,gavg,nout,cout)
c
c...  subroutine to set up i/o and run options
c
      implicit none
      logical ninterp,gavg,nout,cout
c
      integer kti,kto,km,kp,kg,ko
      common/kti,kto,km,kp,kg,ko/
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
      gavg=.false.
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
          j=lnblnk(string)
          mfile=string(3:j)
        else if(index(string,'p=').ne.0) then
          passigned=.true.
          j=lnblnk(string)
          pfile=string(3:j)
        else if(index(string,'g=').ne.0) then
          gassigned=.true.
          j=lnblnk(string)
          gfile=string(3:j)
        else if(index(string,'o=').ne.0) then
          oassigned=.true.
          j=lnblnk(string)
          ofile=string(3:j)
        else if(index(string,'n=').ne.0) then
          j=lnblnk(string)
          read(string(3:j),*) ival
          if(ival.eq.0) then
            ninterp=.false.
          else if(ival.eq.1) then
            ninterp=.true.
          else
            write(kto,800)
            end
          end if
        else if(index(string,'a=').ne.0) then
          j=lnblnk(string)
          read(string(3:j),*) ival
          if(ival.eq.0) then
            gavg=.false.
          else if(ival.eq.1) then
            gavg=.true.
          else
            write(kto,800)
            end
          end if
        else if(index(string,'no=').ne.0) then
          j=lnblnk(string)
          read(string(4:j),*) ival
          if(ival.eq.0) then
            nout=.false.
          else if(ival.eq.1) then
            nout=.true.
          else
            write(kto,800)
            end
          end if
        else if(index(string,'co=').ne.0) then
          j=lnblnk(string)
          read(string(4:j),*) ival
          if(ival.eq.0) then
            cout=.false.
          else if(ival.eq.1) then
            cout=.true.
          else
            write(kto,800)
            end
          end if
        end if
      end do
c
c...  adjust flags for inconsistencies
c
      if(.not.nout) ninterp=.false.
      if(.not.cout) gavg=.false.
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
     & "[a=<gauss_averaging_option>] [no=<nodal_output_option>]",/,
     & "[co=<element_output_option>]",//,
     & "All options have values of 0 (no) or 1 (yes).",/,
     & "Default option values are:",/,
     & "n = 1 (interpolate Gauss values to nodes)",/,
     & "a = 0 (do not average Gauss values to element centroids)",/,
     & "no = 1 (output values at nodes)",/,
     & "co = 0 (do not output values at element centroids)"
c
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
