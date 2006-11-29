      program faultslip
c
c...  This code computes split node boundary conditions due to
c     specified fault slip.
c     The required inputs are:
c     1. A parameter file describing fault information.
c     2. Fault files (one per fault, produced by LaGriT) that
c        describe element/node pairs lying on faults.
c     3. An AVS UCD file that contains nodal fields specifying fault
c        slip.
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
      integer maxnodes,nsd
      double precision eps
      parameter(maxfaults=10000,maxblocks=100,maxentries=10000000,
     & maxattached=2000,maxdefs=10,maxnodes=200000,nsd=3,eps=1.0d-3)
c
c...  global arrays
c
      integer nfault(2,maxentries)
      double precision split(nsd),slip(nsd,maxnodes)
c
c...  info read or deduced from parameter file
c
      integer nfaults,numdefs,faultdef(3,maxdefs)
      double precision cscale,sscale
c
c...  intrinsic functions
c
c
c...  filenames and unit numbers
c
      integer kti,kto,kp,kf,kso
      common/units/kti,kto,kp,ku,kf,kso
c
c...  local variables
c
      integer i,j,k
      integer node,elem
      integer nblocks,nentries,nattached,blk1,blk2,faultnum,ind
      integer numnp,numel,nfnodes,nfelems,nfglob,nfields
      integer elems(maxattached),colors(maxattached),isfield(nsd)
      integer itmp(100)
      double precision fsign
      double precision x(nsd),xtmp(nsd),fnorm(nsd),field(100)
      logical pairused
      character ffile*500,tmp*100
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
      read(kp,*) nblocks,nfaults,cscale,sscale
      if(cscale.eq.0.0d0) cscale=1.0d0
      if(sscale.eq.0.0d0) sscale=1.0d0
c
c...  specify fields in UCD file containing xslip, yslip, zslip, and
c     read them
c
      call pskip(kp)
      read(kp,*) (isfield(j),j=1,nsd)
      read(ku,*) numnp,numel,nfnodes,nfelems,nfglob
      do i=1,numnp
        read(ku,*) tmp
      end do
      do i=1,numel
        read(ku,*) tmp
      end do
      read(ku,*) nfields,(itmp(j),j=1,nfields)
      do i=1,nfields
        read(ku,*) tmp
      end do
      do i=1,numnp
        read(ku,*) node,(field(j),j=1,nfields)
        do j=1,nsd
          slip(j,i)=field(isfield(j))
        end do
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
                  fsign=1.0d0
                  faultnum=k
                  go to 40
                else if(colors(j).eq.faultdef(2,k)) then
                  blk1=colors(j)
                  blk2=faultdef(1,k)
                  fsign=-1.0d0
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
                do k=1,nsd
                  split(k)=0.5d0*fsign*slip(k,node)
                end do
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
      common/units/kti,kto,kp,ku,kf,kso
c
c...  local variables
c
      integer nargs,i,j
      character string*2000
      character pfile*500,ufile*500,ofile*500
      logical fflag(3)
c
      kti=5
      kto=6
      kp=10
      ku=11
      kf=12
      kso=13
c
      do i=1,3
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
        else if(index(string,'o=').ne.0) then
          j=nchar(string)
          ofile=string(3:j)
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
      open(kp,file=pfile,status="old")
      open(ku,file=ufile,status="old")
      open(kso,file=ofile,status="new")
c
 800  format("Usage:",/,
     & "    faultslip p=<parameter_file>",/,
     & "    u=<UCD_slip_specification_file>",/,
     & "    o=<split_node_output_file>")
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
