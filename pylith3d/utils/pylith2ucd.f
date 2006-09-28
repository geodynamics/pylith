c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      program pylith2ucd
c
c...  quick and dirty code to read coordinates and connectivities and
c     output them in UCD format.
c
      implicit none
c
      integer maxnodes,maxelems
      parameter(maxnodes=500000,maxelems=500000)
c
      integer kti,kto,kr,kw
      common/units/kti,kto,kr,kw
      integer numnodes,numels
      common/dimens/numnodes,numels
      integer basebeg,baseend
      common/baseinfo/basebeg,baseend
c
      character basename*500,filenm*500,dummy*80
c
      integer i,j,n,ietype,infin
      integer ien(4,maxelems),mat(maxelems)
      double precision x(3,maxnodes)
c
      call files(basename)
c
c...  read coordinates
c
      filenm=basename(basebeg:baseend)//".coord"
      open(kr,file=filenm,status="old")
      call pskip(kr)
      read(kr,"(a80)") dummy
      call pskip(kr)
      do i=1,numnodes
        read(kr,*) n,(x(j,i),j=1,3)
      end do
      close(kr)
c
c...  read connectivities
c
      filenm=basename(basebeg:baseend)//".connect"
      open(kr,file=filenm,status="old")
      call pskip(kr)
      do i=1,numels
        read(kr,*) n,ietype,mat(i),infin,(ien(j,i),j=1,4)
      end do
      close(kr)
c
c...  output mesh to UCD file
c
      filenm=basename(basebeg:baseend)//".inp"
      open(kw,file=filenm,status="replace")
      call write_ucd_mesh(x,numnodes,
     & ien,numels,
     & ietype,kw)
      close(kw)
      stop
      end
c
c
      subroutine files(basename)
c
c... subroutine to get command-line arguments
c
      implicit none
c
      character basename*(*)
c
      integer kti,kto,kr,kw
      common/units/kti,kto,kr,kw
      integer numnodes,numels
      common/dimens/numnodes,numels
      integer basebeg,baseend
      common/baseinfo/basebeg,baseend
c
      integer nnblnk,nchar
      external nnblnk,nchar
      intrinsic lnblnk,iargc
c
      integer nargs
      logical gotname,gotnodes,gotelems
c
      integer i,j
      character string*500
c
      kti=5
      kto=6
      kr=10
      kw=11
      gotname=.false.
      gotnodes=.false.
      gotelems=.false.
c
      nargs=iargc()
      do i=1,nargs
        call getarg(i,string)
        if(index(string,"b=").ne.0) then
          j=lnblnk(string)
          basename=string(3:j)
          gotname=.true.
          basebeg=nnblnk(basename)
          baseend=nchar(basename)
        else if(index(string,"n=").ne.0) then
          j=lnblnk(string)
          read(string(3:j),*) numnodes
          gotnodes=.true.
        else if(index(string,"e=").ne.0) then
          j=lnblnk(string)
          read(string(3:j),*) numels
          gotelems=.true.
        end if
      end do
c
      if(.not.gotname.or..not.gotnodes.or..not.gotelems) then
        write(kto,700)
        stop
      end if
 700  format("Usage:",/,
     & " pylith2ucd b=<base_name> n=<num_nodes> e=<num_elems>")
      return
      end
c
c
      subroutine write_ucd_mesh(
     & x,numnp,                                                         ! global
     & ien,numelv,                                                      ! elemnt
     & ietype,                                                          ! eltype
     & kucd)                                                            ! ioinfo
c
c...  Routine that outputs nodes and mesh to a UCD file.
c
      implicit none
c
      integer nenmax
      parameter(nenmax=20)
c
c...  subroutine arguments
c
      integer numnp,numelv,ietype,kucd
      integer ien(4,numelv)
      double precision x(3,numnp)
c
c...  local constants
c
      integer inducd(nenmax,10)
      data inducd/
     & 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 5, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 4, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,
     & 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18, 0, 0,
     & 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,14,13,15, 0, 0, 0, 0, 0,
     & 5, 1, 2, 3, 4,10,11,12,13, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, 0,
     & 4, 1, 2, 3, 8, 9,10, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      character eltype(10)*5
      data eltype/"hex","wrick","prism","pyr","tet",
     & "hex","wrick","prism","pyr","tet"/
c
c...  intrinsic functions
c
      intrinsic char,real
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
      integer kti,kto,kr,kw
      common/units/kti,kto,kr,kw
c
c...  local variables
c
      integer nnattr,neattr,nmattr,i,j
      integer nsd,nen,iucd
      integer iel,matmodel
      integer ibyte,intlen,floatlen
      integer itmp(20)
      character magnum*1
c
cdebug      write(6,*) "Hello from write_ucdmesh_f!"
c
      nen=4
      nsd=3
      iucd=1
      magnum=char(7)
      ibyte=1
      intlen=4
      floatlen=4
      nnattr=0
      neattr=0
      nmattr=0
      matmodel=1
c
c...  ascii ucd output
c
      if(iucd.eq.1) then
c
c...  header info
c
        write(kw,"(5i7)") numnp,numelv,nnattr,neattr,nmattr
c
c...  write nodal coordinates
c
        do i=1,numnp
          write(kucd,"(i7,3(2x,1pe15.8))") i,(x(j,i),j=1,nsd)
        end do
c
c...  write element connectivities
c
        do iel=1,numelv
          do j=1,nen
            itmp(j)=ien(j,iel)
          end do
          write(kw,"(2i7,2x,a4,20i7)") iel,matmodel,
     &       eltype(ietype),(itmp(inducd(j,ietype)),j=1,nen)
        end do
        close(kw)
      end if
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
c version
c $Id: pylith2ucd.f,v 1.4 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
