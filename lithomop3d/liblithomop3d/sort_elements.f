c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine sort_elements(ien,mat,infmatmod,ivflist,ivfamily,iens,
     & ivftmp,indxiel,nen,ngauss,maxvfamilies,nvfamilies,nprestrflag,
     & numelv,numnp,nstatesz,nstatesz0,npropsz,ierr,errstrng)
c
c      This routine sorts elements into element families and determines
c      several sizes.
c
c      Infinite element information is not being dealt with in the
c      present version of the code.
c
c      Input values:
c
c        nen:              Number of element nodes.  This must be a
c                          constant for now (no mixing of elements).
c        ngauss:           Number of element gauss points.  This must
c                          also be a constant for now.
c        maxvfamilies:     Maximum possible number of element families.
c                          At present, this corresponds to the number
c                          of material types that have been defined.
c                          In the near future, this will be equal to
c                          the number of available material models.
c        nvfamilies:       Number of element families.
c        nprestrflag:      Prestress flag.  This is equal to 1 if
c                          prestresses are being used.
c        numelv:           Total number of volume elements.
c        numnp:            Total number of nodes.
c        ien(nen,numelv):  Initial (unsorted) element node array.
c                          This array may be deleted after use.
c        mat(numelv):      Material type of each element.
c                          This array may be deleted after use.
c        infmatmod(6,nmatmodmax):  Material model info array:
c                          1:  material model is defined (1) or not (0).
c                          2:  number of state variables for model.
c                          3:  number of material properties for model.
c                          4:  material model is independent (0) or
c                              dependent (1) on state variables.
c                          5:  material response is linear (0) or
c                              nonlinear (1).
c                          6:  number of initial state variables for
c                              model (prestress variables).
c        ivflist(3,maxvfamilies):  Temporary array containing
c                          preliminary element family information:
c                          1:  number of elements in family (can be 0).
c                          2:  material type for family.
c                          3:  material model for family.
c                          This array may be deleted after use.
c        ivftmp(nvfamilies):  Temporary array that is only used locally.
c                          It contains current pointers for the sorted
c                          ien array.
c                          This array may be deleted after use.
c
c      Output values:
c
c        nstatesz:         Size of state variable array.
c        nstatesz0:        Size of initial state variable array.
c        npropsz:          Size of material property array.
c        ivfamily(5,nvfamilies):  Element family info array:
c                          1:  Number of elements in family.
c                          2:  Material model for family.
c                          3:  Pointer into state variable array for
c                              the family.
c                          4:  Pointer into initial state variable array
c                              for the family.
c                          5:  Pointer into material property array for
c                              the family.
c        iens(nen,numelv): Sorted version of ien array.  The ien array
c                          may be deleted after this routine has been
c                          called.
c        indxiel(numelv):  Array containing the original element
c                          numbering.  Once all operations have been
c                          performed that require the original element
c                          numbering, this array may be deleted.
c
c      Error codes:
c          0:  No error
c        101:  Attempt to use undefined material model
c
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nen,ngauss,maxvfamilies,nvfamilies,nprestrflag,numelv
      integer numnp,nstatesz,nstatesz0,npropsz,ierr
      integer ien(nen,numelv),mat(numelv),infmatmod(6,nmatmodmax)
      integer ivflist(3,maxvfamilies),ivfamily(5,nvfamilies)
      integer iens(nen,numelv),ivftmp(nvfamilies),indxiel(numelv)
      character errstrng*(*)
c
c...  local variables
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer ifam,incstate0,nstate0max,nelfamily,j,matmod,nstate
      integer nstate0,nprop,i,iloc
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from read_connect_f!"
      call ifill(iens,izero,nen*numelv)
      call ifill(ivfamily,izero,5*nvfamilies)
      ifam=izero
      nstatesz=izero
      nstatesz0=izero
      npropsz=izero
      incstate0=izero
      nstate0max=izero
      nelfamily=izero
      if(nprestrflag.ne.izero) incstate0=ione
      ierr=izero
c
c...  set up element family information array (ivfamily)
c
      do j=1,maxvfamilies
        if(ivflist(1,j).ne.izero) then
          ifam=ifam+ione
          if(ifam.eq.ione) then
            ivftmp(ifam)=ione
          else
            ivftmp(ifam)=ivftmp(ifam-ione)+nelfamily
          end if
          nelfamily=ivflist(1,j)
          matmod=ivflist(3,j)
          ivfamily(1,ifam)=nelfamily
          ivfamily(2,ifam)=matmod
          nstate=infmatmod(2,matmod)
          nstate0=infmatmod(6,matmod)
          nstate0max=max(nstate0max,nstate0)
          nprop=infmatmod(3,matmod)
          ivfamily(3,ifam)=nstatesz+ione
          ivfamily(4,ifam)=nstatesz0+ione
          ivfamily(5,ifam)=npropsz+ione
          nstatesz=nstatesz+nstate*ngauss*nelfamily
          nstatesz0=nstatesz0+incstate0*nstate0*ngauss*nelfamily
          npropsz=npropsz+nprop
        end if
      end do
c
c...  set size of state0 array if there are no prestresses
c
      if(nprestrflag.eq.izero) nstatesz0=ngauss*nstate0max
c
c...  Loop over elements and sort them into families
c
      do i=1,numelv
        ifam=mat(i)
        matmod=ivfamily(2,ifam)
c
c...  check for illegal material model
c
        if(matmod.gt.nmatmodmax.or.infmatmod(1,matmod).eq.izero) then
          ierr=101
          errstrng="read_connect"
          return
        end if
c
c...  store element nodes in the proper family and update the pointer
c     array for that family
c
        iloc=ivftmp(ifam)
        indxiel(i)=iloc
        call iquate(iens(1,iloc),ien(1,i),nen)
        ivftmp(ifam)=ivftmp(ifam)+ione
      end do
c
      return
      end
c
c version
c $Id: sort_elements.f,v 1.2 2005/04/14 00:55:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
