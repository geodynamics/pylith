c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
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
      subroutine get_units(ifile,ierr,nget,units_defined,units,def)
c
c...  routine to get the variable units for variable(s) 'def'.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ifile,ierr,nget
      character units(nget)*(*),def(nget)*(*)
      logical units_defined(nget)
c
c...  intrinsic functions
c
      intrinsic index
c
c...  user-defined functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer lendef,iopt,i,ii,j,jj,ib,ie
      character string*80,substr1*80,substr2*80
c
      ierr=0
      iopt=1
      do i=1,nget
        units_defined(i)=.false.
      end do
c
      do i=1,nget
        call pskip(ifile)
        read(ifile,"(a80)") string
        ii=index(string,"=")
        if(ii.eq.0) then
          ierr=2
          return
        end if
        substr1=string(1:ii-1)
        call convert_case(substr1,iopt)
        do j=1,nget
          lendef=nchar(def(j))
          jj=index(substr1,def(j)(1:lendef))
          if(jj.ne.0) then
            units_defined(j)=.true.
            substr2=string(ii+1:)
            ib=nnblnk(substr2)
            ie=nchar(substr2)
            if(ie.eq.0) then
              ierr=2
              return
            end if
            units(j)=substr2(ib:ie)
          end if
        end do
      end do
c
      do i=1,nget
        if(.not.units_defined(i)) ierr=2
      end do
c
      return
      end
c
c version
c $Id: get_units.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
