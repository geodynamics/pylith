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
      subroutine scan_prop(nprop,numat,kr,ierr,density_units,
     & young_units,viscosity_coeff_units,cohesion_units,
     & ivisc,iplas,imhist,pfile)
c
c...  subroutine to perform an initial scan of the material properties
c     definition file to determine the number of material types and the
c     units being used for density, volume change, young's modulus,
c     viscosity coefficient, and cohesion.
c     this routine also sets the viscous and plastic solution flags.
c
c     Exceptions should be produced for all error codes since properties
c     must be specified.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         2:  Units not specified
c         3:  Read error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nprop,numat,kr,ierr,ivisc,iplas,imhist
      character density_units*(*),young_units*(*)
      character viscosity_coeff_units*(*),cohesion_units*(*),pfile*(*)
c
c...  defined constants
c
      character def(4)*27
      data def/"density_units","youngs_modulus_units",
     & "viscosity_coefficient_units","cohesion_units"/
c
c...  local variables
c
      integer nget,j,n
      double precision prop(20)
      character units(4)*80
      logical units_defined(4)
c
c...  open input file
c
      ierr=0
      numat=0
      nget=4
      open(kr,file=pfile,status="old",err=20)
c
c...  get units, returning error 2 if they aren't found.
c
      call get_units(kr,ierr,nget,units_defined,units,def)
      if(ierr.eq.2) return
      density_units=units(1)
      young_units=units(2)
      viscosity_coeff_units=units(3)
      cohesion_units=units(4)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      ivisc=0
      iplas=0
      imhist=0
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,(prop(j),j=1,nprop)
        call matflg(prop,nprop,ivisc,iplas,imhist)
        numat=numat+1
        go to 40
c
c...  normal return
c
 10   continue
        close(kr)
        return
c
c...  error opening file
c
 20   continue
        close(kr)
	ierr=1
        return
c
c...  read error
c
 30   continue
        ierr=3
        close(kr)
        return
c
      end
c
c version
c $Id: scan_prop.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
