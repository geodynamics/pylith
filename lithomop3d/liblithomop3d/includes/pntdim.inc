c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2003  All Rights Reserved
c
c <LicenseText>
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
c...  pntdim.inc:  A header that describes the data type and dimensions
c     of arrays containing pointers for the program TECTON.  This
c     header should be included for the main calling routines of the
c     code.
c
c...  arrays containing pointers
c
      integer npepoint(7),npfpoint(11),npgpoint(11),npitpoint(6)
      integer npopoint(4),npppoint(5),npspoint(8),npsppoint(18)
      integer nptpoint(12),npvpoint(3),ippnt(8)
c
c...  scalar pointers associated with arrays
c
      integer ipien,iplm,iplmx,iplmf,ipmat,ipinf,ipprop
      integer ipnflt,ipflt,ipdflt,ipdftn,ipnslp,ipdiff,ipdslp,ipdx,
     & ipdelx,ipdhis,ippslp
      integer ipid,ipidx,ipbdat,ipx,ipd,ipdeld,ipbond,ipskew,ipwink,
     & ipwdat,iphist
      integer ipdcur,ipdxcur,ipdeps,ipdbeta,ipdbetb,ipscur
      integer ipkout,ipsout,ipkflt,ipkx
      integer ipelno,ipside,ipphis,ippres,ippdir
      integer ipbeta,ipbetb,ipdmat,ipstn,ipeps,ipst0,ips,ipstmp
      integer iplnz,ipja,ippcg,ipzcg,ipb,ipbtot,ipbres,ipgvc1,ipgvc2,
     & ippvc,ipdprv,ipeste,ipestf,ipdiag,ipdiagl,ipneex,iplmm,iplmmi
      integer ipmstp,ipalfa,ipdelt,ipmaxi,ipmaxc,iplgdf,ipbbar,iputol,
     & ipftol,ipetol,ipitmx,ipiprint
      integer iptflt,ipwinx,ipxdat
c
c version
c $Id: pntdim.inc,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
