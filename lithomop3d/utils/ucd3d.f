      program ucd3d
c
c     program to produce 3d ucd files from output from the finite
c     element code tecton, and quantities along faults.  the output
c     from this code is meant to be used as input for the avs
c     data visualization system.
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter(npmax=200000,nemax=200000,nmmax=2,nodmx=95,nstmx=48,
     & nfpmx=1000,nefmx=1000,noutmx=20,nlnmx=100)
      real*8 xn(3,npmax),xd(3,npmax),dld(3*npmax),xnt(3,npmax)
      real*8 strn(nemax*9*nstmx),xdim(3),prop(10*nmmax)
      real*8 vnode(nodmx,npmax),smop(npmax),fault(6,nfpmx)
      real*8 strr(6,8,nemax),xlim(3,2)
      integer*4 nsc(nodmx),ide(nefmx),iout(noutmx),ient(8,nemax)
      integer*4 ielq(nodmx),ienp(8,nemax),inodq(nodmx)
      integer*4 isr1(2,npmax),isr2(10,npmax),isr3(8,npmax),ifs(2,nfpmx)
      integer*4 itmp(nfpmx),inq(nodmx,nlnmx)
      integer*4 nnod(nlnmx),inode(nlnmx,nlnmx),nnval(nlnmx),idftn(nfpmx)
      integer*4 ien(8,nemax),mat(nemax),infin(nemax),nfault(3,nfpmx)
      integer*4 nslip(5*nemax),ilock(nfpmx),itmp2(nfpmx),iprint(noutmx)
      character pfil*80,ucdfil*80,titl*70,string*80,unit(8)*10
      character cell(2)*4,ofil*80,pafil*80,labls(nodmx+1)*15
      logical*4 gout,nodout,cstrs,cestn,cpstn,cvstn,cdestn,cdpstn,cdvstn
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      common/dimens/npmx,nemx,nmmx,ndmx,nsmx,nfmx,nefx,nomx,nlmx
      common/scale/cscl,uscl,vscl,wscl,sscl,escl,descl,tscl
      common/units/kti,kto,kplt,kpar,kucd
      common/sflags/cstrs,cestn,cpstn,cvstn,cdestn,cdpstn,cdvstn
      common/vpoint/ivp1,ivp2,ivp3,ivp4,ivp5,ivp6,ivp7,ivp8,ivp9,ivp10,
     & ivp11,ivp12,ivpend
      common/spoint/isp(8)
      data cell/' hex','quad'/
      data nsc/nodmx*1/
      data labls/'X-Displ','Y-Displ','Z-Displ','Tot-Displ',
     & 'X-Vel','Y-Vel','Z-Vel','Tot-Vel',
     & 'Dir-Displ','Dir-Vel','Displ-Phase',
     & 'Sxx','Syy','Szz','Sxy','Syz','Sxz',
     & 'S1','S2','S3',
     & 'Pres','Max-Shear','S-Effective',
     & 'Exx','Eyy','Ezz','Exy','Eyz','Exz',
     & 'E1','E2','E3',
     & 'E-Dilatation','E-Max-Shear','E-Effective',
     & 'Vxx','Vyy','Vzz','Vxy','Vyz','Vxz',
     & 'V1','V2','V3',
     & 'V-Dilatation','V-Max-Shear','V-Effective',
     & 'Pxx','Pyy','Pzz','Pxy','Pyz','Pxz',
     & 'P1','P2','P3',
     & 'P-Dilatation','P-Max-Shear','P-Effective',
     & 'dExx','dEyy','dEzz','dExy','dEyz',
     & 'dExz','dE1','dE2','dE3',
     & 'dE-Dilatation','dE-Max-Shear','dE-Effective',
     & 'dVxx','dVyy','dVzz','dVxy','dVyz',
     & 'dVxz','dV1','dV2','dV3',
     & 'dV-Dilatation','dV-Max-Shear','dV-Effective',
     & 'dPxx','dPyy','dPzz','dPxy','dPyz',
     & 'dPxz','dP1','dP2','dP3',
     & 'dP-Dilatation','dP-Max-Shear','dP-Effective','Zero'/
      npmx=npmax
      nemx=nemax
      nmmx=nmmax
      ndmx=nodmx
      nsmx=nstmx
      nfmx=nfpmx
      nefx=nefmx
      nomx=noutmx
      nlmx=nlnmx
      gout=.false.
      nodout=.false.
      cstrs=.false.
      cestn=.false.
      cpstn=.false.
      cvstn=.false.
      cdestn=.false.
      cdpstn=.false.
      cdvstn=.false.
      ivp1=1
      ivp2=5
      ivp3=9
      ivp4=10
      ivp5=11
      ivp6=12
      ivp7=24
      ivp8=36
      ivp9=48
      ivp10=60
      ivp11=72
      ivp12=84
      ivpend=95
c
c...  set up i/o files.
c
      call files(pfil,pafil,ucdfil)
c
c...  read control parameters from parameters file
c
      call skipa(kpar)
      read(kpar,*) ndval,neval,nline,nout
      if(ndval.ne.0) gout=.true.
      if(neval.ne.0) gout=.true.
      if(nline.ne.0) nodout=.true.
      ntval=0
      call skipa(kpar)
      if(ndval.ne.0) read(kpar,*) (itmp(i),i=1,ndval)
      call iclear(inodq,nodmx)
      do i=1,ndval
        inodq(itmp(i))=1
	if(itmp(i).ge.ivp6.and.itmp(i).le.ivp7-1) cstrs=.true.
	if(itmp(i).ge.ivp7.and.itmp(i).le.ivp8-1) cestn=.true.
	if(itmp(i).ge.ivp8.and.itmp(i).le.ivp9-1) cvstn=.true.
	if(itmp(i).ge.ivp9.and.itmp(i).le.ivp10-1) cpstn=.true.
	if(itmp(i).ge.ivp10.and.itmp(i).le.ivp11-1) cdestn=.true.
	if(itmp(i).ge.ivp11.and.itmp(i).le.ivp12-1) cdvstn=.true.
	if(itmp(i).ge.ivp12.and.itmp(i).le.ivpend) cdpstn=.true.
      end do
      call skipa(kpar)
      if(neval.ne.0) read(kpar,*) (itmp(i),i=1,neval)
      call iclear(ielq,nodmx)
      do i=1,neval
        ielq(itmp(i))=1
      end do
      call skipa(kpar)
      call iclear(nnod,nlnmx)
      call iclear(inode,nlnmx*nlnmx)
      call iclear(nnval,nlnmx)
      call iclear(inq,nodmx*nlnmx)
      if(nline.ne.0) read(kpar,*) (nnod(i),i=1,nline)
      call skipa(kpar)
      do i=1,nline
        read(kpar,*) (inode(j,i),j=1,nnod(i))
      end do
      call skipa(kpar)
      if(nline.ne.0) read(kpar,*) (nnval(i),i=1,nline)
      call skipa(kpar)
      do i=1,nline
        read(kpar,*) (itmp(j),j=1,nnval(i))
        do j=1,nnval(i)
	  if(itmp(j).ge.ivp6.and.itmp(j).le.ivp7-1) cstrs=.true.
	  if(itmp(j).ge.ivp7.and.itmp(j).le.ivp8-1) cestn=.true.
	  if(itmp(j).ge.ivp8.and.itmp(j).le.ivp9-1) cvstn=.true.
	  if(itmp(j).ge.ivp9.and.itmp(j).le.ivp10-1) cpstn=.true.
	  if(itmp(j).ge.ivp10.and.itmp(j).le.ivp11-1) cdestn=.true.
	  if(itmp(j).ge.ivp11.and.itmp(j).le.ivp12-1) cdvstn=.true.
	  if(itmp(j).ge.ivp12.and.itmp(j).le.ivpend) cdpstn=.true.
          inq(itmp(j),i)=1
        end do
      end do
c
c...  set stress pointers
c
      do i=1,7
	isp(i)=1
      end do
      isp(8)=0
      if(cstrs) then
	do i=2,8
	  isp(i)=isp(i)+12
        end do
      end if
      if(cestn) then
	do i=3,8
	  isp(i)=isp(i)+12
        end do
      end if
      if(cvstn) then
	do i=4,8
	  isp(i)=isp(i)+12
        end do
      end if
      if(cpstn) then
	do i=5,8
	  isp(i)=isp(i)+12
        end do
      end if
      if(cdestn) then
	do i=6,8
	  isp(i)=isp(i)+12
        end do
      end if
      if(cdvstn) then
	do i=7,8
	  isp(i)=isp(i)+12
        end do
      end if
      if(cdpstn) then
        isp(8)=isp(8)+12
      end if
      if(isp(8).gt.nstmx) then
	write(kto,*) 'Insufficient storage allocated for stress/strain!'
	write(kto,*) 'Increase nstmx!'
	stop
      end if
      call skipa(kpar)
      if(nout.ne.0) read(kpar,*) (iout(i),i=1,nout)
      call skipa(kpar)
      read(kpar,*) vx,vy,vz,rlam
      call skipa(kpar)
      read(kpar,*) cscl,uscl,vscl,wscl,sscl,escl,descl,tscl,amag,rmax
      if(cscl.eq.0.0) cscl=1.0
      if(uscl.eq.0.0) uscl=1.0
      if(vscl.eq.0.0) vscl=1.0
      if(wscl.eq.0.0) wscl=1.0
      if(sscl.eq.0.0) sscl=1.0
      if(escl.eq.0.0) escl=1.0
      if(descl.eq.0.0) descl=1.0
      if(tscl.eq.0.0) tscl=1.0
      do i=1,3
        xlim(i,1)=1.e20
        xlim(i,2)=-1.e20
      end do
      do i=1,8
	call skipa(kpar)
	read(kpar,810) unit(i)
	i1=nnblnk(unit(i))
	string=unit(i)(i1:10)
	i2=index(' ',string)
	if(i2.eq.0) i2=nchar(string)
	unit(i)=string(1:i2)
      end do
c
c...  read general model information
c
      read(kplt) titl
      write(kto,'(//a70)') titl
      read(kplt) ngem,numnp,nsd,ndof,nstr,nen
      call rread(xn,nsd,numnp,1)
      do i=1,numnp
        do j=1,nsd
          xnt(j,i)=xn(j,i)
          xlim(j,1)=min(xlim(j,1),xn(j,i))
          xlim(j,2)=max(xlim(j,2),xn(j,i))
        end do
      end do
      do i=1,3
        xdim(i)=abs(xlim(i,2)-xlim(i,1))
      end do
      read(kplt) numout
      if(numout.gt.0) read(kplt) (itmp2(i),i=1,numout)
      read(kplt) icontr,ncycle,lastep
      if(icontr.gt.0) read(kplt) (iprint(i),i=1,icontr)
      read(kplt) nsout
      if(nsout.gt.0) read(kplt) (itmp2(i),i=1,nsout)
      read(kplt) ivisc,iplas
      read(kplt) nmat,nprop
      call rread(prop,nprop,nmat,1)
      read(kplt) numel,ngstn
      call iread(mat,numel,1,1)
      call iread(infin,numel,1,1)
      call iread(ien,nen,numel,1)
      do i=1,numel
        do j=1,nen
          ienp(j,i)=ien(j,i)
        end do
      end do
      read(kplt) numfn
      if(numfn.gt.0) call iread(nfault,3,numfn,1)
      if(numfn.gt.0) call rread(fault,ndof,numfn,1)
      read(kplt) numflt
      if(numflt.ne.0) call iread(idftn,numflt,1,1)
      read(kplt) numslp
      if(numslp.gt.0) call iread(nslip,ndof+2,numslp,1)
      read(kplt) numsn
c
c...  adjust arrays for the presence of split and slippery nodes.
c
      call fltadj(xn,fault,ien,nfault,nslip,ifs,ide,numfn,numslp,nf,ns,
     & nfp,nfe)
c
c...  get information on nodes lying on a surface.  this info is used
c     for stress smoothing.
c
      call nsrcnt(xn,ifs,nf,ns,isr1,isr2,isr3,xlim)
c
c...  loop over time steps and create ucd file for each requested step
c
      istep=0
      do i=1,nout
        if(iout(i).gt.icontr) then
          write(kto,*) 'Time step number too large! '
          stop
        end if
        itdif=iout(i)-istep
        if(itdif.gt.0) then
          do j=1,itdif
            read(kplt) itmp2(1)
            nstep=itmp2(1)
            read(kplt) tdum
            call rread(dld,nsd,numnp,1)
            if(nstep.gt.0) then
              read(kplt) deltp
              call rread(dld,nsd,numnp,1)
            end if
            do k=1,numfn
              if(nstep.eq.0) read(kplt) (fault(l,k),l=1,3)
              if(nstep.gt.0) read(kplt) (fault(l,k),l=1,6)
            end do
            do k=1,numsn
              if(nstep.eq.0) read(kplt) m,(fault(l,k),l=1,3)
              if(nstep.gt.0) read(kplt) m,(fault(l,k),l=1,6)
            end do
            do k=1,numsn
              read(kplt) ilock(k+numfn)
            end do
            call rread(strr,nstr,ngstn,numel)
            call rread(strr,nstr,ngstn,numel)
	    if(nstep.gt.0) then
	      if(ivisc.eq.1) call rread(strr,nstr,ngstn,numel)
	      if(iplas.eq.1) call rread(strr,nstr,ngstn,numel)
	      call rread(strr,nstr,ngstn,numel)
	      if(ivisc.eq.1) call rread(strr,nstr,ngstn,numel)
	      if(iplas.eq.1) call rread(strr,nstr,ngstn,numel)
	    end if
            istep=istep+1
          end do
        else if(itdif.lt.0) then
          write(kto,*) 'Invalid time step number requested!'
          stop
        end if
        nstep=0
        jj=iout(i)
        if(jj.gt.0) nstep=iprint(jj)
        write(string,'(i5)') nstep
        i1=nchar(ucdfil)
        is1=nnblnk(string)
        is2=nchar(string)
        ofil=ucdfil(1:i1)//string(is1:is2)//'.inp'
        if(gout) open(unit=kucd,file=ofil,status='new')
        ofil=ucdfil(1:i1)//string(is1:is2)
        if(gout) write(kucd,700) titl,nstep
	nw1=ndval
	if(nw1.ne.0) nw1=nw1+1
	nw2=neval
	if(nw2.ne.0) nw2=nw2+1
        if(gout) write(kucd,710) nnp,numel,nw1,nw2,ntval
	call flush(kucd)
        call stread(xn,xnt,xd,strn,smop,dld,strr,ien,infin,ienp,ient,
     &   isr1,isr2,isr3,vnode,numsn,numfn,numflt,ifs,nfault,fault,ilock,
     &   prop,vx,vy,vz,rlam,itmp,itmp2,idftn,nf,ns,nfe,nprop,time,amag,
     &   xlim,ivisc,iplas)
c
c...  output nodal coordinates
c
        write(kto,*) 'Writing nodal coordinates:'
        if(gout) then
          do j=1,nnp
            write(kucd,720) j,(xd(k,j),k=1,ndof)
          end do
        end if
c
c...  output connectivity and material type info
c
        write(kto,*) 'Writing element connectivities:'
        if(gout) then
          do j=1,numel
            write(kucd,730) j,mat(j),cell(1),(ien(k,j),k=1,8)
          end do
        end if
c
c...  output data values at nodal points
c
        write(kto,*) 'Writing nodal data values:'
        if(ndval.ne.0) call ndout(vnode,inodq,labls,unit,nsc,ndval,kucd)
c
c...  output data values at element centroids
c
        write(kto,*) 'Writing element data values:'
        if(neval.ne.0) call neout(strn,ielq,labls,unit,nsc,neval,kucd)
c
c...  output tabulated data along specified lines of nodes
c
        write(kto,*) 'Writing line data for plotting:'
        if(nodout) call ndlout(xn,vnode,nnod,inode,nnval,inq,labls,
     &   unit,nline,cscl,ofil)
        istep=istep+1
        close(kucd)
      end do
700   format('# ',a60,' step ',i5)
710   format(16i7)
720   format(i7,3(2x,1pe15.8))
730   format(2i7,2x,a4,8i7)
800   format(a60)
810   format(a10)
      stop
      end
c
c
      subroutine centint(xn,ien,str,infin)
c
c...  subroutine to recompute stresses/strains at element centroid
c     rather than 2x2x2 gauss points.
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension xn(3,*),str(numel,ngstn+1,nstmx),ien(8,*),infin(*)
      dimension sh(4,8),rg(8),sg(8),tg(8),det(8),xl(3,8),ia(8)
      data rg/-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0/
      data sg/-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0/
      data tg/ 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0/
      g=1.0/sqrt(3.0)
      do n=1,numel
        do j=1,8
          ia(j)=ien(j,n)
          do k=1,ndof
            xl(k,j)=xn(k,ia(j))
          end do
        end do
        vol=0.0
        do l=1,ngstn
          call shap30(rg(l)*g,sg(l)*g,tg(l)*g,xl,det(l),sh,nen,ndof,ia,
     &     infin(n),n)
          vol=vol+det(l)
        end do
        do j=1,nstr
          str(n,ngstn+1,j)=0.0
          do l=1,ngstn
            rm=det(l)/vol
            str(n,ngstn+1,j)=str(n,ngstn+1,j)+rm*str(n,l,j)
          end do
        end do
      end do
      return
      end
c
c
      subroutine clear(a,na)
c
c...  subroutine to clear a floating-point array
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension a(*)
      do i=1,na
        a(i)=0.0
      end do
      return
      end
c
c
      subroutine covsrt(covar,ncvm,ma,lista,mfit)
c
c...  subroutine to sort covariance matrix
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension covar(ncvm,ncvm),lista(mfit)
      do 12 j=1,ma-1
        do 11 i=j+1,ma
          covar(i,j)=0.0
11      continue
12    continue
      do 14 i=1,mfit-1
        do 13 j=i+1,mfit
          if(lista(j).gt.lista(i)) then
            covar(lista(j),lista(i))=covar(i,j)
          else
            covar(lista(i),lista(j))=covar(i,j)
          endif
13      continue
14    continue
      swap=covar(1,1)
      do 15 j=1,ma
        covar(1,j)=covar(j,j)
        covar(j,j)=0.0
15    continue
      covar(lista(1),lista(1))=swap
      do 16 j=2,mfit
        covar(lista(j),lista(j))=covar(1,j)
16    continue
      do 18 j=2,ma
        do 17 i=1,j-1
          covar(i,j)=covar(j,i)
17      continue
18    continue
      return
      end
c
c
      subroutine cross(vc,v1,v2)
c
c...  subroutine to compute the cross product of two 3-vectors
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension vc(*),v1(*),v2(*)
      vc(1)= v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)= v1(1)*v2(2)-v1(2)*v2(1)
      return
      end
c
c
      subroutine eigsrt(d,v,n,np)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension d(np),v(np,np)
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      end
c
c
      subroutine equate(temp,amat,idim)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
c
c        copies the vector (matrix) amat into temp
c
      dimension temp(*),amat(*)
c
      do i=1,idim
        temp(i)=amat(i)
      end do
      return
      end
c
c
      subroutine extrap(xn,isr,vnode,nintmx,iadd)
c
c...  subroutine to perform linear extrapolation through given nodes
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension xn(3,*),isr(*),vnode(nodmx,*),a(4),cov(4,4),la(4),x(9)
      dimension y(9),sd(9)
      nd=0
      ma=4
      mf=4
      ncvm=4
      i1=isr(1)
      do i=1,ma
        la(i)=i
        a(i)=1.0
      end do
c
c...  assign 'x' values and standard deviations for inversion
c
      do i=2,nintmx+1
        if(isr(i).ne.0) then
          nd=nd+1
          xx=xn(1,isr(i))-xn(1,i1)
          yy=xn(2,isr(i))-xn(2,i1)
          zz=xn(3,isr(i))-xn(3,i1)
          x(nd)=float(isr(i))
          sd(nd)=sqrt(xx*xx+yy*yy+zz*zz)
        end if
      end do
c
c...  loop over stresses/strains and obtain extrapolation coefficients
c
      do i=1,nstr
        do j=1,nd
          y(j)=vnode(i+iadd,isr(1+j))
        end do
        call lfit(x,y,sd,nd,a,ma,la,mf,cov,ncvm,chisq,xn)
        vnode(i+iadd,i1)=a(1)*xn(1,i1)+a(2)*xn(2,i1)+a(3)*xn(3,i1)+a(4)
      end do
      return
      end
c
c
      subroutine files(pfil,pafil,ucdfil)
c
c...  subroutine to define default i/o file names and unit numbers
c     and open necessary files.
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/units/kti,kto,kplt,kpar,kucd
      character pfil*(*),ucdfil*(*),pafil*(*),string*80
      pfil='fedsk.dat'
      pafil='fedsk.par'
      ucdfil='step'
      nargs=iargc()
      do i=1,nargs
        call getarg(i,string)
        if(index(string,'u=').ne.0) then
          j=lnblnk(string)
          ucdfil=string(3:j)
        else if(index(string,'p=').ne.0) then
          j=lnblnk(string)
          pfil=string(3:j)
        else if(index(string,'c=').ne.0) then
          j=lnblnk(string)
          pafil=string(3:j)
        end if
      end do
c
c     define and open i/o units
c
      kti=5
      kto=6
      kplt=12
      kpar=14
      kucd=15
      open(unit=kplt,file=pfil,status='old',form='unformatted')
      open(unit=kpar,file=pafil,status='old')
      return
      end
c
c
      subroutine fltadj(xn,fault,ien,nfault,nslip,ifs,ide,numfn,numslp,
     & nf,ns,nfp,nfe)
c
c...  subroutine to add extra nodes to represent both sides of a
c     discontinuity (slippery nodes or split nodes) and to group
c     discontinuity surfaces together
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      logical*4 nused,eused
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension xn(3,*),fault(ndof,numfn),ifs(2,*),ide(nefmx)
      dimension ien(8,*),nfault(3,numfn),nslip(ndof+2,*)
      nnp=numnp
      nf=0
      ns=0
      nes=0
      nfp=0
      nfe=0
      if(numfn.eq.0.and.numslp.eq.0) return
      call iclear(ifs,2*nfpmx)
c
c...  find split nodes
c
      do i=1,numfn
	eused=.false.
	nused=.false.
	do j=1,nes
	  if(nfault(1,i).eq.ide(j)) eused=.true.
        end do
	do j=1,nf
	  if(nfault(2,i).eq.ifs(1,j)) nused=.true.
        end do
	if(.not.eused.and.(fault(1,i).lt.0.0.or.fault(2,i).lt.0.0.or.
     &   fault(3,i).lt.0.0)) then
	  nes=nes+1
	  ide(nes)=nfault(1,i)
        end if
	if(.not.nused) then
	  nf=nf+1
	  nnp=nnp+1
	  ifs(1,nf)=nfault(2,i)
	  ifs(2,nf)=nnp
	  xn(1,nnp)=xn(1,nfault(2,i))
	  xn(2,nnp)=xn(2,nfault(2,i))
	  xn(3,nnp)=xn(3,nfault(2,i))
        end if
      end do
      do n=1,nes
	nn=ide(n)
	do i=1,8
	  do j=1,nf
	    if(ien(i,nn).eq.ifs(1,j)) ien(i,nn)=ifs(2,j)
	  end do
	end do
	do j=1,numfn
	  if(nfault(1,j).eq.nn) then
	    do k=1,nf
	      if(nfault(2,j).eq.ifs(1,k)) then
		nfault(2,j)=ifs(2,k)
		go to 5
              end if
            end do
5           continue
	  end if
        end do
      end do
      call iclear(ide,nefmx)
      nes=0
c
c...  find slippery nodes
c
      do i=1,numslp
	eused=.false.
	nused=.false.
	do j=1,nes
	  if(nslip(1,i).eq.ide(j)) eused=.true.
        end do
	do j=1,ns
	  if(nslip(2,i).eq.ifs(1,j+nf)) nused=.true.
        end do
	if(.not.eused.and.(nslip(3,i).lt.0.or.nslip(4,i).lt.0.or.
     &   nslip(5,i).lt.0)) then
	  nes=nes+1
	  ide(nes)=nslip(1,i)
        end if
	if(.not.nused) then
	  ns=ns+1
	  nnp=nnp+1
	  ifs(1,ns+nf)=nslip(2,i)
	  ifs(2,ns+nf)=nnp
	  xn(1,nnp)=xn(1,nslip(2,i))
	  xn(2,nnp)=xn(2,nslip(2,i))
	  xn(3,nnp)=xn(3,nslip(2,i))
        end if
      end do
      do n=1,nes
	nn=ide(n)
	do i=1,8
	  do j=1,ns
	    if(ien(i,nn).eq.ifs(1,j+nf)) ien(i,nn)=ifs(2,j+nf)
	  end do
	end do
      end do
      return
      end
c
c
      subroutine funcs(x,xn,afunc,ma)
c
c... subroutine to return the coordinates of point x
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension xn(3,*),afunc(*)
      ii=nint(x)
      afunc(1)=xn(1,ii)
      afunc(2)=xn(2,ii)
      afunc(3)=xn(3,ii)
      afunc(4)=1.0
      return
      end
c
c
      subroutine gaussj(a,n,np,b,m,mp)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (nmax=50)
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix.'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      end
c
c
      subroutine iclear(ia,na)
c
c...  subroutine to clear an integer*4 array
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension ia(*)
      do i=1,na
        ia(i)=0
      end do
      return
      end
c
c
      subroutine infel(v,va,vv,dvv,io)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
c
c...  subroutine to compute local shape functions and derivatives for
c     infinite elements
c
      vv=0.5*(1.0+va*v)
      dvv=0.5*va
      if(io.eq.1) then
        if(va.eq.-1.0) then
          vv=1.0-2.0*v/(1.0+v)
          dvv=-2.0/((1.0+v)*(1.0+v))
        else
          vv=2.0*v/(1.0+v)
          dvv=2.0/((1.0+v)*(1.0+v))
        end if
      else if(io.eq.2) then
        if(va.eq.-1.0) then
          vv=-2.0*v/(1.0-v)
          dvv=-2.0/((1.0-v)*(1.0-v))
        else
          vv=1.0+2.0*v/(1.0-v)
          dvv=2.0/((1.0-v)*(1.0-v))
        end if
      end if
      return
      end
c
c
      subroutine iread(iar,idim1,idim2,idim3)
c
c...  subroutine to read an integer*4 array written in binary format
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/units/kti,kto,kplt,kpar,kucd
      dimension iar(idim1,idim2,idim3)
      read(kplt) iar
      return
      end
c
c
      subroutine jacobi(a,n,np,d,v,nrot)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (nmax=100)
      dimension a(np,np),d(np),v(np,np),b(nmax),z(nmax)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     *         .and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause '50 iterations should never happen'
      return
      end
c
c
      subroutine lfit(x,y,sig,ndata,a,ma,lista,mfit,covar,ncvm,chisq,xn)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (mmax=10)
      dimension x(ndata),y(ndata),sig(ndata),a(ma),lista(ma),
     *    covar(ncvm,ncvm),beta(mmax),afunc(mmax),xn(3,*)
      kk=mfit+1
      do 12 j=1,ma
        ihit=0
        do 11 k=1,mfit
          if (lista(k).eq.j) ihit=ihit+1
11      continue
        if (ihit.eq.0) then
          lista(kk)=j
          kk=kk+1
        else if (ihit.gt.1) then
          pause 'improper set in lista'
        endif
12    continue
      if (kk.ne.(ma+1)) pause 'improper set in lista'
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=0.
13      continue
        beta(j)=0.
14    continue
      do 18 i=1,ndata
        call funcs(x(i),xn,afunc,ma)
        ym=y(i)
        if(mfit.lt.ma) then
          do 15 j=mfit+1,ma
            ym=ym-a(lista(j))*afunc(lista(j))
15        continue
        endif
        sig2i=1./sig(i)**2
        do 17 j=1,mfit
          wt=afunc(lista(j))*sig2i
          do 16 k=1,j
            covar(j,k)=covar(j,k)+wt*afunc(lista(k))
16        continue
          beta(j)=beta(j)+ym*wt
17      continue
18    continue
      if (mfit.gt.1) then
        do 21 j=2,mfit
          do 19 k=1,j-1
            covar(k,j)=covar(j,k)
19        continue
21      continue
      endif
      call gaussj(covar,mfit,ncvm,beta,1,1)
      do 22 j=1,mfit
        a(lista(j))=beta(j)
22    continue
      chisq=0.
      do 24 i=1,ndata
        call funcs(x(i),xn,afunc,ma)
        sum=0.
        do 23 j=1,ma
          sum=sum+a(j)*afunc(j)
23      continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
24    continue
      call covsrt(covar,ncvm,ma,lista,mfit)
      return
      end
c
c
      subroutine lubksb(a,n,np,indx,b)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension a(np,np),indx(n),b(n)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end
c
c
      subroutine ludcmp(a,n,np,indx,d,iel,icase)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (nmax=100,tiny=1.0e-20)
      common/units/kti,kto,kplt,kpar,kucd
      dimension a(np,np),indx(n),vv(nmax)
      d=1.0
      do 12 i=1,n
        aamax=0.0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.0) write(kto,800) iel,icase
        vv(i)=1.0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.0)a(j,j)=tiny
          dum=1.0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.0)a(n,n)=tiny
800   format('singular matrix for element #',i4,'.  case = ',i2)
      return
      end
c
c
      subroutine ludinv(a,ainv,n,np,indx,iel,icase)
c
c...  subroutine to invert a matrix using lu decomposition
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension a(np,*),ainv(np,*),indx(*)
c
c...  create identity matrix
c
      do i=1,n
        do j=1,n
          ainv(i,j)=0.0
        end do
        ainv(i,i)=1.0
      end do
c
c...  decompose matrix
c
      call ludcmp(a,n,np,indx,d,iel,icase)
c
c...  backsubstitute to obtain inverse matrix
c
      do i=1,n
        call lubksb(a,n,np,indx,ainv(1,i))
      end do
      return
      end
c
c
      integer*4 function nchar(string)
c
c       determines the minimum nonblank length of a string
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      character*(*) string
      character blank
      data blank/' '/
      nmax=len(string)
      nchar=0
      do i=1,nmax
        itest=nmax-i+1
        if(string(itest:itest).ne.blank) then
          nchar=itest
          return
        endif
      end do
      return
      end
c
c
      subroutine ndlout(xn,vnode,nnod,inode,nnval,inq,labls,unit,nline,
     & cscl,ofil)
c
c...  subroutine to output data values along specified lines of nodes
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/vpoint/ivp1,ivp2,ivp3,ivp4,ivp5,ivp6,ivp7,ivp8,ivp9,ivp10,
     & ivp11,ivp12,ivpend
      dimension xn(3,*),vnode(nodmx,*),nnod(*),inode(nlnmx,*),nnval(*)
      dimension inq(nodmx,*),vtmp(95),xr(3),xp(3)
      character labls(*)*(*),labll*600,tab*1,ofil*(*),filnm*80,string*5
      character unit(*)*(*),str2*30
      tab='\t'
      knod=20
      do n=1,nline
        write(string,'(i5)') n
        is1=nnblnk(string)
        is2=nchar(string)
        i1=nnblnk(ofil)
        i2=nchar(ofil)
        filnm=ofil(i1:i2)//'.line'//string(is1:is2)//'.tab'
        open(unit=knod,file=filnm,status='new')
        do i=1,600
          labll(i:i)=' '
        end do
	i1=1
	i2=nchar(unit(7))
	str2='('//unit(7)(i1:i2)//')'
	i2=i2+2
	ib=1
	il=nchar(str2)+12
        labll(ib:il)='X-Distance '//str2(i1:i2)//tab
        ib=ib+il
	ie=ib+il-1
        labll(ib:ie)='Y-Distance '//str2(i1:i2)//tab
        ib=ib+il
	ie=ib+il-1
        labll(ib:ie)='Z-Distance '//str2(i1:i2)//tab
        ib=ib+il
	il=nchar(str2)+10
	ie=ib+il-1
        labll(ib:ie)=  'Distance '//str2(i1:i2)
        ib=nchar(labll)+1
        do i=1,nodmx
          ii=inq(i,n)
          if(ii.ne.0) then
	    ind=1
	    if(i.ge.ivp2) ind=2
	    if(i.eq.ivp3) ind=1
	    if(i.eq.ivp4) ind=2
	    if(i.eq.ivp5) ind=3
	    if(i.ge.ivp6) ind=4
	    if(i.ge.ivp7) ind=5
	    if(i.ge.ivp10) ind=6
	    nu=nchar(unit(ind))
            nls=nchar(labls(i))
	    str2=tab//labls(i)(1:nls)//' ('//unit(ind)(1:nu)//')'
            lens=nls+nu+4
            labll(ib:ib+lens)=str2(1:lens)
            ib=ib+lens
          end if
        end do
        write(knod,'(600(:a))') labll(1:ib-1)
        xr(1)=cscl*xn(1,inode(1,n))
        xr(2)=cscl*xn(2,inode(1,n))
        xr(3)=cscl*xn(3,inode(1,n))
        vtmp(4)=0.0
        do i=1,nnod(n)
          xp(1)=cscl*xn(1,inode(i,n))
          xp(2)=cscl*xn(2,inode(i,n))
          xp(3)=cscl*xn(3,inode(i,n))
          d1=xp(1)-xr(1)
          d2=xp(2)-xr(2)
          d3=xp(3)-xr(3)
          vtmp(1)=xp(1)
          vtmp(2)=xp(2)
          vtmp(3)=xp(3)
          vtmp(4)=vtmp(4)+sqrt(d1*d1+d2*d2+d3*d3)
          xr(1)=xp(1)
          xr(2)=xp(2)
          xr(3)=xp(3)
          jj=4
          do j=1,nodmx
            if(inq(j,n).ne.0) then
              jj=jj+1
              vtmp(jj)=vnode(j,inode(i,n))
            end if
          end do
          write(knod,700) vtmp(1),(tab,vtmp(j),j=2,nnval(n)+4)
        end do
        close(knod)
      end do
700   format(1pe15.8,29(a1,1pe15.8))
      return
      end
c
c
      subroutine ndout(vnode,inodq,labls,unit,nsc,ndval,kucd)
c
c...  subroutine to output nodal values to AVS file
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      common/vpoint/ivp1,ivp2,ivp3,ivp4,ivp5,ivp6,ivp7,ivp8,ivp9,ivp10,
     & ivp11,ivp12,ivpend
      dimension vnode(nodmx,*),inodq(*),nsc(*),vtmp(95)
      character labls(*)*(*),unit(*)*(*),string*30
      iz=1
      zero=0.0
      write(kucd,700) ndval+1,(nsc(i),i=1,ndval),iz
      do i=1,nodmx
        if(inodq(i).ne.0) then
	  ind=1
	  if(i.ge.ivp2) ind=2
	  if(i.eq.ivp3) ind=1
	  if(i.eq.ivp4) ind=2
	  if(i.eq.ivp5) ind=3
	  if(i.ge.ivp6) ind=4
	  if(i.ge.ivp7) ind=5
	  if(i.ge.ivp10) ind=6
	  nu=nchar(unit(ind))
	  nls=nchar(labls(i))
	  string=labls(i)(1:nls)//','//unit(ind)(1:nu)
	  lens=nu+nls+1
	  write(kucd,'(30(:a))') string(1:lens)
	end if
      end do
      write(kucd,'(80(:a))') labls(nodmx+1)
      do i=1,nnp
        ii=0
        do j=1,nodmx
          if(inodq(j).ne.0) then
            ii=ii+1
            vtmp(ii)=vnode(j,i)
          end if
        end do
	vtmp(ndval+1)=zero
        write(kucd,710) i,(vtmp(j),j=1,ndval+1)
      end do
700   format(26i7)
710   format(i7,25(2x,1pe15.8))
      return
      end
c
c
      subroutine neout(strn,ielq,labls,unit,nsc,neval,kucd)
c
c...  subroutine to output element values to AVS file
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      common/vpoint/ivp1,ivp2,ivp3,ivp4,ivp5,ivp6,ivp7,ivp8,ivp9,ivp10,
     & ivp11,ivp12,ivpend
      common/spoint/isp(8)
      dimension strn(numel,ngstn+1,nstmx),ielq(*),nsc(*),etmp(95)
      character labls(*)*(*),unit(*)*(*),string*30
      iz=1
      zero=0.0
      write(kucd,700) neval+1,(nsc(i),i=1,neval),iz
      do i=1,nodmx
        if(ielq(i).ne.0) then
	  ind=1
	  if(i.ge.ivp2) ind=2
	  if(i.eq.ivp3) ind=1
	  if(i.eq.ivp4) ind=2
	  if(i.eq.ivp5) ind=3
	  if(i.ge.ivp6) ind=4
	  if(i.ge.ivp7) ind=5
	  if(i.ge.ivp10) ind=6
	  nu=nchar(unit(ind))
	  nls=nchar(labls(i))
	  string=labls(i)(1:nls)//','//unit(ind)(1:nu)
	  lens=nu+nls+1
	  write(kucd,'(30(:a))') string(1:lens)
	end if
      end do
      write(kucd,'(80(:a))') labls(nodmx+1)
      jj=1
      if(ngstn.gt.1) jj=9
      do i=1,numel
        ii=0
        do j=ivp6,nodmx
          if(ielq(j).ne.0) then
            ii=ii+1
	    if(j.ge.ivp6.and.j.lt.ivp7)   ll=j-ivp6 +isp(1)
	    if(j.ge.ivp7.and.j.lt.ivp8)   ll=j-ivp7 +isp(2)
	    if(j.ge.ivp8.and.j.lt.ivp9)   ll=j-ivp8 +isp(3)
	    if(j.ge.ivp9.and.j.lt.ivp10)  ll=j-ivp9 +isp(4)
	    if(j.ge.ivp10.and.j.lt.ivp11) ll=j-ivp10+isp(5)
	    if(j.ge.ivp11.and.j.lt.ivp12) ll=j-ivp11+isp(6)
	    if(j.ge.ivp12)                ll=j-ivp12+isp(7)
	    etmp(ii)=strn(i,jj,ll)
          end if
        end do
	etmp(neval+1)=zero
        write(kucd,710) i,(etmp(j),j=1,neval+1)
      end do
700   format(26i7)
710   format(i7,25(2x,1pe15.8))
      return
      end
c
c
      integer*4 function nnblnk(string)
c
c       determines the position of the first nonblank entry
c       of a string (returns 1 if the first character is
c       not blank)
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      character*(*) string
      character blank
      data blank/' '/
      nmax=len(string)
      nnblnk=nmax
      do i=1,nmax
        if(string(i:i).ne.blank) then
          nnblnk=i
          return
        endif
      end do
      return
      end
c
c
      subroutine nsrcnt(xn,ifs,nf,ns,isr1,isr2,isr3,xlim)
c
c...  subroutine to count the number of surfaces on which a node lies
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/surf/nsr1,nsr2,nsr3,nsrf1,nsrf2,nsrf3
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension xn(3,*),ifs(2,*),isr1(2,*)
      dimension isr2(10,*),isr3(8,*),xlim(3,*)
      logical*4 fault
      nsr1=0
      nsr2=0
      nsr3=0
      eps=1.0
      call iclear(isr1,2*nnp)
      call iclear(isr2,10*nnp)
      call iclear(isr3,8*nnp)
c
c...  loop over nodes and determine how many surfaces each one lies on
c
      do i=1,nnp
        nsurf=0
        fault=.false.
        do j=1,ndof
          if(abs(xn(j,i)-xlim(j,1)).lt.eps) nsurf=nsurf+1
          if(abs(xn(j,i)-xlim(j,2)).lt.eps) nsurf=nsurf+1
        end do
c
c...  split or slippery nodes are re-evaluated at each time step
c
        do j=1,nf+ns
          if(ifs(1,j).eq.i.or.ifs(2,j).eq.i) then
            fault=.true.
            go to 10
          end if
        end do
10      continue
        if(nsurf.eq.1.and.(.not.fault)) then
          nsr1=nsr1+1
          isr1(1,nsr1)=i
        else if(nsurf.eq.2.and.(.not.fault)) then
          nsr2=nsr2+1
          isr2(1,nsr2)=i
        else if(nsurf.eq.3.and.(.not.fault)) then
          nsr3=nsr3+1
          isr3(1,nsr3)=i
        end if
      end do
      nsrf1=nsr1
      nsrf2=nsr2
      nsrf3=nsr3
      return
      end
c
c
      subroutine nsrfnd(xn,ien,isr1,isr2,isr3)
c
c...  subroutine to find nodes for extrapolation to surface nodes
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/surf/nsr1,nsr2,nsr3,nsrf1,nsrf2,nsrf3
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension xn(3,*),ien(8,*),isr1(2,*),isr2(10,*),isr3(8,*),ie(4)
c
c...  if node lies on a single surface, find the node directly beneath
c
      do i=1,nsrf1
        isr1(2,i)=0
        nn=isr1(1,i)
        do n=1,numel
          do j=1,8
            if(ien(j,n).eq.nn) then
              nel=n
              go to 10
            end if
          end do
        end do
10      continue
        dmin=1.e20
        do j=1,8
          ii=ien(j,nel)
          do k=1,nsrf1
            if(isr1(1,k).eq.ii) go to 40
          end do
          do k=1,nsrf2
            if(isr2(1,k).eq.ii) go to 40
          end do
          do k=1,nsrf3
            if(isr3(1,k).eq.ii) go to 40
          end do
          xx=xn(1,ii)-xn(1,nn)
          yy=xn(2,ii)-xn(2,nn)
          zz=xn(3,ii)-xn(3,nn)
          d=sqrt(xx*xx+yy*yy+zz*zz)
          if(d.lt.dmin.and.ii.ne.nn) then
            dmin=d
            id=ii
          end if
40        continue
        end do
        isr1(2,i)=id
      end do
c
c...  if node lies on an edge, find the surrounding nonedge nodes
c
      do i=1,nsrf2
        ne=0
        nn=isr2(1,i)
        call iclear(isr2(2,i),9)
        do n=1,numel
          do j=1,8
            if(ien(j,n).eq.nn) then
              ne=ne+1
              ie(ne)=n
            end if
            if(ne.eq.2) go to 20
          end do
        end do
20      continue
        ll=0
        do j=1,8
          ned1=0
          ned2=0
          do k=1,nsrf2
            if(isr2(1,k).eq.ien(j,ie(1))) ned1=ned1+1
            if(isr2(1,k).eq.ien(j,ie(2))) ned2=ned2+1
          end do
          do k=1,nsrf3
            if(isr3(1,k).eq.ien(j,ie(1))) ned1=ned1+1
            if(isr3(1,k).eq.ien(j,ie(2))) ned2=ned2+1
          end do
          if(ned1.eq.0) then
            nd=0
            do l=1,ll
              if(ien(j,ie(1)).eq.isr2(l+1,i)) nd=nd+1
            end do
            if(nd.eq.0) then
              ll=ll+1
              isr2(ll+1,i)=ien(j,ie(1))
            end if
          end if
          if(ned2.eq.0) then
            nd=0
            do l=1,ll
              if(ien(j,ie(2)).eq.isr2(l+1,i)) nd=nd+1
            end do
            if(nd.eq.0) then
              ll=ll+1
              isr2(ll+1,i)=ien(j,ie(2))
            end if
          end if
        end do
      end do
c
c...  if node lies on a corner, find the surrounding nodes
c
      do i=1,nsrf3
        ne=0
        nn=isr3(1,i)
        call iclear(isr3(2,i),7)
        do n=1,numel
          do j=1,8
            if(ien(j,n).eq.nn) then
              ne=ne+1
              ie(ne)=n
            end if
            if(ne.gt.0) go to 30
          end do
        end do
30      continue
        ll=0
        do j=1,8
          if(ien(j,ie(1)).ne.nn) then
            ll=ll+1
            isr3(ll+1,i)=ien(j,ie(1))
          end if
        end do
      end do
      return
      end
c
c
      subroutine nsrfrm(xn,xnt,xlim,ien,ienp,ient,ifs,ilock,isr1,
     & isr2,isr3,nss,nf,ns)
c
c...  subroutine to form arrays necessary for stress smoothing
c     depending on whether a free surface is locked or unlocked for
c     the current time step
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/surf/nsr1,nsr2,nsr3,nsrf1,nsrf2,nsrf3
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension xn(3,*),xnt(3,*),ienp(8,*),ient(8,*),ifs(2,*),isr1(2,*)
      dimension isr2(10,*),isr3(8,*),xlim(3,*)
      dimension ien(8,*),ilock(*)
      logical*4 inc
      ii=0
      eps=0.001
      nss=numnp
      nsrf1=nsr1
      nsrf2=nsr2
      nsrf3=nsr3
      do i=1,numel
        do j=1,8
          ient(j,i)=ienp(j,i)
        end do
      end do
      do i=1,ns
        inc=.false.
        ii=ii+1
        ind1=ifs(1,i)
        ind2=ifs(2,i)
        nsurf=0
        do n=1,numel
          do k=1,8
            if(ien(k,n).eq.ind2.and.ilock(ii).ne.1) then
              if(.not.inc) nss=nss+1
              xnt(1,nss)=xn(1,ind2)
              xnt(2,nss)=xn(2,ind2)
              xnt(3,nss)=xn(3,ind2)
              ient(k,n)=nss
              inc=.true.
            end if
          end do
        end do
        if(inc) nsurf=nsurf+1
        do k=1,ndof
          if(abs(xn(k,ind1)-xlim(k,1)).lt.eps) nsurf=nsurf+1
          if(abs(xn(k,ind1)-xlim(k,2)).lt.eps) nsurf=nsurf+1
        end do
        if(nsurf.eq.1) then
          nsrf1=nsrf1+1
          isr1(1,nsrf1)=ind1
          if(inc) then
            nsrf1=nsrf1+1
            isr1(1,nsrf1)=nss
          end if
        else if(nsurf.eq.2) then
          nsrf2=nsrf2+1
          isr2(1,nsrf2)=ind1
          if(inc) then
            nsrf2=nsrf2+1
            isr2(1,nsrf2)=nss
          end if
        else if(nsurf.eq.3) then
          nsrf3=nsrf3+1
          isr3(1,nsrf3)=ind1
          if(inc) then
            nsrf3=nsrf3+1
            isr3(1,nsrf3)=nss
          end if
        end if
      end do
      return
      end
c
c
      subroutine pstrs3(sxx,syy,szz,sxy,syz,sxz,sp,ang)
c
c...  subroutine to compute the three principal stresses or strains
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension sp(*),s(3,3),ang(3,*)
      np=3
      s(1,1)=sxx
      s(2,2)=syy
      s(3,3)=szz
      s(1,2)=sxy
      s(2,1)=sxy
      s(1,3)=sxz
      s(3,1)=sxz
      s(2,3)=syz
      s(3,2)=syz
c
c...  perform eigenvalue problem and then sort eigenvalues (principal
c     stresses)
c
      call jacobi(s,np,np,sp,ang,nrot)
c
c...  temporarily change sign of principal stresses so they are sorted
c     properly
c
      sp(1)=-sp(1)
      sp(2)=-sp(2)
      sp(3)=-sp(3)
      call eigsrt(sp,ang,np,np)
      sp(1)=-sp(1)
      sp(2)=-sp(2)
      sp(3)=-sp(3)
      return
      end
c
c
      subroutine reput(vnode,ien,ient,nss,ib)
c
c...  subroutine to put smoothed stresses back onto original nodes
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension vnode(nodmx,*),ient(8,*)
      dimension ien(8,*)
      ie=ib+nstr
      do i=1,nss
        do j=1,6
          vnode(j+ie,i)=vnode(j+ib,i)
        end do
      end do
      do n=1,numel
        do k=1,8
          it=ient(k,n)
          ip=ien(k,n)
          do j=1,6
            vnode(j+ib,ip)=vnode(j+ie,it)
          end do
        end do
      end do
      return
      end
c
c
      subroutine rread(ar,idim1,idim2,idim3)
c
c...  subroutine to read a floating point array written in binary
c     format
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/units/kti,kto,kplt,kpar,kucd
      dimension ar(idim1,idim2,idim3)
      read(kplt) ar
      return
      end
c
c
      subroutine rstrn(strr,str,scl,rdiv,sflg)
c
c...  subroutine to read, scale, and convert to rates stresses and
c     strains
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/units/kti,kto,kplt,kpar,kucd
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension strr(nstr,ngstn,numel),str(numel,ngstn+1,nstr)
      logical*4 sflg
      read(kplt) strr
      if(sflg) then
        do i=1,numel
          do l=1,ngstn
            do k=1,nstr
              str(i,l,k)=scl*strr(k,l,i)/rdiv
             end do
           end do
         end do
       end if
       return
       end
c
c
      subroutine shap30(r,s,t,x,det,sh,nen,ndof,ien,infin,n)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
c
c...program to compute shape functions for a rectangular brick
c
c        r,s,t                         = natural coordinates
c        sh(1,nen),sh(2,nen),sh(3,nen) = x,y,and z derivatives
c                                          of shape functions
c        sht(nsd)                      = intermediate storage for sh
c        sh(4,nen)                     = shape functions
c*       shi(ndof+1,nen)               = infinite element shape
c*                                       functions and derivatives
c*       shtmp(ndof+1,nen)             = local shape functions and
c*                                       derivatives
c        xs(nsd,nsd)                   = jacobian matrix
c        det                           = jacobian matrix
c        x(nsd,nen)                    = local nodal coordinates
c
c         ** note: code only tests for collapse of the rectangular
c                  brick into a prism with a triangular cross section.
c                  the only collapse supported is for the degeneracy of
c                  nodes 3,4 and 7,8.
c
      common/units/kti,kto,kplt,kpar,kucd
      dimension sh(4,*),x(3,*),xs(3,3),sht(3),ien(*),ra(8),sa(8),ta(8)
      dimension shi(4,8),shtmp(4,8),idiv(3),io(3)
      data ra/-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0/
      data sa/-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0/
      data ta/ 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0/
      data idiv/1,10,100/
c
c...compute shape functions and their derivatives in natural
c   coordinates for a brick
c
      do i=1,nen
        rr=1.0+ra(i)*r
        ss=1.0+sa(i)*s
        tt=1.0+ta(i)*t
        sh(4,i)=0.125*rr*ss*tt
        sh(1,i)=0.125*ra(i)*ss*tt
        sh(2,i)=0.125*sa(i)*rr*tt
        sh(3,i)=0.125*ta(i)*rr*ss
      end do
c
c...collapse shape functions to a prism, if requested
c
      if(ien(3).eq.ien(4)) then
        do i=1,4
          sh(i,3)=sh(i,3)+sh(i,4)
          sh(i,7)=sh(i,7)+sh(i,8)
          sh(i,4)=0.0
          sh(i,8)=0.0
        end do
      end if
      call equate(shtmp,sh,nen*(ndof+1))
c
c...if element is an infinite domain element, compute alternate shape
c   functions
c
      if(infin.ne.0) then
        call clear(shi,nen*(ndof+1))
        io(3)=infin/idiv(3)
        io(2)=(infin-io(3)*idiv(3))/idiv(2)
        io(1)=infin-io(3)*idiv(3)-io(2)*idiv(2)
        do i=1,nen
          call infel(r,ra(i),rr,drr,io(1))
          call infel(s,sa(i),ss,dss,io(2))
          call infel(t,ta(i),tt,dtt,io(3))
          shi(4,i)=rr*ss*tt
          shi(1,i)=drr*ss*tt
          shi(2,i)=dss*rr*tt
          shi(3,i)=dtt*rr*ss
        end do
c
c...collapse shape functions to a prism, if requested
c
        if(ien(3).eq.ien(4)) then
          do i=1,4
            shi(i,3)=shi(i,3)+shi(i,4)
            shi(i,7)=shi(i,7)+shi(i,8)
            shi(i,4)=0.0
            shi(i,8)=0.0
          end do
        end if
        call equate(shtmp,shi,nen*(ndof+1))
      end if
c
c...calculate jacobian matrix for (x,y,z) to (r,s,t) transformation
c
      do i=1,3
        do j=1,3
          xs(i,j)=0.0
          do k=1,nen
            xs(i,j)=xs(i,j)+x(j,k)*shtmp(i,k)
          end do
        end do
      end do
c
c...form determinant of jacobean matrix and check for error condition
c
      det=xs(1,1)*xs(2,2)*xs(3,3)+xs(1,2)*xs(2,3)*xs(3,1)+xs(1,3)
     & *xs(2,1)*xs(3,2)-xs(1,3)*xs(2,2)*xs(3,1)-xs(1,2)*xs(2,1)
     & *xs(3,3)-xs(1,1)*xs(2,3)*xs(3,2)
      if(det.le.0.0) then
        write(kto,2000) det,n
	stop
      end if
 2000 format(///' shape function fails!   determinant is ',1pe20.4,
     & '  in element # ',i5)
      return
      end
c
c
      subroutine skipa(iunit)
c
c      routine to skip lines beginning with the string '/*' and blank
c      lines.  this routine ignores leading blanks before the key
c      string.  this routine is for ascii files.
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      character string*80,leader*2
      data leader/'/*'/
   10 continue
        read(iunit,1000,end=20) string
        if(nchar(string).eq.0) goto 10
        inblnk=nnblnk(string)
        if(index(string,leader).eq.inblnk) goto 10
      backspace(iunit)
 1000 format(a80)
   20 return
      end
c
c
      subroutine smooth(vnode,xn,str,smop,ien,infin,isr1,isr2,isr3,nss,
     & iadd)
c
c     subroutine to perform stress smoothing by extrapolating stresses
c     to nodal points
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      common/surf/nsr1,nsr2,nsr3,nsrf1,nsrf2,nsrf3
      dimension vnode(nodmx,*),xn(3,*),smop(*),det(8),rg(8),sg(8)
      dimension str(numel,ngstn+1,*),xl(3,8),sh(4,8),tg(8)
      dimension isr1(2,*),isr2(10,*),isr3(8,*),ien(8,*),infin(*)
      data rg/-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0/
      data sg/-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0/
      data tg/ 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0/
      call clear(smop,nss)
      ns=0
      if(ngstn.eq.8) ns=1
      g=1.0/sqrt(3.0)
      if(ngstn.eq.1) g=0.0
c
c   loop over number of elements
c
      do n=1,numel
c
c     localize nodal coordinates
c
        do i=1,8
          do j=1,3
            xl(j,i)=xn(j,ien(i,n))
          end do
        end do
c
c    compute jacobian at integration points
c
        do l=1,ngstn
          call shap30(g*rg(l),g*sg(l),g*tg(l),xl,det(l),sh,
     &     nen,ndof,ien(1,n),infin(n),n)
        end do
c
c    compute smoothing operator and weighted stresses at nodes
c
        do i=1,8
          m=ien(i,n)
          smop(m)=smop(m)+det(ns*(i-1)+1)
          do j=1,nstr
            vnode(j+iadd,m)=vnode(j+iadd,m)+str(n,ns*(i-1)+1,j)*
     &       det(ns*(i-1)+1)
          end do
        end do
      end do
c
c    invert smoothing operator and compute smoothed stresses
c
      do i=1,nss
        smop(i)=1.0/smop(i)
        do j=1,nstr
          vnode(j+iadd,i)=vnode(j+iadd,i)*smop(i)
        end do
      end do
c
c    apply boundary corrections for nonedge boundary nodes
c
      do i=1,nsrf1
        m=isr1(1,i)
        n=isr1(2,i)
        do j=1,nstr
          vnode(j+iadd,m)=2.0*vnode(j+iadd,m)-vnode(j+iadd,n)
        end do
      end do
c
c    apply boundary corrections for edge boundary nodes
c
      nintmx=9
      do i=1,nsrf2
        call extrap(xn,isr2(1,i),vnode,nintmx,iadd)
      end do
c
c    apply boundary corrections for external corner nodes
c
      nintmx=7
      do i=1,nsrf3
        call extrap(xn,isr3(1,i),vnode,nintmx,iadd)
      end do
      return
      end
c
c
      subroutine ssub(strn,strs)
c
c...  subroutine to subtract one set of stresses/strains from another
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      dimension strn(numel,ngstn+1,nstr),strs(nstr,ngstn,numel)
      do i=1,numel
        do l=1,ngstn
          do k=1,nstr
            strn(i,l,k)=strn(i,l,k)-strs(k,l,i)
          end do
        end do
      end do
      return
      end
c
c
      subroutine stncmp(vnode,ipoint)
c
c...  subroutine to compute stress/strain related quantities
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension vnode(*),ang(3,3),sp(3)
      sxx=vnode(ipoint)
      syy=vnode(ipoint+1)
      szz=vnode(ipoint+2)
      sxy=vnode(ipoint+3)
      syz=vnode(ipoint+4)
      sxz=vnode(ipoint+5)
      call pstrs3(sxx,syy,szz,sxy,syz,sxz,sp,ang)
      vnode(ipoint+6)=sp(1)
      vnode(ipoint+7)=sp(2)
      vnode(ipoint+8)=sp(3)
      vnode(ipoint+9)=-(sp(1)+sp(2)+sp(3))/3.0
      vnode(ipoint+10)=abs(sp(3)-sp(1))*0.5
      vnode(ipoint+11)=sqrt(0.5*(sxx*sxx+syy*syy+szz*szz)+sxy*sxy+
     & syz*syz+sxz*sxz)
      return
      end
c
c
      subroutine stncmpe(strn,iel,igp,ip)
c
c...  wrapper subroutine for stncmp that deals with reversed-order array
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      dimension strn(numel,ngstn+1,nstmx),strs(12)
      do k=1,6
	strs(k)=strn(iel,igp,ip+k-1)
      end do
      ib=1
      call stncmp(strs,ib)
      do k=7,12
	strn(iel,igp,ip+k-1)=strs(k)
      end do
      return
      end
c
c
      subroutine stread(xn,xnt,xd,strn,smop,dld,strr,ien,infin,ienp,
     & ient,isr1,isr2,isr3,vnode,numsn,numfn,numflt,ifs,nfault,fault,
     & ilock,prop,vx,vy,vz,rlam,itmp,itmp2,idftn,nf,ns,nfe,nprop,time,
     & amag,xlim,ivisc,iplas)
c
c...  subroutine to read values for given time step.  all values are
c     extrapolated to nodal points.
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      common/dimens/npmax,nemax,nmmax,nodmx,nstmx,nfpmx,
     & nefmx,noutmx,nlnmx
      common/elmt/numnp,nnp,numel,nmat,nstr,ngstn,ndof,ngem,nen
      common/scale/cscl,uscl,vscl,wscl,sscl,escl,descl,tscl
      common/sflags/cstrs,cestn,cpstn,cvstn,cdestn,cdpstn,cdvstn
      common/vpoint/ivp1,ivp2,ivp3,ivp4,ivp5,ivp6,ivp7,ivp8,ivp9,ivp10,
     & ivp11,ivp12,ivpend
      common/spoint/isp(8)
      common/units/kti,kto,kplt,kpar,kucd
      dimension xn(3,*),xnt(3,*),xd(3,*),strn(numel,9,nstmx),smop(*)
      dimension isr1(2,*),isr2(10,*),isr3(8,*),vnode(nodmx,*),ienp(8,*)
      dimension ifs(2,*),dld(3,numnp),ient(8,*),fault(6,*),xlim(3,*)
      dimension strr(nstr,ngstn,numel),itmp(*),prop(nprop,*)
      dimension ien(8,*),infin(*),nfault(3,*),ilock(*),itmp2(*)
      dimension idftn(numflt)
      logical*4 cstrs,cestn,cpstn,cvstn,cdestn,cdpstn,cdvstn,sflg,smt
      call clear(vnode,nodmx*npmax)
      call clear(fault,6*nfpmx)
      read(kplt) itmp2(1)
      nstep=itmp2(1)
      read(kplt) time
      ia=0
      iaa=0
      iaf=1
      pi=acos(-1.0)
      drad=pi/180.0
      eps=1.0e-6
      smt=cstrs.or.cestn.or.cpstn.or.cvstn.or.cdestn.or.cdpstn.or.cdvstn
c
c...  read current displacements
c
      write(kto,*) 'Reading current nodal displacements:'
      read(kplt) dld
      do i=1,numnp
        vt=0.0
        do j=1,ndof
          vnode(j,i)=dld(j,i)
          vt=vt+vnode(j,i)*vnode(j,i)
        end do
        vnode(4,i)=sqrt(vt)
        vnode(ivp3,i)=vx*vnode(1,i)+vy*vnode(2,i)+vz*vnode(3,i)
      end do
      if(nstep.gt.0) then
        ia=4
	iaa=3
        iaf=4
        read(kplt) deltp
        read(kplt) dld
        do i=1,numnp
          vt=0.0
          do j=1,ndof
            vnode(ivp2-1+j,i)=dld(j,i)/deltp
            vt=vt+vnode(ivp2-1+j,i)*vnode(ivp2-1+j,i)
          end do
          vnode(8,i)=sqrt(vt)
          vnode(ivp4,i)=vx*vnode(5,i)+vy*vnode(6,i)+vz*vnode(7,i)
        end do
      end if
c
c...  add displacements/velocities due to movement on split nodes.
c
      call clear(fault,6*nfpmx)
      do i=1,2*numflt
        if(nstep.eq.0) read(kplt) (fault(j,i),j=1,3)
        if(nstep.gt.0) read(kplt) (fault(j,i),j=1,6)
      end do
      do i=1,nf
	ilock(i)=1
	ip1=ifs(1,i)
	ip2=ifs(2,i)
	do k=1,ivp4-1
	  vnode(k,ip2)=vnode(k,ip1)
	end do
      end do
      do i=1,nf
	ip1=ifs(1,i)
	ip2=ifs(2,i)
	do j=1,numflt
	  if(ip1.eq.idftn(j)) then
	    j1=(j-1)*2+1
	    j2=j1+1
	    ux=fault(1+iaa,j1)
	    uy=fault(2+iaa,j1)
	    uz=fault(3+iaa,j1)
	    v=sqrt(ux*ux+uy*uy+uz*uz)
            if(v.gt.eps) ilock(i)=0
            vnode(1,ip1)=vnode(1,ip1)+fault(1,j1)
            vnode(2,ip1)=vnode(2,ip1)+fault(2,j1)
            vnode(3,ip1)=vnode(3,ip1)+fault(3,j1)
            vnode(5,ip1)=vnode(5,ip1)+fault(4,j1)
            vnode(6,ip1)=vnode(6,ip1)+fault(5,j1)
            vnode(7,ip1)=vnode(7,ip1)+fault(6,j1)
	    vnode(4,ip1)=sqrt(vnode(1,ip1)*vnode(1,ip1)+
     &       vnode(2,ip1)*vnode(2,ip1)+vnode(3,ip1)*vnode(3,ip1))
	    vnode(8,ip1)=sqrt(vnode(5,ip1)*vnode(5,ip1)+
     &       vnode(6,ip1)*vnode(6,ip1)+vnode(7,ip1)*vnode(7,ip1))
            vd=vx*fault(1,j1)+vy*fault(2,j1)+vz*fault(3,j1)
            vv=vx*fault(4,j1)+vy*fault(5,j1)+vz*fault(6,j1)
            vnode(9,ip1)=vnode(9,ip1)+vd
            vnode(10,ip1)=vnode(10,ip1)+vv
            vnode(1,ip2)=vnode(1,ip2)+fault(1,j2)
            vnode(2,ip2)=vnode(2,ip2)+fault(2,j2)
            vnode(3,ip2)=vnode(3,ip2)+fault(3,j2)
            vnode(5,ip2)=vnode(5,ip2)+fault(4,j2)
            vnode(6,ip2)=vnode(6,ip2)+fault(5,j2)
            vnode(7,ip2)=vnode(7,ip2)+fault(6,j2)
	    vnode(4,ip2)=sqrt(vnode(1,ip2)*vnode(1,ip2)+
     &       vnode(2,ip2)*vnode(2,ip2)+vnode(3,ip2)*vnode(3,ip2))
	    vnode(8,ip2)=sqrt(vnode(5,ip2)*vnode(5,ip2)+
     &       vnode(6,ip2)*vnode(6,ip2)+vnode(7,ip2)*vnode(7,ip2))
            vd=vx*fault(1,j2)+vy*fault(2,j2)+vz*fault(3,j2)
            vv=vx*fault(4,j2)+vy*fault(5,j2)+vz*fault(6,j2)
            vnode(9,ip2)=vnode(9,ip2)+vd
            vnode(10,ip2)=vnode(10,ip2)+vv
	    go to 5
          end if
	end do
5       continue
      end do
c
c...  add displacements/velocities due to movement on slippery nodes.
c
      write(kto,*) 'Adding contributions from slippery nodes:'
      call clear(fault,6*nfpmx)
      do i=1,numsn
        if(nstep.eq.0) read(kplt) m,(fault(j,i),j=1,3)
        if(nstep.gt.0) read(kplt) m,(fault(j,i),j=1,6)
        itmp(i)=m
      end do
      do i=1,numsn
        read(kplt) itmp2(i)
      end do
      do i=1,ns
        ip1=ifs(1,i+nf)
        ip2=ifs(2,i+nf)
        do k=1,ivp4-1
          vnode(k,ip2)=vnode(k,ip1)
        end do
      end do
      do i=1,ns
        ip1=ifs(1,i+nf)
        ip2=ifs(2,i+nf)
	ii=0
	do k=1,numsn
	  if(itmp(k).eq.ip1) then
	    ii=k
	    ilock(ii)=itmp2(k)
	  end if
	end do
        vnode(1,ip1)=vnode(1,ip1)+fault(1,ii)
        vnode(2,ip1)=vnode(2,ip1)+fault(2,ii)
        vnode(3,ip1)=vnode(3,ip1)+fault(3,ii)
        vnode(5,ip1)=vnode(5,ip1)+fault(4,ii)
        vnode(6,ip1)=vnode(6,ip1)+fault(5,ii)
        vnode(7,ip1)=vnode(7,ip1)+fault(6,ii)
        vnode(1,ip2)=vnode(1,ip2)-fault(1,ii)
        vnode(2,ip2)=vnode(2,ip2)-fault(2,ii)
        vnode(3,ip2)=vnode(3,ip2)-fault(3,ii)
        vnode(5,ip2)=vnode(6,ip2)-fault(4,ii)
        vnode(6,ip2)=vnode(6,ip2)-fault(5,ii)
        vnode(7,ip2)=vnode(7,ip2)-fault(6,ii)
        vnode(4,ip1)=sqrt(vnode(1,ip1)*vnode(1,ip1)+
     &   vnode(2,ip1)*vnode(2,ip1)+vnode(3,ip1)*vnode(3,ip1))
        vnode(8,ip1)=sqrt(vnode(5,ip1)*vnode(5,ip1)+
     &   vnode(6,ip1)*vnode(6,ip1)+vnode(7,ip1)*vnode(7,ip1))
        vnode(4,ip2)=sqrt(vnode(1,ip2)*vnode(1,ip2)+
     &   vnode(2,ip2)*vnode(2,ip2)+vnode(3,ip2)*vnode(3,ip2))
        vnode(8,ip2)=sqrt(vnode(5,ip2)*vnode(5,ip2)+
     &   vnode(6,ip2)*vnode(6,ip2)+vnode(7,ip2)*vnode(7,ip2))
        vd=vx*fault(1,ii)+vy*fault(2,ii)+vz*fault(3,ii)
        vv=vx*fault(4,ii)+vy*fault(5,ii)+vz*fault(6,ii)
        vnode(9,ip1)=vnode(9,ip1)+vd
        vnode(10,ip1)=vnode(10,ip1)+vv
        vnode(9,ip2)=vnode(9,ip2)-vd
        vnode(10,ip2)=vnode(10,ip2)-vv
      end do
c
c...compute minimum directional displacement as reference for phase
c   values
c
      d0=1.0e20
      do i=1,nnp
        d0=min(d0,vnode(9,i))
      end do
      do i=1,nnp
        vp=2.0*(vnode(9,i)-d0)/rlam
        nw=int(vp)
        vnode(11,i)=rlam*0.5*(vp-float(nw))
      end do
c
c...  update nodal coordinates
c
      write(kto,*) 'Updating nodal coordinates:'
      call update(xn,xd,vnode,nnp,amag,cscl)
c
c...  form stress smoothing arrays
c
      if(smt) then
        write(kto,*) 'Forming stress smoothing arrays:'
        call nsrfrm(xn,xnt,xlim,ien,ienp,ient,ifs,ilock,isr1,
     &   isr2,isr3,nss,nf,ns)
        call nsrfnd(xnt,ient,isr1,isr2,isr3)
      end if
c
c...  read stresses, strains, and strain rates
c
      write(kto,*) 'Reading stresses and strains:'
      rone=1.0
      sflg=cstrs
      call rstrn(strr,strn(1,1,isp(1)),sscl,rone,sflg)
      sflg=cestn
      call rstrn(strr,strn(1,1,isp(2)),escl,rone,sflg)
      if(nstep.gt.0) then
	if(ivisc.eq.1) then
          sflg=cvstn
          call rstrn(strr,strn(1,1,isp(3)),escl,rone,sflg)
	  if(cestn) call ssub(strn(1,1,isp(2)),strr)
        end if
	if(iplas.eq.1) then
          sflg=cpstn
          call rstrn(strr,strn(1,1,isp(4)),escl,rone,sflg)
	  if(cestn) call ssub(strn(1,1,isp(2)),strr)
        end if
        sflg=cdestn
        call rstrn(strr,strn(1,1,isp(5)),descl,deltp,sflg)
	if(ivisc.eq.1) then
          sflg=cdvstn
          call rstrn(strr,strn(1,1,isp(6)),descl,deltp,sflg)
	  if(cdestn) call ssub(strn(1,1,isp(5)),strr)
        end if
	if(iplas.eq.1) then
          sflg=cdpstn
          call rstrn(strr,strn(1,1,isp(7)),descl,deltp,sflg)
	  if(cdestn) call ssub(strn(1,1,isp(5)),strr)
        end if
      end if
c
c...  smooth desired stresses and strains and extrapolate them to nodes
c
      write(kto,*) 'Smoothing stresses and strains:'
      if(cstrs) then
	ib=ivp6-1
        call smooth(vnode,xnt,strn(1,1,isp(1)),smop,ient,infin,isr1,
     &   isr2,isr3,nss,ib)
	call reput(vnode,ien,ient,nss,ib)
      end if
      if(cestn) then
	ib=ivp7-1
        call smooth(vnode,xnt,strn(1,1,isp(2)),smop,ient,infin,isr1,
     &   isr2,isr3,nss,ib)
	call reput(vnode,ien,ient,nss,ib)
      end if
      if(nstep.gt.0) then
        if(cvstn) then
	  ib=ivp8-1
          call smooth(vnode,xnt,strn(1,1,isp(3)),smop,ient,infin,isr1,
     &     isr2,isr3,nss,ib)
	  call reput(vnode,ien,ient,nss,ib)
        end if
        if(cpstn) then
	  ib=ivp9-1
          call smooth(vnode,xnt,strn(1,1,isp(4)),smop,ient,infin,isr1,
     &     isr2,isr3,nss,ib)
	  call reput(vnode,ien,ient,nss,ib)
        end if
        if(cdestn) then
	  ib=ivp10-1
          call smooth(vnode,xnt,strn(1,1,isp(5)),smop,ient,infin,isr1,
     &     isr2,isr3,nss,ib)
	  call reput(vnode,ien,ient,nss,ib)
        end if
        if(cdvstn) then
	  ib=ivp11-1
          call smooth(vnode,xnt,strn(1,1,isp(6)),smop,ient,infin,isr1,
     &     isr2,isr3,nss,ib)
	  call reput(vnode,ien,ient,nss,ib)
        end if
        if(cdpstn) then
	  ib=ivp12-1
          call smooth(vnode,xnt,strn(1,1,isp(7)),smop,ient,infin,isr1,
     &     isr2,isr3,nss,ib)
	  call reput(vnode,ien,ient,nss,ib)
        end if
      end if
c
c...  compute stress/strain-related quantities at nodes
c
      write(kto,*) 'Computing stress/strain quantities at nodes:'
      do i=1,nnp
	if(cstrs) call stncmp(vnode(1,i),ivp6)
	if(cestn) call stncmp(vnode(1,i),ivp7)
	if(nstep.gt.0) then
	  if(cvstn) call stncmp(vnode(1,i),ivp8)
	  if(cpstn) call stncmp(vnode(1,i),ivp9)
	  if(cdestn) call stncmp(vnode(1,i),ivp10)
	  if(cdvstn) call stncmp(vnode(1,i),ivp11)
	  if(cdpstn) call stncmp(vnode(1,i),ivp12)
	 end if
      end do
c
c...  compute and save stresses/strains at element centroids if they
c     are not already computed there
c
      write(kto,*) 'Computing stress/strain quantities at centroids:'
      if(ngstn.gt.1) then
	if(cstrs) call centint(xn,ien,strn(1,1,isp(1)),infin)
	if(cestn) call centint(xn,ien,strn(1,1,isp(2)),infin)
	if(nstep.gt.0) then
	  if(cvstn) call centint(xn,ien,strn(1,1,isp(3)),infin)
	  if(cpstn) call centint(xn,ien,strn(1,1,isp(4)),infin)
	  if(cdestn) call centint(xn,ien,strn(1,1,isp(5)),infin)
	  if(cdvstn) call centint(xn,ien,strn(1,1,isp(6)),infin)
	  if(cdpstn) call centint(xn,ien,strn(1,1,isp(7)),infin)
        end if
      end if
c
c...  compute stress-related quantities at element centroids
c
      ll=1
      if(ngstn.gt.1) ll=9
      ib=1
      do i=1,numel
	if(cstrs) call stncmpe(strn,i,ll,isp(1))
	if(cestn) call stncmpe(strn,i,ll,isp(2))
	if(nstep.gt.0) then
	  if(cvstn) call stncmpe(strn,i,ll,isp(3))
	  if(cpstn) call stncmpe(strn,i,ll,isp(4))
	  if(cdestn) call stncmpe(strn,i,ll,isp(5))
	  if(cdvstn) call stncmpe(strn,i,ll,isp(6))
	  if(cdpstn) call stncmpe(strn,i,ll,isp(7))
        end if
      end do
c
c...  scale results
c
      write(kto,*) 'Scaling results:'
      do i=1,nnp
        do j=1,ivp6-1
          scl=uscl
          if(j.ge.ivp2) scl=vscl
          if(j.eq.ivp3) scl=uscl
          if(j.eq.ivp4) scl=vscl
          if(j.eq.ivp5) scl=wscl
          vnode(j,i)=scl*vnode(j,i)
        end do
      end do
      return
      end
c
c
      subroutine update(xn,xd,vnode,nnp,amag,cscl)
c
c...  subroutine to add current displacements to nodal coordinates
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      dimension xn(3,*),xd(3,*),vnode(6,*)
      do i=1,nnp
        do j=1,3
          xd(j,i)=cscl*(xn(j,i)+amag*vnode(j,i))
        end do
      end do
      return
      end
