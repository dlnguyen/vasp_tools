***********************    program   optic     ************************
*                                                                     *
*    Function    : Calculate optical constants of the material        *
*    Available   : Dielectric constants, SHG and LEO coefficients     *
*    Author      : Chun-gang Duan                                     * 
*    Last updates: Mar 6th, 2002 : dielectric constants               * 
*                  Apr 30th,2001 : SHG coefficients (method 2)        *  
*                  Mar 28th,2001 : SHG coefficients (method 1)        *  
*                  Jan 2003 (gyg, ntu):  modified to take inputs from *
*                                        a PAW(VASP) or FLAPW(WIEN2k) *
*                                        calculation.                 *
*                  Mar 2003 (gyg, ntu):  modified for metals and      *
*                                                     nonmetals.      *
*                  Sep 2003 (gyg, ntu):  modified for magnetic cases. *
*                                        optical conductivity and     *
*                                        EELS spectra now printed out *
*    Reference   : 1. Claudio Aversa and J.E. Sipe,                   *
*                     PRB, Vol. 52, 14636 (1995)                      *
*                  2. James L.P. Hughes and J.E. Sipe,                *
*                     PRB, Vol. 53, 10751 (1996)                      *
*                  3. S.N. Rashkeev, W.R.L. Lambrecht and B. Segall   *
*                     PRB, Vol. 57, 3905 (1998)                       *
*                  4. C.G. Duan, J. Li, Z.Q. Gu, D.S. Wang            *
*                     PRB, Vol. 60, 9435 (1999)                       *
*                  5. J. Lin et. al.                                  *
*                     PRB, Vol. 60, 13380 (1999)                      *
***********************************************************************
      program OPTIC

c     include 'param.h'
      implicit real*8 (a-h,o-z)
      parameter (ND=7)
      parameter (NBD=700)
c     Define Constants
      parameter (Pi=3.14159265359d0)
      parameter (Hartree=27.2116d0)
      parameter (Hbar=1.0545887d0)
      parameter (CC=2.99793d0)
      parameter (CCau=1.37037d2)
      parameter (EE=1.6021892d0)
      parameter (Bohr=5.2917706d0)
      parameter (ESU=5.830345d0)

      logical ortho
c
c     FILE  1: define case name, 'opticpack.def' 
c     FILE  5: input,  'case.opticin'
c     FILE  6: output, 'case.opticout'
c     FILE  7: momentum matrix elements, 'case.mme'
c              P_mn=<Psi_m|P|Psi_n>, Psi: Bloch (eigen) Function
c                                    P  : momentum operator
c              This MME file is actually basis independent, i.e,
c              theoretically it should be the same for different
c              energy band methods. Details about the format please
c              refer to file 'readme'. 
c     FILE 11: 1)imaginary part of dielectric constant, 'case.epsim'
c              2)imaginary part of Chi_2, 'case.chiim'
c     FILE 12: 1)real part of dielectric constant, 'case.epsre'
c              2)real part of Chi_2, 'case.chire'
c     FILE 13: 1)module of dielectric constant, 'case.epsabs'
c              2)module of Chi_2, 'case.chiabs'
c     FILE 14: optical conductivity (real part), 'case.sigre'
c     FILE 15: optical conductivity (imaginary part), 'case.sigim'
c     FILE 16: refractive index (n), 'case.refra'
c     FILE 17: linear absorption coefficients, 'case.alpha1'
c     FILE 18: joint density of states, 'case.jdos'
c     FILE 20
c       ~ 100: reserved
c

c     timer variables
c      real*8 cp,cp0,cp1
c     arrays
      double complex, allocatable ::  mme(:,:,:)    !(3,nb,nb)
c     real*8, allocatable, dimension (:) :: en,ev,jdos,jdos1,jdos2 ! ne, nb
      real*8, allocatable, dimension (:) :: en,ev,jdos ! ne, nb
      real*8, allocatable, dimension (:,:) :: chi_im, chi_re, chi_abs !ne, l
      real*8, allocatable, dimension (:,:) :: refra, optc, alpha !ne, l
      real*8, allocatable, dimension (:,:) :: eels,  optcim      !ne, l
      real*8  chi_zero(18), pcom(18),pcom2(18)
      real*8  sym(3,3,48), trx(3,3,48)
      real*8  g(3,3), rg(3,3),qk(3), ptensor(3,3),ptensor2(3,3,3)
      double complex ptmp,ppp
      dimension icom(18),jcom(18)
      dimension imat(3,3,48),qkp(3),symp(3,3)
      dimension nbcal(2),nbmin(2),nbmax(2)

c     nb     : number of bands
c     ne     : number of energy sampling points

c     start timer
c      call timer(cp0)
c      cp0 = 0.d0

c     open file
      call start_optic(ioption)

c     output title
	write(6,1001) 

	read(5,*) ifstatic
	read(5,*) iprint
	if(iprint.gt.3.or.iprint.lt.1) iprint=1
        read(5,*) ival,ef,imetal,ispin,iso,itr
        nspin = 1
        if(ispin.eq.2) nspin = 2
        read(5,*) emin, erange, de, sigma, nbcal0
        do is = 1,nspin
        nbcal(is) = nbcal0
        enddo
c     ifstatic:  if calculate static value only
c     iprint  :  print level, 1:low 2:mid 3:high
cgyg  ival    :  number of valence bands (not used for a metal)
cgyg  ef      :  fermi energy (eV) (not used for a non-metal)
cgyg  imetal  :  1 for a metal; otherwise for a non-metal
cgyg  a non-metal can also be treated as a metal.
c     emin    :  starting energy (eV)
c     erange  :  energy range (eV)
c     de      :  energy increment (eV)
c     sigma   :  Gaussian factor (delta function) (eV)
c     nbcal   :  number of bands to be considered, could be changed in
c                later calculations due to matrix size

      ne=erange/de+1

      read(5,*) isci
c     choice of using scissors operator
c     isci   :  0     do not use
c               1     naive  form, only shift energy
c               2     also renormalize momentum matrix       
      eshift=0.d0
      if(isci.gt.0) read(5,*) eshift
c     eshift :  energy shift (eV)

c     choose components to be calculated
c     (1) for dielectric constants: chi_(icom)     
c         icom :  
c              1......xx
c              2......yy
c              3......zz
c              4......yz(zy)
c              5......zx(xz)
c              6......xy(yx)
c     (2) for SHG coefficients :    chi_(icom,jcom)
c         icom:             jcom:
c              1......x          1......xx
c              2......y          2......yy
c              3......z          3......zz
c                                4......yz(zy)
c                                5......zx(xz)
c                                6......xy(yx)
c     (3) same as (2), different formula, static value only

      read(5,*) itot

      if(ioption.eq.1) then
         write(6, 1005) ioption
         if(itot.gt.6) 
     &   stop 'Error: At most 6 components to be calculated!'
         do i=1, itot
            read(5,*) icom(i)
            if(icom(i).gt.6.or.icom(i).lt.1) 
     &      stop 'Error: Wrong components'
         enddo
         if(itot.eq.0) then 
            itot=6
            do j=1,6
               icom(j)=j
            enddo
         endif
      elseif(ioption.eq.2.or.ioption.eq.3) then
         write(6, 1006) ioption
	   if(itot.gt.18) 
     &   stop 'Error: At most 18 components to be calculated!'
         do i=1, itot
            read(5,*) icom(i),jcom(i)
            if(icom(i).gt.3.or.icom(i).lt.1.or.jcom(i).lt.1.or.
     &      jcom(i).gt.6) stop 'Error: Wrong components'
         enddo
         if(itot.eq.0) then 
            itot=18
            jcount=0
            do i=1,3
               do j=1,6
                  jcount=jcount+1
                  icom(jcount)=i
	            jcom(jcount)=j
	         enddo
            enddo
         endif
      else 
         stop 'No this option so far!'
      endif

c     display input values
	if(ifstatic.ne.0) write(6,1050)
      write(6,1100) emin, erange, de, sigma 
      write(6,1105) ne, iso, itr
      write(6,1110) isci, eshift
	if(ioption.eq.1) write(6,1111) itot,(icom(i),i=1,itot)
	if(ioption.eq.2) write(6,1112) itot,(icom(i),jcom(i),i=1,itot)
	if(ioption.eq.3) write(6,1112) itot,(icom(i),jcom(i),i=1,itot)
	write(6,*) 

      nunit=7
c     read structure information 
cgyg  read(nunit) nk, omega, ival, nsym
cgyg  read(nunit) g, rg
      read(2,*) omega
      read(2,*) alat
      read(2,*) ortho
      do j = 1,3
      read(2,*) (g(i,j),i=1,3)
      enddo     
c     change into a.u.
c     alat = alat/0.52918d0
c     do j = 1,3
c       do i = 1,3
c         g(i,j) = g(i,j)/0.52918d0
c       enddo
c     enddo
      call recip(g,rg,omega0)    
c
      icon=ival+1
c     nk     : number of k-points
c     omega  : cell volume (a.u.)
c     ival   : number of valence bands
c     nsym   : number of symmetric operations
c     g(rg)  : (reciprocal) lattice vectors (a.u.)
c              g1: (g(1,1),g(2,1),g(3,1))
c              g2: (g(1,2),g(2,2),g(3,2))
c              g3: (g(1,3),g(2,3),g(3,3))
c              note: g*(rg~)=E

	call getnb(nbmin,nbmax,nunit,nk,nspin)
	if(iprint.gt.1) then
        do is = 1,nspin
      write(6,'(" spin:",i5)') is
      write(6,1114) nbmin(is),nbmax(is)
        enddo
	endif

c     nbmin     : minimum number of bands of one k point
c     nbmax     : maximum number of bands of one k point

      do is = 1,nspin
      if(nbcal(is).gt.nbmin(is).or.nbcal(is).lt.icon) then
         write(6,'(" spin:",i5)') is
         write(6,1115) nbcal(is),nbmin(is)
         nbcal(is)=nbmin(is)
      endif
      enddo
      nbmaxx = nbmax(1)
      if (nspin.eq.2) nbmaxx = max0(nbmax(1),nbmax(2))

c     read space group information
cgyg  read(nunit) (((sym(i,j,k),i=1,3),j=1,3),k=1,nsym)
      read(2,*) nsym
      do k = 1,nsym
        do j = 1,3
          read(2,'(3i2)') (imat(i,j,k),i=1,3)
          if( .not. ortho ) then
          do i = 1,3
           trx(i,j,k) = dble(imat(i,j,k))
          enddo
          endif
        enddo
          read(2,*) ndum
      enddo
cgyg
      if (ortho) then
        do k = 1,nsym
          do i = 1,3
          do j = 1,3
            sym(i,j,k) = dble(imat(i,j,k))
          enddo
          enddo
        enddo
      endif

      if(iprint.gt.1) write(6,1120) g,rg
      if(imetal.eq.1) then
      write(6,1121) omega, nbcal0, ef, nsym
      else
      write(6,1122) omega, nbcal0, ival, nsym
      endif

      if (ortho) then
      call transgroup1(nsym,sym,trx,g,rg)
      else
      call transgroup(nsym,sym,trx,g,rg)
      endif

        if(iprint.gt.1) then
      do k=1,nsym
         write(6,1200) k, ((trx(i,j,k),j=1,3),i=1,3)
         write(6,1201) ((sym(i,j,k),j=1,3),i=1,3)
      enddo
        endif

c     allocate memory
      ALLOCATE (en(ne), chi_im(ne,itot), chi_re(ne,itot), 
     &          chi_abs(ne,itot))
      ALLOCATE (ev(nbmaxx), mme(3,nbmaxx,nbmaxx))
      ALLOCATE (refra(ne,itot),optc(ne,itot),optcim(ne,itot),
     &          eels(ne,3),alpha(ne,itot))

c     initialize arrays
      do i=1,ne
         en(i)=emin+de*(i-1)
      enddo
      nsize=ne*itot
      call inidouble(chi_im,  nsize)
      call inidouble(chi_re,  nsize)
      call inidouble(chi_abs, nsize)
      call inidouble(chi_zero,itot )
      call inidouble(pcom,itot )
      
c     determine unit and factors
      gaussfac=sqrt(Pi)*sigma/Hartree
cgyg  change to Hartree unit now to save time
      sigma = sigma/Hartree
      do i=1,ne 
         en(i)=en(i)/Hartree
      enddo
      if(ioption.eq.1) then
c        ALLOCATE (jdos(ne),jdos1(ne),jdos2(ne))
         ALLOCATE (jdos(ne))
         do is = 1,nspin
         do ie = 1,ne
         jdos(ie) = 0.d0
c        jdos1(ie) = 0.d0
c        jdos2(ie) = 0.d0
         enddo
         enddo
         unit=1.d0
         factor=16.d0*Pi/omega
         if(iso.eq.1) factor=0.5d0*factor
      elseif(ioption.eq.2.or.ioption.eq.3) then
         unit=1.0d-2
	   fac=ESU*Pi/omega
           if(iso.eq.1) fac=0.5d0*fac
      endif

c     set k points counter
      ktot=0

      nunit=7
      read(nunit,'(I6)') nk

      do 900 ik=1,nk

c     ini mme
      call inicomplex(mme,nbmaxx*nbmaxx*3)
      
c     read k-point info (qk:coordinates, wfk:weight)
      read(nunit,'(4F14.8)') qk(1),qk(2),qk(3),wfk
cgyg  qk must be in units of (primitive vectors)/2pi of reciprocal lattice. 
c
      write(6,1300) ik,qk,wfk

      do 900 is = 1,nspin

c     read eigen values (eV)
      read(nunit,'(I6)') nb
      read(nunit,'(5E13.5)') (ev(j),j=1,nb)
	if(iprint.gt.2) then
      write(6,'(" spin:",i5)') is
      write(6,1310) (ev(j),j=1,nb)
      endif

cgyg  read momentum matrix elements (in units of eV/bohr radius)
      do j=1,nb
         do i=j,nb
            read(nunit,'(6E13.5)') (mme(ntm,i,j),ntm=1,3)
c          do ntm=1,3
c           mme(ntm,i,j)=mme(ntm,i,j)*(0.d0,-1.d0)
c          enddo
         enddo
      enddo

c     unfold mme
      do j=1,nb
         do i=1,j-1
            mme(1,i,j)=dconjg(mme(1,j,i))
            mme(2,i,j)=dconjg(mme(2,j,i))
            mme(3,i,j)=dconjg(mme(3,j,i))
         enddo
      enddo
      write(6,*) 'read mme'
cgyg for a metal, determine the top valence band (ival) and bottom
cgyg conduction band (icon) using the fermi energy (ef)

      if (imetal.eq.1) then
          ival = nb
        do i = 1,nb
          if (ev(i).gt.ef) then
          ival = i - 1
          goto 10
          endif
        enddo
        write(6,'(" warning: no conduction band",2i5)')
     &           ival,nb
 10     if(ival.eq.0) 
     &  write(6,'(" warning: no valence band",i5)') ival
        icon = ival + 1
        write(6,'(" ival =",i5," icon =",i5)') ival,icon
      endif

c     scissors operator 
      if(isci.gt.0) then
c     momentum renormalization
         if(isci.eq.2) then
            do j=icon,nbcal(is)
               do i=1,ival
                  do k=1,3
                     shiftfac=eshift/(ev(j)-ev(i))
                     mme(k,i,j)=mme(k,i,j)*(1+shiftfac)
                     mme(k,j,i)=dconjg(mme(k,i,j))
                  enddo
               enddo
            enddo
         endif
c        energy shift 
         do i=icon,nbcal(is)
            ev(i)=ev(i)+eshift
         enddo
      endif

c     find the group of k and its star 
      call subgroup(qk,nsym,trx,nsubt)
      
c     count total k points
cgyg  nsubt: rank of the group of k; ndivt: number of the k stars.
      ndivt=nsym/nsubt
      if(itr.eq.1) ndivt = 1
      ktot=ktot+ndivt
      write(6,'(" nsubt =",i5," ndivt =",i5)') nsubt,ndivt


100   continue
c     calculate the linear susceptibility
      do 110 i=1,ival
      do 110 j=icon,nbcal(is)
      eji=ev(j)-ev(i)
      eji=eji/Hartree

c     construct P tensor
      do jj=1,3
         do ii=1,3
            ptmp=mme(ii,i,j)*mme(jj,j,i)
            preal=dreal(ptmp)
            ptensor(ii,jj)=preal
         enddo
      enddo

      if(itr.eq.1) then
      do jj=1,3
         do ii=1,3
            ptensor(ii,jj)=ndivt*ptensor(ii,jj)
         enddo
      enddo
      else
c     symmetrization 
      call symptensor(ptensor,sym,nsym,nsubt)
      endif
      
      wdeno=factor/(eji*eji*eji)
      wdeno2=factor*Pi*0.5d0/(eji*eji)

c     find corresponding components
      call findcom(pcom,icom,ptensor,itot)

c     calculate the imaginary part
      if(ifstatic.eq.0) then
      do ie=1,ne
c        ehar=en(ie)/Hartree
c        edi=eji-ehar
         edi=eji-en(ie)
c        expa=edi*Hartree/sigma
         expa=edi/sigma
         expa=expa*expa
c        delta function
         if(expa.lt.50.0) then
c           ex=exp(-expa)/gaussfac
            ex=exp(-expa)
            ex=ex*wdeno2
            do jc=1,itot
               chi_im(ie,jc)=chi_im(ie,jc)+pcom(jc)*ex
            enddo
            jdos(ie) = jdos(ie)+ex
         endif
       enddo
       endif
      
c     calculate chi(0) directly 
      do jc=1,itot
         chi_zero(jc)=chi_zero(jc)+pcom(jc)*wdeno
      enddo

110   continue
900   continue
 
c     output totoal k points
      write(6,1350) ktot

      if(ioption.eq.3) then
      do j=1,itot
         chi_zero(j)=chi_zero(j)/ktot
      enddo
      goto 910
      endif

c     normalize results
      do j=1,itot
         chi_zero(j)=chi_zero(j)/ktot
         do ie=1,ne
            chi_im(ie,j)=chi_im(ie,j)/gaussfac/ktot
         enddo
      enddo
      if(ioption.eq.1) then
         do ie=1,ne
            jdos(ie)=jdos(ie)/gaussfac/ktot
         enddo
c        do ie=1,ne
c           jdos1(ie)=jdos(ie)/en(ie)
c           jdos2(ie)=jdos(ie)/en(ie)**2
c        enddo
      endif 

cgyg change back to eV unit.
      sigma = sigma*Hartree
      do ie =1,ne
         en(ie)=en(ie)*Hartree
      enddo
 
c     using Kramers-Kronig relation to obtain the real part of
c     the dielectric function
      if(ifstatic.eq.0) then
         call kk_trans(en,ne,de,chi_im,chi_re,itot)
	endif

	if(ioption.eq.1) then
        do j=1, itot
         if(icom(j).le.3) then
         chi_zero(j)=chi_zero(j)+1.d0
	   do ie=1,ne
	      chi_re(ie,j)=chi_re(ie,j)+1.d0
	   enddo
         endif
	enddo
        endif

 910  continue
c     output static value(s)
      write(6,1380)
      if(ioption.eq.1) then
         write(6,1400) (chi_zero(i),i=1,itot)
         write(6,1401) (chi_re(1,i),i=1,itot)
	else if (ioption.eq.2.or.ioption.eq.3) then
         write(6,1409)
	   write(6,1410) (chi_zero(i),i=1,itot)
         write(6,1411) (chi_re(1,i),i=1,itot)
	endif

      if(ioption.eq.3) goto 920

c	calculate the module of chi
      do j=1, itot
	   do ie=1,ne
	      aa=chi_im(ie,j)
	      bb=chi_re(ie,j)
	      chi_abs(ie,j)=sqrt(aa*aa+bb*bb)
	   enddo
	enddo

      if(ioption.eq.1) then

c gyg calculate optical conductivity
      do j=1, itot
         del = 0.d0
         if(icom(j).le.3)  del = 1.d0
	   do ie=1,ne
	      optc(ie,j)=chi_im(ie,j)*en(ie)/(Hartree*4*Pi)
              optcim(ie,j)=(del-chi_re(ie,j)*en(ie)/Hartree)/(4*Pi)
	   enddo
	enddo

c gyg calculate electron energy loss spectrum
      do j=1, 3   
           do ie=1,ne
              eels(ie,j)=chi_im(ie,j)/(chi_re(ie,j)**2+chi_im(ie,j)**2)
           enddo
        enddo

c     calculate the linear refractive index
      do j=1, itot
	   do ie=1,ne
	      aa=chi_abs(ie,j)
	      bb=aa+chi_re(ie,j)
	      refra(ie,j)=sqrt(bb*0.5d0)
	   enddo
	enddo

c     calculate the linear absorption coefficients
c               (in units of 10^5 cm-1)
c      facabs=1.d4/(Hartree*CCau*Bohr) 
      facabs=EE/(Hbar*CC) 
      do j=1, itot
	   do ie=1,ne
	      aa=chi_im(ie,j)*en(ie)/refra(ie,j)
	      alpha(ie,j)=aa*facabs
	   enddo
	enddo

c     f sum rule check
      write(6,1450)
      do i=1,itot
	effn=0
      do ie=1,ne-1
         effn=effn+0.5*(optc(ie,i)+optc(ie+1,i))*de
      enddo
      effn=effn*2*Omega/Hartree/Pi
      write(6,1455) i, effn
      enddo
      
      endif
               
c     generate output files for plotting 
c     write titles first
      write(11,1481)
      write(12,1482)
      write(13,1483)

	if(ioption.eq.1) then
      write(14,1484)
      write(15,1487)
      write(16,1485)
c     write(17,1486)
      write(18,1488)
      write(19,1489)
c gyg convert optical conductivity in units 10^15 sec^-1
      do j=1,itot
      do ie=1,ne
         optc(ie,j)  =1.51925d0*27.2116d0*optc(ie,j)
         optcim(ie,j)=1.51925d0*27.2116d0*optcim(ie,j)
      enddo
      enddo
	endif

      do ie=1,ne
         write(11,1500) en(ie),(chi_im(ie,j),j=1,itot)
         write(12,1500) en(ie),(chi_re(ie,j),j=1,itot)
         write(13,1500) en(ie),(chi_abs(ie,j),j=1,itot)
         if(ioption.eq.1) then
         write(14,1500) en(ie),(optc(ie,j),j=1,itot)
         write(15,1500) en(ie),(optcim(ie,j),j=1,itot)
         write(16,1500) en(ie),(refra(ie,j),j=1,itot)
         write(17,1500) en(ie),(alpha(ie,j),j=1,itot)
         if (ie .ge. 2) write(18,1510) en(ie),jdos(ie),jdos(ie)/en(ie),
     +                         jdos(ie)/en(ie)**2
         write(19,1500) en(ie),(eels(ie,j),j=1,3)
	 endif
      enddo
c     if(ioption.eq.1) deallocate(jdos,jdos1,jdos2)
      if(ioption.eq.1) deallocate(jdos)

 920  continue

c      call timer(cp1)
c     cp1 = 0.d0
c      cp=cp1-cp0

c     output time usage
c      call output_time(cp, ihour, imin, isec)
c      write(6,2001) ihour, imin, isec

c     release memory
      DEALLOCATE (en, chi_im, chi_re, chi_abs)
      DEALLOCATE (ev, mme)

      write(6,2002)   
      call finish_optic(ioption)

1001  format('********** PROGRAM OPTIC **********')
1005  format(/, i2, ': Calculate dielectric constant')
1006  format(/, i2, ': Calculate SHG coefficients')
1050  format(/, ' Only calculate static value(s)')
1100  format(/, ' Emin=',f10.3,', erange=',f10.3, /,
     &          '   de=',f10.3,', sigma =',f10.3)
1105  format(' ne =', i6, ' iso =', i6, ' itr =', i6)
1110  format(/,' isci=',i3, ', eshift =', f10.3,/)
1111  format(i3, ' components to be calculated:', 6i4,/)
1112  format(i3, ' components to be calculated:',/ 
     &       3(2x,6(i1,',',i1,3x),/))
1114  format(' nbmin=', i6,' nbmax=',i6)
1115  format(/,' *** Warning: number of total bands changed from',
     &       i6,' to',i6,' ***') 
1120  format(/,' Lattice vectors:', /,3f12.6, /, 
     &         3f12.6,/,3f12.6,/,
     &       ' Reciprocal lattice vectors:',/,3f12.6,/,
     &         3f12.6,/,3f12.6,/)
1121  format(' Vol  =', f12.6, / 
     &       ' Number of bands =', i10, /,
     &       ' A metal: Fermi energy (Ef) =',f10.6,/,
     &       ' Number of symmetric operations =', i8,/)
1122  format(' Vol  =', f12.6, /
     &       ' Number of bands =', i10, /,
     &       ' A non-metal: Number of valence bands =', i8, /,
     &       ' Number of symmetric operations =', i8,/)
1200  format(' Operation', i3,/, ' R(trx):', 3f10.4,/,8x,3f10.4,
     &         /,8x,3f10.4)
1201  format(' R(sym):', 3f10.4,/,8x,3f10.4,/,8x,3f10.4)
1300  format(/,' k-point',i6,' : (',f8.3,',',f8.3,',',f8.3,')',2x,f10.6)
1310  format(1x,5f12.4)
1350  format(/,' Total effective k-points:',i6)
1380  format(/,' Static value(s):') 
1400  format(/,' Epsl   =', 6f10.4)
1401  format(' Epsl_kk=', 6f10.4)
1409  format(/,' In units of 10^-6 esu:') 
1410  format(/,' Chi_direct=', 3(6f10.6,/,12x))
1411  format(' Chi_kk    =', 3(6f10.6,/,12x))
1450  format(/,' Check f sum rule:')
1455  format(i4, ' Neff=', f10.5)
1481  format('### imaginary part of chi or epsilon')
1482  format('### real part of chi or epsilon')
1483  format('### absolute value of chi or epsilon')
1484  format('### real part of conductivity (in units of 10^15 sec^-1)')
1487  format('### imag part of conductivity (in units of 10^15 sec^-1)')
1485  format('### linear refractive index')
1486  format('### linear absorption coefficients',
     &       '(in units of 10^5 cm-1)')
1488  format('### joint density of states')
1489  format('### electron energy loss spectra')
1500  format(f8.3, 18f10.4)
1510  format(f8.3, 6(1pe14.5))
2001  format(/' Total time used:', I10,'H ',I2,'M ',I2,'S')
2002  format('*********  END OF OPTIC  **********')
      end
!-----------------------------------------------------------------------
      subroutine start_optic(ioption)
      character* 20 casename
      character* 30 filename 
      open(1, file='opticpack.def', status='old', err=100) 
      read(1,*,err=101) casename
      nlen=lenstr(casename)
      close (1)
      filename=casename(1:nlen)//'.strucin'
      open(2, file=filename, status='old', err =104)
      filename=casename(1:nlen)//'.opticin'
      open(5, file=filename, status='old', err =102)
      read(5,*) ioption
c     ioption:  1     dielectric constants
c               2     SHG coefficients
c               3     SHG (method 2, satisfy the Kleinman Symmetry)
c               4     LEO coefficients

      filename=casename(1:nlen)//'.opticout'
      open(6, file=filename)
      filename=casename(1:nlen)//'.mme'
      open(7, file=filename, status='old', form='formatted', err=103)
      if(ioption.eq.1) then
         filename=casename(1:nlen)//'.epsim'
         open(11,file=filename)
         filename=casename(1:nlen)//'.epsre'
         open(12,file=filename)
         filename=casename(1:nlen)//'.epsabs'
         open(13,file=filename)
      else if(ioption.eq.2) then
	 filename=casename(1:nlen)//'.chiim'
         open(11,file=filename)
         filename=casename(1:nlen)//'.chire'
         open(12,file=filename)
         filename=casename(1:nlen)//'.chiabs'
         open(13,file=filename)
      endif
      
	if(ioption.eq.1) then
      filename=casename(1:nlen)//'.sigre'
      open(14,file=filename)
      filename=casename(1:nlen)//'.sigim'
      open(15,file=filename)
      filename=casename(1:nlen)//'.refra'
      open(16,file=filename)
      filename=casename(1:nlen)//'.alpha1'
      open(17,file=filename)
      filename=casename(1:nlen)//'.jdos'
      open(18,file=filename)
      filename=casename(1:nlen)//'.eels'
      open(19,file=filename)
      endif

      return
 100  stop 'error in opening file :opticpack.def'
 101  stop 'error in reading file :opticpack.def'
 102  stop 'error in opening file :*.opticin'
 103  stop 'error in opening file :*.mme'
 104  stop 'error in opening file :*.strucin'
      end
!-----------------------------------------------------------------------
      subroutine finish_optic(ioption)
      close (5)
      close (6)
      close (7)
      close (11)
      close (12)
      close (13)
      
      if(ioption.eq.1) then
      close (14)
      close (15)
      close (16)
      endif

      return
      end
!-----------------------------------------------------------------------
      function lenstr (string)
      integer lenstr, i
      character string*(*)
      lenstr = 0
      do i = 1, len (string)
         if (string (i:i) .ne. ' ') lenstr = i
      enddo
      end
!-----------------------------------------------------------------------
      subroutine transgroup1(nsym,sym,trx,g,rg)
      implicit real*8(a-h,o-z)
      dimension sym(3,3,48),trx(3,3,48)
      dimension g(3,3), rg(3,3)
      do k=1,nsym
         do j=1,3
            do i=1,3
               trx(i,j,k)=0.d0
               do l=1,3
                  do m=1,3
                  trx(i,j,k)=trx(i,j,k)+rg(m,i)*
     *                       sym(m,l,k)*g(l,j)
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine transgroup(nsym,sym,trx,g,rg)
      implicit real*8(a-h,o-z)
      dimension sym(3,3,48),trx(3,3,48)
      dimension g(3,3), rg(3,3) 
      do k=1,nsym
         do j=1,3
            do i=1,3
cgyg           trx(i,j,k)=0.d0
               sym(i,j,k)=0.d0
               do l=1,3
                  do m=1,3
                  sym(i,j,k)=sym(i,j,k)+g(i,m)*trx(m,l,k)*rg(j,l)
c                 sym(i,j,k)=sym(i,j,k)+rg(i,m)*trx(m,l,k)*g(j,l)
cgyg              trx(i,j,k)=trx(i,j,k)+g(i,m)*
cgyg *                       sym(m,l,k)*rg(j,l)
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine subgroup(rk,nsym,trx,nsubt)
      implicit real*8(a-h,o-z)
      dimension trx(3,3,48),rktran(3),rk(3)

      nsubt=0
      do 10 j=1,nsym
         do k=1,3
            rktran(k)=0.d0
            do l=1,3
               rktran(k)=rktran(k)+trx(k,l,j)*rk(l)
            enddo
         enddo
         isame=0
         do k=1,3
            dif=dabs(rktran(k)-rk(k))
            idif=int(dif)
            dif=dif-real(idif)
            if(dif.lt.1.d-5) isame=isame+1
         enddo
         if(isame.eq.3) then
            nsubt=nsubt+1
         endif
10    continue 

	if(nsubt.eq.0) stop 'Something wrong with the basis or group'

      return
      end
!-----------------------------------------------------------------------
      subroutine symptensor(ptensor,trx,nt,nsubt)
      implicit real*8(a-h,o-z)
      dimension ptensor(3,3),trx(3,3,48),ptmp(3,3),pnew(3,3),
     &          pfin(3,3)
      
      call inidouble(pfin,9)
      do mtrans=1,nt
         do j=1,3
            do i=1,3
	         pnew(i,j)=0.d0
               do ii=1,3
                   pnew(i,j)=pnew(i,j)+
     &             ptensor(i,ii)*trx(j,ii,mtrans)
               enddo
            enddo
         enddo
      
   	   do j=1,3
            do i=1,3
               ptmp(i,j)=0.d0
               do ii=1,3
                   ptmp(i,j)=ptmp(i,j)+
     &             pnew(ii,j)*trx(i,ii,mtrans)
               enddo
           enddo
         enddo

         do j=1,3
	      do i=1,3
	         pfin(i,j)=pfin(i,j)+ptmp(i,j)
	      enddo
	   enddo

      enddo

      do j=1,3
         do i=1,3
            ptensor(i,j)=pfin(i,j)/nsubt
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine symptensor2(ptensor,trx,nt,nsubt)
      implicit real*8(a-h,o-z)
      dimension ptensor(3,3,3),trx(3,3,48),ptmp(3,3,3),
     &          pnew(3,3,3),pfin(3,3,3)
      
	call inidouble(pfin,27)
      do mtrans=1,nt

        call inidouble(pnew,27)
         do k=1,3
            do j=1,3
               do i=1,3
                  do kk=1,3
                     pnew(i,j,k)=pnew(i,j,k)+
     &               ptensor(i,j,kk)*trx(k,kk,mtrans)
                  enddo
               enddo
	      enddo
	   enddo

        call inidouble(ptmp,27)
         do k=1,3
            do j=1,3
               do i=1,3
                  do jj=1,3
                     ptmp(i,j,k)=ptmp(i,j,k)+
     &               pnew(i,jj,k)*trx(j,jj,mtrans)
                  enddo
               enddo
	      enddo
	   enddo

        call inidouble(pnew,27)
         do k=1,3
            do j=1,3
               do i=1,3
                  do ii=1,3
                     pnew(i,j,k)=pnew(i,j,k)+
     &               ptmp(ii,j,k)*trx(i,ii,mtrans)
                  enddo
               enddo
	      enddo
	   enddo

         do k=1,3 
            do j=1,3
	         do i=1,3
	             pfin(i,j,k)=pfin(i,j,k)+pnew(i,j,k)
	         enddo
            enddo
	   enddo
    
      enddo
    
      do i=1,3
         do j=1,3
            do k=1,3
               ptensor(i,j,k)=pfin(i,j,k)/nsubt
            enddo
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine kk_trans(en,ne,de,chi_im,chi_re,itot)

c     include 'param.h'     
      implicit real*8 (a-h,o-z)
c     parameter (ND=7)
c     parameter (NBD=700)
c     Define Constants
      parameter (Pi=3.14159265359d0)
c     parameter (Hartree=27.2116d0)
c     parameter (Hbar=1.0545887d0)
c     parameter (CC=2.99793d0)
c     parameter (CCau=1.37037d2)
c     parameter (EE=1.6021892d0)
c     parameter (Bohr=5.2917706d0)
c     parameter (ESU=5.830345d0)

      dimension en(ne),chi_im(ne,itot),chi_re(ne,itot)

      do ic=1,itot
         do l=1,ne
            dsum=0.d0   
            do i=1,ne-1
               if(abs(en(i)-en(l)).lt.0.0001) goto 500
               if(abs(en(i+1)-en(l)).lt.0.0001) goto 500
               a=chi_im(i,ic)
               b=chi_im(i+1,ic)
               dsum=dsum+0.5*(a*en(i)/(en(i)*en(i)-en(l)*en(l))
     &              +b*en(i+1)/(en(i+1)*en(i+1)-en(l)*en(l)))
500            continue
            enddo
            chi_re(l,ic)=(2/Pi)*de*dsum
         enddo
      enddo
      
      return
      end
!-----------------------------------------------------------------------
      subroutine findcom(pcom,icom,ptensor,itot)
      implicit real*8(a-h,o-z)
      dimension pcom(itot),icom(itot),ptensor(9)
      
      do i=1,itot
         if(icom(i).eq.1) pcom(i)=ptensor(1)
         if(icom(i).eq.2) pcom(i)=ptensor(5)
         if(icom(i).eq.3) pcom(i)=ptensor(9)
         if(icom(i).eq.4) pcom(i)=ptensor(8)
         if(icom(i).eq.5) pcom(i)=ptensor(3)
         if(icom(i).eq.6) pcom(i)=ptensor(4)
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine findcom2(pcom,icom,jcom,ptensor,itot)
      implicit real*8(a-h,o-z)
      dimension pcom(itot),icom(itot),jcom(itot),ptensor(3,3,3)
      
      do i=1,itot
	 if(jcom(i).eq.1) pcom(i)=ptensor(icom(i),1,1)
         if(jcom(i).eq.2) pcom(i)=ptensor(icom(i),2,2)
         if(jcom(i).eq.3) pcom(i)=ptensor(icom(i),3,3)
         if(jcom(i).eq.4) pcom(i)=ptensor(icom(i),2,3)
         if(jcom(i).eq.5) pcom(i)=ptensor(icom(i),3,1)
         if(jcom(i).eq.6) pcom(i)=ptensor(icom(i),1,2)
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine getnb(nbmin,nbmax,nunit,nk,nspin)
      implicit real*8(a-h,o-z)
      complex*16 mme(3)
      dimension ev(4000)
      dimension nbmin(2),nbmax(2)
c	dimension g(3,3),rg(3,3)
c       dimension imat(3,3,48)
      do is = 1,nspin
      nbmin(is)=1000000
      nbmax(is)=0
      enddo
      read(nunit,'(I6)') nk
      write(80,'(i5)') nk
cgyg  read(nunit) (((sym,i=1,3),j=1,3),k=1,nsym)
      do 360 ik=1,nk
      read(nunit,'(4F14.8)') qk1,qk2,qk3,wfk
      write(80,'(i4,3f10.6,2x,f10.6)') ik,qk1,qk2,qk3,wfk
      do is = 1,nspin
      read(nunit,'(I6)') nb
      read(nunit,'(5E13.5)') (ev(j),j=1,nb)
      write(80,'(i5)') nb
      write(80,'(5f10.5)') (ev(j),j=1,nb)
      do j=1,nb
         do i=j,nb
            read(nunit,'(6E13.5)') (mme(k),k=1,3)
c           write(80,'(2i4,3(1x,2e11.4)') j,i,(mme(k),k=1,3)
         enddo
      enddo
      if(nbmin(is).gt.nb) nbmin(is)=nb
      if(nbmax(is).lt.nb) nbmax(is)=nb
      enddo
360   continue

      rewind nunit 
cgyg  read(nunit) nk, omega, ival, nsym
cgyg  read(nunit) g, rg
c     read(nunit) nk, omega, nsym
      return
      end
!-----------------------------------------------------------------------
	subroutine inidouble(a,n)
	real*8 a(n)
	do i=1,n
	   a(i)=0.d0
	enddo
	end
!-----------------------------------------------------------------------
	subroutine iniint(a,n)
	integer a(n)
	do i=1,n
	   a(i)=0
	enddo
	end
!-----------------------------------------------------------------------
	subroutine inicomplex(a,n)
	double complex a(n)
	do i=1,n
	   a(i)=(0.d0,0.d0)
	enddo
	end
      SUBROUTINE RECIP(RBAS,GBAS,VC)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EPS(3,3,3),RBAS(3,3),GBAS(3,3),V(3),BET(3),GBET(3)
      DATA EPS/27*0.D0/
      EPS(1,2,3)=1.D0
      EPS(2,3,1)=1.D0
      EPS(3,1,2)=1.D0
      EPS(1,3,2)=-1.D0
      EPS(3,2,1)=-1.D0
      EPS(2,1,3)=-1.D0
      DO  3  I=1,3
      I1=I+1
      IF(I1.GT.3) I1=I1-3
      I2=I+2
      IF(I2.GT.3) I2=I2-3
      V(I)=0.D0
      DO 4  J=1,3
      SIM=0.D0
      DO  5  K=1,3
      DO  5  L=1,3
    5 SIM=SIM+EPS(J,K,L)*RBAS(K,I1)*RBAS(L,I2)
      GBAS(J,I)=SIM
    4 V(I)=V(I)+SIM*RBAS(J,I)
      BET(I)=0.D0
      GBET(I)=0.D0
      DO  6 J=1,3
      GBAS(J,I)=GBAS(J,I)/V(I)
      BET(I)=BET(I)+RBAS(J,I)**2
    6 GBET(I)=GBET(I)+GBAS(J,I)**2
      BET(I)=DSQRT(BET(I))
      GBET(I)=DSQRT(GBET(I))
    3 CONTINUE
      WRITE(6,1061) RBAS
 1061 FORMAT(//,'  BASIS OF DIRECT LATTICE (A.U.)',//,
     *  3(3F10.6//))
 1060 FORMAT(//,'  BASIS OF RECIPROCAL LATTICE (WITHOUT 2*PI)',//,
     *  3(3F10.6//))
      VC=DABS(V(1))
      WRITE(6,1070) VC
 1070 FORMAT(' WIGNER-SEITZ CELL VOLUME : ',F15.8)
      WRITE(6,1060) GBAS
      RETURN
      END
c      subroutine timer(time)
c      real*8 time                                                  
c      time=second()
c      end

c      function second()
c      real t(2)
c      second=etime(t)
c      end

***********************************************************************
      subroutine output_time(cp, i, j, k)
      real*8 cp
c     hours
      i=cp/3600.0
c     minutes
      j=(cp-i*3600.0)/60.0
c     seconds
      k=cp-i*3600.0-60.0*j
      return
	end
