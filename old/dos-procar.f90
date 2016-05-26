      implicit real*8(a-h,o-z)
      parameter (nkd =900)
      parameter (nbd = 84 )
      parameter (natmd = 20)
      parameter (ned = 1001)

! F90
      dimension dump(20)
      real(kind=8), allocatable :: oc(:,:,:,:), eig_ev(:,:), wt(:)
      real(kind=8), allocatable :: ee(:), gpdos(:), gpdost(:)
      real(kind=8), allocatable :: gpdos_s(:), gpdos_p(:), gpdos_d(:)

!     dimension dump(20),oc(nkd,nbd,natmd,4),eig_ev(nkd,nbd),wt(nkd)
!     dimension ee(ned),gpdos(ned),gpdost(ned)
!     dimension gpdos_s(ned),gpdos_p(ned),gpdos_d(ned)
! F90
      open(7,file='PROCAR') 
      pi = 3.141592654

      write(*,*) 'Spin polarized calculation? (no=1,yes=2):'
      read (*,*) ispin
      if ((ispin.ne.1).and.(ispin.ne.2)) then
      write(*,*) ' INPUT ERROR, ispin must equal to 1 or 2 '
      stop
      endif
      if (ispin.eq.1) then
      open(31,file='ldos.dat') 
      open(32,file='dos-tot.dat') 
      elseif (ispin.eq.2) then
      open(31,file='ldos-up.dat') 
      open(32,file='dos-tot-up.dat') 
      open(41,file='ldos-dn.dat') 
      open(42,file='dos-tot-dn.dat') 
      endif

      write(6,*) ' Enter fermi energy: '
      read(5,*) fermi
!c    write(6,102) fermi
      read(7,103) dump                    
!c    write(6,103) dump                    

      write(6,*) 'which atom (na) you want to plot LDOS:'
      read(5,*) na
      write(6,*) 'enter the gaussian smearing factor:'
      read(5,*) gaussian

      do 9000 isp=1,ispin
        read(7,104) nk,nband,nion
!c      write(6,104) nk,nband,nion
        if (nk .gt. nkd) stop ' nk too large '
        emin = 1000.0
        emax = -1000.0

! F90
        allocate(wt(nk))
        allocate(eig_ev(nk,nband))
        allocate(oc(nk,nband,nion+1,4))
! F90

        do 1000 k = 1,nk
          read(7,103) dump
          read(7,105) kp,pt1,pt2,pt3,wt(k)
!c        write(6,105) kp,pt1,pt2,pt3,wt(k)
          read(7,103) dump
          do  nb = 1,nband
            read(7,106) nb1,eig_ev(k,nb),occ
!c          write(6,106) nb1,eig_ev(k,nb),occ
            eig_ev(k,nb) = eig_ev(k,nb)-fermi 
            if (eig_ev(k,nb) .gt. emax) emax = eig_ev(k,nb)
            if (eig_ev(k,nb) .lt. emin) emin = eig_ev(k,nb)
            read(7,103) dump
            read(7,103) dump
!c          write(6,*) 'nion=',nion
            niont = nion +1
            if (nion .eq. 1) niont = 1
            do  ion = 1,niont
              read(7,107) (oc(k,nb,ion,j),j=1,4)
!c            write(6,107) (oc(k,nb,ion,j),j=1,4)
            enddo
            read(7,103) dump
!c          write(6,103) dump
          enddo
1000    continue
        weight = 0.0
        do k = 1, nk
          weight = weight + wt(k)
        enddo
        do k = 1,nk
          wt(k) =  wt(k) / weight
        enddo
!c      estart_ev = - 5.0
!c      eend_ev = 2.0
        estart_ev = int(emin -1.0)
        eend_ev = int(emax +1.0)
        netot =  701
        de_ev = (eend_ev - estart_ev)/(netot-1)

! F90
        allocate(ee(netot))
        allocate(gpdos(netot))
        allocate(gpdost(netot))
        allocate(gpdos_s(netot))
        allocate(gpdos_p(netot))
        allocate(gpdos_d(netot))
! F90

        do 2000 ne = 1, netot
          ee(ne) = estart_ev + (ne-1) * de_ev
          gpdos(ne) = 0.0d0
          gpdos_s(ne) = 0.0d0
          gpdos_p(ne) = 0.0d0
          gpdos_d(ne) = 0.0d0
          gpdost(ne) = 0.0d0
 2000   continue
        ascal = de_ev   / ( gaussian * sqrt(pi) )
        ascal = 1.0     / ( gaussian * sqrt(pi) )
        scal_spin = 2.0
        do 5000 k=1,nk
          do 4000 nb=1,nband
            ddos = scal_spin * oc(k,nb,na,4)*wt(k) 
            ddos_s = scal_spin * oc(k,nb,na,1)*wt(k) 
            ddos_p = scal_spin * oc(k,nb,na,2)*wt(k) 
            ddos_d = scal_spin * oc(k,nb,na,3)*wt(k) 
            ddost = scal_spin *wt(k) 
              do ne = 1,netot
                dij = ( eig_ev(k,nb) - ee(ne) )**2 / (gaussian**2)
                gpdos(ne) = gpdos(ne) + ascal * ddos * exp(-dij)
                gpdos_s(ne) = gpdos_s(ne) + ascal * ddos_s * exp(-dij)
                gpdos_p(ne) = gpdos_p(ne) + ascal * ddos_p * exp(-dij)
                gpdos_d(ne) = gpdos_d(ne) + ascal * ddos_d * exp(-dij)
                gpdost(ne) = gpdost(ne) + ascal * ddost * exp(-dij)
              enddo
 4000     continue
 5000   continue
        do 6000 ne = 1,netot
          if (gpdos(ne) .lt. 1.0E-29) gpdos(ne)=0.0
          if (gpdos_s(ne) .lt. 1.0E-29) gpdos_s(ne)=0.0
          if (gpdos_p(ne) .lt. 1.0E-29) gpdos_p(ne)=0.0
          if (gpdos_d(ne) .lt. 1.0E-29) gpdos_d(ne)=0.0
          if (gpdost(ne) .lt. 1.0E-29) gpdost(ne)=0.0
          write(21+10*isp,108) ee(ne),gpdos(ne),gpdos_s(ne),gpdos_p(ne) &
     &                        ,gpdos_d(ne)  
           write(22+10*isp,108) ee(ne),gpdost(ne)
 6000   continue

! F90
        deallocate(wt)
        deallocate(eig_ev)
        deallocate(oc)
        deallocate(ee)
        deallocate(gpdos)
        deallocate(gpdost)
        deallocate(gpdos_s)
        deallocate(gpdos_p)
        deallocate(gpdos_d)
! F90

 9000 continue
 101  format(10x,f9.5) 
 102  format(f10.5) 
 103  format(20a4)
 104  format(16x,i3,20x,i5,19x,i4)
 105  format(10x,i3,5x,3f11.8,13x,f11.8)
 106  format(4x,i4,9x,f14.8,7x,f12.8)
 107  format(3x,4f7.3)
 108  format(f10.5,4e12.4)
      stop
      end
