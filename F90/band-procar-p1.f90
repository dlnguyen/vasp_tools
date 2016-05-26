      implicit real*8(a-h,o-z)
      parameter (nbd = 200)
      parameter (nkd = 500)
      parameter (nxd = 300)
      parameter (natmd = 150)

! F90
      real(kind=8) :: a(3,3), b(3,3), c(3)
      real(kind=8), allocatable :: e(:,:), sk(:,:)
      real(kind=8), allocatable :: xx(:), wei(:)
      real(kind=8), allocatable :: oc(:,:,:,:)
      dimension dump(20)
!     dimension a(3,3),b(3,3),c(3),e(nkd,nbd),sk(nkd,3)
!     dimension xx(nxd) ,wei(nkd)
!     dimension dump(20),oc(nkd,nbd,natmd,4)
! F90

!c     open(5,file='band.inp')
      open(7,file='PROCAR',form='FORMATTED',status='OLD') 
      pi = 3.141592654
      read(7,103) dump                    

      write(*,*) 'Spin polarized calculation? (no=1,yes=2):'
      read (5,*) ispin
      if ((ispin.ne.1).and.(ispin.ne.2)) then
      write(*,*) ' INPUT ERROR, ispin must equal to 1 or 2 '
      stop
      endif

!c     write(*,*) 'Enter # of interval (npoints) and division (ndiv):'
!c     read (*,*) npoints,ndiv
      open(9,file='KPOINTS',form='FORMATTED',status='OLD') 
      read(9,100) temp
      read(9,*) ndiv

      write(*,*) 'Enter the range of energy to plot:'
      read (5,*) er1,er2
      emin=min(er1,er2)
      emax=max(er1,er2)
      write(*,*) 'Enter the value of fermi energy:'
      read(5,*)  ef
      write(6,*) ' which atom you want in your projection band '
      read(5,*) natomp
      write(6,*) ' which orbital you want s,p,d,total  1/2/3/4  '
      read(5,*) norb
      write(6,*) ' enter the scaling factor in your projection band '
      read(5,*) scaling

      if (ispin.eq.1) then
      open(11,file='band-p.dat')
      elseif (ispin.eq.2) then
      open(11,file='band-p-up.dat')
      open(12,file='band-p-dn.dat')
      endif

      open(8,file='POSCAR',form='FORMATTED',status='OLD')
      read(8,100) temp
!c      write(6,100) temp
 100   format(20a4)
       read (8,*) aa
!c      WRITE(6,*) aa
!c   
!c      *** read lattice constant from POSCAR**
!c            
      do i=1,3
         read (8,*) (a(i,j),j=1,3)
!c       WRITE(6,500) (a(i,j),j=1,3)
 500   format (3f12.8)
      enddo
      do i=1,3
         do j=1,3
       a(i,j)=aa*a(i,j) 
         enddo
!c       WRITE(6,500) (a(i,j),j=1,3)
      enddo
!c
!c     *** read lattice vector from POSCAR***
!c
      volume=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
&           +a(1,3)*a(2,1)*a(3,2)-a(1,1)*a(2,3)*a(3,2) &
&           -a(1,2)*a(2,1)*a(3,3)-a(1,3)*a(2,2)*a(3,1)
      do i=1,3
          if (i .eq. 1) then
            j=2
            k=3
          else if (i .eq. 2) then
            j=3
            k=1
          else
            j=1
            k=2
          endif 
        c(1)=a(j,2)*a(k,3)-a(j,3)*a(k,2)
        c(2)=a(j,3)*a(k,1)-a(j,1)*a(k,3) 
        c(3)=a(j,1)*a(k,2)-a(j,2)*a(k,1)
        do j=1,3
           b(i,j)=2*pi*c(j)/volume
!c         WRITE (6,*) b(i,j)
        enddo
       enddo


      do 9000 isp=1,ispin

        read(7,104) nk,nband,nion

! F90      
        allocate(sk(1:nk,1:3))
        allocate(wei(1:nk))
        allocate(e(1:nk,1:nband))
        allocate(oc(nk,nband,nion,4))
! F90 
        do 1000 k = 1,nk
          read(7,103) dump
!c        write(6,103) dump
          read(7,105) kp,(sk(k,j),j=1,3),wei(k)
!c        write(6,1105) kp,(sk(k,j),j=1,3),wei(k)
          read(7,103) dump
!c        write(6,103) dump
          do  nb = 1,nband
            read(7,106) nb1,e(k,nb),occ
!c          write(6,1106) nb1,e(k,nb),occ
            read(7,103) dump
!c          write(6,103) dump
            read(7,103) dump
!c          write(6,103) dump
!c          write(6,*) 'nion=',nion
            niont = nion +1
            if (nion .eq. 1) niont = 1
!c          write(6,*) 'niont=',niont
            do  ion = 1,niont
              read(7,107) (oc(k,nb,ion,j),j=1,4)
!c            write(6,107) (oc(k,nb,ion,j),j=1,4)
            enddo
            read(7,103) dump
!c          write(6,103) dump
          enddo
 1000   continue

        weight = 0.0
        do k = 1, nk
          weight = weight + wei(k)
        enddo

        do k = 1,nk
          wei(k) =  wei(k) / weight
        enddo

 101  format(10x,f9.5) 
 102  format(f10.5) 
 103  format(20a4)
 104  format(16x,i3,20x,i5,19x,i4)
 105  format(10x,i3,5x,3f11.8,13x,f11.8)
1105  format(' nkpoint ',10x,i3,5x,3f11.8,13x,f11.8)
 106  format(4x,i4,9x,f14.8,7x,f12.8)
1106  format(' nband ',4x,i4,9x,f14.8,7x,f12.8)
 107  format(3x,4f7.3)


!c
!c     *** find reciprocal lattice vector ***
! F90
      allocate(xx(nk))
! F90
      xx(1) = 0.0
      nn = 1
      do k = 1,nk-1
        dkx=(sk(k+1,1)-sk(k,1))*b(1,1) + (sk(k+1,2)-sk(k,2))*b(2,1) &
       &   + (sk(k+1,3)- sk(k,3))*b(3,1)
        dky=(sk(k+1,1)-sk(k,1))*b(1,2) + (sk(k+1,2)-sk(k,2))*b(2,2) &
       &   + (sk(k+1,3)- sk(k,3))*b(3,2)
        dkz=(sk(k+1,1)-sk(k,1))*b(1,3) + (sk(k+1,2)-sk(k,2))*b(2,3) &
       &   + (sk(k+1,3)- sk(k,3))*b(3,3)
        del =  sqrt ( dkx**2 + dky**2 + dkz**2 )
        nn = nn +1
        xx(nn) = xx(nn-1) + del
      enddo
      scale = xx(nn)
      do i = 1,nn
      xx(i) = xx(i)/scale
      enddo

      do n=1,nband
        if (mod(n,2).ne.0) then
         do k=1,nk
          ee = e(k,n) - ef
          if ( ee .gt. emax ) ee = emax
          if ( ee .lt. emin ) ee = emin
!c     write(6,107) (oc(k,nb,ion,j),j=1,4)
         write (10+isp,300) xx(k),ee,scaling*oc(k,n,natomp,norb)
         enddo
        elseif (mod(n,2).eq.0) then
         do i=nk,1,-1
          ee = e(i,n) - ef
          if ( ee .gt. emax ) ee = emax
          if ( ee .lt. emin ) ee = emin
          write (10+isp,300) xx(i),ee,scaling*oc(i,n,natomp,norb)
         enddo
        endif
      enddo
300   format (f12.8,2x,f12.8,f10.5)

      if (mod(nband,2) .ne. 0) then
          write (10+isp,300) xx(nk),emin,0.0
          write (10+isp,300) xx(1),emin,0.0
      else
          write (10+isp,300) xx(1),emin,0.0
      endif
!c
!c     *** write xx-ee ***
!c
        npoints=nk/ndiv
       do n=2,npoints
         kk=(n-1)*ndiv
        write (10+isp,300) xx(kk),emin,0.0
        write (10+isp,300) xx(kk),emax,0.0
        write (10+isp,300) xx(kk),emin,0.0
      enddo
        write (10+isp,300) xx(nk),emin,0.0
        write (10+isp,300) xx(nk),emax,0.0
        write (10+isp,300) xx(1),emax,0.0
        write (10+isp,300) xx(1),emin,0.0
        zero=0.0
        write (10+isp,300) xx(1),zero,0.0
        write (10+isp,300) xx(nk),zero,0.0

! F90      
        deallocate(sk)
        deallocate(wei)
        deallocate(e)
        deallocate(xx)
        deallocate(oc)
! F90 
 9000  continue

      stop
      end
