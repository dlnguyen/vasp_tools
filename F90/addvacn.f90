!
! f90 version          
! For slab calculation, change the thickness of the vaccum in the slab 
!     to the value we wanted.
!
! modified by T.C. Leung  May 2016
!
      implicit real*8 (a-h,o-z)
      character(len=80) :: dump, dump1,dump2 
      real(kind=8), allocatable :: x(:), y(:), z(:), zz(:), dz(:) 
      real(kind=8) :: a(3,3)
      integer, allocatable :: nat(:), natn(:) 
      character(len=8), allocatable :: dum(:,:)

      open(8,file='poscar')
      write(6,*) ' Enter your input file POSCAR|CONTCAR (0|1) '
      read(5,*) inp
      
      if (inp .eq. 0) open(7,file='POSCAR')
      if (inp .eq. 1) open(7,file='CONTCAR')
       write(6,*) ' Enter your VASP version 4.x|5.x (0|1) '
       read(5,*) iver
!c
!c     vac is the thickness of the vaccum we wanted
!c
!c     vac = 10.0
      write(6,*) ' Enter the number of atomic type (nty) '
      read(5,*) nty
      write(6,*) ' Enter the thickness of the vaccum     '
      read(5,*) vac
      write(6,*) ' Do you use the selective dynamics yes|no (1|0)? '
      read(5,*) sd        
      read(7,1) dump
      write(8,1) dump
  1   format(a80)
      read(7,*) aa
      write(8,2) aa
  2   format(f20.17)
      do i=1,3
      read(7,*) (a(i,j),j=1,3)
!c     write(8,3) (a(i,j),j=1,3)
  3   format(1x,3f22.16)
      enddo
      if ( iver .eq. 1) then
      read(7,1) dump2
      write(6,1) dump2
!c     write(8,1) dump
      endif

! f90
      allocate(nat(1:nty))
! f90

      read(7,*) (nat(i),i=1,nty)
  4   format(20i4)
      natom=0
      do i=1,nty
      natom=natom+nat(i)
      enddo  
      read(7,1) dump
      if (sd .eq. 1) read(7,1) dump1

! f90
      allocate(x(1:natom))
      allocate(y(1:natom))
      allocate(z(1:natom))
      allocate(zz(1:natom))
      allocate(dz(1:natom))
      allocate(dum(1:natom,1:3))
! f90

      do i=1,natom
      if (sd .eq. 0) read(7,*) x(i),y(i),z(i)
      if (sd .eq. 1) read(7,*) x(i),y(i),z(i),(dum(i,j),j=1,3)
  5   format(3f20.16,3x,3a4)
      zz(i) = z(i)                        
      if (zz(i) .lt. 0.0) zz(i)=zz(i)+1.0
      enddo
      do i = 1, natom-1
        do j = i+1,natom
       if ( zz(j) .lt. zz(i) ) then
       zdump = zz(i)
       zz(i) = zz(j)
       zz(j) = zdump
       endif
         enddo
       enddo
       dzmax = 0.0
       do i = 1, natom
       dz(i) = zz(i+1)-zz(i)
       if ( i .eq. natom ) dz(i) = 1.0 + zz(1) - zz(i) 
       if (dz(i) .gt. dzmax) then
       dzmax = dz(i)
       ii = i
       endif
       write(6,6) i,zz(i),dz(i),ii,dzmax
  6    format(i5,2f10.5,i5,f10.5)
       enddo
       vacz = ( zz(ii) + zz(ii+1) )/ 2.0
       write(6,7) vacz
  7    format(' vacz =',f10.5)
!c     find the thickness of the vaccum
      del = dzmax
      c=sqrt(a(3,1)**2 +a(3,2)**2 + a(3,3)**2)
      d0=aa*c      
      vac0=del*d0
      d=d0+(vac-vac0)
!c     we want the vaccum to have 10 A, therefore we scale the a(3,3) by
!c     the zscale
      do i = 1,natom
      if ( z(i) .gt. vacz ) z(i) = z(i) - 1.0
      write(6,6) i,z(i)
      enddo
      zscale = d/d0
      do i=1,3
      a(3,i)=a(3,i)*zscale
      enddo
      do i=1,3
      write(8,3) (a(i,j),j=1,3)
      enddo
      if ( iver .eq. 1) then
      write(6,1) dump2
      write(8,1) dump2
      endif
      write(8,4) (nat(i),i=1,nty)
      write(8,1) dump
      if (sd .eq. 1) write(8,1) dump1
      
      do i=1,natom 
      if (sd .eq. 0) write(8,5) x(i),y(i),z(i)/zscale
      if (sd .eq. 1) write(8,5) x(i),y(i),z(i)/zscale,(dum(i,j),j=1,3)
      enddo
      stop
      end
