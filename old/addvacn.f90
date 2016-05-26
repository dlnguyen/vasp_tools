!c
!c We will change the thickness of the vaccum to the value we wanted.
!c
      implicit real*8 (a-h,o-z)

 
! f90
      character(len=80) :: dump, dump1 
      real(kind=8), allocatable :: x(:), y(:), z(:) 
      real(kind=8) :: a(3,3)
      integer, allocatable :: nat(:), natn(:) 
!      dimension dump(20),x(1000),y(1000),z(1000),a(3,3),nat(20),
!     &          natn(20),dump1(20)
! f90


!c    &          dum(1000,3),natn(20),dump1(20)

! f90
      character(len=8), allocatable :: dum(:,:)
!      character*8  dum(1000,3)
! f90

      open(8,file='poscar')
      write(6,*) ' Enter your input file POSCAR|CONTCAR (0|1) '
      read(5,*) inp
      
      if (inp .eq. 0) open(7,file='POSCAR')
      if (inp .eq. 1) then
       write(6,*) ' Enter your VASP version 4.x|5.x (0|1) '
       read(5,*) ver
       open(7,file='CONTCAR')
      end if
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
      if (inp .eq. 1 .and. ver .eq. 1) read(7,1) dump

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
      allocate(dum(1:natom,1:3))
! f90

      do i=1,natom
      if (sd .eq. 0) read(7,*) x(i),y(i),z(i)
      if (sd .eq. 1) read(7,*) x(i),y(i),z(i),(dum(i,j),j=1,3)
  5   format(3f20.16,3x,3a4)
      if (z(i) .gt. 0.5) z(i)=z(i)-1.0
      enddo
!c     find the max of z
      zmax=0.0
      zmin=0.0
      do i=1,natom
      if (z(i) .gt. zmax)  zmax=z(i)
      if (z(i) .lt. zmin)  zmin=z(i)
      enddo
!c     find the thickness of the vaccum
      del=(0.5-zmax)+(0.5+zmin)
      c=sqrt(a(3,1)**2 +a(3,2)**2 + a(3,3)**2)
      d0=aa*c      
      vac0=del*d0
      d=d0+(vac-vac0)
!c     we want the vaccum to have 10 A, therefore we scale the a(3,3) by
!c     the zscale
      zscale = d/d0
      do i=1,3
      a(3,i)=a(3,i)*zscale
      enddo
      do i=1,3
      write(8,3) (a(i,j),j=1,3)
      enddo
      write(8,4) (nat(i),i=1,nty)
      write(8,1) dump
      if (sd .eq. 1) write(8,1) dump1
      
      do i=1,natom 
      if (sd .eq. 0) write(8,5) x(i),y(i),z(i)/zscale
      if (sd .eq. 1) write(8,5) x(i),y(i),z(i)/zscale,(dum(i,j),j=1,3)
      enddo
      stop
      end
