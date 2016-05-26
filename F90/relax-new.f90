! The program is the Fortran 90 version of relax-new.f
! Date: 2012/4/19
! Author: T.C. Leung           
! Function:
! Find the increase or decrease percent between layers of the relax slab and those of the bulk. 

! f90
      real(kind=8), allocatable :: z(:)
      integer, allocatable :: nat(:)
!     dimension z(500),nat(5)
! f90

      dimension dump(20)
      write(6,*) ' Enter your VASP version 4.x|5.x (0|1) '
      read(5,*) ver
      open(4,file='CONTCAR')
      open(7,file='relax.out')
      read(4,1) dump
 1    format(20a4)
      read(4,*) aa
      read(4,1) dump
      read(4,1) dump
      read(4,*) a1,a2,a3
      cc = aa* sqrt( a1**2 +a2**2 + a3**2 )
      if (ver .eq. 1) read(4,1) dump
      write(6,*) ' Enter no of type '
      read(5,*) ntyp

! f90   
      allocate(nat(1:ntyp))
! f90

      read(4,*) (nat(i),i=1,ntyp)
      natom = 0
      do i=1,ntyp
      natom = natom + nat(i)
      enddo
      write(6,*) ' Selective dynamics yes|no (1|0) ? '
      read(5,*) isel
      if (isel .eq. 1) read(4,1) dump
      read(4,1) dump
      write(6,*) ' Enter del0 '
      read(5,*) del0

! f90
      allocate(z(1:natom))
! f90

      do i=1,natom
      read(4,*) x,y,z(i)
      if (z(i) .gt. 0.5) z(i) = z(i) - 1.0
      enddo

      do i = 1,natom
      do j = i+1,natom
      if (z(j) .gt. z(i) ) then
      temp = z(i)
      z(i) = z(j)
      z(j) = temp
      endif
      enddo
      enddo

      do i =1,natom
!c     write(6,*) i,z(i)
      enddo
      do i = 1,natom-1
      dels= (z(i) - z(i+1))*cc
      if (dels .gt. 5.0) go to 100
      enddo
 100  nn = i+1
!c     write(6,*) nn
      do i = nn, natom
      z(i) = z(i) + 1.0
      enddo
      do i = 1,natom
      do j = i+1,natom
      if (z(j) .gt. z(i) ) then
      temp = z(i)
      z(i) = z(j)
      z(j) = temp
      endif
      enddo
      enddo
      do i =1,natom
!c     write(6,*) i,z(i)
      enddo

!c     write(6,*) '    del    del(angstrom)  del(%)   '
      write(7,*) '    del    del(angstrom)  del(%)   '
      do i = 1,natom-1
      del = z(i) - z(i+1)
      relax =( (del-del0)/del0 ) *100.0
!c     write(6,2) del,del*cc,relax
      write(7,2) del,del*cc,relax
 2    format(3(f10.5,2x))
      enddo
      stop 
      end
