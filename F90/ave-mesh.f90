      implicit real*8 (a-h,o-z)
!     dimension dump(20),num(2),x(100),y(100),z(100),chg(100,100,800)
      dimension dump(20),num(100)

      real(kind=8), allocatable :: x(:), y(:), z(:)
      real(kind=8), allocatable :: chg(:,:,:)
      real(kind=8), allocatable :: chgt(:)
      dimension a(3,3)

      write(*,*) 'Choose input file (LOCPOT=1, CHGCAR=2):'
      read(*,*) input
      if (input.eq.1) then
      open(15,file='LOCPOT')
      write(6,*) ' Enter fermi energy: '
      read(*,*) ef
      elseif (input.eq.2) then
      open(15,file='CHGCAR')
      ef=0.0
      elseif ((input.ne.1).and.(input.ne.2)) then
      write(*,*) ' INPUT ERROR, input must equal to 1 or 2 '
      stop
      endif


      write(*,*) 'Spin polarized calculation? (no=1, yes=2):'
      read (*,*) ispin
      if (ispin.eq.1) then
       if (input.eq.1) open(16,file='workfn.dat')
       if (input.eq.2) open(16,file='chgave.dat')
      elseif (ispin.eq.2) then
       if (input.eq.1) then
       open(16,file='workfn-up.dat')
       open(26,file='workfn-dn.dat')
       elseif (input.eq.2) then
       open(16,file='chgave-tot.dat')
       open(26,file='chgave-spn.dat')
       endif
      elseif ((ispin.ne.1).and.(ispin.ne.2)) then
      write(*,*) ' INPUT ERROR, ispin must equal to 1 or 2 '
      stop
      endif


      open(19,file='atom.dat')

       write(6,*) ' Enter your VASP version 4.x|5.x (0|1) '
       read(5,*) iver
       write(6,*) ' Enter the number of atomic type (nty) '
       read(5,*) nty


      read(15,1) dump
      read(15,*) scale
      do i=1,3
      read(15,*) (a(i,j),j=1,3)
      enddo
 1    format(20a4)
      aa=sqrt(a(3,1)**2+a(3,2)**2+a(3,3)**2)
      if (iver .eq. 1) then
      read(15,1) dump
!     write(6,1) dump
      endif
      read(15,*) (num(i),i=1,nty)
!     write(6,*) (num(i),i=1,nty)
 2    format(10i4)
      natm=0
      do i=1,nty
      natm=natm+num(i)
      enddo
      read(15,1) dump
!     write(6,1) dump
      zero=0.0

      write(6,*) ' natm = ',natm
! Huang    
      allocate(x(1:natm))
      allocate(y(1:natm))
      allocate(z(1:natm))
! Huang

      do i=1,natm
      read(15,*) x(i),y(i),z(i)
      enddo
      do i=1,natm-1
      do j=i+1,natm
      if (z(j) .lt. z(i))  then
      zz=z(j)
      z(j)=z(i)
      z(i)=zz    
      endif
      enddo
 30   format(2f10.5)
      enddo
      do i = 1,natom 
      write(19,30) z(i)*aa*scale,zero
      enddo
      nat=natm/2+1
      do i=1,nat
!     write(6,4) i,z(i)
   4  format(i5,f10.5)
      enddo
      read(15,1) dump

      do 600 is = 1,ispin
      read(15,*) nx,ny,nz
!     write(6,*) 'is, nx,ny,nz ',is,nx,ny,nz
      nn=nx*ny
!Huang 
      if (is .eq. 1) allocate(chg(1:nx,1:ny,1:nz))
!Huang
      read(15,*) (((chg(j,k,iz),j=1,nx),k=1,ny),iz=1,nz)
!     write(6,5) (((chg(j,k,iz),j=1,nx),k=1,ny),iz=1,nz)
   5  format(5(e18.10,1x))
!Huang
      if (is .eq. 1) allocate(chgt(1:nz))
!Huang
      chgscal = 1.0
      if ( input .eq. 2 ) chgscal = aa*scale
      do i=1,nz
      chgt(i)=0.0
      do j=1,nx
      do k=1,ny
      chgt(i)=chgt(i)+chg(j,k,i)
      enddo
      enddo
      x1=float(i-1)/float(nz)
      chgt(i)=chgt(i)/float(nn)-ef
      write(6+10*is,10) x1*aa*scale,chgt(i)/chgscal               
 10   format(2f12.5)
       enddo

      if (input.eq.1)then
      emax=-9999.0
      do i=1,nz
      if (chgt(i) .gt. emax) emax=chgt(i)
      enddo
      workfn=emax
!     write(6,201) ef
      write(6,301) workfn
      endif
!201  format(' fermi energy =',f10.5)
 301  format(' workfunction =',f10.5)

      if ((is.eq.1).and.(ispin.eq.2))then
       if (input.eq.1) then
       read(15,*) (ttt,i=1,natm)

       elseif (input.eq.2) then
       do nn=1,natm
       read(15,'(24x,2i4)') nn1,nn2
       read(15,*) (ttt,i=1,nn2)
       enddo
       read(15,*) (ttt,i=1,natm)
       endif
      endif
 600  continue
      stop  
      end
