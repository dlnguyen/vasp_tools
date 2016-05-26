      implicit real*8(a-h,o-z)
      parameter (nat=200)
!      parameter (ntypd = 4)
      parameter (nshelld=10)
!      parameter (natmd=80)
      dimension num(5)
      dimension a(3,3)
!      dimension neigh(ntypd,natmd,ntypd,natmd)
!      dimension natom(3)
!      dimension r(ntypd,natmd,ntypd,natmd,nat)
      dimension rr(nshelld)
!      dimension no(ntypd,natmd,ntypd,natmd,nat)
      dimension mo(nshelld)
!      dimension itype(ntypd,natmd,ntypd,natmd,nat)
!      dimension t(ntypd,natmd,3)
!      character*4 fname(ntypd)
      integer(kind=4), allocatable :: natom(:)
      real(kind=8), allocatable :: t(:,:,:)
      real(kind=8), allocatable :: r(:,:,:,:,:)
      real(kind=8), allocatable :: no(:,:,:,:,:)
      real(kind=8), allocatable :: itype(:,:,:,:,:)
      real(kind=8), allocatable :: neigh(:,:,:,:)
      character(len=4), allocatable :: fname(:)
      
      open(22,file='neig.dat')

      write(6,*) ' Enter the input file POSCAR|CONTCAR (1|2) '
      read(5,*) ifile
      if (ifile .eq. 1) open(8,file='POSCAR')
      if (ifile .eq. 2) open(8,file='CONTCAR')
      write(6,*) ' Enter the type of atom '
      read(5,*) ntype
!      if (ntype .gt. ntypd) stop ' ntype too large '
!yhh
      allocate(fname(1:ntype))
!yhh
      write(6,*) ' Selective dynamics yes|no (1|0) '
      read(5,*) isel
      write(6,*) ' Enter the name of each atom (10a4) '
      read(5,1) (fname(j),j=1,ntype)
      write(6,*) ' Enter the Rcut (angstrom) '
      read(5,*) rcut
      write(6,*) ' Enter your VASP version 4.x|5.x (0|1) '
      read(5,*) ver
      read(8,600) temp
!     write(6,600) temp
 600  format(20a4)
      read(8,*) aa
!     write(6,*) aa
 
!      *** read lattice constant from POSCAR**
 
      do i=1,3
      read(8,*) (a(i,j),j=1,3)
!     write(6,*) (a(i,j),j=1,3)
 500  format (3f12.8)
      enddo
      do i=1,3
       do j=1,3
       a(i,j)=aa*a(i,j)
       enddo
      enddo
!yhh

      allocate(natom(1:ntype))
!yhh
      if (ifile .eq. 2 .and. ver .eq. 1) read(8,600) temp 
      read(8,*) (natom(i),i=1,ntype)
!yhh
      nmax=0
      do i=1,ntype
       if (natom(i) > nmax) nmax=natom(i)
      enddo

      allocate(t(1:ntype,1:nmax,1:3))
      allocate(r(1:ntype,1:nmax,1:ntype,1:nmax,1:nat))
      allocate(no(1:ntype,1:nmax,1:ntype,1:nmax,1:nat))
      allocate(itype(1:ntype,1:nmax,1:ntype,1:nmax,1:nat))
      allocate(neigh(1:ntype,1:nmax,1:ntype,1:nmax))
!yhh
!     read(8,1) (fname(j),j=1,ntype)
!     write(6,1) (fname(j),j=1,ntype)
 1    format(10a4)
      if (isel .eq. 1) read(8,600) temp
      do i = 1,ntype
!      if (natom(i) .gt. natmd) stop ' natom too large '
      enddo
!     write(6,*) ' natom(i) =',(natom(i),i=1,ntype)
      read(8,600) temp
      do i=1,ntype
       do j=1,natom(i)
       read(8,*) t(i,j,1),t(i,j,2),t(i,j,3)
       enddo
      enddo

!     *** read lattice vector from POSCAR***

      do 200 i = 1,ntype
       do 200 j = 1,natom(i)
      
        do 100 ip = 1,ntype
         do 100 jp = 1,natom(ip)
         nn = 0
          do 50 n1=-5,5
           do 50 n2=-5,5
            do 50 n3=-5,5
            dtij1 = t(i,j,1) - t(ip,jp,1)
            dtij2 = t(i,j,2) - t(ip,jp,2)
            dtij3 = t(i,j,3) - t(ip,jp,3)
            xp = (n1+dtij1)*a(1,1) + (n2+dtij2)*a(2,1) +&
                                       (n3+dtij3)*a(3,1)
            yp = (n1+dtij1)*a(1,2) + (n2+dtij2)*a(2,2) +&
                                       (n3+dtij3)*a(3,2)
            zp = (n1+dtij1)*a(1,3) + (n2+dtij2)*a(2,3) +&
                                       (n3+dtij3)*a(3,3)
!           write(6,*) ' test2 ',xp,yp,zp
            dr = sqrt(xp**2+yp**2+zp**2)
! 	    write(6,1234) i,j,ip,jp,dr,n1,n2,n3,xp,yp,zp
1234	 format(4i3,f10.3,3i3,3f10.3)
!	 write(6,'(4f10.3)') xp,yp,zp,dr
             if ((dr .gt. rcut) .or. (dr .eq. 0.0)) go to 50
             nn = nn+ 1
             if (nn .gt. nat) stop ' nn too large '
             r(i,j,ip,jp,nn) = dr
             itype(i,j,ip,jp,nn) = ip
  50  continue

      neigh(i,j,ip,jp) = nn
!     do n = 1,neigh(i,j,ip,jp)
!     write(6,*) i,j,ip,jp,r(i,j,ip,jp,n)
!     enddo

 100  continue
 200  continue

      do 1000 i = 1,ntype
       do 1000 j = 1,natom(i)
        do 1000 ip = 1,ntype
         do 1000 jp = 1,natom(ip)
           do 60 m = 1,neigh(i,j,ip,jp)-1
            do 60 n = m+1,neigh(i,j,ip,jp)
             if (r(i,j,ip,jp,n) .lt. r(i,j,ip,jp,m)) then
             dump = r(i,j,ip,jp,n)
             r(i,j,ip,jp,n) = r(i,j,ip,jp,m)
             r(i,j,ip,jp,m) = dump
             idump = itype(i,j,ip,jp,n)
             itype(i,j,ip,jp,n) = itype(i,j,ip,jp,m)
             itype(i,j,ip,jp,m) = idump
             endif
  60  continue        
 1000 continue
!     crite(6,*) neigh(1,1,1,2)
!     do n=1,neigh(1,1,1,2)
!     write(6,*) r(1,1,1,2,n)
!     enddo

      do 2200 i = 1,ntype
       do 2200 j = 1,natom(i)
      write(22,111) fname(i),j
      write(6,111) fname(i),j
        do 2000 ip = 1,ntype
         do 2000 jp = 1,natom(ip)
      
!        write(6,*) ' neigh =',neigh(i,j,ip,jp)
         nshell = 1
          do n= 2,neigh(i,j,ip,jp)
           if (abs(r(i,j,ip,jp,n)-r(i,j,ip,jp,n-1)) .gt. 0.0001) then
           nshell = nshell + 1
            if(nshell .gt. nshelld) stop ' nshell too large '
            rr(nshell) = r(i,j,ip,jp,n)
            mo(nshell) = n
!           write(6,*) n,nshell,rr(nshell),mo(nshell)
           endif
          enddo
          rr(1) = r(i,j,ip,jp,1)
          num(1) = neigh(i,j,ip,jp)
          if (neigh(i,j,ip,jp).eq. 0.0) go to 2000
          if (nshell .eq. 1) go to 1999
!         mo(1) = mo(2) - 1
          mo(1) = 1
          mo(nshell+1) = neigh(i,j,ip,jp)
!         write(6,*) rr(1),mo(1),nshell,mo(nshell+1)
         
      do ii = 1,nshell
!     write(6,*) ii,i,j,ip,jp,mo(ii),rr(ii)
      enddo

!     num(1) = mo(2)-mo(1)
      do m = 1,nshell   
      num(m) = mo(m+1)-mo(m)
      num(nshell) = mo(nshell+1)-mo(nshell)+1
!     write(6,*) m,num(m),rr(m),num(nshell),rr(nshell) 
 111  format('------------------------------------',/,'[',a4,'-',i3,']')
 112  format(2x,a4,'-',i1,5(f7.3,i3))
 113  format(2x,a4,'-',i2,5(f7.3,i3))
 114  format(2x,a4,'-',i3,5(f7.3,i3))
      enddo
      if (nshell .gt. 5) nshell=5
1999  continue
!         write(6,112) fname(ip),jp,(rr(ii),num(ii),ii=1,nshell)
      if (jp .le. 9) then 
        write(22,112) fname(ip),jp,(rr(ii),num(ii),ii=1,nshell)
        write(6,112) fname(ip),jp,(rr(ii),num(ii),ii=1,nshell)
      end if
      if (jp .gt. 9 .and. jp .le. 99) then
        write(22,113) fname(ip),jp,(rr(ii),num(ii),ii=1,nshell)
        write(6,113) fname(ip),jp,(rr(ii),num(ii),ii=1,nshell)
      end if
      if (jp .gt. 99 ) then
        write(22,114) fname(ip),jp,(rr(ii),num(ii),ii=1,nshell)
        write(6,114) fname(ip),jp,(rr(ii),num(ii),ii=1,nshell)
      end if
2000  continue
2200  continue

      stop
      end
