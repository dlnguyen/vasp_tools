!=========================================================================
! Program output-pdb-cif.f :
!
! File transfer program.  POSCAR -->  poscar.pdb,  poscar.xsf and  poscar.cif
!                        CONTCAR --> contcar.pdb, contcar.xsf and contcar.cif
!
! Jyh-Pin Chou,      June  2006
!
! This program is freely distributed and comes with no warranty. Please
! send any comments to the author (cinnamon@alumni.ccu.edu.tw)
!=========================================================================
      implicit real*8(a-h,o-z)
      character*64 temp,tmp,ext,output,chgtmp
      character*2 A_name(10),an(88)
      parameter(nptype=100)
      parameter(ndimx=100)
      parameter(ndimy=100)
      parameter(ndimz=500)
      dimension nt(nptype),iA_an(nptype),chg(ndimx,ndimy,ndimz)
      logical alive1,alive2,alive3,alive4,alive5,alive6
      dimension a(3,3)
      inquire(file='POSCAR' , exist=alive1)   
      inquire(file='POTCAR' , exist=alive4)
      inquire(file='CONTCAR', exist=alive2)  
      inquire(file='OUTCAR' , exist=alive3)  
      inquire(file='CHGCAR' , exist=alive5)
      inquire(file='PARCHG' , exist=alive6)
      pi=4.0d0*atan(1.0d0) ;   pi2=2.0d0*pi  ;  nchoice=1   ;   ntype=0
         an( 1)='H' ;an( 2)='He';an( 3)='Li';an( 4)='Be';an( 5)='B' 
         an( 6)='C' ;an( 7)='N' ;an( 8)='O' ;an( 9)='F' ;an(10)='Ne'
         an(11)='Na';an(12)='Mg';an(13)='Al';an(14)='Si';an(15)='P' 
         an(16)='S' ;an(17)='Cl';an(18)='Ar';an(19)='K' ;an(20)='Ca'
         an(21)='Sc';an(22)='Ti';an(23)='V' ;an(24)='Cr';an(25)='Mn'
         an(26)='Fe';an(27)='Co';an(28)='Ni';an(29)='Cu';an(30)='Zn'
         an(31)='Ga';an(32)='Ge';an(33)='As';an(34)='Se';an(35)='Br'
         an(36)='Kr';an(37)='Rb';an(38)='Sr';an(39)='Y' ;an(40)='Zr'
         an(41)='Nb';an(42)='Mo';an(43)='Tc';an(44)='Ru';an(45)='Rh'
         an(46)='Pd';an(47)='Ag';an(48)='Cd';an(49)='In';an(50)='Sn'
         an(51)='Sb';an(52)='Te';an(53)='I' ;an(54)='Xe';an(55)='Cs'
         an(56)='Ba';an(57)='La';an(58)='Ce';an(59)='Pr';an(60)='Nd'
         an(61)='Pm';an(62)='Sm';an(63)='Eu';an(64)='Gd';an(65)='Tb'
         an(66)='Dy';an(67)='Ho';an(68)='Er';an(69)='Tm';an(70)='Yb'
         an(71)='Lu';an(72)='Hf';an(73)='Ta';an(74)='W' ;an(75)='Re'
         an(76)='Os';an(77)='Ir';an(78)='Pt';an(79)='Au';an(80)='Hg'
         an(81)='Tl';an(82)='Pb';an(83)='Bi';an(84)='Po';an(85)='At'
         an(86)='Rn';an(87)='Fr';an(88)='Ra'
!==================================================================
      write(6,*) ' Enter your VASP version 4.x|5.x (0|1) '
      read(5,*) ver
      write(6,*) ' Selective dynamic yes|no (1|0) '
      read(5,*) isel
      do i=1,10
       A_name(i)='  '
      end do
      if(alive4)print *,'POTCAR Exist !'
      if(alive4)print *,'Open and read the kinds and the names of atoms'
      if(alive4) open(3,file='POTCAR')
      if(alive4) then
       open(99,status='scratch')
       do  ;  read(3,'(2a64)',iostat=istat) temp,tmp
        if(temp(4:8) .eq. 'TITEL') then
         ntype=ntype+1  ;  write(99,*) temp                               !    detect the kinds of atoms in POTCAR
        end if
        if(istat .ne. 0) exit
       end do
       rewind(99)  
       do i=1,ntype
        read(99,*) temp,temp,temp,A_name(i)                               !    detect the name of atoms in POTCAR
       end do  ;  close(99)
      end if

      if( alive3 .or. alive4) write(6,20) ntype,A_name
  20  format(2x,i3,1x,13hkinds of atom,2x,1h:,10a3)
!==================================================================
      if(alive3 .and. .not. alive4) print *,'OUTCAR Exist !'
      if(alive3 .and. .not. alive4) print *,'open and read the kinds and
     &the names of all atoms'
      if(alive3 .and. .not. alive4) open(3,file='OUTCAR')
      if(alive3 .and. .not. alive4) then
        do  ;  read(3,'(2a64)',iostat=istat) temp,tmp
          if(temp(4:8) .eq. 'TITEL') ntype=ntype+1                       !    detect the kinds of atoms in OUTCAR
          if(istat .ne. 0) exit
        end do
        rewind(3)  
        do  ;  read(3,'(a64)') temp
        if(temp(1:6) .eq. ' INCAR') exit
        end do  
        open(99,status='scratch')
        do i=1,ntype
          read(3,*) temp,temp,tmp  ;  write(99,*) tmp                    !     detect the name of atoms in OUTCAR
        end do  ;  rewind(99)
        do i=1,ntype
          read(99,*) A_name(i)
        end do
        close(99) 
      end if
!==================================================================
      if( .not. alive3 .and. .not. alive4) then
        print  *,'How many kinds of atoms in your system ?'
        read(5,*) ntype  
        print  *,'Input the names of the atoms : '
        print  *,'   Notice : there are 2 characters. (include space) '
        read(5,'(a2)') (A_name(i),i=1,ntype)
      end if  
      do
        print *,'Extend your system ? (y/n) :  (make supercell)'
        read(5,*) ext
!     ext='n'
        if(ext .eq. 'y') then
          call extender(ntype,alive1,alive2,isel)
          exit
        end if
        if(ext .eq. 'n') exit
        print *,'You must input y/n :'
      end do
!==================================================================
      if(alive1 .and. alive2) then
        if(ext .eq. 'y') open( 1,file="poscar-ext" )
        if(ext .eq. 'y') open(21,file="contcar-ext")
        if(ext .eq. 'n') open( 1,file="POSCAR")
        if(ext .eq. 'n') open(21,file="CONTCAR")
      else if(alive1 .and. .not. alive2) then
        if(ext .eq. 'y') open( 1,file="poscar-ext")   
        if(ext .eq. 'n') open( 1,file="POSCAR")
      else if(alive2 .and. .not. alive1) then
        if(ext .eq. 'y') open(21,file="contcar-ext")
        if(ext .eq. 'n') open(21,file="CONTCAR") 
      end if
      nfi =  1   ;  if(.not. alive1) nfi=21
      nff =  1   ;  if(      alive2) nff=21
!==================================================================
      do nf=nfi,nff,20
       read(nf,'(A1)') temp          ;      read(nf,*) aa
       read(nf,*) a(1,1),a(2,1),a(3,1)
       read(nf,*) a(1,2),a(2,2),a(3,2)
       read(nf,*) a(1,3),a(2,3),a(3,3)
c      if(ver .eq. 1 .and. alive2 .and. nf .eq. 21) read(nf,*)
c      if(ver.eq.1.and.alive2.and.nf.eq.21.and.ext.eq.'n') read(nf,*)
       if (ver .eq. 1) read(nf,*)
       read(nf,*,err=99999) (nt(i),i=1,ntype)
c      write(6,*) ' testing ',(nt(i),i=1,ntype)
c       write(6,*)'nf=',nf
       read(nf,'(A20)') temp                
c      write(6,'(A20)') temp
       if (isel .eq. 1) then
       read(nf,'(A20)') temp
c      write(6,'(A20)') temp
       endif
       natom = 0 
       do i=1,ntype
        natom = natom + nt(i)
       end do
c      write(6,*) ' natom =',natom
!      ##################################################################
         if(nf .eq.  1) open(2,file="poscar.pdb")
         if(nf .gt. 20) open(22,file="contcar.pdb")
         write(nf+1,'(6hCOMPND)')
         do j=1,ntype
          do i=1,nt(j)
           read(nf,*) x,y,z
c          write(6,*) x,y,z
           if(x > 0.5) x = x - 1.0
           if(y > 0.5) y = y - 1.0
           if(z > 0.5) z = z - 1.0
           xx = x*a(1,1) + y*a(1,2) + z*a(1,3)
           yy = x*a(2,1) + y*a(2,2) + z*a(2,3)
           zz = x*a(3,1) + y*a(3,2) + z*a(3,3)
           write(nf+1,3) i,A_name(j),xx*aa,yy*aa,zz*aa
          end do
         end do
         write(nf+1,4) natom
         write(nf+1,'(3hEND)') 
         if(nf .eq.  1) close(2)
         if(nf .gt. 20) close(22)
         rewind(nf)
c        if (ver .eq. 1 .and. alive2 .and. nf .eq. 21) then
c        if (ver.eq.1.and.alive2.and.nf.eq.21.and.ext.eq.'n')then
         if (ver .eq. 1) then
          ii=8
         else
          ii=7
         end if
         do i=1,ii
          read(nf,'(A20)') temp
c         write(6,'(A20)') temp
         end do
         if (isel .eq. 1) then 
         read(nf,'(A20)') temp
         endif
!      ##################################################################
         if(nf .eq.  1) open(2,file='poscar.cif')
         if(nf .gt. 20) open(22,file='contcar.cif')
              va = dsqrt( a(1,1)*a(1,1)+a(2,1)*a(2,1)+a(3,1)*a(3,1) )
              vb = dsqrt( a(1,2)*a(1,2)+a(2,2)*a(2,2)+a(3,2)*a(3,2) )
              vc = dsqrt( a(1,3)*a(1,3)+a(2,3)*a(2,3)+a(3,3)*a(3,3) )
         theta_a = acos ((a(1,3)*a(1,1)+a(2,3)*a(2,1)+a(3,3)*a(3,1) ) 
     &/ (vc*va) )
         theta_b = acos ((a(1,2)*a(1,3)+a(2,2)*a(2,3)+a(3,2)*a(3,3) ) 
     &/ (vb*vc) )
         theta_c = acos ((a(1,1)*a(1,2)+a(2,1)*a(2,2)+a(3,1)*a(3,2) ) 
     &/ (va*vb) )
         write(nf+1,'(5hdata_)')
         write(nf+1,'(5hloop_)')
         write(nf+1,'(26h_symmetry_equiv_pos_as_xyz)') 
         write(nf+1,'(7h  x,y,z)')
         write(nf+1,'(14h_cell_Length_a,f15.4)') aa*va
         write(nf+1,'(14h_cell_Length_b,f15.4)') aa*vb
         write(nf+1,'(14h_cell_Length_c,f15.4)') aa*vc
         write(nf+1,'(17h_cell_angle_alpha,f15.4)') theta_a*360/pi2
         write(nf+1,'(17h_cell_angle_beta ,f15.4)') theta_b*360/pi2
         write(nf+1,'(17h_cell_angle_gamma,f15.4)') theta_c*360/pi2
         write(nf+1,'(5hloop_)')
         write(nf+1,'(16h_atom_site_label)')
         write(nf+1,'(22h_atom_site_type_symbol)')
         write(nf+1,'(18h_atom_site_fract_x)')
         write(nf+1,'(18h_atom_site_fract_y)')
         write(nf+1,'(18h_atom_site_fract_z)')
         do j=1,ntype
          do i=1,nt(j)
     	   read(nf,*) x,y,z
!          if(x>0.5) x = x - 1 ; if(y>0.5) y = y-1 ; if(z>0.5) z = z - 1
           if(i .lt. 10   .and.                 
     &len_trim(A_name(j)) .eq. 1) write(nf+1,111) 
     &A_name(j),i,A_name(j),x,y,z
           if(i .ge. 10   .and. i .lt. 100   .and. 
     &len_trim(A_name(j)) .eq. 1) write(nf+1,112) 
     &A_name(j),i,A_name(j),x,y,z
           if(i .ge. 100  .and. i .lt. 1000  .and. 
     &len_trim(A_name(j)) .eq. 1) write(nf+1,113) 
     &A_name(j),i,A_name(j),x,y,z
           if(i .ge. 1000 .and. i .lt. 10000 .and. 
     &len_trim(A_name(j)) .eq. 1) write(nf+1,114) 
     &A_name(j),i,A_name(j),x,y,z
           if(i .lt. 10   .and.                    
     &len_trim(A_name(j)) .eq. 2) write(nf+1,211) 
     &A_name(j),i,A_name(j),x,y,z
           if(i .ge. 10   .and. i .lt. 100   .and. 
     &len_trim(A_name(j)) .eq. 2) write(nf+1,212) 
     &A_name(j),i,A_name(j),x,y,z
           if(i .ge. 100  .and. i .lt. 1000  .and. 
     &len_trim(A_name(j)) .eq. 2) write(nf+1,213) 
     &A_name(j),i,A_name(j),x,y,z
           if(i .ge. 1000 .and. i .lt. 10000 .and.
     &len_trim(A_name(j)) .eq. 2) write(nf+1,214) 
     &A_name(j),i,A_name(j),x,y,z
          end do
         end do
         if(nf .eq.  1) close(2)
         if(nf .gt. 20) close(22)
!      ##################################################################
         if(nf .eq.  1) open(2,file='poscar.xsf')
         if(nf .gt. 20) open(22,file='contcar.xsf')
!        if(alive5 .and. chgtmp .ne. 'n') then
         if(alive5 .or. alive6) then
         if( nf .eq. 1)    print *,' Do you want to add the information
     &of CHGCAR/PARCHG in poscar.xsf (y|n)?'
         if( nf .eq. 21)   print *,' Do you want to add the information
     &of CHGCAR/PARCHG in contcar.xsf (y|n)?'
           read(5,*) chgtmp
           if(chgtmp .eq. 'y') then
           print *,' Which file do you want to add CHGCAR|PARCHG (1|2) '
           read(5,*) ichg
             if (ichg .eq. 1) open(10,file='CHGCAR')
             if (ichg .eq. 2) open(10,file='PARCHG')
	       read(10,'(a64)') temp
c              write(6,'(a64)') temp
               read(10,*) aa
c              write(6,*) ' a= ',aa
               do i=1,3
               read(10,*) (a(i,j),j=1,3)
c              write(6,*) (a(i,j),j=1,3)
               enddo
               if (ver .eq. 1) read(10,'(a64)') temp
               read(10,*) (nt(j),j=1,ntype)
               nat=0
               do i=1,ntype
               nat=nat+nt(i)
               enddo
	       read(10,'(a64)') temp
c              write(6,'(a64)') temp

          
         write(nf+1,'("# XSF file, first line : MOLECULE 
     &0-D, POLYMER 1-D, SLAB 2-D, CRYSTAL 3-D")')
         write(nf+1,'("# or you can use DIM-GROUP and give the 
     &      number of dimension and group under the key word")')
         write(6,*) ' Enter the dimension of the system (0|1|2|3) '
         read(5,*) idimension
     
         if (idimension .eq. 0 ) write(nf+1,'(9h MOLECULE)')
         if (idimension .eq. 1 ) write(nf+1,'(8h POLYMER)')
         if (idimension .eq. 2 ) write(nf+1,'(5h SLAB)')
         if (idimension .eq. 3 ) write(nf+1,'(8h CRYSTAL)')
         write(nf+1,'(8h PRIMVEC)')
         write(nf+1,'(3f16.10)') a(1,1)*aa,a(2,1)*aa,a(3,1)*aa
         write(nf+1,'(3f16.10)') a(1,2)*aa,a(2,2)*aa,a(3,2)*aa
         write(nf+1,'(3f16.10)') a(1,3)*aa,a(2,3)*aa,a(3,3)*aa
         write(nf+1,'(10h PRIMCOORD)')
         write(nf+1,'(i3,i2)') nat,1
         write(6,*) ' Enter the shift of atoms along z-axis (angstrom)'
         read(5,*) zshift
         do j=1,ntype
          do i=1,88
           if(A_name(j) .eq. an(i)) iA_an(j)=i
          end do
         end do
         do j=1,ntype
          do i=1,nt(j)
           read(10,*) x,y,z
           if (z .gt. 0.5) z = z - 1.0
           xx = x*a(1,1) + y*a(1,2) + z*a(1,3)
           yy = x*a(2,1) + y*a(2,2) + z*a(2,3)
           zz = x*a(3,1) + y*a(3,2) + z*a(3,3)
           write(nf+1,'(i3,1x,3f17.10)') iA_an(j),xx*aa,yy*aa,zz*aa+zshift
          end do
         end do
             read(10,*)
	     read(10,*) nx,ny,nz  
c            write(6,*) nx,ny,nz
         if (nx .gt. ndimx) stop ' nx too large '
         if (ny .gt. ndimy) stop ' ny too large '
         if (nz .gt. ndimz) stop ' nz too large '
             write(6,*) ' nx,ny,nz ',nx,ny,nz
             del = aa*a(3,3)/float(nz)
             write(6,*) ' Spacing for grids along z-direction ',del
             read(10,*) (((chg(jj,kk,iz),jj=1,nx),kk=1,ny),iz=1,nz)
             close(10)
             write(nf+1,'(" ")')
             write(nf+1,'("BEGIN_BLOCK_DATAGRID_3D")')
             write(nf+1,'("CHGCAR_from_VASP      ")')
             write(nf+1,'("DATAGRID_3D_VASP      ")')
             write(nf+1,'(3i5)') nx,ny,nz
             write(nf+1,'(3i5)') 0,0,0
             write(nf+1,'(3f10.5)') a(1,1)*aa,a(2,1)*aa,a(3,1)*aa
             write(nf+1,'(3f10.5)') a(1,2)*aa,a(2,2)*aa,a(3,2)*aa
             write(nf+1,'(3f10.5)') a(1,3)*aa,a(2,3)*aa,a(3,3)*aa
             write(nf+1,9) (((chg(jj,kk,iz),jj=1,nx),kk=1,ny),iz=1,nz)
             write(nf+1,'("END_DATAGRID_3D_VASP      ")')
             write(nf+1,'("END_BLOCK_DATAGRID_3D_VASP      ")')
           end if
	 end if
         if(nf .eq.  1) close(2)
         if(nf .gt. 20) close(22)
      end do

    3   format('ATOM    ',i3,1x,a2,'           1    ',3f8.3,
     &         '  1.00  0.00') 
    4   format('TER      ',i3)
    7   format(a32)
    8   format(a32,f20.10)
    9   format(6e15.5)
  211   format(   a2,3h000,i1,3x,a2,1x,3f15.9,3x)
  212   format(   a2,2h00, i2,3x,a2,1x,3f15.9,3x)
  213   format(   a2,1h0,  i3,3x,a2,1x,3f15.9,3x)
  214   format(   a2,      i4,3x,a2,1x,3f15.9,3x)
  111   format(1x,a1,3h000,i1,4x,a1,1x,3f15.9,3x)
  112   format(1x,a1,2h00, i2,4x,a1,1x,3f15.9,3x)
  113   format(1x,a1,1h0,  i3,4x,a1,1x,3f15.9,3x)
  114   format(1x,a1,      i4,4x,a1,1x,3f15.9,3x)
        stop

99999 continue
      print   *,'ERROR, check program or input files'

      end 
!==================================================================
!==================================================================
       subroutine extender(nkatom,alive1,alive2,isel)
       implicit real*8(a-h,o-z)
       dimension av(3,3)
       logical alive1,alive2
       character*32 temp
       character*2 tx,ty,tz
       parameter(nka=80)
       dimension n(nka),x0(nka,nka),y0(nka,nka),z0(nka,nka)
       dimension tx(nka,nka),ty(nka,nka),tz(nka,nka)
       if(alive1) open( 1,file='POSCAR')   
       if(alive1) open( 7,file='poscar-ext')  
       if(alive2) open(21,file='CONTCAR')  
       if(alive2) open(27,file='contcar-ext')  

       print *,'Enter n1,n2,n3'   ;  read(5,*) n1,n2,n3
       nfi=1 ; if(.not. alive1) nfi=21
       nff=1 ; if(      alive2) nff=21
!     ##################################################################
       do nf=nfi,nff,20
        read (nf,*) temp         ;       read (nf,*) aa  
        write(nf+6,*) temp       ;       write(nf+6,'(f21.17)') aa
        read(nf,*) av(1,1),av(1,2),av(1,3)
        read(nf,*) av(2,1),av(2,2),av(2,3)
        read(nf,*) av(3,1),av(3,2),av(3,3)
        if (alive2 .and. nf .eq. 21) read(nf,*) temp
        read(nf,*) (n(i),i=1,nkatom)
        read(nf,*) temp          
        if ( isel .eq. 1)  read(nf,*) temp
        write(nf+6,'(3(2x,f20.16))') (n1*av(1,j),j=1,3)
        write(nf+6,'(3(2x,f20.16))') (n2*av(2,j),j=1,3)
        write(nf+6,'(3(2x,f20.16))') (n3*av(3,j),j=1,3)
        write(nf+6,'(8i6)') (n1*n2*n3*n(i),i=1,nkatom)
        if (isel .eq. 1) write(nf+6,'(18hSelective Dynamics)')
        write(nf+6,'(6hDirect)')
!     ----------------------------------------------------------------------------
        do 100 j=1,nkatom
         do i = 1,n(j)
          if (isel .eq. 0) 
     &     read(nf,*) x0(i,j),y0(i,j),z0(i,j)                        
          if (isel .eq. 1) 
     &     read(nf,*) x0(i,j),y0(i,j),z0(i,j),tx(i,j),ty(i,j),tz(i,j)
         end do
         do ni = 1,n(j)
          do 1000 ii=1,n1
           do 1000 jj=1,n2
            do 1000 kk=1,n3
             x  = ( x0(ni,j) + (ii-1) )/float(n1)
             y  = ( y0(ni,j) + (jj-1) )/float(n2)
             z  = ( z0(ni,j) + (kk-1) )/float(n3)
        if (isel .eq. 0)  write(nf+6,'(3f20.13)') x,y,z
        if (isel .eq. 1)  
     &  write(nf+6,'(3f20.13,3a2)') x,y,z,tx(ni,j),ty(ni,j),tz(ni,j)
 1000      continue
          end do 
 100    continue
        close(nf);close(nf+6)
       end do
!     ############################################################################
       return
       end
