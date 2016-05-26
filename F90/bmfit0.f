c     
c     Use the information in SUM and POSCAR to 
c     generate the input file for bmfit.f
c
c     modified by T.C. Leung May 2016
c
      implicit real*8 (a-h,o-z)
      parameter ( natmd = 100 )
      dimension a(3,3),aa(natmd),ee(natmd)
      character*64 temp1,temp2
      open(1,file='POSCAR')
      read(1,*) temp1
      read(1,*) acel  
      do i=1,3
      read(1,*) (a(i,j),j=1,3)
      enddo
      vol=a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)
     &     + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)
     &     - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2)
      write(6,*) ' vol of the cell is ',vol * acel**3
      open(5,file='SUM')
      read(5,*) temp1
      emin = 1000.0
      nn = 0
      do i = 1,100
      nn = nn + 1
      read(5,*,end=99) aa(i),temp1,temp2,ee(i)
      aa(i) = vol * aa(i)**3 
      if (ee(i) .lt. emin) then 
      emin = ee(i)
      inum = i
      endif
      enddo
  99  continue
      nn = nn - 1
c     write(6,*) ' mini at ',aa(inum),ee(inum)
c     write(6,*) ' num of data ',nn
      if ( nn .gt. natmd ) stop ' natmd too small '
      open(8,file='bmfit.dat')
      t=0.001
      t1=0.000001
      t2=0.1
      del = 0.05
      ndiv = 9
      shift = 0.0
      ift = 1
      nplot = 0
      write(8,11) ee(inum),t,aa(inum)
      write(8,11) ee(inum)-0.1, ee(inum)+0.1
      write(8,11) t1,t2
      write(8,11) aa(1),aa(nn)
      write(8,11) del
      write(8,12) ndiv
      write(8,12) nn,ift,nplot
      do i=1,nn
      write(8,11) aa(i),ee(i)
      enddo
      write(8,11) shift
 11   format(3f10.5)
 12   format(3i4)
      stop
      end

