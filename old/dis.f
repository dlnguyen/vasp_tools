      dimension a1(3),a2(3),a3(3)
      dimension xv(2),yv(2),zv(2)
      a1(1)=8.d0
      a1(2)=0.d0
      a1(3)=0.d0
      a2(1)=0.d0
      a2(2)=8.d0
      a2(3)=0.d0
      a3(1)=0.d0
      a3(2)=0.d0
      a3(3)=8.d0
      au=1.0d0
      read(5,*) 
      read(5,*) 
      read(5,*) 
      read(5,*)
      read(5,*)
c
      natom=2
      do imove=1,200
      mmod=mod(imove,3)
      read(*,*)
c      if(imove.eq.1600) then
c      write(1,3) natom
c      write(1,*)
c      endif
c      if(mmod.eq.1) then
      write(*,3) natom
      write(*,*)
c      endif
   3  format(i5)
c
      do i=1,natom
      read(5,*) x0,y0,z0
c      if(z0.gt.0.60) z0=z0-1.d0
c      if(y0.gt.0.80) y0=y0-1.d0
c      if(x0.gt.0.99) x0=x0-1.d0
c      if(z0.lt.0.05) z0=z0+1.d0
c      if(y0.lt.0.05) y0=y0+1.d0
c      if(x0.lt.0.01) x0=x0+1.d0
c      x0=dmod(x0+0.5d0,1.d0)
c      y0=dmod(y0+0.5d0,1.d0)
c      z0=dmod(z0+0.5d0,1.d0)
      xv(i)=x0*a1(1)
      yv(i)=y0*a2(2)
      zv(i)=z0*a3(3)
      xv(i)=au*xv(i)
      yv(i)=au*yv(i)
      zv(i)=au*zv(i)
      xvp=xv(i)
      yvp=yv(i)
      zvp=zv(i)
      if(zvp.gt.4.00) zvp=zvp-8.d0
      if(yvp.gt.4.00) yvp=yvp-8.d0
      if(xvp.gt.4.00) xvp=xvp-8.d0
c
c      if(imove.eq.1600) then
c      if(i.le.2) write(1,1) xv,yv,zv
c      if(i.gt.2) write(1,2) xv,yv,zv
c      endif
c      if(mmod.eq.1) then
      if(i.eq.1) write(6,1) xvp,yvp,zvp
      if(i.eq.2) write(6,2) xvp,yvp,zvp
c      endif
c      write(6,1) xv(i),yv(i),zv(i)
   
c    1 format(4HATOM,i7,3H  H,11x,1H1,3x,3f8.2)
c    2 format(4HATOM,i7,3H Pd,11x,1H1,3x,3f8.2)
    1 format(3H  C,3x,3f8.2)
    2 format(3H  O,3x,3f8.2)
c
c
      enddo
      do iii=1,natom
      read(5,*)
      enddo
      xd=xv(1)-xv(2)
      if(abs(xd).gt.4.0) xd=8.0-abs(xd)
      yd=yv(1)-yv(2)
      if(abs(yd).gt.4.0) yd=8.0-abs(yd)
      zd=zv(1)-zv(2)
      if(abs(zd).gt.4.0) zd=8.0-abs(zd)
      dis=SQRT(xd*xd+yd*yd+zd*zd)
      if(dis.gt.4) write(8,*) imove,xv(1),xv(2),yv(1),yv(2),zv(1),zv(2)
      write(9,4) imove,dis
      enddo
   4  format(i5,f12.7)
c
      stop
      end
