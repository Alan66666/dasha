c    *************************************************
c             3-D COMPRESSIBLE HALL MHD CODE
c    *************************************************
c
c 1. Scheme  --- 4 order Runge-Kutta scheme. 
c 2. History --- completed in April, 1997.
c 3. Author  --- Z.W. Ma (Department of Physics and Astronomy, UI)
c 4. Project ---- To study the substorm dynamics with
c                empirical magnetotail geometry.
c 5. Parallelized by Lu Xingqiang in July, 2009
c                  Introduction
c The present version of this code simulate the whole simulation box
c in x-y-z domain. Periodic (or free) boundary conditions are used
c at the outgoin boundary. (Symmetric B.C. are used on one of the
c boundaries for symmetry case). The incoming boundary is a driven
c (or free) B.C.
c     **********************************************
c-----------------
c     Basic variable definitions:
c     x(mx,my,mz,1) to x(mx,my,mz,8) represents, respectively,
c     rho, rho*vx, rho*vy, rho*vz, bx, by, bz and energy.
c-----------------
c
c
c
      program main
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      logical :: exist      
	  
c*****Parallelized by Lu Xingqiang in July, 2009************
c Initiate MPI
      call mpi_init(ierror)    
      call mpi_comm_size(mpi_comm_world,nsize,ierror)
      call mpi_comm_rank(mpi_comm_world,nrank,ierror)
c
      if (nsize.ne.nprt) then
      print*,'The number of processors is not equal to', nprt
c     call mpi_abort(mpi_comm_world,1)
      stop
      endif
      if(mod(mx,nprx).gt.1.or.mod(mz,nprz).gt.1)then
      print*,'(mod(mx,nprx) or mod(mz,nprz) must be less than 1'
      stop
      endif
      call mpi_comm_rank(mpi_comm_world,nrank,ierr)
      myleft=nrank-1
      if((mod(myleft,nprx).eq.nprx-1)
     1   .or.(myleft.lt.0)) myleft=mpi_proc_null
      myright=nrank+1
      if(mod(myright,nprx).eq.0) myright=mpi_proc_null
      mylow=nrank-nprx
      if(mylow.lt.0) mylow=mpi_proc_null
      myupper=nrank+nprx
      if(myupper.gt.nsize-1) myupper=mpi_proc_null
      mepx=mod(nrank,nprx)
      mepz=int(nrank/nprx)
c
c锟斤拷锟姐范围锟斤拷锟斤拷锟斤拷锟竭界但锟斤拷锟斤拷锟斤拷锟节边界）
      xfirst=3
      xlast=mxpr-2
      zfirst=3
      zlast=mzpr-2
      if(mepx.eq.(nprx-1))xlast=mxpr-2+mod(mx,nprx)
      if(mepz.eq.(nprz-1))zlast=mzpr-2+mod(mz,nprz)
c     print*,nrank,xlast,zlast
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c
c锟斤拷锟姐范围锟斤拷锟斤拷锟斤拷锟斤拷锟竭界）
      x_first=3
      x_last=mxpr-2
      z_first=3
      z_last=mzpr-2
      if(mepx.eq.0)x_first=5
      if(mepx.eq.(nprx-1))x_last=mxpr-4+mod(mx,nprx)
      if(mepz.eq.0)z_first=5
      if(mepz.eq.(nprz-1))z_last=mzpr-4+mod(mz,nprz)
c     print*,nrank,x_first,x_last,z_first,z_last
c     print*,nrank,mepx,mepz
c     call mpi_type_contiguous(mxpr,mpi_double_precision,htype,ierr)
c     call mpi_type_vector(mzpr,1,(mxpr+1)+1,
c    1     mpi_double_precision,vtype,ierr)
c     call mpi_type_vector(mzpr,1,(mxpr+1)*(my+1)+1,
c    1     mpi_double_precision,3vtype,ierr)
c*****Parallelized by Lu Xingqiang in July, 2009************
      dtime=1.
c      dtime=5.
c     dnstep=5.
      dnstep=10.
      nstop=900000!1100
c      nstop=3
      nst=1
      nstnew=1
      dtstatus=1
      call input
c     if(nrank.eq.2.or.nrank.eq.5)then
c     print*,nrank,'xx=',xx
c     print*,nrank,'zz=',zz
c     endif
      call initia
c      if(nrank .eq. 24) print*,'main0 x(5)',(x(5,13,101,di8),di8=1,8)
c      if(nrank .eq. 24) print*,'main0 x(4)',(x(4,13,101,di8),di8=1,8)
c
c     call recrd2
c     print*,nst
c     stop
c
      call setdt
c      if(nrank.eq.0)print*,'main',nstep,time,dt
c
      dt0=dt
c
      if(.not.lrstrt)then
c     call recrd
c     if(uniformx.and.uniformz) then
c     call recrd1
c     else
c     call recrd2
c     call recrd3
c     endif
      endif
      if(lrstrt) then
   20 call readin
c      call recrdint
      call recrd2
      call recrd3  
c     nst=nst+nint
c     if(nst.le.nend) goto 20
c     stop
      endif
      timein=time
      nstin=nst
c   
  200 continue
c      if(nrank .eq. 24) print*,'main0 x(5)',(x(5,13,101,di8),di8=1,8)
c      if(nrank .eq. 24) print*,'main0 x(4)',(x(4,13,101,di8),di8=1,8)
      call setdt
      if(nrank.eq.0)print*,'main',nstep,time,dt
c     stop
      if(dt.lt.(0.1*dt0)) then
      dtstatus=0
      call setdt
      nst=nst+nint
c      call recrdpe
      call recrd2
      if(nrank.eq.0)print*,'recrd2 finished'
c      call recrd3
c      if(nrank.eq.0)print*,'recrd3 finished'
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      stop
      endif
c     dt=0.01
c
      if (nrank.eq.0) then
c  
        inquire(file='./data/stepnm', exist=exist)
        if (exist) then
          open(16, file='./data/stepnm', status='old', position='append'
     1     , action='write')
        else
          open(16,file='./data/stepnm', status='new', action='write')
        end if
        write(16,*) nstep, time, dt
        close(16)
c  
      endif
	  
      call stepon

c     print*,nst
c     stop
      nstep=nstep+1
      time=time+dt
c      if(nstep .le. 7) call recrdpe
      
      if(nrank.eq.0)then
      write(20,10)dt,time,nstep
   10 format(2(1x,f10.5),1x,i5)
      endif
c
c     if(mod(nstep,50).eq.0) then
      if(abs(time-timein-(nst-nstin+1)*dtime).le.dt.and.
     1  (nstep-nstp(nst)).ge.dnstep) then
      nstp(nst)=nstep
      nst=nst+nint
      if(mod(nst,3).eq.1) then
      call recrd
      endif
c     if(uniformx.and.uniformz) then
c     call recrd1
c     else
c      call recrdpe
      call recrd2
      call recrd3
c     endif
      endif
c
      if(nstep.gt.nstop) goto 300
      if(nst.lt.nend) goto 200
  300 close(11)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
c     stop
      end
c
      subroutine input
c --------------------------
c  This routine inputs parameters and define basic variables
c   LRSTRT: =.f., starting from t=0; =.t., continueing from
c           steps as given by NST.
c   NEND:   the final steps intended for the current run ,
c           including the steps from the previous run.
c   NSP:    time step interval for diagnostic plots.
c   DELTAV: thickness of the initial velocity shear.
c   DELTAC: thickness of the initial current sheet (for
c           normalization coventions , see routine INITIA.
c   RATIO:  ratio of dy to dx, or the ratio of the width to
c           length of the box, if the number of grid points is such
c           that my=mx.
c   NSMTH:  the number of rows starting from incoming
c           boundary for the smoothing region in x-direction.
c   ETA:    exact inverse of magnetic Renolds number.
c   GAMMA:  adiabatic constant.
c   BETAS:   ratio of kinetic to magnetic pressure at magnetosheath side
c   MACHS:   Alfven mach number of incoming flow  at magnetosheath side.
c   MACHM:   Alfven mach number of incoming flow at magnetopause side.
c   ANGS:    bz=b*cos(ang) and by=b*sin(ang), at magnetosheath side.
c   ANGM:    bz=b*cos(ang) and by=b*sin(ang), at magnetopause side.
c   TS0:   initia magnetosheath temperature
c   TM0:   initia magnetopause temperature.
c   BS0:   initia magnetosheath magnetic field strength
c   BM0:   initia magnetopause magnetic field strength.
c   NCASE:  case number of the run.
c   NST:    beginning data file number for the current run.
c   NINT:   data file number increment.
c --------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      include 'mhd03.for'
c
c      call gridpnt
      call gridpntnew2
c     do jx=1,mx1/2-44
c     xxt(jx)=xmin0+(jx-1)*(xmax-xmin0-6.)/(mx1-91)
c     enddo
c     do jx=mx1/2-44,mx1/2+46
c     xxt(jx)=3.*cos((jx-mx1/2-46)*pi/90.)
c     enddo
c     do jx=mx1/2+46,mx1
c     xxt(jx)=3.+(jx-mx1/2-46)*(xmax-xmin0-6.)/(mx1-91)
c     enddo
c     do jz=1,mz1/2-44
c     zzt(jz)=zmin0+(jz-1)*(zmax-zmin0-6.)/(mz1-91)
c     enddo
c     do jz=mz1/2-44,mz1/2+1
c     zzt(jz)=-sqrt(9.-xxt(mx1/2+1-abs(jz-1-mz1/2+45))**2)
c     enddo
c     do jz=mz1/2+1,mz1/2+46
c     zzt(jz)=sqrt(9.-xxt(mx1/2+1-abs(jz-1-mz1/2-45))**2)
c     enddo
c     do jz=mz1/2+46,mz1
c     zzt(jz)=3.+(jz-mz1/2-46)*(zmax-zmin0-6.)/(mz1-91)
c     enddo
c     do jx=1,mx1-1
c     do jz=1,mz1-1
c     dxt(jx)=xxt(jx+1)-xxt(jx)
c     dzt(jz)=zzt(jz+1)-zzt(jz)
c     enddo
c     enddo
c     dxt(mx1)=dxt(mx1-1)
c     dzt(mz1)=dzt(mz1-1)
      if(periodx) then
      xxt(mx)=xxt(mx1)+(xxt(2)-xxt(1))
      endif
      if(periody) then
      yy(my)=yy(my1)+(yy(2)-yy(1))
      endif
      if(periodz) then
      zzt(mz)=zzt(mz1)+(zzt(2)-zzt(1))
      endif
      if(nrank.eq.0)then
      open(unit=11,file='grid.dat',status='unknown',form='formatted')
      write(11,99)(xxt(jx),dxt(jx),jx=1,mx),
     1(zzt(jz),dzt(jz),jz=1,mz),(yy(jy),dy(jy),jy=1,my)
   99 format(6(1x,f13.5))
      close(11)
      endif
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c
      do jx=xfirst,xlast
      xx(jx)=xxt(mepx*(mxpr-4)+jx-2)
      dx(jx)=dxt(mepx*(mxpr-4)+jx-2)
      enddo
c     do jy=1,my
c     yy(jy)=yy(jy)
c     dy(jy)=dy(jy)
c     enddo
      do jz=zfirst,zlast
      zz(jz)=zzt(mepz*(mzpr-4)+jz-2)
      dz(jz)=dzt(mepz*(mzpr-4)+jz-2)
      enddo
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c
      call message1(xx,mxpr,myleft,myright)
      call message1(dx,mxpr,myleft,myright)
      call message1(zz,mzpr,mylow,myupper)
      call message1(dz,mzpr,mylow,myupper)
c
c 1. coefficient ax0,ay0,az0....dx0,dy0,dz0 related to the five-point second difference 
c with a fourth order accuracz and multiply by a factor of '2*h**4', 
c where h is the grid size in a single direction
c
c 2. coefficient ax1,ay1,az1....dx1,dy1,dz1 related to the five-point first differential 
c with a fourth order accuracz
c
c 3. coefficient ax2,ay2,az2....dx2,dy2,dz2 related to for the five-point second differential 
c with a fourth order accuracz for the uniform meshing and 
c third-order accuracz for the non-uniform grid
c
c 4. coefficient ax3,ay3,az3....dx3,dy3,dz3 related to for the fourth differential 
c with a fourth order accuracz for the uniform meshing and 
c third-order accuracz for the non-uniform grid
c
c 5. coefficient axp,axm....cxp,cxm work with left- or right-side first
c differential or boundary conditions with a third order accuracz
c
c 6. coefficient azp,azm....czp,czm work with low- or upper-side first 
c differential or boundary conditions with a third order accuracz
c
      do 1 jx=x_first,x_last
      dxp1=xx(jx+1)-xx(jx)
      dxm1=xx(jx)-xx(jx-1)
      dxp2=xx(jx+2)-xx(jx)
      dxm2=xx(jx)-xx(jx-2)
      f1= dxp1**2+dxp1*dxm2
      f2= dxp1**3-dxp1*dxm2**2
      f3= dxp1**4+dxp1*dxm2**3
      g1=-dxm1**2+dxm1*dxm2
      g2= dxm1**3-dxm1*dxm2**2
      g3=-dxm1**4+dxm1*dxm2**3
      h1= dxp2**2+dxp2*dxm2
      h2= dxp2**3-dxp2*dxm2**2
      h3= dxp2**4+dxp2*dxm2**3
      ca= (f1*h3-f3*h1)*(g2*h3-g3*h2)
     1   -(g1*h3-g3*h1)*(f2*h3-f3*h2)
      
c      ax0(jx)=dxp1*dxm1*dxp2*dxm2*h3*(g2*h3-g3*h2)/ca
c      bx0(jx)=dxp1*dxm1*dxp2*dxm2*h3*(h2*f3-h3*f2)/ca
c      cx0(jx)=dxp1*dxm1*dxp2*dxm2*h3*(f2*g3-f3*g2)/ca
c      dx0(jx)=-(dxp1*ax0(jx)+dxm1*bx0(jx)
c     1   +dxp2*cx0(jx))/dxm2
c     2*h**4*f''(x)=ax0*[f(x+h)-f(x)]+bx0*[f(x)-f(x-h)]+cx0*[f(x+2h)-f(x)]+dx0*[f(x)-f(x-2h)]+o(h^5)
      
      ax2(jx)=2.*h3*(g2*h3-g3*h2)/ca
      bx2(jx)=2.*h3*(h2*f3-h3*f2)/ca
      cx2(jx)=2.*h3*(f2*g3-f3*g2)/ca
      dx2(jx)=-(dxp1*ax2(jx)+dxm1*bx2(jx)
     1   +dxp2*cx2(jx))/dxm2
c     f''(x)=ax2*[f(x+h)-f(x)]+bx2*[f(x)-f(x-h)]+cx2*[f(x+2h)-f(x)]+dx2*[f(x)-f(x-2h)]+o(h^5)
      
c      ax3(jx)=24.*h3*(g1*h2-g2*h1)/ca
c      bx3(jx)=24.*h3*(h1*f2-h2*f1)/ca
c      cx3(jx)=24.*h3*(f1*g2-f2*g1)/ca
c      dx3(jx)=-(dxp1*ax3(jx)+dxm1*bx3(jx)
c     1   +dxp2*cx3(jx))/dxm2
c     f''''(x)=ax3*[f(x+h)-f(x)]+bx3*[f(x)-f(x-h)]+cx3*[f(x+2h)-f(x)]+dx3*[f(x)-f(x-2h)]+o(h^5)
    1 continue
      do 12 jy=3,my-2
      dyp1=yy(jy+1)-yy(jy)
      dym1=yy(jy)-yy(jy-1)
      dyp2=yy(jy+2)-yy(jy)
      dym2=yy(jy)-yy(jy-2)
      f1= dyp1**2+dyp1*dym2
      f2= dyp1**3-dyp1*dym2**2
      f3= dyp1**4+dyp1*dym2**3
      g1=-dym1**2+dym1*dym2
      g2= dym1**3-dym1*dym2**2
      g3=-dym1**4+dym1*dym2**3
      h1= dyp2**2+dyp2*dym2
      h2= dyp2**3-dyp2*dym2**2
      h3= dyp2**4+dyp2*dym2**3
      ca= (f1*h3-f3*h1)*(g2*h3-g3*h2)
     1   -(g1*h3-g3*h1)*(f2*h3-f3*h2)
c      ay0(jy)=dyp1*dym1*dyp2*dym2*h3*(g2*h3-g3*h2)/ca
c      by0(jy)=dyp1*dym1*dyp2*dym2*h3*(h2*f3-h3*f2)/ca
c      cy0(jy)=dyp1*dym1*dyp2*dym2*h3*(f2*g3-f3*g2)/ca
c      dy0(jy)=-(dyp1*ay0(jy)+dym1*by0(jy)
c     1   +dyp2*cy0(jy))/dym2
      ay2(jy)=2.*h3*(g2*h3-g3*h2)/ca
      by2(jy)=2.*h3*(h2*f3-h3*f2)/ca
      cy2(jy)=2.*h3*(f2*g3-f3*g2)/ca
      dy2(jy)=-(dyp1*ay2(jy)+dym1*by2(jy)
     1   +dyp2*cy2(jy))/dym2
c      ay3(jy)=24.*h3*(g1*h2-g2*h1)/ca
c      by3(jy)=24.*h3*(h1*f2-h2*f1)/ca
c      cy3(jy)=24.*h3*(f1*g2-f2*g1)/ca
c      dy3(jy)=-(dyp1*ay3(jy)+dym1*by3(jy)
c     1   +dyp2*cy3(jy))/dym2
   12 continue
      do 2 jz=z_first,z_last
      dzp1=zz(jz+1)-zz(jz)
      dzm1=zz(jz)-zz(jz-1)
      dzp2=zz(jz+2)-zz(jz)
      dzm2=zz(jz)-zz(jz-2)
      f1= dzp1**2+dzp1*dzm2
      f2= dzp1**3-dzp1*dzm2**2
      f3= dzp1**4+dzp1*dzm2**3
      g1=-dzm1**2+dzm1*dzm2
      g2= dzm1**3-dzm1*dzm2**2
      g3=-dzm1**4+dzm1*dzm2**3
      h1= dzp2**2+dzp2*dzm2
      h2= dzp2**3-dzp2*dzm2**2
      h3= dzp2**4+dzp2*dzm2**3
      ca= (f1*h3-f3*h1)*(g2*h3-g3*h2)
     1   -(g1*h3-g3*h1)*(f2*h3-f3*h2)
c      az0(jz)=dzp1*dzm1*dzp2*dzm2*h3*(g2*h3-g3*h2)/ca
c      bz0(jz)=dzp1*dzm1*dzp2*dzm2*h3*(h2*f3-h3*f2)/ca
c      cz0(jz)=dzp1*dzm1*dzp2*dzm2*h3*(f2*g3-f3*g2)/ca
c      dz0(jz)=-(dzp1*az0(jz)+dzm1*bz0(jz)
c     1   +dzp2*cz0(jz))/dzm2
      az2(jz)=2.*h3*(g2*h3-g3*h2)/ca
      bz2(jz)=2.*h3*(h2*f3-h3*f2)/ca
      cz2(jz)=2.*h3*(f2*g3-f3*g2)/ca
      dz2(jz)=-(dzp1*az2(jz)+dzm1*bz2(jz)
     1   +dzp2*cz2(jz))/dzm2
c      az3(jz)=24.*h3*(g1*h2-g2*h1)/ca
c      bz3(jz)=24.*h3*(h1*f2-h2*f1)/ca
c      cz3(jz)=24.*h3*(f1*g2-f2*g1)/ca
c      dz3(jz)=-(dzp1*az3(jz)+dzm1*bz3(jz)
c     1   +dzp2*cz3(jz))/dzm2
    2 continue
      
      do 3 jx=x_first,x_last
      dxp1=xx(jx+1)-xx(jx)
      dxm1=xx(jx)-xx(jx-1)
      dxp2=xx(jx+2)-xx(jx)
      dxm2=xx(jx)-xx(jx-2)
      f1= dxp1   +dxp1**4/dxm2**3
      f2= dxp1**2-dxp1**4/dxm2**2
      f3= dxp1**3+dxp1**4/dxm2
      g1= dxm1   -dxm1**4/dxm2**3
      g2=-dxm1**2+dxm1**4/dxm2**2
      g3= dxm1**3-dxm1**4/dxm2
      h1= dxp2   +dxp2**4/dxm2**3
      h2= dxp2**2-dxp2**4/dxm2**2
      h3= dxp2**3+dxp2**4/dxm2
      ca= (f1*h3-f3*h1)*(g2*h3-g3*h2)
     1   -(g1*h3-g3*h1)*(f2*h3-f3*h2)
      ax1(jx)=h3*(g2*h3-g3*h2)/ca
      bx1(jx)=h3*(h2*f3-h3*f2)/ca
      cx1(jx)=h3*(f2*g3-f3*g2)/ca
      dx1(jx)=(dxp1**2*ax1(jx)-dxm1**2*bx1(jx)
     1   +dxp2**2*cx1(jx))/dxm2**2
    3 continue
      do 9 jy=3,my-2
      dyp1=yy(jy+1)-yy(jy)
      dym1=yy(jy)-yy(jy-1)
      dyp2=yy(jy+2)-yy(jy)
      dym2=yy(jy)-yy(jy-2)
      f1= dyp1   +dyp1**4/dym2**3
      f2= dyp1**2-dyp1**4/dym2**2
      f3= dyp1**3+dyp1**4/dym2
      g1= dym1   -dym1**4/dym2**3
      g2=-dym1**2+dym1**4/dym2**2
      g3= dym1**3-dym1**4/dym2
      h1= dyp2   +dyp2**4/dym2**3
      h2= dyp2**2-dyp2**4/dym2**2
      h3= dyp2**3+dyp2**4/dym2
      ca= (f1*h3-f3*h1)*(g2*h3-g3*h2)
     1   -(g1*h3-g3*h1)*(f2*h3-f3*h2)
      ay1(jy)=h3*(g2*h3-g3*h2)/ca
      by1(jy)=h3*(h2*f3-h3*f2)/ca
      cy1(jy)=h3*(f2*g3-f3*g2)/ca
      dy1(jy)=(dyp1**2*ay1(jy)-dym1**2*by1(jy)
     1   +dyp2**2*cy1(jy))/dym2**2
    9 continue
      do 4 jz=z_first,z_last
      dzp1=zz(jz+1)-zz(jz)
      dzm1=zz(jz)-zz(jz-1)
      dzp2=zz(jz+2)-zz(jz)
      dzm2=zz(jz)-zz(jz-2)
      f1= dzp1   +dzp1**4/dzm2**3
      f2= dzp1**2-dzp1**4/dzm2**2
      f3= dzp1**3+dzp1**4/dzm2
      g1= dzm1   -dzm1**4/dzm2**3
      g2=-dzm1**2+dzm1**4/dzm2**2
      g3= dzm1**3-dzm1**4/dzm2
      h1= dzp2   +dzp2**4/dzm2**3
      h2= dzp2**2-dzp2**4/dzm2**2
      h3= dzp2**3+dzp2**4/dzm2
      ca= (f1*h3-f3*h1)*(g2*h3-g3*h2)
     1   -(g1*h3-g3*h1)*(f2*h3-f3*h2)
      az1(jz)=h3*(g2*h3-g3*h2)/ca
      bz1(jz)=h3*(h2*f3-h3*f2)/ca
      cz1(jz)=h3*(f2*g3-f3*g2)/ca
      dz1(jz)=(dzp1**2*az1(jz)-dzm1**2*bz1(jz)
     1   +dzp2**2*cz1(jz))/dzm2**2
    4 continue
      do 5 jx=xfirst,xlast-3
      dxp1=xx(jx+1)-xx(jx)
      dxp2=xx(jx+2)-xx(jx)
      dxp3=xx(jx+3)-xx(jx)
      f1=dxp1-dxp1**3/dxp3**2
      f2=dxp1**2-dxp1**3/dxp3
      g1=dxp2-dxp2**3/dxp3**2
      g2=dxp2**2-dxp2**3/dxp3
      ca=f1*g2-f2*g1
      axp(jx)=g2/ca
      bxp(jx)=-f2/ca
      cxp(jx)=(1-axp(jx)*dxp1-bxp(jx)*dxp2)/dxp3
    5 continue
      do 6 jx=xfirst+3,xlast
      dxm1=xx(jx)-xx(jx-1)
      dxm2=xx(jx)-xx(jx-2)
      dxm3=xx(jx)-xx(jx-3)
      f1=dxm1-dxm1**3/dxm3**2
      f2=dxm1**2-dxm1**3/dxm3
      g1=dxm2-dxm2**3/dxm3**2
      g2=dxm2**2-dxm2**3/dxm3
      ca=f1*g2-f2*g1
      axm(jx)=g2/ca
      bxm(jx)=-f2/ca
      cxm(jx)=(1-axm(jx)*dxm1-bxm(jx)*dxm2)/dxm3
    6 continue
      do 10 jy=1,my-3
      dyp1=yy(jy+1)-yy(jy)
      dyp2=yy(jy+2)-yy(jy)
      dyp3=yy(jy+3)-yy(jy)
      f1=dyp1-dyp1**3/dyp3**2
      f2=dyp1**2-dyp1**3/dyp3
      g1=dyp2-dyp2**3/dyp3**2
      g2=dyp2**2-dyp2**3/dyp3
      ca=f1*g2-f2*g1
      ayp(jy)=g2/ca
      byp(jy)=-f2/ca
      cyp(jy)=(1-ayp(jy)*dyp1-byp(jy)*dyp2)/dyp3
   10 continue
      do 11 jy=4,my
      dym1=yy(jy)-yy(jy-1)
      dym2=yy(jy)-yy(jy-2)
      dym3=yy(jy)-yy(jy-3)
      f1=dym1-dym1**3/dym3**2
      f2=dym1**2-dym1**3/dym3
      g1=dym2-dym2**3/dym3**2
      g2=dym2**2-dym2**3/dym3
      ca=f1*g2-f2*g1
      aym(jy)=g2/ca
      bym(jy)=-f2/ca
      cym(jy)=(1-aym(jy)*dym1-bym(jy)*dym2)/dym3
   11 continue
      do 7 jz=zfirst,zlast-3
      dzp1=zz(jz+1)-zz(jz)
      dzp2=zz(jz+2)-zz(jz)
      dzp3=zz(jz+3)-zz(jz)
      f1=dzp1-dzp1**3/dzp3**2
      f2=dzp1**2-dzp1**3/dzp3
      g1=dzp2-dzp2**3/dzp3**2
      g2=dzp2**2-dzp2**3/dzp3
      ca=f1*g2-f2*g1
      azp(jz)=g2/ca
      bzp(jz)=-f2/ca
      czp(jz)=(1-azp(jz)*dzp1-bzp(jz)*dzp2)/dzp3
    7 continue
      do 8 jz=zfirst+3,zlast
      dzm1=zz(jz)-zz(jz-1)
      dzm2=zz(jz)-zz(jz-2)
      dzm3=zz(jz)-zz(jz-3)
      f1=dzm1-dzm1**3/dzm3**2
      f2=dzm1**2-dzm1**3/dzm3
      g1=dzm2-dzm2**3/dzm3**2
      g2=dzm2**2-dzm2**3/dzm3
      ca=f1*g2-f2*g1
      azm(jz)=g2/ca
      bzm(jz)=-f2/ca
      czm(jz)=(1-azm(jz)*dzm1-bzm(jz)*dzm2)/dzm3
    8 continue
      if(myleft.eq.mpi_proc_null)then
      cxp1=(xx(xfirst+1)-xx(xfirst))**2
      cxp2=(xx(xfirst+2)-xx(xfirst))**2
      cxs2=(xx(xfirst+1)-xx(xfirst))/(xx(xfirst+2)-xx(xfirst+1))
      endif
      if(myright.eq.mpi_proc_null)then
      cxm1=(xx(xlast)-xx(xlast-1))**2
      cxm2=(xx(xlast)-xx(xlast-2))**2
      cxe2=(xx(xlast)-xx(xlast-1))/(xx(xlast-1)-xx(xlast-2))
      endif
      cyp1=(yy(2)-yy(1))**2
      cyp2=(yy(3)-yy(1))**2
      cym1=(yy(my)-yy(my-1))**2
      cym2=(yy(my)-yy(my-2))**2
      if(mylow.eq.mpi_proc_null)then
      czp1=(zz(zfirst+1)-zz(zfirst))**2
      czp2=(zz(zfirst+2)-zz(zfirst))**2
      czs2=(zz(zfirst+1)-zz(zfirst))/(zz(zfirst+2)-zz(zfirst+1))
      endif
      if(myupper.eq.mpi_proc_null)then
      czm1=(zz(zlast)-zz(zlast-1))**2
      czm2=(zz(zlast)-zz(zlast-2))**2
      cze2=(zz(zlast)-zz(zlast-1))/(zz(zlast-1)-zz(zlast-2))
      endif
      time=0.
      nstep=0
      return
      end
c
      subroutine initia
c----------------
c Defines coordinates system and specifies initial configuration.
c Normlization convention:
c   1. Density --- normalised to asymtotic value, i.e., rho=1
c   2. Magnetic field --- normalised to asymtotic value, i.e.,
c                         b0=1.
c   3. Velocity --- normalised to asymtotic Alfven speed, VA=1, a
c                   natural result of 1. and 2.
c   4. Length --- normalised to a=10*dx, i.e., dx=0.1
c   5. Time --- normalised to a/VA.
c   6. Parallelized by Lu Xingqiang in July, 2009
c---------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      dimension sxx(mx),szz(mz2),spsi(mx,mz2),scur(mx,mz2),spr(mx,mz2)
c     common /cdata/ xl(20)
      dimension work1(mxpr),psi00(mxpr,mzpr,nsize)

      dimension bxo(mxpr*my*mzpr),byo(mxpr*my*mzpr),bzo(mxpr*my*mzpr)
     1      ,bxin(mxpr,my,mzpr),byin(mxpr,my,mzpr),bzin(mxpr,my,mzpr)
	 
      dimension bxofc(my*mzpr),byofc(my*mzpr),bzofc(my*mzpr)
	  
      integer n_x,n_z
      character*13 cntb
      character*15 cntbfc	  
      character*3 cn	  
	  
c Coordinates system
c                A x
c                l
c                l
c                l
c   ----------------------------> y
c                l
c                l
c                l
c                l
c
c assign asymmetric quantities:
      iopt=5
      do jx=1,mx
      sxx(jx)=xxt(jx)
      enddo
      if(symmetryz) then
      do jz=1,mz
      szz(jz)=zzt(jz)
      szz(mz2-jz+1)=-szz(jz)
      enddo
      else
      do jz=1,mz
      szz(jz)=zzt(jz)
      enddo
      endif
c useless:sxx,szz

      call mpi_comm_rank(mpi_comm_world,nrank,ierror)

      if (mxdat) then
        cntb='./inp/binx'//cn(nrank)      
        open(881,file=cntb)
        DO
          read(881,*,iostat=ierr) bxo
          if (ierr/=0) Exit
        ENDDO
        close(881) 
             
        cntb='./inp/biny'//cn(nrank)       
        open(882,file=cntb)
        DO
          read(882,*,iostat=ierr) byo         
          if (ierr/=0) Exit
        ENDDO
        close(882) 
        
        cntb='./inp/binz'//cn(nrank)       
        open(883,file=cntb)
        DO
          read(883,*,iostat=ierr) bzo         
          if (ierr/=0) Exit
        ENDDO
        close(883)                
              
        do 353 jy=1,my
          do 353 jx=x_first-2,x_last+2
            do 353 jz=z_first-2,z_last+2
            n_x=x_last-x_first+5;
            n_z=z_last-z_first+5;
            if((mepx.eq.0).or.(mepz.eq.0)) then
             if (mepx.eq.0) then
              if (mepz.eq.0) then 
               bxin(jx,jy,jz)=bxo(n_x*n_z*(jy-1)+n_z*(jx-1-2)+jz-2)
               byin(jx,jy,jz)=byo(n_x*n_z*(jy-1)+n_z*(jx-1-2)+jz-2)
               bzin(jx,jy,jz)=bzo(n_x*n_z*(jy-1)+n_z*(jx-1-2)+jz-2) 
              else
               bxin(jx,jy,jz)=bxo(n_x*n_z*(jy-1)+n_z*(jx-1-2)+jz)
               byin(jx,jy,jz)=byo(n_x*n_z*(jy-1)+n_z*(jx-1-2)+jz)
               bzin(jx,jy,jz)=bzo(n_x*n_z*(jy-1)+n_z*(jx-1-2)+jz)
              endif
            else
               bxin(jx,jy,jz)=bxo(n_x*n_z*(jy-1)+n_z*(jx-1)+jz-2)
               byin(jx,jy,jz)=byo(n_x*n_z*(jy-1)+n_z*(jx-1)+jz-2)
               bzin(jx,jy,jz)=bzo(n_x*n_z*(jy-1)+n_z*(jx-1)+jz-2)
            endif
            else
             bxin(jx,jy,jz)=bxo(n_x*n_z*(jy-1)+n_z*(jx-1)+jz)
             byin(jx,jy,jz)=byo(n_x*n_z*(jy-1)+n_z*(jx-1)+jz)
             bzin(jx,jy,jz)=bzo(n_x*n_z*(jy-1)+n_z*(jx-1)+jz) 			 
            endif
  353   continue
      endif
	  
	  if (mxdatfc) then
	  if(myleft.eq.mpi_proc_null)then
        cntbfc='./fcinp/binx'//cn(nrank)      
        open(991,file=cntbfc)
        DO
          read(991,*,iostat=ierr) bxofc
          if (ierr/=0) Exit
        ENDDO
        close(991) 
             
        cntbfc='./fcinp/biny'//cn(nrank)       
        open(992,file=cntbfc)
        DO
          read(992,*,iostat=ierr) byofc         
          if (ierr/=0) Exit
        ENDDO
        close(992) 
        
        cntbfc='./fcinp/binz'//cn(nrank)       
        open(993,file=cntbfc)
        DO
          read(993,*,iostat=ierr) bzofc         
          if (ierr/=0) Exit
        ENDDO
        close(993)                
              
        do 354 jy=1,my
          do 354 jz=z_first-2,z_last+2
            n_z=z_last-z_first+5;
              if (mepz.eq.0) then 
               bxinfc(jy,jz)=bxofc(n_z*(jy-1)+jz-2)*(-50.)
               byinfc(jy,jz)=byofc(n_z*(jy-1)+jz-2)*(-50.)
               bzinfc(jy,jz)=bzofc(n_z*(jy-1)+jz-2)*(-50.)
              else
               bxinfc(jy,jz)=bxofc(n_z*(jy-1)+jz)*(-50.)
               byinfc(jy,jz)=byofc(n_z*(jy-1)+jz)*(-50.)
               bzinfc(jy,jz)=bzofc(n_z*(jy-1)+jz)*(-50.)
              endif
  354   continue
      endif
	  endif
	  	  
      do 101 jz=z_first-2,z_last+2
      do 101 jy=1,my
      do 101 jx=x_first-2,x_last+2
      rr0(jx,jy,jz)=sqrt(xx(jx)**2+yy(jy)**2+zz(jz)**2)
      rr1(jx,jy,jz)=sqrt((xx(jx)-0.0)**2+yy(jy)**2+zz(jz)**2)
  101 continue
c      
      do 10 jz=z_first-2,z_last+2
      do 10 jy=1,my
      do 10 jx=x_first-2,x_last+2
c rho	  n  
      x(jx,jy,jz,1)  = 1.d-3*exp((-(rr0(jx,jy,jz)-1.04)*
     1    (rr0(jx,jy,jz)-1.04)/(2.*0.15*0.15)))+0.1*1.d-3
c vx,vy,vz	  
      x(jx,jy,jz,2)  = 0.
      x(jx,jy,jz,3)  = 0.
      x(jx,jy,jz,4)  = 0.
c bx,by,bz	  
      if (mxdat) then
      x(jx,jy,jz,5)  = bxin(jx,jy,jz)*(-25./3.)
      x(jx,jy,jz,6)  = byin(jx,jy,jz)*(-25./3.)
      x(jx,jy,jz,7)  = bzin(jx,jy,jz)*(-25./3.)
      else
      if(rr1(jx,jy,jz).ge.0.1) then
      x(jx,jy,jz,5)  = 14.0*(-3.)*xx(jx)*zz(jz)/rr0(jx,jy,jz)**5
      x(jx,jy,jz,6)  = 14.0*(-3.)*yy(jy)*zz(jz)/rr0(jx,jy,jz)**5
      x(jx,jy,jz,7)  = 14.0*(xx(jx)**2+yy(jy)**2-2.*zz(jz)**2)
     1                 /rr0(jx,jy,jz)**5
      else
      x(jx,jy,jz,5)  = 0.0
      x(jx,jy,jz,6)  = 0.0
      x(jx,jy,jz,7)  = 0.0
      endif
      endif 
c eq	  
      rhoi(jx,jy,jz)=x(jx,jy,jz,1)
      bxi(jx,jy,jz)=x(jx,jy,jz,5)
      byi(jx,jy,jz)=x(jx,jy,jz,6)
      bzi(jx,jy,jz)=x(jx,jy,jz,7)
c      pri=0.02
      pri(jx,jy,jz)=0.2*1.d-3*exp((-(rr0(jx,jy,jz)-1.04)*
     1        (rr0(jx,jy,jz)-1.04)/(2.*0.15*0.15)))+0.1*0.2*1.d-3
c      pri(jx,jy,jz)=0.3465*0.5
      x(jx,jy,jz,8) = pri(jx,jy,jz)
   10 continue
c
      do 102 m=1,3
      do 102 jz=z_first-2,z_last+2
      do 102 jy=1,my
      do 102 jx=x_first-2,x_last+2
        xresd(jx,jy,jz,m)=0.
  102 continue
      do 104 m=1,8
      do 104 jz=z_first-2,z_last+2
      do 104 jy=1,my
      do 104 jx=x_first-2,x_last+2
        xfold(jx,jy,jz,m)=x(jx,jy,jz,m)
  104 continue
      call right(0d0,1.d-9,1)
      do 103 m=1,3
      do 103 jz=z_first-2,z_last+2
      do 103 jy=1,my
      do 103 jx=x_first-2,x_last+2
        xresd(jx,jy,jz,m)=xdif(jx,jy,jz,m+1)
  103 continue
c
cmao      call diagn_res
c     call recrd3
c     print*,nst
c     stop
      call foreta(time)
      call formu(time)
      call current
c
c      call energy(time)
c      call recrdpe
c      call recrdint
      call recrd2
      call recrd3  
      return
      end
c
c
c
      subroutine extrap(bin,xi,zi,nx,nz,bl,xl,zl,kx,kz)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bin(nx,nz)
      dimension xl(kx),zl(kz),bl(kx,kz)
      do 30 jz=1,kz
      do 30 jx=1,kx
      do 10 ix=1,nx-1
      if((xi(ix)-1.e-5).le.xl(jx).and.(xi(ix+1)+1.e-5).ge.xl(jx)) then 
      beta1=(xl(jx)-xi(ix))/(xi(ix+1)-xi(ix))
      lx=ix
      goto 15
      else
      endif
   10 continue
   15 do 20 iz=1,nz-1
      if((zi(iz)-1.e-5).le.zl(jz).and.(zi(iz+1)+1.e-5).ge.zl(jz)) then
      beta2=(zl(jz)-zi(iz))/(zi(iz+1)-zi(iz))
      lz=iz
      goto 25
      else
      endif
   20 continue
   25 bl(jx,jz)=bin(lx,lz)*(1-beta1)*(1-beta2)
     1           +bin(lx+1,lz)*beta1*(1-beta2)
     1           +bin(lx,lz+1)*(1-beta1)*beta2
     1           +bin(lx+1,lz+1)*beta1*beta2
   30 continue
c
      return
      end
c
c
c
      subroutine extrap3(bin,xi,yi,zi,nx,ny,nz,bl,xl,yl,zl,kx,ky,kz)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),yi(ny),zi(nz),bin(nx,ny,nz)
      dimension xl(kx),yl(ky),zl(kz),bl(kx,ky,kz)
      do 40 jz=1,kz
      do 40 jy=1,ky
      do 40 jx=1,kx
      do 10 ix=1,nx-1
      if((xi(ix)-1.e-5).le.xl(jx).and.(xi(ix+1)+1.e-5).ge.xl(jx)) then 
      beta1=(xl(jx)-xi(ix))/(xi(ix+1)-xi(ix))
      lx=ix
      goto 15
      else
      endif
   10 continue
   15 do 20 iy=1,ny-1
      if((yi(iy)-1.e-5).le.yl(jy).and.(yi(iy+1)+1.e-5).ge.yl(jy)) then
      beta2=(yl(jy)-yi(iy))/(yi(iy+1)-yi(iy))
      ly=iy
      goto 25
      else
      endif
   20 continue
   25 do 30 iz=1,nz-1
      if((zi(iz)-1.e-5).le.zl(jz).and.(zi(iz+1)+1.e-5).ge.zl(jz)) then
      beta3=(zl(jz)-zi(iz))/(zi(iz+1)-zi(iz))
      lz=iz
      goto 35
      else
      endif
   30 continue
   35 bl(jx,jy,jz)=bin(lx,ly,lz)*(1-beta1)*(1-beta2)*(1-beta3)
     1            +bin(lx+1,ly,lz)*beta1*(1-beta2)*(1-beta3)
     1            +bin(lx,ly+1,lz)*(1-beta1)*beta2*(1-beta3)
     1            +bin(lx+1,ly+1,lz)*beta1*beta2*(1-beta3)
     1            +bin(lx,ly,lz+1)*(1-beta1)*(1-beta2)*beta3
     1            +bin(lx+1,ly,lz+1)*beta1*(1-beta2)*beta3
     1            +bin(lx,ly+1,lz+1)*(1-beta1)*beta2*beta3
     1            +bin(lx+1,ly+1,lz+1)*beta1*beta2*beta3
   40 continue
c
      return
      end
c
c
c
      subroutine extrho(bin,xi,zi,nx,nz,bl,xx,zz,mx,mz,err)
      implicit real*8 (a-h,o-z)
      parameter(kk=5,nx0=2*kk,nz0=2*kk,nxp=nx0+1,nzp=nz0+1)
      dimension xi(nx),zi(nz),bin(nx,nz)
      dimension xx(mx),zz(mz),bl(mx,mz)
      dimension xc(nxp),zc(nzp),bm(nxp,nzp)
      nxh=nxp/2
      nzh=nzp/2
      do 5 kz=1,mz
      do 5 kx=1,mx
      call sindex(xi,zi,nx,nz,xx(kx),zz(kz),lx,lz)
      if(lx.le.nxh) then
        mkx1=1
        mkx2=nxp
      else
        if((nx-lx).le.nxh) then
        mkx1=nx-nx0
        mkx2=nx
        else
        mkx1=lx-nxh
        mkx2=lx+nxh
        endif
      endif
      if(lz.le.nzh) then
        mkz1=1
        mkz2=nzp
      else
        if((nz-lz).le.nzh) then
        mkz1=nz-nz0
        mkz2=nz
        else
        mkz1=lz-nzh
        mkz2=lz+nzh
        endif
      endif
      do 4 jz=mkz1,mkz2
      do 4 jx=mkx1,mkx2
      xc(jx-mkx1+1)=xi(jx)
      zc(jz-mkz1+1)=zi(jz)
      bm(jx-mkx1+1,jz-mkz1+1)=bin(jx,jz)
    4 continue
      call polin2(xc,zc,bm,nxp,nzp,xx(kx),zz(kz),aa,err)
      bl(kx,kz)=aa
    5 continue
      return
      end
c
c
c
      SUBROUTINE POLIN2(X1A,X2A,YA,M,N,X1,X2,Y,DY)
      implicit real*8 (a-h,o-z)
      PARAMETER (MMAX=20,NMAX=20)
      DIMENSION X1A(M),X2A(N),YA(M,N),YNTMP(NMAX),YMTMP(MMAX)
      DO 12 J=1,M
        DO 11 K=1,N
          YNTMP(K)=YA(J,K)
   11   CONTINUE
        CALL POLINT(X2A,YNTMP,N,X2,YMTMP(J),DY)
   12 CONTINUE
      CALL POLINT(X1A,YMTMP,M,X1,Y,DY)
      RETURN
      END
c
c
c
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      implicit real*8 (a-h,o-z)
      PARAMETER (NMAX=20) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF(DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
   11 CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
c         IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
   12   CONTINUE
        IF(2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
   13 CONTINUE
      RETURN
      END
c
c
c
      subroutine sindex(xi,zi,nx,nz,xl,zl,lx,lz)
      implicit real*8 (a-h,o-z)
      dimension xi(nx),zi(nz)
      do 10 ix=1,nx-1
      if((xi(ix)-1.d-4).le.xl.and.(xi(ix+1)+1.d-4).ge.xl) then
       lx=ix
      goto 15
      else
      endif
   10 continue
   15 do 20 iz=1,nz-1
      if((zi(iz)-1.d-4).le.zl.and.(zi(iz+1)+1.d-4).ge.zl) then
      lz=iz
      goto 30
      else
      endif
   20 continue
c
   30 continue
      return
      end
c
c
c
      subroutine stepon
c
c     This routine time-advances X's by fourth order in time and second
c     order in space Runge-Kotta differential scheme.
c     note: X is always the up-to-date value while Xm being the
c           intermediate value, and Xdif is increment
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c
c     ml=1
      xwfront=-10.+0.08*time
c     if(xx(x_last).gt.xwfront) goto 100
c     revised by jiliang in Sep.,2019
      do 11 m=1,8
      do 11 jz=z_first-2,z_last+2
      do 11 jy=1,my
      do 11 jx=x_first-2,x_last+2
      xfold(jx,jy,jz,m)=x(jx,jy,jz,m)
   11 continue
      
      tt=time+dt/2
      ddt=dt/2.
c
      call right(tt,ddt,1)
      do 1 m=1,8
      do 1 jz=z_first-2,z_last+2
      do 1 jy=1,my
      do 1 jx=x_first-2,x_last+2
c       xfold(jx,jy,jz,m)=x(jx,jy,jz,m)
        xm(jx,jy,jz,m)=xfold(jx,jy,jz,m)+xdif(jx,jy,jz,m)*dt/6.
        x(jx,jy,jz,m)=xfold(jx,jy,jz,m)+xdif(jx,jy,jz,m)*dt/2.
    1 continue
c     call recrd2
c     print*,nst
c     stop
c
      tt=time+dt/2
      ddt=dt/2.
      call right(tt,ddt,1)
      do 2 m=1,8
      do 2 jz=z_first-2,z_last+2
      do 2 jy=1,my
      do 2 jx=x_first-2,x_last+2
        xm(jx,jy,jz,m)=xm(jx,jy,jz,m)+xdif(jx,jy,jz,m)*dt/3.
        x(jx,jy,jz,m)=xfold(jx,jy,jz,m)+xdif(jx,jy,jz,m)*dt/2.
    2 continue
c
      tt=time+dt
      ddt=dt
      call right(tt,ddt,1)
      do 3 m=1,8
      do 3 jz=z_first-2,z_last+2
      do 3 jy=1,my
      do 3 jx=x_first-2,x_last+2
        xm(jx,jy,jz,m)=xm(jx,jy,jz,m)+xdif(jx,jy,jz,m)*dt/3.
        x(jx,jy,jz,m)=xfold(jx,jy,jz,m)+xdif(jx,jy,jz,m)*dt
    3 continue
c
      tt=time+dt
      ddt=dt/6.
      call right(tt,ddt,2)
      do 4 m=1,8
      do 4 jz=z_first-2,z_last+2
      do 4 jy=1,my
      do 4 jx=x_first-2,x_last+2
        x(jx,jy,jz,m)=xm(jx,jy,jz,m)+xdif(jx,jy,jz,m)*dt/6.
    4 continue
c
c     if(mod(nstep,4).eq.0) then
      do 5 jz=z_first-2,z_last+2
      do 5 jy=1,my
      do 5 jx=x_first-2,x_last+2
        xm(jx,jy,jz,1)=x(jx,jy,jz,1)-rhoi(jx,jy,jz)
        xm(jx,jy,jz,2)=x(jx,jy,jz,2)
        xm(jx,jy,jz,3)=x(jx,jy,jz,3)
        xm(jx,jy,jz,4)=x(jx,jy,jz,4)
        xm(jx,jy,jz,5)=x(jx,jy,jz,5)-bxi(jx,jy,jz)
        xm(jx,jy,jz,6)=x(jx,jy,jz,6)-byi(jx,jy,jz)
        xm(jx,jy,jz,7)=x(jx,jy,jz,7)-bzi(jx,jy,jz)
        xm(jx,jy,jz,8)=x(jx,jy,jz,8)-pri(jx,jy,jz)
    5 continue
c
      caf0=0.d0
c      caf0=0.95d0
      call avrgtt(tt,1,8,caf0,2)
c     call avrgrtt(1,8,caf0,1)
c     call avrgtt(tt,1,8,caf0,1)
c     call avrgtt(tt,8,8,caf0,3)
      do 6 jz=z_first-2,z_last+2
      do 6 jy=1,my
      do 6 jx=x_first-2,x_last+2
        x(jx,jy,jz,1)=xm(jx,jy,jz,1)+rhoi(jx,jy,jz)
        x(jx,jy,jz,2)=xm(jx,jy,jz,2)
        x(jx,jy,jz,3)=xm(jx,jy,jz,3)
        x(jx,jy,jz,4)=xm(jx,jy,jz,4)
        x(jx,jy,jz,5)=xm(jx,jy,jz,5)+bxi(jx,jy,jz)
        x(jx,jy,jz,6)=xm(jx,jy,jz,6)+byi(jx,jy,jz)
        x(jx,jy,jz,7)=xm(jx,jy,jz,7)+bzi(jx,jy,jz)
        x(jx,jy,jz,8)=xm(jx,jy,jz,8)+pri(jx,jy,jz)
    6 continue
c     endif
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c
c     call message4(x,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)

      call current
c
  100 continue
      if(mod(nstep,10).eq.0) then
c     if(myleft.eq.mpi_proc_null)then
c     open(unit=17,file='velo.dat',status='unknown',form='formatted')
c     write(17,15)(x(jx,my/2,21,2),jx=2,17),zz(21),time
c  15 format(8(1x,e13.5))
c     endif
      call diagn
c     call energy(time)
      endif
c     if(mod(nstep,50).eq.0.and.time.gt.0.) call satellite
      return
      end
c
c
c
      subroutine diagn
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      dimension w0(mxpr,mzpr)
c  define statement functions
c  d2fc= d2 f / dx2   with central difference
      d2fc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*((fp-f0)/(xp1-x0)-(f0-fm)/(x0-xm1))/(xp1-xm1)
      call current
      do 10 jz=2,mzpr
      do 10 jx=2,mxpr
        w0(jx,jz)=cur(jx,my/2+1,jz,1)
   10 continue
      call funmax(w0,cxmax0,cxmin0,xx,zz,mxpr,mzpr)
      do 20 jz=2,mzpr
      do 20 jx=2,mxpr
        w0(jx,jz)=cur(jx,my/2+1,jz,2)
   20 continue
      call funmax(w0,cymax0,cymin0,xx,zz,mxpr,mzpr)
      do 30 jz=2,mzpr
      do 30 jx=2,mxpr
        w0(jx,jz)=cur(jx,my/2+1,jz,3)
   30 continue
      call funmax(w0,czmax0,czmin0,xx,zz,mxpr,mzpr)
      call mpi_reduce(cxmax0,cxmax,1,mpi_double_precision,mpi_max,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(cymax0,cymax,1,mpi_double_precision,mpi_max,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(czmax0,czmax,1,mpi_double_precision,mpi_max,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(cxmin0,cxmin,1,mpi_double_precision,mpi_min,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(cymin0,cymin,1,mpi_double_precision,mpi_min,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(czmin0,czmin,1,mpi_double_precision,mpi_min,0,
     1     mpi_comm_world,ierr)
      if(nrank.eq.0)then
      open(unit=9,file='current.dat',status='unknown',form='formatted')
      write(9,15)cxmax,cxmin,cymax,cymin,czmax,czmin,time
   15 format(7(1x,e13.4))
      endif
      if(abs(cymax).ge.1.d7) stop
c
      return
      end
c
c
c
      subroutine diagn_res
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      dimension w0(mxpr,mzpr)
c  define statement functions
c  d2fc= d2 f / dx2   with central difference
      d2fc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*((fp-f0)/(xp1-x0)-(f0-fm)/(x0-xm1))/(xp1-xm1)
      call current
      do 10 jz=2,mzpr
      do 10 jx=2,mxpr
        w0(jx,jz)=xresd(jx,my/2+1,jz,1)
   10 continue
      call funmax(w0,cxmax0,cxmin0,xx,zz,mxpr,mzpr)
      do 20 jz=2,mzpr
      do 20 jx=2,mxpr
        w0(jx,jz)=xresd(jx,my/2+1,jz,2)
   20 continue
      call funmax(w0,cymax0,cymin0,xx,zz,mxpr,mzpr)
      do 30 jz=2,mzpr
      do 30 jx=2,mxpr
        w0(jx,jz)=xresd(jx,my/2+1,jz,3)
   30 continue
      call funmax(w0,czmax0,czmin0,xx,zz,mxpr,mzpr)
      call mpi_reduce(cxmax0,cxmax,1,mpi_double_precision,mpi_max,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(cymax0,cymax,1,mpi_double_precision,mpi_max,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(czmax0,czmax,1,mpi_double_precision,mpi_max,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(cxmin0,cxmin,1,mpi_double_precision,mpi_min,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(cymin0,cymin,1,mpi_double_precision,mpi_min,0,
     1     mpi_comm_world,ierr)
      call mpi_reduce(czmin0,czmin,1,mpi_double_precision,mpi_min,0,
     1     mpi_comm_world,ierr)
      if(nrank.eq.0)then
      open(unit=9,file='residue.dat',status='unknown',form='formatted')
      write(9,15)cxmax,cxmin,cymax,cymin,czmax,czmin,time
   15 format(7(1x,e13.4))
      close(9)
      endif
c      if(nrank.eq.0)then
c      print*,(xresd(jx,my/2+1,zlast/2,2),jx=xfirst,xlast)
c      close(9)
c      endif
c      if(nrank.eq.120)then
c      print*,(xresd(jx,my/2+1,zlast/2,2),jx=xfirst,xlast)
c      close(9)
c      endif
c
      return
      end
c
c
c
      subroutine right(t,ddt,km)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c  define statement functions
c  d2fc= d2 f / dx2  with three-point second-order(uniform grids)/first-order(non-uniform grids) accuracz central difference
c      d2fc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*((fp-f0)/(xp1-x0)-(f0-fm)/(x0-xm1))/(xp1-xm1)
      
c  d2fc= d2 f / dx2  with five-point forth-order(uniform grids)/third-order(non-uniform grids) accuracz central difference
      d2fc(fm2,fm,f0,fp,fp2,a,b,c,d)=
     1 a*(fp-f0)+b*(f0-fm)+c*(fp2-f0)+d*(f0-fm2)
      
c  d1fc= d f / dx  with five-point fourth-order(uniform grids) accuracz central difference
      d1fc(fm2,fm,f0,fp,fp2,a,b,c,d)=
     1   a*(fp-f0)+b*(f0-fm)+c*(fp2-f0)+d*(f0-fm2)
      
c  d1fc= d f / dx  with three-point second-order(uniform grids) accuracz central difference
c     d1fc(fm,f0,fp,xm1,x0,xp1)=
c    1   ((xm1-x0)/(xp1-x0)*(fp-f0)
c    1    -(xp1-x0)/(xm1-x0)*(fm-f0))/(xm1-xp1)
      
c  d1fm= d f / dx  with  one-sided difference involving -2 -1 and 0
c  points
c      d1fm(fm2,fm1,f0,xm2,xm1,x0)=
c     1   ( (xm2-x0)/(xm1-x0)*(fm1-f0)
c     1    -(xm1-x0)/(xm2-x0)*(fm2-f0) ) / (xm2-xm1)
      
c  d1fp= d f / dx  with one-sided difference involving 0  1 and 2
c  points
c      d1fp(f0,fp1,fp2,x0,xp1,xp2)=
c     1   ( (xp2-x0)/(xp1-x0)*(fp1-f0)
c     1    -(xp1-x0)/(xp2-x0)*(fp2-f0) ) / (xp2-xp1)
c
      call current
c
      do 2 m=1,4
      call flux(m)
      do 1 jz=z_first,z_last
      do 1 jy=3,my-2
      do 1 jx=x_first,x_last
      xdif(jx,jy,jz,m)=-d1fc(fs(jx-2,jy,jz),fs(jx-1,jy,jz),fs(jx,jy,jz)
     1  ,fs(jx+1,jy,jz),fs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     2                 -d1fc(gs(jx,jy-2,jz),gs(jx,jy-1,jz),gs(jx,jy,jz)
     3  ,gs(jx,jy+1,jz),gs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
     4                 -d1fc(hs(jx,jy,jz-2),hs(jx,jy,jz-1),hs(jx,jy,jz)
     5  ,hs(jx,jy,jz+1),hs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
    1 continue
    2 continue
c
      do 102 jz=z_first-2,z_last+2
      do 102 jy=1,my
      do 102 jx=x_first-2,x_last+2
      fs(jx,jy,jz)=x(jx,jy,jz,2)/x(jx,jy,jz,1)
      gs(jx,jy,jz)=x(jx,jy,jz,3)/x(jx,jy,jz,1)
      hs(jx,jy,jz)=x(jx,jy,jz,4)/x(jx,jy,jz,1)
c
      xdif(jx,jy,jz,2)=xdif(jx,jy,jz,2)-xresd(jx,jy,jz,1)
     1   +cur(jx,jy,jz,2)*x(jx,jy,jz,7)-cur(jx,jy,jz,3)*x(jx,jy,jz,6)
      xdif(jx,jy,jz,3)=xdif(jx,jy,jz,3)-xresd(jx,jy,jz,2)
     1   +cur(jx,jy,jz,3)*x(jx,jy,jz,5)-cur(jx,jy,jz,1)*x(jx,jy,jz,7)
      xdif(jx,jy,jz,4)=xdif(jx,jy,jz,4)-xresd(jx,jy,jz,3)
     1   +cur(jx,jy,jz,1)*x(jx,jy,jz,6)-cur(jx,jy,jz,2)*x(jx,jy,jz,5)
c
c     revised by ji liang in Sep.,2019
      if (totalE) then
c      if((nrank.eq.0) .and. (jx.eq.7) .and. (jy.eq.7) .and. (jz.eq.7))
c     1 write(*,*)'totalE true1'
c for Etotal=-u1 curl B +eta*J1
      efld(jx,jy,jz,1)=(-x(jx,jy,jz,3)*x(jx,jy,jz,7)+x(jx,jy,jz,4)
     1    *x(jx,jy,jz,6))/x(jx,jy,jz,1)+etan(jx,jy,jz)*cur(jx,jy,jz,1)
      efld(jx,jy,jz,2)=(-x(jx,jy,jz,4)*x(jx,jy,jz,5)+x(jx,jy,jz,2)
     1    *x(jx,jy,jz,7))/x(jx,jy,jz,1)+etan(jx,jy,jz)*cur(jx,jy,jz,2)
      efld(jx,jy,jz,3)=(-x(jx,jy,jz,2)*x(jx,jy,jz,6)+x(jx,jy,jz,3)
     1    *x(jx,jy,jz,5))/x(jx,jy,jz,1)+etan(jx,jy,jz)*cur(jx,jy,jz,3)
      else
      if((nrank.eq.0) .and. (jx.eq.7) .and. (jy.eq.7) .and. (jz.eq.7))
     1 write(*,*)'totalE false1'
c for Epart=-u1 clurl B1 + eta*J1
      efld(jx,jy,jz,1)=(-x(jx,jy,jz,3)*(x(jx,jy,jz,7)-bzi(jx,jy,jz))
     1                  +x(jx,jy,jz,4)*(x(jx,jy,jz,6)-byi(jx,jy,jz)))
     1                 /x(jx,jy,jz,1)+etan(jx,jy,jz)*cur(jx,jy,jz,1)
      efld(jx,jy,jz,2)=(-x(jx,jy,jz,4)*(x(jx,jy,jz,5)-bxi(jx,jy,jz))
     1                  +x(jx,jy,jz,2)*(x(jx,jy,jz,7)-bzi(jx,jy,jz)))
     1                 /x(jx,jy,jz,1)+etan(jx,jy,jz)*cur(jx,jy,jz,2)
      efld(jx,jy,jz,3)=(-x(jx,jy,jz,2)*(x(jx,jy,jz,6)-byi(jx,jy,jz))
     1                  +x(jx,jy,jz,3)*(x(jx,jy,jz,5)-bxi(jx,jy,jz)))
     1                 /x(jx,jy,jz,1)+etan(jx,jy,jz)*cur(jx,jy,jz,3)
      endif
  102 continue
c
      call avrgE(t,1,3,4)
c
      if (hall) then
      do 122 jz=z_first-2,z_last+2
      do 122 jy=1,my
      do 122 jx=x_first-2,x_last+2
c     if (t.le.1000) then
c     cht(jx,jy,jz)=0.
c     else
      cht(jx,jy,jz)=0.5*chall*(tanh((rr0(jx,jy,jz)-1.1)/8.)
     1           -tanh((rr0(jx,jy,jz)-2.)/10.))
c    1   *(tanh((t-1000)/100.))*exp(-16.*zz(jz)**2)
c    1   *exp(-(xx(jx)+10)**2/4.)*exp(-(yy(jy)/0.01)**2)
c     endif
  122 continue
      do 132 jz=z_first,z_last
      do 132 jy=3,my-2
      do 132 jx=x_first,x_last
      xdh(jx,jy,jz,1)=cht(jx,jy,jz)*(cur(jx,jy,jz,2)*x(jx,jy,jz,7)
     1   -cur(jx,jy,jz,3)*x(jx,jy,jz,6)
     1   -d1fc(x(jx-2,jy,jz,8),x(jx-1,jy,jz,8),x(jx,jy,jz,8)
     2   ,x(jx+1,jy,jz,8),x(jx+2,jy,jz,8),ax1(jx),bx1(jx)
     3   ,cx1(jx),dx1(jx)))/x(jx,jy,jz,1)
      xdh(jx,jy,jz,2)=cht(jx,jy,jz)*(cur(jx,jy,jz,3)*x(jx,jy,jz,5)
     1   -cur(jx,jy,jz,1)*x(jx,jy,jz,7)
     1   -d1fc(x(jx,jy-2,jz,8),x(jx,jy-1,jz,8),x(jx,jy,jz,8)
     2   ,x(jx,jy+1,jz,8),x(jx,jy+2,jz,8),ay1(jy),by1(jy)
     3   ,cy1(jy),dy1(jy)))/x(jx,jy,jz,1)
      xdh(jx,jy,jz,3)=cht(jx,jy,jz)*(cur(jx,jy,jz,1)*x(jx,jy,jz,6)
     1   -cur(jx,jy,jz,2)*x(jx,jy,jz,5)
     1   -d1fc(x(jx,jy,jz-2,8),x(jx,jy,jz-1,8),x(jx,jy,jz,8)
     2   ,x(jx,jy,jz+1,8),x(jx,jy,jz+2,8),az1(jz),bz1(jz)
     3   ,cz1(jz),dz1(jz)))/x(jx,jy,jz,1)
  132 continue
      call bndryh
      call avrgh(1,3,5)
c      call avrgh(1,3,2)
      do 112 jz=z_first-2,z_last+2
      do 112 jy=1,my
      do 112 jx=x_first-2,x_last+2
      efld(jx,jy,jz,1)=efld(jx,jy,jz,1)+xdh(jx,jy,jz,1)
      efld(jx,jy,jz,2)=efld(jx,jy,jz,2)+xdh(jx,jy,jz,2)
      efld(jx,jy,jz,3)=efld(jx,jy,jz,3)+xdh(jx,jy,jz,3)
  112 continue
      endif
c
      do 101 jz=z_first,z_last
      do 101 jy=3,my-2
      do 101 jx=x_first,x_last
      xdif(jx,jy,jz,5)=-d1fc(efld(jx,jy-2,jz,3),efld(jx,jy-1,jz,3)
     1    ,efld(jx,jy,jz,3),efld(jx,jy+1,jz,3),efld(jx,jy+2,jz,3)
     1    ,ay1(jy),by1(jy),cy1(jy),dy1(jy))+d1fc(efld(jx,jy,jz-2,2)
     1    ,efld(jx,jy,jz-1,2),efld(jx,jy,jz,2),efld(jx,jy,jz+1,2)
     1    ,efld(jx,jy,jz+2,2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
      xdif(jx,jy,jz,6)=-d1fc(efld(jx,jy,jz-2,1),efld(jx,jy,jz-1,1)
     1    ,efld(jx,jy,jz,1),efld(jx,jy,jz+1,1),efld(jx,jy,jz+2,1)
     1    ,az1(jz),bz1(jz),cz1(jz),dz1(jz))+d1fc(efld(jx-2,jy,jz,3)
     1    ,efld(jx-1,jy,jz,3),efld(jx,jy,jz,3),efld(jx+1,jy,jz,3)
     1    ,efld(jx+2,jy,jz,3),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
      xdif(jx,jy,jz,7)=-d1fc(efld(jx-2,jy,jz,2),efld(jx-1,jy,jz,2)
     1    ,efld(jx,jy,jz,2),efld(jx+1,jy,jz,2),efld(jx+2,jy,jz,2)
     1    ,ax1(jx),bx1(jx),cx1(jx),dx1(jx))+d1fc(efld(jx,jy-2,jz,1)
     1    ,efld(jx,jy-1,jz,1),efld(jx,jy,jz,1),efld(jx,jy+1,jz,1)
     1    ,efld(jx,jy+2,jz,1),ay1(jy),by1(jy),cy1(jy),dy1(jy))
      if (totalE) then
c      if((nrank.eq.0) .and. (jx.eq.7) .and. (jy.eq.7) .and. (jz.eq.7))
c     1 write(*,*)'totalE true2'
c dB/dt=-curl(Etotal)
      xdif(jx,jy,jz,5)=xdif(jx,jy,jz,5)
      xdif(jx,jy,jz,6)=xdif(jx,jy,jz,6)
      xdif(jx,jy,jz,7)=xdif(jx,jy,jz,7)
      else
      if((nrank.eq.0) .and. (jx.eq.7) .and. (jy.eq.7) .and. (jz.eq.7))
     1 write(*,*)'totalE false2'
c dB/dt=-curl(Epart)+curl(u1 curl B0)
      xdif(jx,jy,jz,5)=xdif(jx,jy,jz,5)+(
     1-bxi(jx,jy,jz)*
     1 (d1fc(fs(jx-2,jy,jz),fs(jx-1,jy,jz),fs(jx,jy,jz)
     1 ,fs(jx+1,jy,jz),fs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     1 +d1fc(gs(jx,jy-2,jz),gs(jx,jy-1,jz),gs(jx,jy,jz)
     1 ,gs(jx,jy+1,jz),gs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
     1 +d1fc(hs(jx,jy,jz-2),hs(jx,jy,jz-1),hs(jx,jy,jz)
     1 ,hs(jx,jy,jz+1),hs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
     2+(bxi(jx,jy,jz)*d1fc(fs(jx-2,jy,jz),fs(jx-1,jy,jz),fs(jx,jy,jz)
     2 ,fs(jx+1,jy,jz),fs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))+
     2  byi(jx,jy,jz)*d1fc(fs(jx,jy-2,jz),fs(jx,jy-1,jz),fs(jx,jy,jz)
     2 ,fs(jx,jy+1,jz),fs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))+
     2  bzi(jx,jy,jz)*d1fc(fs(jx,jy,jz-2),fs(jx,jy,jz-1),fs(jx,jy,jz)
     2 ,fs(jx,jy,jz+1),fs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
c     3-(fs(jx,jy,jz)*d1fc(bxi(jx-2,jy,jz),bxi(jx-1,jy,jz),bxi(jx,jy,jz),
c     3 bxi(jx+1,jy,jz),bxi(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
c     3 +gs(jx,jy,jz)*d1fc(bxi(jx,jy-2,jz),bxi(jx,jy-1,jz),bxi(jx,jy,jz),
c     3 bxi(jx,jy+1,jz),bxi(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
c     3 +hs(jx,jy,jz)*d1fc(bxi(jx,jy,jz-2),bxi(jx,jy,jz-1),bxi(jx,jy,jz),
c     3 bxi(jx,jy,jz+1),bxi(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
     3-(fs(jx,jy,jz)*(-3.*zz(jz)*(rr0(jx,jy,jz)**2-5.*xx(jx)**2)
     3                /rr0(jx,jy,jz)**7)
     3 +gs(jx,jy,jz)*(3.*xx(jx)*zz(jz)*5.*yy(jy)
     3                /rr0(jx,jy,jz)**7)
     3 +hs(jx,jy,jz)*(-3.*xx(jx)*(rr0(jx,jy,jz)**2-5.*zz(jz)**2)
     3                /rr0(jx,jy,jz)**7))
     4 )
      xdif(jx,jy,jz,6)=xdif(jx,jy,jz,6)+(
     1-byi(jx,jy,jz)*
     1 (d1fc(fs(jx-2,jy,jz),fs(jx-1,jy,jz),fs(jx,jy,jz)
     1 ,fs(jx+1,jy,jz),fs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     1 +d1fc(gs(jx,jy-2,jz),gs(jx,jy-1,jz),gs(jx,jy,jz)
     1 ,gs(jx,jy+1,jz),gs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
     1 +d1fc(hs(jx,jy,jz-2),hs(jx,jy,jz-1),hs(jx,jy,jz)
     1 ,hs(jx,jy,jz+1),hs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
     2+(bxi(jx,jy,jz)*d1fc(gs(jx-2,jy,jz),gs(jx-1,jy,jz),gs(jx,jy,jz)
     2 ,gs(jx+1,jy,jz),gs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))+
     2  byi(jx,jy,jz)*d1fc(gs(jx,jy-2,jz),gs(jx,jy-1,jz),gs(jx,jy,jz)
     2 ,gs(jx,jy+1,jz),gs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))+
     2  bzi(jx,jy,jz)*d1fc(gs(jx,jy,jz-2),gs(jx,jy,jz-1),gs(jx,jy,jz)
     2 ,gs(jx,jy,jz+1),gs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
c     3-(fs(jx,jy,jz)*d1fc(byi(jx-2,jy,jz),byi(jx-1,jy,jz),byi(jx,jy,jz),
c     3 byi(jx+1,jy,jz),byi(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
c     3 +gs(jx,jy,jz)*d1fc(byi(jx,jy-2,jz),byi(jx,jy-1,jz),byi(jx,jy,jz),
c     3 byi(jx,jy+1,jz),byi(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
c     3 +hs(jx,jy,jz)*d1fc(byi(jx,jy,jz-2),byi(jx,jy,jz-1),byi(jx,jy,jz),
c     3 byi(jx,jy,jz+1),byi(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
     3-(fs(jx,jy,jz)*(3.*yy(jy)*zz(jz)*5.*xx(jx)
     3                /rr0(jx,jy,jz)**7)
     3 +gs(jx,jy,jz)*(-3.*zz(jz)*(rr0(jx,jy,jz)**2-5.*yy(jy)**2)
     3                /rr0(jx,jy,jz)**7)
     3 +hs(jx,jy,jz)*(-3.*yy(jy)*(rr0(jx,jy,jz)**2-5.*zz(jz)**2)
     3                /rr0(jx,jy,jz)**7))
     4 )
      xdif(jx,jy,jz,7)=xdif(jx,jy,jz,7)+(
     1-bzi(jx,jy,jz)*
     1 (d1fc(fs(jx-2,jy,jz),fs(jx-1,jy,jz),fs(jx,jy,jz)
     1 ,fs(jx+1,jy,jz),fs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     1 +d1fc(gs(jx,jy-2,jz),gs(jx,jy-1,jz),gs(jx,jy,jz)
     1 ,gs(jx,jy+1,jz),gs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
     1 +d1fc(hs(jx,jy,jz-2),hs(jx,jy,jz-1),hs(jx,jy,jz)
     1 ,hs(jx,jy,jz+1),hs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
     2+(bxi(jx,jy,jz)*d1fc(hs(jx-2,jy,jz),hs(jx-1,jy,jz),hs(jx,jy,jz)
     2 ,hs(jx+1,jy,jz),hs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))+
     2  byi(jx,jy,jz)*d1fc(hs(jx,jy-2,jz),hs(jx,jy-1,jz),hs(jx,jy,jz)
     2 ,hs(jx,jy+1,jz),hs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))+
     2  bzi(jx,jy,jz)*d1fc(hs(jx,jy,jz-2),hs(jx,jy,jz-1),hs(jx,jy,jz)
     2 ,hs(jx,jy,jz+1),hs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
c     3-(fs(jx,jy,jz)*d1fc(bzi(jx-2,jy,jz),bzi(jx-1,jy,jz),bzi(jx,jy,jz),
c     3 bzi(jx+1,jy,jz),bzi(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
c     3 +gs(jx,jy,jz)*d1fc(bzi(jx,jy-2,jz),bzi(jx,jy-1,jz),bzi(jx,jy,jz),
c     3 bzi(jx,jy+1,jz),bzi(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
c     3 +hs(jx,jy,jz)*d1fc(bzi(jx,jy,jz-2),bzi(jx,jy,jz-1),bzi(jx,jy,jz),
c     3 bzi(jx,jy,jz+1),bzi(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
     3-(fs(jx,jy,jz)*((2.*xx(jx)*rr0(jx,jy,jz)**2-5.*xx(jx)*
     3                (xx(jx)**2+yy(jy)**2-2.*zz(jz)**2))
     3                /rr0(jx,jy,jz)**7)
     3 +gs(jx,jy,jz)*((2.*yy(jy)*rr0(jx,jy,jz)**2-5.*yy(jy)*
     3                (xx(jx)**2+yy(jy)**2-2.*zz(jz)**2))
     3                /rr0(jx,jy,jz)**7)
     3 +hs(jx,jy,jz)*((-4.*zz(jz)*rr0(jx,jy,jz)**2-5.*zz(jz)*
     3                (xx(jx)**2+yy(jy)**2-2.*zz(jz)**2))
     3                /rr0(jx,jy,jz)**7))
     4 )
      endif
c
      xdif(jx,jy,jz,8)=-d1fc(x(jx-2,jy,jz,8),x(jx-1,jy,jz,8)
     1    ,x(jx,jy,jz,8),x(jx+1,jy,jz,8),x(jx+2,jy,jz,8)
     1    ,ax1(jx),bx1(jx),cx1(jx),dx1(jx))*fs(jx,jy,jz)
     2    -d1fc(x(jx,jy-2,jz,8),x(jx,jy-1,jz,8),x(jx,jy,jz,8)
     2    ,x(jx,jy+1,jz,8),x(jx,jy+2,jz,8),ay1(jy),by1(jy)
     2    ,cy1(jy),dy1(jy))*gs(jx,jy,jz)
     3    -d1fc(x(jx,jy,jz-2,8),x(jx,jy,jz-1,8),x(jx,jy,jz,8)
     3    ,x(jx,jy,jz+1,8),x(jx,jy,jz+2,8),az1(jz),bz1(jz)
     3    ,cz1(jz),dz1(jz))*hs(jx,jy,jz)
     4    +(-d1fc(fs(jx-2,jy,jz),fs(jx-1,jy,jz),fs(jx,jy,jz)
     4  ,fs(jx+1,jy,jz),fs(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     4      -d1fc(gs(jx,jy-2,jz),gs(jx,jy-1,jz),gs(jx,jy,jz)
     4  ,gs(jx,jy+1,jz),gs(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
     4      -d1fc(hs(jx,jy,jz-2),hs(jx,jy,jz-1),hs(jx,jy,jz)
     4  ,hs(jx,jy,jz+1),hs(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz)))
     5    *gamma*x(jx,jy,jz,8)+(gamma-1.)*etan(jx,jy,jz)
     5    *(cur(jx,jy,jz,1)**2+cur(jx,jy,jz,2)**2+cur(jx,jy,jz,3)**2)
  101 continue
c
      coef=0.05
      if(damp) then
      do 5 m=2,4
      do 5 jz=z_first,z_last
      do 5 jy=2,my-1
      do 5 jx=x_first,x_last
      xdif(jx,jy,jz,m)=xdif(jx,jy,jz,m)-coef*x(jx,jy,jz,m)
    5 continue
      do 6 jz=z_first,z_last
      do 6 jy=2,my-1
      do 6 jx=x_first,x_last
      xdif(jx,jy,jz,8)=xdif(jx,jy,jz,8)-coef*x(jx,jy,jz,1)*
     1     (fs(jx,jy,jz)**2+gs(jx,jy,jz)**2+hs(jx,jy,jz)**2)
    6 continue
      endif
c
      if(viscos) then
      do 7 jz=z_first,z_last
      do 7 jy=3,my-2
      do 7 jx=x_first,x_last
c      af=         (d2fc(fs(jx-1,jy,jz),fs(jx,jy,jz),fs(jx+1,jy,jz)
c     1               ,xx(jx-1),xx(jx),xx(jx+1))
c     2            +d2fc(fs(jx,jy-1,jz),fs(jx,jy,jz),fs(jx,jy+1,jz)
c     3               ,yy(jy-1),yy(jy),yy(jy+1))    
c     4            +d2fc(fs(jx,jy,jz-1),fs(jx,jy,jz),fs(jx,jy,jz+1)
c     5               ,zz(jz-1),zz(jz),zz(jz+1)))
c      ag=         (d2fc(gs(jx-1,jy,jz),gs(jx,jy,jz),gs(jx+1,jy,jz)
c     1               ,xx(jx-1),xx(jx),xx(jx+1))
c     2            +d2fc(gs(jx,jy-1,jz),gs(jx,jy,jz),gs(jx,jy+1,jz)
c     3               ,yy(jy-1),yy(jy),yy(jy+1))    
c     4            +d2fc(gs(jx,jy,jz-1),gs(jx,jy,jz),gs(jx,jy,jz+1)
c     5               ,zz(jz-1),zz(jz),zz(jz+1)))
c      ah=         (d2fc(hs(jx-1,jy,jz),hs(jx,jy,jz),hs(jx+1,jy,jz)
c     1               ,xx(jx-1),xx(jx),xx(jx+1))
c     2            +d2fc(hs(jx,jy-1,jz),hs(jx,jy,jz),hs(jx,jy+1,jz)
c     3               ,yy(jy-1),yy(jy),yy(jy+1))    
c     4            +d2fc(hs(jx,jy,jz-1),hs(jx,jy,jz),hs(jx,jy,jz+1)
c     5               ,zz(jz-1),zz(jz),zz(jz+1)))
c     xdif(jx,jy,jz,2)=xdif(jx,jy,jz,2)+fmu0*af
c     xdif(jx,jy,jz,3)=xdif(jx,jy,jz,3)+fmu0*ag
c     xdif(jx,jy,jz,4)=xdif(jx,jy,jz,4)+fmu0*ah
c     xdif(jx,jy,jz,8)=xdif(jx,jy,jz,8)+fmu0*
c    1    (fs(jx,jy,jz)*af+gs(jx,jy,jz)*ag+hs(jx,jy,jz)*ah)
      af=(d2fc(fs(jx-2,jy,jz),fs(jx-1,jy,jz),fs(jx,jy,jz),fs(jx+1,jy,jz)
     1        ,fs(jx+2,jy,jz),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
     2   +d2fc(fs(jx,jy-2,jz),fs(jx,jy-1,jz),fs(jx,jy,jz),fs(jx,jy+1,jz)
     3        ,fs(jx,jy+2,jz),ay2(jy),by2(jy),cy2(jy),dy2(jy))
     4   +d2fc(fs(jx,jy,jz-2),fs(jx,jy,jz-1),fs(jx,jy,jz),fs(jx,jy,jz+1)
     5        ,fs(jx,jy,jz+2),az2(jz),bz2(jz),cz2(jz),dz2(jz)))
      ag=(d2fc(gs(jx-2,jy,jz),gs(jx-1,jy,jz),gs(jx,jy,jz),gs(jx+1,jy,jz)
     1        ,gs(jx+2,jy,jz),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
     2   +d2fc(gs(jx,jy-2,jz),gs(jx,jy-1,jz),gs(jx,jy,jz),gs(jx,jy+1,jz)
     3        ,gs(jx,jy+2,jz),ay2(jy),by2(jy),cy2(jy),dy2(jy))
     4   +d2fc(gs(jx,jy,jz-2),gs(jx,jy,jz-1),gs(jx,jy,jz),gs(jx,jy,jz+1)
     5        ,gs(jx,jy,jz+2),az2(jz),bz2(jz),cz2(jz),dz2(jz)))
      ah=(d2fc(hs(jx-2,jy,jz),hs(jx-1,jy,jz),hs(jx,jy,jz),hs(jx+1,jy,jz)
     1        ,hs(jx+2,jy,jz),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
     2   +d2fc(hs(jx,jy-2,jz),hs(jx,jy-1,jz),hs(jx,jy,jz),hs(jx,jy+1,jz)
     3        ,hs(jx,jy+2,jz),ay2(jy),by2(jy),cy2(jy),dy2(jy))
     4   +d2fc(hs(jx,jy,jz-2),hs(jx,jy,jz-1),hs(jx,jy,jz),hs(jx,jy,jz+1)
     5        ,hs(jx,jy,jz+2),az2(jz),bz2(jz),cz2(jz),dz2(jz)))
      xdif(jx,jy,jz,2)=xdif(jx,jy,jz,2)+fmu(jx,jy,jz)*af
      xdif(jx,jy,jz,3)=xdif(jx,jy,jz,3)+fmu(jx,jy,jz)*ag
      xdif(jx,jy,jz,4)=xdif(jx,jy,jz,4)+fmu(jx,jy,jz)*ah
      xdif(jx,jy,jz,8)=xdif(jx,jy,jz,8)+fmu(jx,jy,jz)*
     1    (fs(jx,jy,jz)*af+gs(jx,jy,jz)*ag+hs(jx,jy,jz)*ah)
    7 continue
      endif
c
      if(diffu) then
      do 8 jz=z_first,z_last
      do 8 jy=3,my-2
      do 8 jx=x_first,x_last
      xdif(jx,jy,jz,1)=xdif(jx,jy,jz,1)+frho(jx,jy,jz)
     1  *(d2fc(x(jx-2,jy,jz,1),x(jx-1,jy,jz,1),x(jx,jy,jz,1)
     1        ,x(jx+1,jy,jz,1),x(jx+2,jy,jz,1)
     1        ,ax2(jx),bx2(jx),cx2(jx),dx2(jx))
     2   +d2fc(x(jx,jy-2,jz,1),x(jx,jy-1,jz,1),x(jx,jy,jz,1)
     2        ,x(jx,jy+1,jz,1),x(jx,jy+2,jz,1)
     2        ,ay2(jy),by2(jy),cy2(jy),dy2(jy))
     3   +d2fc(x(jx,jy,jz-2,1),x(jx,jy,jz-1,1),x(jx,jy,jz,1)
     3        ,x(jx,jy,jz+1,1),x(jx,jy,jz+2,1)
     3        ,az2(jz),bz2(jz),cz2(jz),dz2(jz)))
      xdif(jx,jy,jz,8)=xdif(jx,jy,jz,8)+fpr(jx,jy,jz)
     1  *(d2fc(x(jx-2,jy,jz,8),x(jx-1,jy,jz,8),x(jx,jy,jz,8)
     1        ,x(jx+1,jy,jz,8),x(jx+2,jy,jz,8)
     1        ,ax2(jx),bx2(jx),cx2(jx),dx2(jx))
     2   +d2fc(x(jx,jy-2,jz,8),x(jx,jy-1,jz,8),x(jx,jy,jz,8)
     2        ,x(jx,jy+1,jz,8),x(jx,jy+2,jz,8)
     2        ,ay2(jy),by2(jy),cy2(jy),dy2(jy))
     3   +d2fc(x(jx,jy,jz-2,8),x(jx,jy,jz-1,8),x(jx,jy,jz,8)
     3        ,x(jx,jy,jz+1,8),x(jx,jy,jz+2,8)
     3        ,az2(jz),bz2(jz),cz2(jz),dz2(jz)))
    8 continue
      endif
c
      call bndry(t,ddt,km)
c     if(t.ne.0.) call smthxy(t,ddt,km,1,8,1)
      caf0=0.d0
c      caf0=0.6d0
      if(t.ne.0.) call avrg(t,ddt,km,1,8,caf0,2)
c     if(t.ne.0.) call avrgr(t,ddt,km,1,8,caf0,1)
c     if(t.ne.0.) call avrg(t,ddt,km,1,8,caf0,1)
      return
      end
c
c
c
      subroutine bndry(t,ddt,km)
c
c------------------------------
c Set the boundaries for X's
c------------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c     dimension rr0(mxpr,my,mzpr)
      dimension zbound(mxpr,my,2)
c
c     revised by jiliang in 05 to 09, 2019
c
      do 100 m=1,8
      do 100 jz=z_first-2,z_last+2
      do 100 jy=1,my
      do 100 jx=x_first-2,x_last+2
      if(rr1(jx,jy,jz).le.rcr) then
      xdif(jx,jy,jz,m)  = 0.
c      else
c      if ((m.ge.5).and.(m.le.7)) then
c      xdif(jx,jy,jz,m)=xdif(jx,jy,jz,m)
c     1    *tanh((rr0(jx,jy,jz)/3.-1.)**2/2.)
c      else
c      xdif(jx,jy,jz,m)=xdif(jx,jy,jz,m)
c     1    *tanh((rr0(jx,jy,jz)/3.-1.)**2/4.)
c      endif
      endif
  100 continue
      call message4(xdif,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c
      if(.not.periodx) then
      if(myleft.eq.mpi_proc_null)then
      do 1 jz=z_first-2,z_last+2
      do 1 jy=1,my
      x(4,jy,jz,1)=rhoi(4,jz,jz)+1.*0.5*(tanh((t-9.)/1.)-
     1        tanh((t-39.)/1.5))
     3       *0.5*(tanh((zz(jz)+0.6)/0.15)-tanh((zz(jz)-0.6)/0.15))
     4       *0.5*(tanh((yy(jy)+0.4)/0.15)-tanh((yy(jy)-0.4)/0.15))
      x(4,jy,jz,2)=0.8*x(4,jy,jz,1)*0.5*(tanh((t-9.)/1.)-
     1        tanh((t-39.)/1.5))
     2       *0.5*(tanh((zz(jz)+0.2)/0.15)-tanh((zz(jz)-0.2)/0.15))
     3       *0.5*(tanh((yy(jy)+0.3)/0.15)-tanh((yy(jy)-0.3)/0.15))
      x(4,jy,jz,3)=0.
      x(4,jy,jz,4)=0.
      x(4,jy,jz,5)=bxi(4,jy,jz)
     1    +bxinfc(jy,jz)*0.0017*dsqrt(t*t*t*7.25*7.25*7.25)
     2    /(dexp(t*7.25*0.0069)+1.)
      x(4,jy,jz,6)=byi(4,jy,jz)
     1    +byinfc(jy,jz)*0.0017*dsqrt(t*t*t*7.25*7.25*7.25)
     2    /(dexp(t*7.25*0.0069)+1.)
      x(4,jy,jz,7)=bzi(4,jy,jz)
     1    +bzinfc(jy,jz)*0.0017*dsqrt(t*t*t*7.25*7.25*7.25)
     2    /(dexp(t*7.25*0.0069)+1.)
      x(4,jy,jz,8)=pri(4,jy,jz)
     1         +0.2*0.5*(tanh((t-9.)/1.)-
     1        tanh((t-39.)/1.5))
     3       *0.5*(tanh((zz(jz)+0.6)/0.15)-tanh((zz(jz)-0.6)/0.15))
     4       *0.5*(tanh((yy(jy)+0.4)/0.15)-tanh((yy(jy)-0.4)/0.15))

      x(3,jy,jz,1)=rhoi(3,jz,jz)+1.*0.5*(tanh((t-9.)/1.)-
     1        tanh((t-39.)/1.5))
     3       *0.5*(tanh((zz(jz)+0.6)/0.15)-tanh((zz(jz)-0.6)/0.15))
     4       *0.5*(tanh((yy(jy)+0.4)/0.15)-tanh((yy(jy)-0.4)/0.15))
      x(3,jy,jz,2)=0.8*x(3,jy,jz,1)*0.5*(tanh((t-9.)/1.)-
     1        tanh((t-39.)/1.5))
     2       *0.5*(tanh((zz(jz)+0.2)/0.15)-tanh((zz(jz)-0.2)/0.15))
     3       *0.5*(tanh((yy(jy)+0.3)/0.15)-tanh((yy(jy)-0.3)/0.15))
      x(3,jy,jz,3)=0.
      x(3,jy,jz,4)=0.
      x(3,jy,jz,5)=bxi(3,jy,jz)
     1    +bxinfc(jy,jz)*0.0017*dsqrt(t*t*t*7.25*7.25*7.25)
     2    /(dexp(t*7.25*0.0069)+1.)
      x(3,jy,jz,6)=byi(3,jy,jz)
     1    +byinfc(jy,jz)*0.0017*dsqrt(t*t*t*7.25*7.25*7.25)
     2    /(dexp(t*7.25*0.0069)+1.)
      x(3,jy,jz,7)=bzi(3,jy,jz)
     1    +bzinfc(jy,jz)*0.0017*dsqrt(t*t*t*7.25*7.25*7.25)
     2    /(dexp(t*7.25*0.0069)+1.)
      x(3,jy,jz,8)=pri(3,jy,jz)
     1         +0.2*0.5*(tanh((t-9.)/1.)-
     1        tanh((t-39.)/1.5))
     3       *0.5*(tanh((zz(jz)+0.6)/0.15)-tanh((zz(jz)-0.6)/0.15))
     4       *0.5*(tanh((yy(jy)+0.4)/0.15)-tanh((yy(jy)-0.4)/0.15))
      do 11 m=1,8
      if(km.eq.1) then
      xdif(3,jy,jz,m)=(x(3,jy,jz,m)-xfold(3,jy,jz,m))/(ddt+1.d-9)
      xdif(4,jy,jz,m)=(x(4,jy,jz,m)-xfold(4,jy,jz,m))/(ddt+1.d-9)
      else
      xdif(3,jy,jz,m)=(x(3,jy,jz,m)-xm(3,jy,jz,m))/(ddt+1.d-9)
      xdif(4,jy,jz,m)=(x(4,jy,jz,m)-xm(4,jy,jz,m))/(ddt+1.d-9)
      endif
   11 continue
    1 continue
      endif
      if(myright.eq.mpi_proc_null)then
      do 2 m=1,8
      do 2 jz=z_first-2,z_last+2
      do 2 jy=1,my
c-----f1=0 3p-----

      xdif(xlast,jy,jz,m)=0
      xdif(xlast-1,jy,jz,m)=0.5*xdif(xlast-2,jy,jz,m)
c
    2 continue
      endif
      else
c
      do 10 m=1,8
      do 10 jz=z_first-1,z_last+1
      do 10 jy=1,my
      if(myleft.eq.mpi_proc_null)then
      xdif(2,jy,jz,m)=xdif(x_last-1,jy,jz,m)
      endif
      if(myright.eq.mpi_proc_null)then
      xdif(xlast,jy,jz,m)=xdif(3,jy,jz,m)
      endif
   10 continue
      endif
c
      if(.not.periody) then
      do 3 m=1,8
      do 3 jz=z_first-2,z_last+2
      do 3 jx=x_first-2,x_last+2

c-----f2=0 3p-----
      xdif(jx,1,jz,m)=0.
      xdif(jx,2,jz,m)=0.5*xdif(jx,3,jz,m)

      xdif(jx,my,jz,m)=0.
      xdif(jx,my-1,jz,m)=0.5*xdif(jx,my-2,jz,m)

    3 continue
      else
      do 30 m=1,8
      do 30 jz=z_first-1,z_last+1
      do 30 jx=x_first-1,x_last+1
      xdif(jx,1,jz,m)=xdif(jx,my-1,jz,m)
      xdif(jx,my,jz,m)=xdif(jx,2,jz,m)
   30 continue
      endif
c
      if(.not.periodz) then
      if(mylow.eq.mpi_proc_null)then
      do 5 m=1,8
      do 5 jy=1,my
      do 5 jx=x_first-2,x_last+2
      xdif(jx,jy,3,m)=0.
      xdif(jx,jy,4,m)=0.5*xdif(jx,jy,5,m)

    5 continue

      endif
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c     call message4(xdif,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c
      if(myupper.eq.mpi_proc_null)then
      if(symmetryz) then
      do 6 jy=1,my
      do 6 jx=x_first-1,x_last+1
      xdif(jx,jy,zlast,1)=xdif(jx,my-jy+1,zlast-2,1)
      xdif(jx,jy,zlast,2)=xdif(jx,my-jy+1,zlast-2,2)
      xdif(jx,jy,zlast,3)=-xdif(jx,my-jy+1,zlast-2,3)
      xdif(jx,jy,zlast,4)=-xdif(jx,my-jy+1,zlast-2,4)
      xdif(jx,jy,zlast,5)=-xdif(jx,my-jy+1,zlast-2,5)
      xdif(jx,jy,zlast,6)=xdif(jx,my-jy+1,zlast-2,6)
      xdif(jx,jy,zlast,7)=xdif(jx,my-jy+1,zlast-2,7)
      xdif(jx,jy,zlast,8)=xdif(jx,my-jy+1,zlast-2,8)
    6 continue
      else
      do 60 m=1,8
      do 60 jy=1,my
      do 60 jx=x_first-2,x_last+2

c-----f2=0 3p-----
      xdif(jx,jy,zlast,m)=0.
      xdif(jx,jy,zlast-1,m)=0.5*xdif(jx,jy,zlast-2,m)

   60 continue

      endif
      endif
      else
c
      if (nrank.eq.(nsize-1)) then
      do  m=1,8
      do  jx=xfirst,xlast
      do  jy=1,my
      wf3(jx,jy,m)=xdif(jx,jy,z_last-1,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
      call MPI_SEND( wf3, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 2,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.0) then
      call MPI_RECV( wf3,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,2,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,8
      do  jx=xfirst,xlast
      do  jy=1,my
      xdif(jx,jy,2,m)=wf3(jx,jy,m)
      wf3(jx,jy,m)=xdif(jx,jy,3,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      call MPI_SEND( wf3,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,3,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.(nsize-1)) then
      call MPI_RECV( wf3, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 3,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,8
      do  jx=xfirst,xlast
      do  jy=1,my
      xdif(jx,jy,zlast,m)=wf3(jx,jy,m)
      enddo
      enddo
      enddo
      endif
c     do 101 m=1,8
c     do 101 jx=xfirst,xlast
c     do 101 jy=1,my
c     if(mylow.eq.mpi_proc_null)then
c     xdif(jx,jy,2,m)=xdif(jx,jy,zlast-1,m)
c     endif
c     if(myupper.eq.mpi_proc_null)then
c     xdif(jx,jy,zlast,m)=xdif(jx,jy,3,m)
c     endif
c 101 continue
      endif



c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c     call message4(xdif,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c
c     call message4(xdif,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
cccccccccccccccccccccccccccccccccccccccccccc uy uz limiter
      if(uyzlimiter) then
      if(mylow.eq.mpi_proc_null)then
      do 51 jz=3,4
      do 51 jy=1,my
      do 51 jx=x_first,x_last+2
      if(xdif(jx,jy,jz,4) .ge. 0.) then
      xdif(jx,jy,jz,4)=0.
      endif
c      if(xdif(jx,jy,jz,2) .le. 0.) then
c      xdif(jx,jy,jz,2)=0.
c      endif
   51 continue
      endif
      if(myupper.eq.mpi_proc_null)then
      do 61 jz=zlast-1,zlast
      do 61 jy=1,my
      do 61 jx=x_first,x_last+2
      if(xdif(jx,jy,jz,4) .le. 0.) then
      xdif(jx,jy,jz,4)=0.
      endif
c      if(xdif(jx,jy,jz,2) .le. 0.) then
c      xdif(jx,jy,jz,2)=0.
c      endif
   61 continue
      endif
      
      do 31 jy=1,2
      do 31 jz=z_first-2,z_last+2
      do 31 jx=x_first,x_last+2
      if(xdif(jx,jy,jz,3) .ge. 0.) then
      xdif(jx,jy,jz,3)=0.
      endif
c      if(xdif(jx,jy,jz,2) .le. 0.) then
c      xdif(jx,jy,jz,2)=0.
c      endif
   31 continue
      do 32 jy=my-1,my
      do 32 jz=z_first-2,z_last+2
      do 32 jx=x_first,x_last+2
      if(xdif(jx,jy,jz,3) .le. 0.) then
      xdif(jx,jy,jz,3)=0.
      endif
c      if(xdif(jx,jy,jz,2) .le. 0.) then
c      xdif(jx,jy,jz,2)=0.
c      endif
   32 continue
      endif
cccccccccccccccccccccccccccccccccccccccccccc
c
      call message4(xdif,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c     do 101 jz=z_first-1,z_last+1
c     do 101 jy=1,my
c     do 101 jx=x_first-1,x_last+1
c     rr0(jx,jy,jz)=sqrt(xx(jx)**2+yy(jy)**2+zz(jz)**2)
c 101 continue
      return
      end
c
c
c
      subroutine bndrytt(t,ms,me)
c
c------------------------------
c Set the boundaries for X's
c------------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      dimension zbound(mxpr,my,2)
c
c     revised by jiliang in 05 to 09, 2019
c
      do 100 m=1,8
      do 100 jz=z_first-2,z_last+2
      do 100 jy=1,my
      do 100 jx=x_first-2,x_last+2
      if(rr1(jx,jy,jz).le.rcr) then
      xm(jx,jy,jz,m)  = 0.
c      else
c      if ((m.ge.5).and.(m.le.7)) then
c      xm(jx,jy,jz,m)=xm(jx,jy,jz,m)*tanh(4.*(rr0(jx,jy,jz)/3.-1.)**2)
c      else
c      xm(jx,jy,jz,m)=xm(jx,jy,jz,m)*tanh((rr0(jx,jy,jz)/3.-1.)**2/4.)
c      endif
      endif
  100 continue
      call message4(xm,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c
      if(.not.periodx) then
      do 1 m=ms,me
      do 1 jz=z_first-2,z_last+2
      do 1 jy=1,my
      if(myleft.eq.mpi_proc_null)then
c     if (m.ne.2) then
c     xm(2,jy,jz,m)=(cxp2*xm(3,jy,jz,m)
c    1        -cxp1*xm(4,jy,jz,m))/(cxp2-cxp1)
c     else
c      xm(4,jy,jz,m)=xm(4,jy,jz,m)
c      xm(3,jy,jz,m)=xm(3,jy,jz,m)

      xm(4,jy,jz,m)=(axp(4)*xm(5,jy,jz,m)+bxp(4)*xm(6,jy,jz,m)+
     1      cxp(4)*xm(7,jy,jz,m))/(axp(4)+bxp(4)+cxp(4))
      xm(3,jy,jz,m)=(axp(3)*xm(4,jy,jz,m)+bxp(3)*xm(5,jy,jz,m)+
     1      cxp(3)*xm(6,jy,jz,m))/(axp(3)+bxp(3)+cxp(3))

c     endif
      endif
      if(myright.eq.mpi_proc_null)then
c-----f1=0 3p-----

      xm(xlast,jy,jz,m)=0
      xm(xlast-1,jy,jz,m)=0.5*xm(xlast-2,jy,jz,m)

      endif
    1 continue
      else
c
      do 11 i=0,nprz-1
      if (nrank.eq.((i+1)*nprx-1)) then
      do  m=1,8
      do  jz=z_first-1,z_last+1
      do  jy=1,my
      wf1(jy,jz-1,m)=xm(x_last-1,jy,jz,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
      call MPI_SEND( wf1, (mzpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     i*nprx, 2,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.(i*nprx)) then
      call MPI_RECV( wf1,(mzpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     (i+1)*nprx-1,2,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,8
      do  jz=z_first-1,z_last+1
      do  jy=1,my
      xm(2,jy,jz,m)=wf1(jy,jz-1,m)
      wf1(jy,jz-1,m)=xm(3,jy,jz,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      call MPI_SEND( wf1,(mzpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     (i+1)*nprx-1,3,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.((i+1)*nprx-1)) then
      call MPI_RECV( wf1, (mzpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     i*nprx, 3,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,8
      do  jz=z_first-1,z_last+1
      do  jy=1,my
      xm(x_last,jy,jz,m)=wf1(jy,jz-1,m)
      enddo
      enddo
      enddo
      endif
   11 continue
c
c     do 10 m=ms,me
c     do 10 jz=z_first-1,z_last+1
c     do 10 jy=2,my-1
c     if(myleft.eq.mpi_proc_null)then
c     xm(2,jy,jz,m)=xm(xlast-1,jy,jz,m)
c     endif
c     if(myright.eq.mpi_proc_null)then
c     xm(xlast,jy,jz,m)=xm(3,jy,jz,m)
c     endif
c  10 continue
      
      endif
c
      if(.not.periody) then
      do 3 m=ms,me
      do 3 jz=z_first-2,z_last+2
      do 3 jx=x_first,x_last+2

c-----f2=0 3p-----

      xm(jx,1,jz,m)=0
      xm(jx,2,jz,m)=0.5*xm(jx,3,jz,m)

      xm(jx,my,jz,m)=0
      xm(jx,my-1,jz,m)=0.5*xm(jx,my-2,jz,m)

    3 continue
      else
      do 30 m=ms,me
      do 30 jz=z_first-1,z_last+1
      do 30 jx=x_first-1,x_last+1
      xm(jx,1,jz,m)=xm(jx,my-1,jz,m)
      xm(jx,my,jz,m)=xm(jx,2,jz,m)
   30 continue
      endif
c
      if(.not.periodz) then
      if(mylow.eq.mpi_proc_null)then
      do 5 m=ms,me
      do 5 jy=1,my
      do 5 jx=x_first,x_last+2

c-----f2=0 3p-----

      xm(jx,jy,3,m)=0
      xm(jx,jy,4,m)=0.5*xm(jx,jy,5,m)

    5 continue

      endif
      if(myupper.eq.mpi_proc_null)then
      if(symmetryz) then
      do 6 jy=1,my
      do 6 jx=x_first-1,x_last+1
      xm(jx,jy,zlast,1)=xm(jx,my-jy+1,zlast-2,1)
      xm(jx,jy,zlast,2)=xm(jx,my-jy+1,zlast-2,2)
      xm(jx,jy,zlast,3)=-xm(jx,my-jy+1,zlast-2,3)
      xm(jx,jy,zlast,4)=-xm(jx,my-jy+1,zlast-2,4)
      xm(jx,jy,zlast,5)=-xm(jx,my-jy+1,zlast-2,5)
      xm(jx,jy,zlast,6)=xm(jx,my-jy+1,zlast-2,6)
      xm(jx,jy,zlast,7)=xm(jx,my-jy+1,zlast-2,7)
      xm(jx,jy,zlast,8)=xm(jx,my-jy+1,zlast-2,8)
    6 continue
      else
      do 60 m=ms,me
      do 60 jy=1,my
      do 60 jx=x_first,x_last+2

c-----f2=0 3p-----

      xm(jx,jy,zlast,m)=0
      xm(jx,jy,zlast-1,m)=0.5*xm(jx,jy,zlast-2,m)

   60 continue

      endif
      endif
      else
c
c     do 11 i=0,nprz-1
      if (nrank.eq.(nsize-1)) then
      do  m=1,8
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      wf3(jx,jy,m)=xm(jx,jy,z_last-1,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
      call MPI_SEND( wf3, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 2,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.0) then
      call MPI_RECV( wf3,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,2,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,8
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      xm(jx,jy,2,m)=wf3(jx,jy,m)
      wf3(jx,jy,m)=xm(jx,jy,3,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      call MPI_SEND( wf3,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,3,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.(nsize-1)) then
      call MPI_RECV( wf3, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 3,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,8
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      xm(jx,jy,z_last,m)=wf3(jx,jy,m)
      enddo
      enddo
      enddo
      endif
c  11 continue
c     do 101 m=ms,me
c     do 101 jx=x_first-1,x_last+1
c     do 101 jy=1,my
c     if(mylow.eq.mpi_proc_null)then
c     xm(jx,jy,2,m)=xm(jx,jy,zlast-1,m)
c     endif
c     if(myupper.eq.mpi_proc_null)then
c     xm(jx,jy,zlast,m)=xm(jx,jy,3,m)
c     endif
c 101 continue
      endif
c
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
cccccccccccccccccccccccccccccccccccccccccccc uy uz limiter
      if(uyzlimiter) then
      if(mylow.eq.mpi_proc_null)then
      do 51 jz=3,4
      do 51 jy=1,my
      do 51 jx=x_first,x_last+2
      if(xm(jx,jy,jz,4) .ge. 0.) then
      xm(jx,jy,jz,4)=0.
      endif
c      if(xm(jx,jy,jz,2) .le. 0.) then
c      xm(jx,jy,jz,2)=0.
c      endif
   51 continue
      endif
      if(myupper.eq.mpi_proc_null)then
      do 61 jz=zlast-1,zlast
      do 61 jy=1,my
      do 61 jx=x_first,x_last+2
      if(xm(jx,jy,jz,4) .le. 0.) then
      xm(jx,jy,jz,4)=0.
      endif
c      if(xm(jx,jy,jz,2) .le. 0.) then
c      xm(jx,jy,jz,2)=0.
c      endif
   61 continue
      endif
      
      do 31 jy=1,2
      do 31 jz=z_first-2,z_last+2
      do 31 jx=x_first,x_last+2
      if(xm(jx,jy,jz,3) .ge. 0.) then
      xm(jx,jy,jz,3)=0.
      endif
c      if(xm(jx,jy,jz,2) .le. 0.) then
c      xm(jx,jy,jz,2)=0.
c      endif
   31 continue
      do 32 jy=my-1,my
      do 32 jz=z_first-2,z_last+2
      do 32 jx=x_first,x_last+2
      if(xm(jx,jy,jz,3) .le. 0.) then
      xm(jx,jy,jz,3)=0.
      endif
c      if(xm(jx,jy,jz,2) .le. 0.) then
c      xm(jx,jy,jz,2)=0.
c      endif
   32 continue
      endif
cccccccccccccccccccccccccccccccccccccccccccc
c
      call message4(xm,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c
      return
      end
c
      subroutine flux(m)
c
c-------------------------------
c  Calculate fluxes
c  Notations: X1    X2     X3     X4     X5   X6   X7  X8
c             rho   rhovx  rhovy  rhovz  bx   by   bz  e
c-------------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c
c [0] Preparation
c Find current stored in the working array, the three components of array w
c corresponds to jx, jy, jz.
c
c     call foreta(time)
      if(m.eq.1) then
      do 1 jz=z_first-2,z_last+2
      do 1 jy=1,my
      do 1 jx=x_first-2,x_last+2
c
c [1] Continuity eq.
c
      fs(jx,jy,jz)=x(jx,jy,jz,2)
      gs(jx,jy,jz)=x(jx,jy,jz,3)
      hs(jx,jy,jz)=x(jx,jy,jz,4)
    1 continue
c
      else
      if(m.eq.2) then
c [2] Momentum eq.
      do 2 jz=z_first-2,z_last+2
      do 2 jy=1,my
      do 2 jx=x_first-2,x_last+2
      fs(jx,jy,jz)=x(jx,jy,jz,2)**2/x(jx,jy,jz,1)
     1            +(x(jx,jy,jz,8)-pri(jx,jy,jz))
      gs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,3)/x(jx,jy,jz,1)
      hs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,4)/x(jx,jy,jz,1)
    2 continue
      else
      if(m.eq.3) then
      do 3 jz=z_first-2,z_last+2
      do 3 jy=1,my
      do 3 jx=x_first-2,x_last+2
      fs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,3)/x(jx,jy,jz,1)
      gs(jx,jy,jz)=x(jx,jy,jz,3)**2/x(jx,jy,jz,1)
     1            +(x(jx,jy,jz,8)-pri(jx,jy,jz))
      hs(jx,jy,jz)=x(jx,jy,jz,4)*x(jx,jy,jz,3)/x(jx,jy,jz,1)
c
    3 continue
      else
      if(m.eq.4) then
      do 4 jz=z_first-2,z_last+2
      do 4 jy=1,my
      do 4 jx=x_first-2,x_last+2
      fs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,4)/x(jx,jy,jz,1)
      gs(jx,jy,jz)=x(jx,jy,jz,3)*x(jx,jy,jz,4)/x(jx,jy,jz,1)
      hs(jx,jy,jz)=x(jx,jy,jz,4)**2/x(jx,jy,jz,1)
     1            +(x(jx,jy,jz,8)-pri(jx,jy,jz))
    4 continue
      else
      if(m.eq.5) then
      do 5 jz=z_first-2,z_last+2
      do 5 jy=1,my
      do 5 jx=x_first-2,x_last+2
c
c [3] Magnetic induction eq.
c
      vcrbz=(x(jx,jy,jz,2)*x(jx,jy,jz,6)-x(jx,jy,jz,3)
     1        *x(jx,jy,jz,5))/x(jx,jy,jz,1)
      vcrby=(x(jx,jy,jz,4)*x(jx,jy,jz,5)-x(jx,jy,jz,2)
     1        *x(jx,jy,jz,7))/x(jx,jy,jz,1)
c
      fs(jx,jy,jz)=0.
      gs(jx,jy,jz)=-vcrbz+etan(jx,jy,jz)*cur(jx,jy,jz,3)
      hs(jx,jy,jz)=vcrby-(etan(jx,jy,jz)*cur(jx,jy,jz,2))
c    1             -eta0*cj(jx,jz))
    5 continue
      else
      if(m.eq.6) then
      do 6 jz=z_first-2,z_last+2
      do 6 jy=1,my
      do 6 jx=x_first-2,x_last+2
c
c [3] Magnetic induction eq.
c
      vcrbz=(x(jx,jy,jz,2)*x(jx,jy,jz,6)-x(jx,jy,jz,3)
     1        *x(jx,jy,jz,5))/x(jx,jy,jz,1)
      vcrbx=(x(jx,jy,jz,3)*x(jx,jy,jz,7)-x(jx,jy,jz,6)
     1        *x(jx,jy,jz,4))/x(jx,jy,jz,1)
c
      fs(jx,jy,jz)=vcrbz-etan(jx,jy,jz)*cur(jx,jy,jz,3)
      gs(jx,jy,jz)=0.
      hs(jx,jy,jz)=-vcrbx+etan(jx,jy,jz)*cur(jx,jy,jz,1)
    6 continue
      else
      if(m.eq.7) then
      do 7 jz=z_first-2,z_last+2
      do 7 jy=1,my
      do 7 jx=x_first-2,x_last+2
c
c [3] Magnetic induction eq.
c
      vcrby=(x(jx,jy,jz,4)*x(jx,jy,jz,5)-x(jx,jy,jz,2)
     1        *x(jx,jy,jz,7))/x(jx,jy,jz,1)
      vcrbx=(x(jx,jy,jz,3)*x(jx,jy,jz,7)-x(jx,jy,jz,6)
     1        *x(jx,jy,jz,4))/x(jx,jy,jz,1)
c
      fs(jx,jy,jz)=-vcrby+(etan(jx,jy,jz)*cur(jx,jy,jz,2))
c    1             -eta0*cj(jx,jz))
      gs(jx,jy,jz)=vcrbx-etan(jx,jy,jz)*cur(jx,jy,jz,1)
      hs(jx,jy,jz)=0.
    7 continue
      else
      if(m.eq.8) then
      do 8 jz=z_first-2,z_last+2
      do 8 jy=1,my
      do 8 jx=x_first-2,x_last+2
c
c [4] Energy eq.
c
      b2=x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2
      bdotv=(x(jx,jy,jz,5)*x(jx,jy,jz,2)+x(jx,jy,jz,6)*x(jx,jy,jz,3)
     1      +x(jx,jy,jz,7)*x(jx,jy,jz,4))/x(jx,jy,jz,1)
      eng=x(jx,jy,jz,8)+.5*b2
c
      fs(jx,jy,jz)=eng*x(jx,jy,jz,2)/x(jx,jy,jz,1)
     1            -bdotv*x(jx,jy,jz,5)
     2            +etan(jx,jy,jz)*(cur(jx,jy,jz,2)*x(jx,jy,jz,7)
     3            -cur(jx,jy,jz,3)*x(jx,jy,jz,6))
      gs(jx,jy,jz)=eng*x(jx,jy,jz,3)/x(jx,jy,jz,1)
     1            -bdotv*x(jx,jy,jz,6)
     2            +etan(jx,jy,jz)*(cur(jx,jy,jz,3)*x(jx,jy,jz,5)
     3            -cur(jx,jy,jz,1)*x(jx,jy,jz,7))
      hs(jx,jy,jz)=eng*x(jx,jy,jz,4)/x(jx,jy,jz,1)
     1            -bdotv*x(jx,jy,jz,7)
     2            +etan(jx,jy,jz)*(cur(jx,jy,jz,1)*x(jx,jy,jz,6)
     3            -cur(jx,jy,jz,2)*x(jx,jy,jz,5))
c
    8 continue
      else
      endif
      endif
      endif
      endif
      endif
      endif
      endif
      endif
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c      if(nrank .eq. 24) print*,'flux fs_x',m,fs(3,13,101),fs(4,13,101)
c     1                       ,fs(5,13,101),fs(6,13,101),fs(7,13,101)
c      if(nrank .eq. 24) print*,'flux gs_y',m,gs(5,11,101),gs(5,12,101)
c     1                       ,gs(5,13,101),gs(5,14,101),gs(5,15,101)
c      if(nrank .eq. 24) print*,'flux hs_z',m,hs(5,13,99),hs(5,13,100)
c     1                       ,hs(5,13,101),hs(5,13,102),hs(5,13,103)
c     call message3(fs,mxpr,my,mzpr,myleft,myright,mylow,myupper)
c     call message3(gs,mxpr,my,mzpr,myleft,myright,mylow,myupper)
c     call message3(hs,mxpr,my,mzpr,myleft,myright,mylow,myupper)
      return
      end
c
c
c
      subroutine current
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      dimension v(mxpr,my,mzpr),u(mxpr,my,mzpr)
c
c  define statement functions
c  d1fc= df/dx with fourth-order accuracz central difference
      d1fc(fm2,fm,f0,fp,fp2,a,b,c,d)=
     1   a*(fp-f0)+b*(f0-fm)+c*(fp2-f0)+d*(f0-fm2)
c  d1fc= df/dx with central difference
c      d1fc(fm,f0,fp,xm1,x0,xp1)=
c     1  ((xm1-x0)/(xp1-x0)*(fp-f0)
c     1   -(xp1-x0)/(xm1-x0)*(fm-f0))/(xm1-xp1)
c
      do 10 jz=z_first-2,z_last+2
      do 10 jy=1,my
      do 10 jx=x_first-2,x_last+2
      u(jx,jy,jz)=x(jx,jy,jz,5)-bxi(jx,jy,jz)
      v(jx,jy,jz)=x(jx,jy,jz,6)-byi(jx,jy,jz)
      w(jx,jy,jz)=x(jx,jy,jz,7)-bzi(jx,jy,jz)
   10 continue
c
      do 1 jz=z_first,z_last
      do 1 jy=3,my-2
      do 1 jx=x_first,x_last
      cur(jx,jy,jz,1)=d1fc(w(jx,jy-2,jz),w(jx,jy-1,jz),w(jx,jy,jz)
     1    ,w(jx,jy+1,jz),w(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
     2               -d1fc(v(jx,jy,jz-2),v(jx,jy,jz-1),v(jx,jy,jz)
     3    ,v(jx,jy,jz+1),v(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
      cur(jx,jy,jz,2)=d1fc(u(jx,jy,jz-2),u(jx,jy,jz-1),u(jx,jy,jz)
     1    ,u(jx,jy,jz+1),u(jx,jy,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     2               -d1fc(w(jx-2,jy,jz),w(jx-1,jy,jz),w(jx,jy,jz)
     3    ,w(jx+1,jy,jz),w(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
      cur(jx,jy,jz,3)=d1fc(v(jx-2,jy,jz),v(jx-1,jy,jz),v(jx,jy,jz)
     1    ,v(jx+1,jy,jz),v(jx+2,jy,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     2               -d1fc(u(jx,jy-2,jz),u(jx,jy-1,jz),u(jx,jy,jz)
     3    ,u(jx,jy+1,jz),u(jx,jy+2,jz),ay1(jy),by1(jy),cy1(jy),dy1(jy))
    1 continue
c
      call bndryc
c
      call avrgc(1,3,4)
c     call avrgc(1,3,2)
c     call avrgre(1,3,1)
c
      return
      end
c
c
c
      subroutine diverB
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c
c  define statement functions
c  d1fc= df/dx with fourth-order accuracz central difference
      d1fc(fm2,fm,f0,fp,fp2,a,b,c,d)=
     1   a*(fp-f0)+b*(f0-fm)+c*(fp2-f0)+d*(f0-fm2)
c  d1fc= df/dx with central difference
c      d1fc(fm,f0,fp,xm1,x0,xp1)=
c     1  ((xm1-x0)/(xp1-x0)*(fp-f0)
c     1   -(xp1-x0)/(xm1-x0)*(fm-f0))/(xm1-xp1)
c
      do 1 jz=z_first,z_last
      do 1 jy=3,my-2
      do 1 jx=x_first,x_last
      divB(jx,jy,jz)=d1fc(x(jx-2,jy,jz,5),x(jx-1,jy,jz,5),x(jx,jy,jz,5)
     1 ,x(jx+1,jy,jz,5),x(jx+2,jy,jz,5),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     2              +d1fc(x(jx,jy-2,jz,6),x(jx,jy-1,jz,6),x(jx,jy,jz,6)
     3 ,x(jx,jy+1,jz,6),x(jx,jy+2,jz,6),ay1(jy),by1(jy),cy1(jy),dy1(jy))
     4              +d1fc(x(jx,jy,jz-2,7),x(jx,jy,jz-1,7),x(jx,jy,jz,7)
     5 ,x(jx,jy,jz+1,7),x(jx,jy,jz+2,7),az1(jz),bz1(jz),cz1(jz),dz1(jz))
    1 continue
c
      call bndrydiv
c
c      call avrgdiv(1,1,4)
c
      return
      end
c
c
c
      subroutine bndryc
c
c------------------------------
c Set the boundaries for current's
c------------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      do 100 m=1,3
      do 100 jz=z_first-2,z_last+2
      do 100 jy=1,my
      do 100 jx=x_first-2,x_last+2
      if(rr1(jx,jy,jz).le.rcr) then
      cur(jx,jy,jz,m)  = 0.
c      else
c      cur(jx,jy,jz,m)=cur(jx,jy,jz,m)*tanh((rr0(jx,jy,jz)/3.-1.)**2/4.)
      endif
  100 continue
      if(.not.periodx) then
      do 2 m=1,3
      do 2 jz=z_first-2,z_last+2
      do 2 jy=1,my
      if(myleft.eq.mpi_proc_null)then
c-----f1=0 3p-----
      cur(4,jy,jz,m)=(axp(4)*cur(5,jy,jz,m)+bxp(4)*cur(6,jy,jz,m)+
     1      cxp(4)*cur(7,jy,jz,m))/(axp(4)+bxp(4)+cxp(4))
      cur(3,jy,jz,m)=(axp(3)*cur(4,jy,jz,m)+bxp(3)*cur(5,jy,jz,m)+
     1      cxp(3)*cur(6,jy,jz,m))/(axp(3)+bxp(3)+cxp(3))
c
      endif
      if(myright.eq.mpi_proc_null)then
c-----f1=0 3p-----

      cur(xlast,jy,jz,m)=0.
      cur(xlast-1,jy,jz,m)=0.5*cur(xlast-2,jy,jz,m) 
c
      endif
    2 continue
      else
c
      do 21 m=1,3
      do 21 jz=z_first-1,z_last+1
      do 21 jy=2,my-1
      if(myleft.eq.mpi_proc_null)then
      cur(2,jy,jz,m)=cur(x_last-1,jy,jz,m)
      endif
      if(myright.eq.mpi_proc_null)then
      cur(x_last,jy,jz,m)=cur(3,jy,jz,m)
      endif
   21 continue
      endif
c
      if(.not.periody) then
      do 3 m=1,3
      do 3 jz=z_first-2,z_last+2
      do 3 jx=x_first-2,x_last+2

      cur(jx,1,jz,m)=0
      cur(jx,2,jz,m)=0.5*cur(jx,3,jz,m)

      cur(jx,my,jz,m)=0
      cur(jx,my-1,jz,m)=0.5*cur(jx,my-2,jz,m)

    3 continue
      else
      do 30 m=1,3
      do 30 jz=z_first-1,z_last+1
      do 30 jx=x_first-1,x_last+1
      cur(jx,1,jz,m)=cur(jx,my-1,jz,m)
      cur(jx,my,jz,m)=cur(jx,2,jz,m)
   30 continue
      endif
c
      if(.not.periodz) then
      if(mylow.eq.mpi_proc_null)then
      do 4 m=1,3
      do 4 jy=1,my
      do 4 jx=x_first-2,x_last+2

      cur(jx,jy,3,m)=0
      cur(jx,jy,4,m)=0.5*cur(jx,jy,5,m)

    4 continue
      endif
c
      if(myupper.eq.mpi_proc_null)then
      if(symmetryz) then
      do 5 jy=1,my
      do 5 jx=x_first-1,x_last+1
      cur(jx,jy,zlast,1)=-cur(jx,my-jy+1,zlast-2,1)
      cur(jx,jy,zlast,2)=cur(jx,my-jy+1,zlast-2,2)
      cur(jx,jy,zlast,3)=cur(jx,my-jy+1,zlast-2,3)
    5 continue
      else
      do 6 m=1,3
      do 6 jy=1,my
      do 6 jx=x_first-2,x_last+2

      cur(jx,jy,zlast,m)=0
      cur(jx,jy,zlast-1,m)=0.5*cur(jx,jy,zlast-2,m)

    6 continue
      endif
      endif
      else
      if (nrank.eq.(nsize-1)) then
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      wf4(jx,jy,m)=cur(jx,jy,z_last-1,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
      call MPI_SEND( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 2,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.0) then
      call MPI_RECV( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,2,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
c
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      cur(jx,jy,2,m)=wf4(jx,jy,m)
      wf4(jx,jy,m)=cur(jx,jy,3,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      call MPI_SEND( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,3,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.(nsize-1)) then
      call MPI_RECV( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 3,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      cur(jx,jy,z_last,m)=wf4(jx,jy,m)
      enddo
      enddo
      enddo
      endif
      endif
c
      call message4(cur,mxpr,my,mzpr,3,myleft,myright,mylow,myupper)
c
      return
      end
c
c
c
      subroutine bndryh
c
c------------------------------
c Set the boundaries for current's
c------------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      do 100 m=1,3
      do 100 jz=z_first-2,z_last+2
      do 100 jy=1,my
      do 100 jx=x_first-2,x_last+2
      if(rr1(jx,jy,jz).le.rcr) then
      xdh(jx,jy,jz,m)  = 0.
c      else
c      xdh(jx,jy,jz,m)=xdh(jx,jy,jz,m)*tanh((rr0(jx,jy,jz)/3.-1.)**2/4.)
      endif
  100 continue
c      if(nrank .eq. 0) print*,'1',time,nstep
c      if(nrank .eq. 0) print*,nrank,x_first-1,x_first-2
c      if(nrank .eq. 0) print*,'2',time,nstep
c      if(nrank .eq. 0) print*,'3',time,nstep
c      if(nrank .eq. 0) print*,'4',time,nstep
c      if(nrank .eq. 0) print*,'5',time,nstep
      if(.not.periodx) then
      do 2 m=1,3
      do 2 jz=z_first-2,z_last+2
      do 2 jy=1,my
      if(myleft.eq.mpi_proc_null)then
c-----f1=0 3p-----
      xdh(4,jy,jz,m)=(axp(4)*xdh(5,jy,jz,m)+bxp(4)*xdh(6,jy,jz,m)+
     1      cxp(4)*xdh(7,jy,jz,m))/(axp(4)+bxp(4)+cxp(4))
      xdh(3,jy,jz,m)=(axp(3)*xdh(4,jy,jz,m)+bxp(3)*xdh(5,jy,jz,m)+
     1      cxp(3)*xdh(6,jy,jz,m))/(axp(3)+bxp(3)+cxp(3))

      endif
      if(myright.eq.mpi_proc_null)then
c-----f1=0 3p-----

      xdh(xlast,jy,jz,m)=0
      xdh(xlast-1,jy,jz,m)=0.5*xdh(xlast-2,jy,jz,m)
c
      endif
    2 continue
      else
c
      do 21 m=1,3
      do 21 jz=z_first-1,z_last+1
      do 21 jy=2,my-1
      if(myleft.eq.mpi_proc_null)then
      xdh(2,jy,jz,m)=xdh(x_last-1,jy,jz,m)
      endif
      if(myright.eq.mpi_proc_null)then
      xdh(x_last,jy,jz,m)=xdh(3,jy,jz,m)
      endif
   21 continue
      endif
c
      if(.not.periody) then
      do 3 m=1,3
      do 3 jz=z_first-2,z_last+2
      do 3 jx=x_first-2,x_last+2

c-----f2=0 3p-----

      xdh(jx,1,jz,m)=0
      xdh(jx,2,jz,m)=0.5*xdh(jx,3,jz,m)

      xdh(jx,my,jz,m)=0
      xdh(jx,my-1,jz,m)=0.5*xdh(jx,my-2,jz,m)

    3 continue
      else
      do 30 m=1,3
      do 30 jz=z_first-1,z_last+1
      do 30 jx=x_first-1,x_last+1
      xdh(jx,1,jz,m)=xdh(jx,my-1,jz,m)
      xdh(jx,my,jz,m)=xdh(jx,2,jz,m)
   30 continue
      endif
c
      if(.not.periodz) then
      if(mylow.eq.mpi_proc_null)then
      do 4 m=1,3
      do 4 jy=1,my
      do 4 jx=x_first-2,x_last+2

c-----f2=0 3p-----

      xdh(jx,jy,3,m)=0
      xdh(jx,jy,4,m)=0.5*xdh(jx,jy,5,m)

    4 continue
      endif
c
      if(myupper.eq.mpi_proc_null)then
      if(symmetryz) then
      do 5 jy=1,my
      do 5 jx=x_first-1,x_last+1
      xdh(jx,jy,zlast,1)=-xdh(jx,my-jy+1,zlast-2,1)
      xdh(jx,jy,zlast,2)=xdh(jx,my-jy+1,zlast-2,2)
      xdh(jx,jy,zlast,3)=xdh(jx,my-jy+1,zlast-2,3)
    5 continue
      else
      do 6 m=1,3
      do 6 jy=1,my
      do 6 jx=x_first-2,x_last+2

c-----f2=0 3p-----

      xdh(jx,jy,zlast,m)=0
      xdh(jx,jy,zlast-1,m)=0.5*xdh(jx,jy,zlast-2,m)

    6 continue
      endif
      endif
      else
      if (nrank.eq.(nsize-1)) then
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      wf4(jx,jy,m)=xdh(jx,jy,z_last-1,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
      call MPI_SEND( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 2,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.0) then
      call MPI_RECV( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,2,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
c
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      xdh(jx,jy,2,m)=wf4(jx,jy,m)
      wf4(jx,jy,m)=xdh(jx,jy,3,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      call MPI_SEND( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,3,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.(nsize-1)) then
      call MPI_RECV( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 3,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      xdh(jx,jy,z_last,m)=wf4(jx,jy,m)
      enddo
      enddo
      enddo
      endif
      endif
c
      call message4(xdh,mxpr,my,mzpr,3,myleft,myright,mylow,myupper)
c
      return
      end
c
c
c
      subroutine bndryE
c
c------------------------------
c Set the boundaries for current's
c------------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      do 100 m=1,3
      do 100 jz=z_first-2,z_last+2
      do 100 jy=1,my
      do 100 jx=x_first-2,x_last+2
      if(rr1(jx,jy,jz).le.rcr) then
      efld(jx,jy,jz,m)  = 0.
c      else
c      efld(jx,jy,jz,m)=efld(jx,jy,jz,m)*tanh((rr0(jx,jy,jz)/3.-1.)**2/4.)
      endif
  100 continue
      if(.not.periodx) then
      do 2 m=1,3
      do 2 jz=z_first-2,z_last+2
      do 2 jy=1,my
      if(myleft.eq.mpi_proc_null)then
c-----f1=0 3p-----
      efld(4,jy,jz,m)=(axp(4)*efld(5,jy,jz,m)+bxp(4)*efld(6,jy,jz,m)+
     1      cxp(4)*efld(7,jy,jz,m))/(axp(4)+bxp(4)+cxp(4))
      efld(3,jy,jz,m)=(axp(3)*efld(4,jy,jz,m)+bxp(3)*efld(5,jy,jz,m)+
     1      cxp(3)*efld(6,jy,jz,m))/(axp(3)+bxp(3)+cxp(3))
c
      endif
      if(myright.eq.mpi_proc_null)then
c-----f1=0 3p-----
      efld(xlast,jy,jz,m)=0
      efld(xlast-1,jy,jz,m)=0.5*efld(xlast-2,jy,jz,m)
c
      endif
    2 continue
      else
c
      do 21 m=1,3
      do 21 jz=z_first-1,z_last+1
      do 21 jy=2,my-1
      if(myleft.eq.mpi_proc_null)then
      efld(2,jy,jz,m)=efld(x_last-1,jy,jz,m)
      endif
      if(myright.eq.mpi_proc_null)then
      efld(x_last,jy,jz,m)=efld(3,jy,jz,m)
      endif
   21 continue
      endif
c
      if(.not.periody) then
      do 3 m=1,3
      do 3 jz=z_first-2,z_last+2
      do 3 jx=x_first-2,x_last+2

c-----f2=0 3p-----
      efld(jx,1,jz,m)=0
      efld(jx,2,jz,m)=0.5*efld(jx,3,jz,m)

      efld(jx,my,jz,m)=0
      efld(jx,my-1,jz,m)=0.5*efld(jx,my-2,jz,m)

    3 continue
      else
      do 30 m=1,3
      do 30 jz=z_first-1,z_last+1
      do 30 jx=x_first-1,x_last+1
      efld(jx,1,jz,m)=efld(jx,my-1,jz,m)
      efld(jx,my,jz,m)=efld(jx,2,jz,m)
   30 continue
      endif
c
      if(.not.periodz) then
      if(mylow.eq.mpi_proc_null)then
      do 4 m=1,3
      do 4 jy=1,my
      do 4 jx=x_first-2,x_last+2

c-----f2=0 3p-----
      efld(jx,jy,3,m)=0
      efld(jx,jy,4,m)=0.5*efld(jx,jy,5,m)

    4 continue
      endif
c
      if(myupper.eq.mpi_proc_null)then
      if(symmetryz) then
      do 5 jy=1,my
      do 5 jx=x_first-1,x_last+1
      efld(jx,jy,zlast,1)=-efld(jx,my-jy+1,zlast-2,1)
      efld(jx,jy,zlast,2)=efld(jx,my-jy+1,zlast-2,2)
      efld(jx,jy,zlast,3)=efld(jx,my-jy+1,zlast-2,3)
    5 continue
      else
      do 6 m=1,3
      do 6 jy=1,my
      do 6 jx=x_first-2,x_last+2
c-----f2=0 3p-----
      efld(jx,jy,zlast,m)=0
      efld(jx,jy,zlast-1,m)=0.5*efld(jx,jy,zlast-2,m)
    6 continue
      endif
      endif
      else
      if (nrank.eq.(nsize-1)) then
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      wf4(jx,jy,m)=efld(jx,jy,z_last-1,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
      call MPI_SEND( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 2,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.0) then
      call MPI_RECV( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,2,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
c
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      efld(jx,jy,2,m)=wf4(jx,jy,m)
      wf4(jx,jy,m)=efld(jx,jy,3,m)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      call MPI_SEND( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,3,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.(nsize-1)) then
      call MPI_RECV( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 3,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,3
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      efld(jx,jy,z_last,m)=wf4(jx,jy,m)
      enddo
      enddo
      enddo
      endif
      endif
c
      call message4(efld,mxpr,my,mzpr,3,myleft,myright,mylow,myupper)
c
      return
      end
c
c
c
      subroutine bndrydiv
c
c------------------------------
c Set the boundaries for current's
c------------------------------
c
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      do 100 jz=z_first-2,z_last+2
      do 100 jy=1,my
      do 100 jx=x_first-2,x_last+2
      if(rr0(jx,jy,jz).le.3.) then
c      divB(jx,jy,jz)  = 0.
c      else
c      divB(jx,jy,jz)=divB(jx,jy,jz)*tanh((rr0(jx,jy,jz)/3.-1.)**2/4.)
      endif
  100 continue
c
      if(.not.periodx) then
      do 2 jz=z_first-2,z_last+2
      do 2 jy=1,my
      if(myleft.eq.mpi_proc_null)then
c-----f1=0 3p-----
      divB(4,jy,jz)=(axp(4)*divB(5,jy,jz)+bxp(4)*divB(6,jy,jz)+
     1      cxp(4)*divB(7,jy,jz))/(axp(4)+bxp(4)+cxp(4))
      divB(3,jy,jz)=(axp(3)*divB(4,jy,jz)+bxp(3)*divB(5,jy,jz)+
     1      cxp(3)*divB(6,jy,jz))/(axp(3)+bxp(3)+cxp(3))
c
c-----f2=0 2p-----
c      divB(4,jy,jz)=divB(5,jy,jz)
c     1   +1*(divB(5,jy,jz)-divB(6,jy,jz))
c      divB(3,jy,jz)=divB(4,jy,jz)
c     1   +1*(divB(4,jy,jz)-divB(5,jy,jz))
c-----f2=0 3p-----
c      aa1=5./2.
c      bb1=-2.
c      cc1=1./2.
c      divB(4,jy,jz)=aa1*divB(5,jy,jz)+bb1*divB(6,jy,jz)
c     1              +cc1*divB(7,jy,jz)
c      divB(3,jy,jz)=aa1*divB(4,jy,jz)+bb1*divB(5,jy,jz)
c     1              +cc1*divB(6,jy,jz)
c-----f2=0 4p-----

c-----f2=0 5p-----
c
c-----f3=0 3p-----
c      divB(4,jy,jz)=3*divB(5,jy,jz)-3*divB(6,jy,jz)+
c     1      1*divB(7,jy,jz)
c      divB(3,jy,jz)=3*divB(4,jy,jz)-3*divB(5,jy,jz)+
c     1      1*divB(6,jy,jz)
c-----f3=0 4p-----

c-----f3=0 5p-----

c-----f3=0 6p-----

      endif
      if(myright.eq.mpi_proc_null)then
c-----f1=0 3p-----
      divB(xlast-1,jy,jz)=(axm(xlast-1)*divB(xlast-2,jy,jz)
     1  +bxm(xlast-1)*divB(xlast-3,jy,jz)+cxm(xlast-1)
     1  *divB(xlast-4,jy,jz))/(axm(xlast-1)+bxm(xlast-1)+cxm(xlast-1))
      divB(xlast,jy,jz)=(axm(xlast)*divB(xlast-1,jy,jz)+bxm(xlast)
     1  *divB(xlast-2,jy,jz)+cxm(xlast)*divB(xlast-3,jy,jz))
     1  /(axm(xlast)+bxm(xlast)+cxm(xlast))
c
c-----f2=0 2p-----
c      divB(xlast-1,jy,jz)=divB(xlast-2,jy,jz)
c     1   +1*(divB(xlast-2,jy,jz)-divB(xlast-3,jy,jz))
c      divB(xlast,jy,jz)=divB(xlast-1,jy,jz)
c     1   +1*(divB(xlast-1,jy,jz)-divB(xlast-2,jy,jz))
c-----f2=0 3p-----
c      aa1=5./2.
c      bb1=-2.
c      cc1=1./2.
c      divB(xlast-1,jy,jz)=aa1*divB(xlast-2,jy,jz)
c     1                    +bb1*divB(xlast-3,jy,jz)
c     1                    +cc1*divB(xlast-4,jy,jz)
c      divB(xlast,jy,jz)=aa1*divB(xlast-1,jy,jz)
c     1                  +bb1*divB(xlast-2,jy,jz)
c     1                  +cc1*divB(xlast-3,jy,jz)
c-----f2=0 4p-----

c-----f2=0 5p-----
c
c-----f3=0 3p-----
c      divB(xlast-1,jy,jz)=3*divB(xlast-2,jy,jz)
c     1  -3*divB(xlast-3,jy,jz)+1
c     1  *divB(xlast-4,jy,jz)
c      divB(xlast,jy,jz)=3*divB(xlast-1,jy,jz)-3
c     1  *divB(xlast-2,jy,jz)+1*divB(xlast-3,jy,jz)
c-----f3=0 4p-----

c-----f3=0 5p-----

c-----f3=0 6p-----

      endif
    2 continue
      else
c
      do 21 jz=z_first-1,z_last+1
      do 21 jy=2,my-1
      if(myleft.eq.mpi_proc_null)then
      divB(2,jy,jz)=divB(x_last-1,jy,jz)
      endif
      if(myright.eq.mpi_proc_null)then
      divB(x_last,jy,jz)=divB(3,jy,jz)
      endif
   21 continue
      endif
c
      if(.not.periody) then
      do 3 jz=z_first-2,z_last+2
      do 3 jx=x_first-2,x_last+2
c-----f1=0 3p-----
c      divB(jx,2,jz)=(ayp(2)*divB(jx,3,jz)+byp(2)*divB(jx,4,jz)+
c     1      cyp(2)*divB(jx,5,jz))/(ayp(2)+byp(2)+cyp(2))
c      divB(jx,my-1,jz)=(aym(my-1)*divB(jx,my-2,jz)+bym(my-1)
c     1      *divB(jx,my-3,jz)+cym(my-1)*divB(jx,my-4,jz))
c     1      /(aym(my-1)+bym(my-1)+cym(my-1))
c      divB(jx,1,jz)=(ayp(1)*divB(jx,2,jz)+byp(1)*divB(jx,3,jz)+
c     1      cyp(1)*divB(jx,4,jz))/(ayp(1)+byp(1)+cyp(1))
c      divB(jx,my,jz)=(aym(my)*divB(jx,my-1,jz)+bym(my)
c     1      *divB(jx,my-2,jz)+cym(my)*divB(jx,my-3,jz))
c     1      /(aym(my)+bym(my)+cym(my))
c-----f1=0 4p-----
c      aa1=48./25.
c      bb1=-36./25.
c      cc1=16./25.
c      dd1=-3./25.
c      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
c     1              +cc1*divB(jx,5,jz)+dd1*divB(jx,6,jz)
c      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
c     1                 +bb1*divB(jx,my-3,jz)
c     1                 +cc1*divB(jx,my-4,jz)
c     1                 +dd1*divB(jx,my-5,jz)
c      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
c     1              +cc1*divB(jx,4,jz)+dd1*divB(jx,5,jz)
c      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
c     1               +bb1*divB(jx,my-2,jz)
c     1               +cc1*divB(jx,my-3,jz)
c     1               +dd1*divB(jx,my-4,jz)
c-----f1=0 5p-----
c      aa1=300./137.
c      bb1=-300./137.
c      cc1=200./137.
c      dd1=-75./137.
c      ee1=12./137.
c      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
c     1              +cc1*divB(jx,5,jz)+dd1*divB(jx,6,jz)
c     1              +ee1*divB(jx,7,jz)
c      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
c     1                 +bb1*divB(jx,my-3,jz)
c     1                 +cc1*divB(jx,my-4,jz)
c     1                 +dd1*divB(jx,my-5,jz)
c     1                 +ee1*divB(jx,my-6,jz)
c      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
c     1              +cc1*divB(jx,4,jz)+dd1*divB(jx,5,jz)
c     1              +ee1*divB(jx,6,jz)
c      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
c     1               +bb1*divB(jx,my-2,jz)
c     1               +cc1*divB(jx,my-3,jz)
c     1               +dd1*divB(jx,my-4,jz)
c     1               +ee1*divB(jx,my-5,jz)
c-----f2=0 2p-----
c      divB(jx,2,jz)=divB(jx,3,jz)
c     1   +1*(divB(jx,3,jz)-divB(jx,4,jz))
c      divB(jx,my-1,jz)=divB(jx,my-2,jz)
c     1   +1*(divB(jx,my-2,jz)-divB(jx,my-3,jz))
c      divB(jx,1,jz)=divB(jx,2,jz)
c     1   +1*(divB(jx,2,jz)-divB(jx,3,jz))
c      divB(jx,my,jz)=divB(jx,my-1,jz)
c     1   +1*(divB(jx,my-1,jz)-divB(jx,my-2,jz))
c-----f2=0 3p-----
      aa1=5./2.
      bb1=-2.
      cc1=1./2.
      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
     1              +cc1*divB(jx,5,jz)
      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
     1                 +bb1*divB(jx,my-3,jz)
     1                 +cc1*divB(jx,my-4,jz)
      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
     1              +cc1*divB(jx,4,jz)
      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
     1               +bb1*divB(jx,my-2,jz)
     1               +cc1*divB(jx,my-3,jz)
c-----f2=0 4p-----
c      aa1=104./35.
c      bb1=-114./35.
c      cc1=56./35.
c      dd1=-11./35.
c      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
c     1              +cc1*divB(jx,5,jz)+dd1*divB(jx,6,jz)
c      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
c     1                 +bb1*divB(jx,my-3,jz)
c     1                 +cc1*divB(jx,my-4,jz)
c     1                 +dd1*divB(jx,my-5,jz)
c      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
c     1              +cc1*divB(jx,4,jz)+dd1*divB(jx,5,jz)
c      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
c     1               +bb1*divB(jx,my-2,jz)
c     1               +cc1*divB(jx,my-3,jz)
c     1               +dd1*divB(jx,my-4,jz)
c-----f2=0 5p-----
c      aa1=154./45.
c      bb1=-214./45.
c      cc1=156./45.
c      dd1=-61./45.
c      ee1=10./45.
c      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
c     1              +cc1*divB(jx,5,jz)+dd1*divB(jx,6,jz)
c     1              +ee1*divB(jx,7,jz)
c      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
c     1                 +bb1*divB(jx,my-3,jz)
c     1                 +cc1*divB(jx,my-4,jz)
c     1                 +dd1*divB(jx,my-5,jz)
c     1                 +ee1*divB(jx,my-6,jz)
c      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
c     1              +cc1*divB(jx,4,jz)+dd1*divB(jx,5,jz)
c     1              +ee1*divB(jx,6,jz)
c      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
c     1               +bb1*divB(jx,my-2,jz)
c     1               +cc1*divB(jx,my-3,jz)
c     1               +dd1*divB(jx,my-4,jz)
c     1               +ee1*divB(jx,my-5,jz)
c
c-----f3=0 3p-----
c      divB(jx,2,jz)=3*divB(jx,3,jz)-3*divB(jx,4,jz)+
c     1      1*divB(jx,5,jz)
c      divB(jx,my-1,jz)=3*divB(jx,my-2,jz)-3
c     1      *divB(jx,my-3,jz)+1*divB(jx,my-4,jz)
c      divB(jx,1,jz)=3*divB(jx,2,jz)-3*divB(jx,3,jz)+
c     1      1*divB(jx,4,jz)
c      divB(jx,my,jz)=3*divB(jx,my-1,jz)-3
c     1      *divB(jx,my-2,jz)+1*divB(jx,my-3,jz)
c-----f3=0 4p-----
c      aa1=18./5.
c      bb1=-24./5.
c      cc1=14./5.
c      dd1=-3./5.
c      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
c     1              +cc1*divB(jx,5,jz)+dd1*divB(jx,6,jz)
c      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
c     1                 +bb1*divB(jx,my-3,jz)
c     1                 +cc1*divB(jx,my-4,jz)
c     1                 +dd1*divB(jx,my-5,jz)
c      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
c     1              +cc1*divB(jx,4,jz)+dd1*divB(jx,5,jz)
c      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
c     1               +bb1*divB(jx,my-2,jz)
c     1               +cc1*divB(jx,my-3,jz)
c     1               +dd1*divB(jx,my-4,jz)
c-----f3=0 5p-----
c      aa1=71./17.
c      bb1=-118./17.
c      cc1=98./17.
c      dd1=-41./17.
c      ee1=7./17.
c      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
c     1              +cc1*divB(jx,5,jz)+dd1*divB(jx,6,jz)
c     1              +ee1*divB(jx,7,jz)
c      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
c     1                 +bb1*divB(jx,my-3,jz)
c     1                 +cc1*divB(jx,my-4,jz)
c     1                 +dd1*divB(jx,my-5,jz)
c     1                 +ee1*divB(jx,my-6,jz)
c      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
c     1              +cc1*divB(jx,4,jz)+dd1*divB(jx,5,jz)
c     1              +ee1*divB(jx,6,jz)
c      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
c     1               +bb1*divB(jx,my-2,jz)
c     1               +cc1*divB(jx,my-3,jz)
c     1               +dd1*divB(jx,my-4,jz)
c     1               +ee1*divB(jx,my-5,jz)
c-----f3=0 6p-----
c      aa1=232./49.
c      bb1=-461./49.
c      cc1=496./49.
c      dd1=-307./49.
c      ee1=104./49.
c      ff1=-15./49.
c      divB(jx,2,jz)=aa1*divB(jx,3,jz)+bb1*divB(jx,4,jz)
c     1              +cc1*divB(jx,5,jz)+dd1*divB(jx,6,jz)
c     1              +ee1*divB(jx,7,jz)+ff1*divB(jx,8,jz)
c      divB(jx,my-1,jz)=aa1*divB(jx,my-2,jz)
c     1                 +bb1*divB(jx,my-3,jz)
c     1                 +cc1*divB(jx,my-4,jz)
c     1                 +dd1*divB(jx,my-5,jz)
c     1                 +ee1*divB(jx,my-6,jz)
c     1                 +ff1*divB(jx,my-7,jz)
c      divB(jx,1,jz)=aa1*divB(jx,2,jz)+bb1*divB(jx,3,jz)
c     1              +cc1*divB(jx,4,jz)+dd1*divB(jx,5,jz)
c     1              +ee1*divB(jx,6,jz)+ff1*divB(jx,7,jz)
c      divB(jx,my,jz)=aa1*divB(jx,my-1,jz)
c     1               +bb1*divB(jx,my-2,jz)
c     1               +cc1*divB(jx,my-3,jz)
c     1               +dd1*divB(jx,my-4,jz)
c     1               +ee1*divB(jx,my-5,jz)
c     1               +ff1*divB(jx,my-6,jz)
c
c-----f4=0-----
c      divB(jx,2,jz)=4*divB(jx,3,jz)-6*divB(jx,4,jz)+
c     1              4*divB(jx,5,jz)-1*divB(jx,6,jz)
c      divB(jx,my-1,jz)=4*divB(jx,my-2,jz)
c     1                -6*divB(jx,my-3,jz)
c     1                +4*divB(jx,my-4,jz)
c     1                -1*divB(jx,my-5,jz)
c      divB(jx,1,jz)=4*divB(jx,2,jz)-6*divB(jx,3,jz)+
c     1              4*divB(jx,4,jz)-1*divB(jx,5,jz)
c      divB(jx,my,jz)=4*divB(jx,my-1,jz)
c     1              -6*divB(jx,my-2,jz)
c     1              +4*divB(jx,my-3,jz)
c     1              -1*divB(jx,my-4,jz)
c-----f5=0-----
c      divB(jx,2,jz)=5*divB(jx,3,jz)-10*divB(jx,4,jz)+
c     1              10*divB(jx,5,jz)-5*divB(jx,6,jz)+
c     1              1*divB(jx,7,jz)
c      divB(jx,my-1,jz)=5*divB(jx,my-2,jz)
c     1                -10*divB(jx,my-3,jz)
c     1                +10*divB(jx,my-4,jz)
c     1                -5*divB(jx,my-5,jz)
c     1                +1*divB(jx,my-6,jz)
c      divB(jx,1,jz)=5*divB(jx,2,jz)-10*divB(jx,3,jz)+
c     1              10*divB(jx,4,jz)-5*divB(jx,5,jz)+
c     1              1*divB(jx,6,jz)
c      divB(jx,my,jz)=5*divB(jx,my-1,jz)
c     1              -10*divB(jx,my-2,jz)
c     1              +10*divB(jx,my-3,jz)
c     1              -5*divB(jx,my-4,jz)
c     1              +1*divB(jx,my-5,jz)
    3 continue
      else
      do 30 jz=z_first-1,z_last+1
      do 30 jx=x_first-1,x_last+1
      divB(jx,1,jz)=divB(jx,my-1,jz)
      divB(jx,my,jz)=divB(jx,2,jz)
   30 continue
      endif
c
      if(.not.periodz) then
      if(mylow.eq.mpi_proc_null)then
      do 4 jy=1,my
      do 4 jx=x_first-2,x_last+2
c-----f1=0 3p-----
c      divB(jx,jy,4)=(azp(4)*divB(jx,jy,5)+bzp(4)*divB(jx,jy,6)+
c     1      czp(4)*divB(jx,jy,7))/(azp(4)+bzp(4)+czp(4))
c      divB(jx,jy,3)=(azp(3)*divB(jx,jy,4)+bzp(3)*divB(jx,jy,5)+
c     1      czp(3)*divB(jx,jy,6))/(azp(3)+bzp(3)+czp(3))
c-----f1=0 4p-----
c      aa1=48./25.
c      bb1=-36./25.
c      cc1=16./25.
c      dd1=-3./25.
c      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
c     1              +cc1*divB(jx,jy,7)+dd1*divB(jx,jy,8)
c      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
c     1              +cc1*divB(jx,jy,6)+dd1*divB(jx,jy,7)
c-----f1=0 5p-----
c      aa1=300./137.
c      bb1=-300./137.
c      cc1=200./137.
c      dd1=-75./137.
c      ee1=12./137.
c      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
c     1              +cc1*divB(jx,jy,7)+dd1*divB(jx,jy,8)
c     1              +ee1*divB(jx,jy,9)
c      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
c     1              +cc1*divB(jx,jy,6)+dd1*divB(jx,jy,7)
c     1              +ee1*divB(jx,jy,8)
c-----f2=0 2p-----
c      divB(jx,jy,4)=divB(jx,jy,5)
c     1   +1*(divB(jx,jy,5)-divB(jx,jy,6))
c      divB(jx,jy,3)=divB(jx,jy,4)
c     1   +1*(divB(jx,jy,4)-divB(jx,jy,5))
c-----f2=0 3p-----
      aa1=5./2.
      bb1=-2.
      cc1=1./2.
      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
     1              +cc1*divB(jx,jy,7)
      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
     1              +cc1*divB(jx,jy,6)
c-----f2=0 4p-----
c      aa1=104./35.
c      bb1=-114./35.
c      cc1=56./35.
c      dd1=-11./35.
c      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
c     1              +cc1*divB(jx,jy,7)+dd1*divB(jx,jy,8)
c      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
c     1              +cc1*divB(jx,jy,6)+dd1*divB(jx,jy,7)
c-----f2=0 5p-----
c      aa1=154./45.
c      bb1=-214./45.
c      cc1=156./45.
c      dd1=-61./45.
c      ee1=10./45.
c      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
c     1              +cc1*divB(jx,jy,7)+dd1*divB(jx,jy,8)
c     1              +ee1*divB(jx,jy,9)
c      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
c     1              +cc1*divB(jx,jy,6)+dd1*divB(jx,jy,7)
c     1              +ee1*divB(jx,jy,8)
c
c-----f3=0 3p-----
c      divB(jx,jy,4)=3*divB(jx,jy,5)-3*divB(jx,jy,6)+
c     1      1*divB(jx,jy,7)
c      divB(jx,jy,3)=3*divB(jx,jy,4)-3*divB(jx,jy,5)+
c     1      1*divB(jx,jy,6)
c-----f3=0 4p-----
c      aa1=18./5.
c      bb1=-24./5.
c      cc1=14./5.
c      dd1=-3./5.
c      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
c     1              +cc1*divB(jx,jy,7)+dd1*divB(jx,jy,8)
c      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
c     1              +cc1*divB(jx,jy,6)+dd1*divB(jx,jy,7)
c-----f3=0 5p-----
c      aa1=71./17.
c      bb1=-118./17.
c      cc1=98./17.
c      dd1=-41./17.
c      ee1=7./17.
c      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
c     1              +cc1*divB(jx,jy,7)+dd1*divB(jx,jy,8)
c     1              +ee1*divB(jx,jy,9)
c      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
c     1              +cc1*divB(jx,jy,6)+dd1*divB(jx,jy,7)
c     1              +ee1*divB(jx,jy,8)
c-----f3=0 6p-----
c      aa1=232./49.
c      bb1=-461./49.
c      cc1=496./49.
c      dd1=-307./49.
c      ee1=104./49.
c      ff1=-15./49.
c      divB(jx,jy,4)=aa1*divB(jx,jy,5)+bb1*divB(jx,jy,6)
c     1              +cc1*divB(jx,jy,7)+dd1*divB(jx,jy,8)
c     1              +ee1*divB(jx,jy,9)+ff1*divB(jx,jy,10)
c      divB(jx,jy,3)=aa1*divB(jx,jy,4)+bb1*divB(jx,jy,5)
c     1              +cc1*divB(jx,jy,6)+dd1*divB(jx,jy,7)
c     1              +ee1*divB(jx,jy,8)+ff1*divB(jx,jy,9)
c
c-----f4=0-----
c      divB(jx,jy,4)=4*divB(jx,jy,5)-6*divB(jx,jy,6)+
c     1               4*divB(jx,jy,7)-1*divB(jx,jy,8)
c      divB(jx,jy,3)=4*divB(jx,jy,4)-6*divB(jx,jy,5)+
c     1               4*divB(jx,jy,6)-1*divB(jx,jy,7)
c-----f5=0-----
c      divB(jx,jy,4)=5*divB(jx,jy,5)-10*divB(jx,jy,6)+
c     1               10*divB(jx,jy,7)-5*divB(jx,jy,8)+
c     1               1*divB(jx,jy,9)
c      divB(jx,jy,3)=5*divB(jx,jy,4)-10*divB(jx,jy,5)+
c     1               10*divB(jx,jy,6)-5*divB(jx,jy,7)+
c     1               1*divB(jx,jy,8)
    4 continue
      endif
c
      if(myupper.eq.mpi_proc_null)then
      if(symmetryz) then
      do 5 jy=1,my
      do 5 jx=x_first-1,x_last+1
      divB(jx,jy,zlast)=divB(jx,my-jy+1,zlast-2)
    5 continue
      else
      do 6 jy=1,my
      do 6 jx=x_first-2,x_last+2
c-----f1=0 3p-----
c      divB(jx,jy,zlast-1)=(azm(zlast-1)*divB(jx,jy,zlast-2)
c     1  +bzm(zlast-1)*divB(jx,jy,zlast-3)+czm(zlast-1)
c     1  *divB(jx,jy,zlast-4))/(azm(zlast-1)+bzm(zlast-1)+czm(zlast-1))
c      divB(jx,jy,zlast)=(azm(zlast)*divB(jx,jy,zlast-1)+bzm(zlast)
c     1  *divB(jx,jy,zlast-2)+czm(zlast)*divB(jx,jy,zlast-3))
c     1  /(azm(zlast)+bzm(zlast)+czm(zlast))
c-----f1=0 4p-----
c      aa1=48./25.
c      bb1=-36./25.
c      cc1=16./25.
c      dd1=-3./25.
c      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
c     1                    +bb1*divB(jx,jy,zlast-3)
c     1                    +cc1*divB(jx,jy,zlast-4)
c     1                    +dd1*divB(jx,jy,zlast-5)
c      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
c     1                  +bb1*divB(jx,jy,zlast-2)
c     1                  +cc1*divB(jx,jy,zlast-3)
c     1                  +dd1*divB(jx,jy,zlast-4)
c-----f1=0 5p-----
c      aa1=300./137.
c      bb1=-300./137.
c      cc1=200./137.
c      dd1=-75./137.
c      ee1=12./137.
c      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
c     1                    +bb1*divB(jx,jy,zlast-3)
c     1                    +cc1*divB(jx,jy,zlast-4)
c     1                    +dd1*divB(jx,jy,zlast-5)
c     1                    +ee1*divB(jx,jy,zlast-6)
c      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
c     1                  +bb1*divB(jx,jy,zlast-2)
c     1                  +cc1*divB(jx,jy,zlast-3)
c     1                  +dd1*divB(jx,jy,zlast-4)
c     1                  +ee1*divB(jx,jy,zlast-5)
c-----f2=0 2p-----
c      divB(jx,jy,zlast-1)=divB(jx,jy,zlast-2)
c     1   +1*(divB(jx,jy,zlast-2)-divB(jx,jy,zlast-3))
c      divB(jx,jy,zlast)=divB(jx,jy,zlast-1)
c     1   +1*(divB(jx,jy,zlast-1)-divB(jx,jy,zlast-2))
c-----f2=0 3p-----
      aa1=5./2.
      bb1=-2.
      cc1=1./2.
      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
     1                    +bb1*divB(jx,jy,zlast-3)
     1                    +cc1*divB(jx,jy,zlast-4)
      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
     1                  +bb1*divB(jx,jy,zlast-2)
     1                  +cc1*divB(jx,jy,zlast-3)
c-----f2=0 4p-----
c      aa1=104./35.
c      bb1=-114./35.
c      cc1=56./35.
c      dd1=-11./35.
c      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
c     1                    +bb1*divB(jx,jy,zlast-3)
c     1                    +cc1*divB(jx,jy,zlast-4)
c     1                    +dd1*divB(jx,jy,zlast-5)
c      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
c     1                  +bb1*divB(jx,jy,zlast-2)
c     1                  +cc1*divB(jx,jy,zlast-3)
c     1                  +dd1*divB(jx,jy,zlast-4)
c-----f2=0 5p-----
c      aa1=154./45.
c      bb1=-214./45.
c      cc1=156./45.
c      dd1=-61./45.
c      ee1=10./45.
c      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
c     1                    +bb1*divB(jx,jy,zlast-3)
c     1                    +cc1*divB(jx,jy,zlast-4)
c     1                    +dd1*divB(jx,jy,zlast-5)
c     1                    +ee1*divB(jx,jy,zlast-6)
c      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
c     1                  +bb1*divB(jx,jy,zlast-2)
c     1                  +cc1*divB(jx,jy,zlast-3)
c     1                  +dd1*divB(jx,jy,zlast-4)
c     1                  +ee1*divB(jx,jy,zlast-5)
c
c-----f3=0 3p-----
c      divB(jx,jy,zlast-1)=3*divB(jx,jy,zlast-2)
c     1  -3*divB(jx,jy,zlast-3)+1
c     1  *divB(jx,jy,zlast-4)
c      divB(jx,jy,zlast)=3*divB(jx,jy,zlast-1)-3
c     1  *divB(jx,jy,zlast-2)+1*divB(jx,jy,zlast-3)
c-----f3=0 4p-----
c      aa1=18./5.
c      bb1=-24./5.
c      cc1=14./5.
c      dd1=-3./5.
c      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
c     1                    +bb1*divB(jx,jy,zlast-3)
c     1                    +cc1*divB(jx,jy,zlast-4)
c     1                    +dd1*divB(jx,jy,zlast-5)
c      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
c     1                  +bb1*divB(jx,jy,zlast-2)
c     1                  +cc1*divB(jx,jy,zlast-3)
c     1                  +dd1*divB(jx,jy,zlast-4)
c-----f3=0 5p-----
c      aa1=71./17.
c      bb1=-118./17.
c      cc1=98./17.
c      dd1=-41./17.
c      ee1=7./17.
c      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
c     1                    +bb1*divB(jx,jy,zlast-3)
c     1                    +cc1*divB(jx,jy,zlast-4)
c     1                    +dd1*divB(jx,jy,zlast-5)
c     1                    +ee1*divB(jx,jy,zlast-6)
c      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
c     1                  +bb1*divB(jx,jy,zlast-2)
c     1                  +cc1*divB(jx,jy,zlast-3)
c     1                  +dd1*divB(jx,jy,zlast-4)
c     1                  +ee1*divB(jx,jy,zlast-5)
c-----f3=0 6p-----
c      aa1=232./49.
c      bb1=-461./49.
c      cc1=496./49.
c      dd1=-307./49.
c      ee1=104./49.
c      ff1=-15./49.
c      divB(jx,jy,zlast-1)=aa1*divB(jx,jy,zlast-2)
c     1                    +bb1*divB(jx,jy,zlast-3)
c     1                    +cc1*divB(jx,jy,zlast-4)
c     1                    +dd1*divB(jx,jy,zlast-5)
c     1                    +ee1*divB(jx,jy,zlast-6)
c     1                    +ff1*divB(jx,jy,zlast-7)
c      divB(jx,jy,zlast)=aa1*divB(jx,jy,zlast-1)
c     1                  +bb1*divB(jx,jy,zlast-2)
c     1                  +cc1*divB(jx,jy,zlast-3)
c     1                  +dd1*divB(jx,jy,zlast-4)
c     1                  +ee1*divB(jx,jy,zlast-5)
c     1                  +ff1*divB(jx,jy,zlast-6)
c
c-----f4=0-----
c      divB(jx,jy,zlast-1)=4*divB(jx,jy,zlast-2)
c     1                    -6*divB(jx,jy,zlast-3)
c     1                    +4*divB(jx,jy,zlast-4)
c     1                    -1*divB(jx,jy,zlast-5)
c      divB(jx,jy,zlast)=4*divB(jx,jy,zlast-1)
c     1                  -6*divB(jx,jy,zlast-2)
c     1                  +4*divB(jx,jy,zlast-3)
c     1                  -1*divB(jx,jy,zlast-4)
c-----f5=0-----
c      divB(jx,jy,zlast-1)=5*divB(jx,jy,zlast-2)
c     1                    -10*divB(jx,jy,zlast-3)
c     1                    +10*divB(jx,jy,zlast-4)
c     1                    -5*divB(jx,jy,zlast-5)
c     1                    +1*divB(jx,jy,zlast-6)
c      divB(jx,jy,zlast)=5*divB(jx,jy,zlast-1)
c     1                  -10*divB(jx,jy,zlast-2)
c     1                  +10*divB(jx,jy,zlast-3)
c     1                  -5*divB(jx,jy,zlast-4)
c     1                  +1*divB(jx,jy,zlast-5)
    6 continue
      endif
      endif
      else
      if (nrank.eq.(nsize-1)) then
      do  m=1,1
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      wf4(jx,jy,m)=divB(jx,jy,z_last-1)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
      call MPI_SEND( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 2,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.0) then
      call MPI_RECV( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,2,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
c
      do  m=1,1
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      divB(jx,jy,2)=wf4(jx,jy,m)
      wf4(jx,jy,m)=divB(jx,jy,3)
      enddo
      enddo
      enddo
c mpi-------------------------------------------------------------------
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      call MPI_SEND( wf4,(mxpr-2)*my*8,MPI_DOUBLE_PRECISION,
     1     nsize-1,3,MPI_COMM_WORLD, ierror)
      endif
c     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
      if (nrank.eq.(nsize-1)) then
      call MPI_RECV( wf4, (mxpr-2)*my*8, MPI_DOUBLE_PRECISION,
     1     0, 3,MPI_COMM_WORLD, status, ierror)
c mpi-------------------------------------------------------------------
      do  m=1,1
      do  jx=x_first-1,x_last+1
      do  jy=1,my
      divB(jx,jy,z_last)=wf4(jx,jy,m)
      enddo
      enddo
      enddo
      endif
      endif
c
      call message5(divB,mxpr,my,mzpr,myleft,myright,mylow,myupper)
c
      return
      end
c
c
c
      subroutine gridpnt
      include 'mhd01.for'
      dimension  help(2000),xxh(2000),dxh(2000)
c     if(.not.periodx) then
c------The mesh grid divide into 4 parts: [-30,-10],[-10,0],[0,10],[10,80] for the dayside portion
      nnxh=(mx1-1)/2
      nnx=2*nnxh
      mmx=nnx+1
      xxt(1)=xmin0
      dxt(mx1)=0.
      if(.not.halfx) nnx=2*nnxh
      dxi=(xmax-xmin)/float(nnx)
      if(uniformx) then
      nnxh=(mx1-1)/2
      nnx=2*nnxh
      if(.not.halfx) xmin=xmin0
      do 11 jx=1,mmx
      dxt(jx)=dxi
      xxt(jx)=xmin+(jx-1)*dxi
   11 continue
      else
      nnxh=(mx1-1)/8
      nnx1=400
      xxmin=0.
      xxmax=20.
      dxmin=0.02
      width=12.
      dxmax=((xxmax-xxmin)-nnx1*dxmin)/(0.33*nnx1)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      do 12 jx=1,nnx1
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx1*2/3.)/float(nnx1))
      dxt(nnx1-jx+1)=dxh(jx)
   12 continue
      xxmin=0.
      xxmax=10.
      nnx2=300
      dxmin=0.02
      width=8.
      dxmax=((xxmax-xxmin)-nnx2*dxmin)/(0.33*nnx2)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      deltax=dxt(nnx1)-(dxp+dxm*tanh((width-width*nnx2*2/3.)
     1          /float(nnx2)))
      do 13 jx=1,nnx2
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx2*2/3.)/float(nnx2))
      dxt(jx+nnx1)=deltax+dxh(jx)
      dxt(nnx1+2*nnx2-jx+1)=deltax+dxh(jx)
   13 continue
      xxmin=0.
      xxmax=70.
      nnx3=440
      dxmin=0.02
      width=15.
      dxmax=((xxmax-xxmin)-nnx3*dxmin)/(0.33*nnx3)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      deltax=dxt(nnx1+2*nnx2)-(dxp+dxm*tanh((width-width*nnx3*2/3.)
     1          /float(nnx3)))
      do 14 jx=1,nnx3
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx3*2/3.)/float(nnx3))
      dxt(jx+nnx1+2*nnx2)=deltax+dxh(jx)
   14 continue
      xxt(1)=-30.
      do 15 jx=1,mx1-1
      xxt(jx+1)=xxt(jx)+dxt(jx)
   15 continue
      delta=(80.-xxt(mx1))/float(mx1-1)
      do 16 jx=1,mx1-1
      dxt(jx)=dxt(jx)+delta
      xxt(jx+1)=xxt(jx)+dxt(jx)
   16 continue
      endif
      if(nrank.eq.0)then
      open(unit=11,file='gridx.dat',status='unknown',form='formatted')
      write(11,99)(jx,xxt(jx),dxt(jx),jx=1,mx)
   99 format(2(i5,2(1x,f13.5)))
      close(11)
      endif
c
      if(uniformx.and..not.halfx) goto 40
c     if(halfx) then
c     if(mx.ne.mmx) then
c     do 14 jx=1,mmx
c     help(jx)=-xxt(mmx-jx+1)
c  14 continue
c     do 15 jx=1,mmx
c     xxt(jx)=help(jx)
c  15 continue
c     xxt(mx)=-xxt(nnx)
c     do 114 jx=1,mmx
c     help(jx)=dxt(mmx-jx+1)
c 114 continue
c     do 115 jx=1,mmx
c     dxt(jx)=help(jx)
c 115 continue
c     dxt(mx)=dxt(nnx)
c     endif
c     else
c     do 16 jx=1,nnx+1
c     help(jx)=-xxt(nnx-jx+2)
c     help(jx+nnx)=xxt(jx)
c  16 continue
c     do 17 jx=1,2*nnx+1
c     xxt(jx)=help(jx)
c  17 continue
c     do 18 jx=1,nnx+1
c     help(jx)=dxt(nnx-jx+2)
c     help(jx+nnx)=dxt(jx)
c  18 continue
c     do 19 jx=1,2*nnx+1
c     dxt(jx)=help(jx)
c  19 continue
c     endif
c     else
c     dxi=(xmax-xmin)/float(nnx)
c     do 300 jx=1,mx
c     dx(jx)=dxi
c     xxt(jx)=xmin+(jx-1)*dxi
c 300 continue
c     endif
c
   40 nnyh=(my1-1)/2
      nny=2*nnyh
c     if(.not.periody) then
      mmy=nny+1
      yy(1)=ymin
      dy(my1)=0.
      if(.not.halfy) nny=nnyh
      dyi=(ymax-ymin)/float(nny)
      if(uniformy) then
      if(.not.halfy) ymin=-ymax
      do 20 jy=1,mmy
      dy(jy)=dyi
      yy(jy)=ymin+(jy-1)*dyi
   20 continue
      else
      width=5.
      if(dymax.eq.dymin) then
      nyhr=(yhr-ymin)/dymin+1
      nny1=nny+1-nyhr
      dymax=(2.*(ymax-yhr)-(nny+1-nyhr)*dymin)
     1        /float(nny+1-nyhr)
      dyp=0.5*(dymax+dymin)
      dym=0.5*(dymax-dymin)
      do 21 jy=1,nyhr-1
      dy(jy)=dyp+dym*tanh((-width*nny1/2.)/nny1)
      yy(jy+1)=yy(jy)+dy(jy)
   21 continue
      else
      nyhr=2*(yhr-ymin)/(dymax+dymin)+1
      nny1=nny+1-nyhr
      dyp1=0.5*(dymax+dymin)
      dym1=0.5*(dymax-dymin)
      dymax=(2.*(ymax-yhr)-(nny+1-nyhr)*dymin)
     1        /float(nny+1-nyhr)
      dyp=0.5*(dymax+dymin)
      dym=0.5*(dymax-dymin)
      ddy=dyp+dym*tanh((-width*nny1/2.)/nny1)
     1   -(dyp1-dym1*tanh((width*nyhr/2.)/nyhr))
      do 121 jy=1,nyhr-1
      dy(jy)=ddy+dyp1-dym1*tanh((width*(jy+1)-width*nyhr/2.)/nyhr)
      yy(jy+1)=yy(jy)+dy(jy)
  121 continue
      endif
      do 22 jy=nyhr,nny
      dy(jy)=dyp+dym*tanh((width*(jy-nyhr)-width*nny1/2.)/nny1)
      yy(jy+1)=yy(jy)+dy(jy)
   22 continue
      delta=(ymax-yy(nny+1))/float(nny)
      do 23 jy=1,nny
      dy(jy)=dy(jy)+delta
      yy(jy+1)=yy(jy)+dy(jy)
   23 continue
      dy(nny+1)=dy(nny)
      endif
      if(uniformy.and..not.halfy) goto 50
      if(halfy) then
      if(my1.ne.mmy) then
      do 24 jy=1,mmy
      help(jy)=-yy(mmy-jy+1)
   24 continue
      do 25 jy=1,mmy
      yy(jy)=help(jy)
   25 continue
      yy(my1)=-yy(nny)
      do 124 jy=1,mmy
      help(jy)=dy(mmy-jy+1)
  124 continue
      do 125 jy=1,mmy
      dy(jy)=help(jy)
  125 continue
      dy(my1)=dy(nny)
      endif
      else
      do 26 jy=1,nny+1
      help(jy)=-yy(nny-jy+2)
      help(jy+nny)=yy(jy)
   26 continue
      do 27 jy=1,2*nny+1
      yy(jy)=help(jy)
   27 continue
      do 28 jy=1,nny+1
      help(jy)=dy(nny-jy+2)
      help(jy+nny)=dy(jy)
   28 continue
      do 29 jy=1,2*nny+1
      dy(jy)=help(jy)
   29 continue
      endif
c     else
c     dyi=(ymax-ymin)/float(nny)
c     do 30 jy=1,my
c     dy(jy)=dyi
c     yy(jy)=ymin+(jy-1)*dyi
c  30 continue
c     endif
c
   50 nnzh=(mz-1)/2
      nnz=2*nnzh
      mmz=nnz+1
      zzt(1)=zmin
      dzt(mz)=0.
      if(.not.halfz) nnz=nnzh
      dzi=(zmax-zmin)/float(nnz)
      if(uniformz) then
      if(.not.halfz) zmin=zmin0
      do 31 jz=1,mmz
      dzt(jz)=dzi
      zzt(jz)=zmin+(jz-1)*dzi
   31 continue
      else
      width=20.
      if(dzmax.eq.dzmin) then
      nzhr=(zhr-zmin)/dzmin+1
      nnz1=nnz+1-nzhr
      dzmax=(2.*(zmax-zhr)-(nnz+1-nzhr)*dzmin)
     1        /float(nnz+1-nzhr)
      dzp=0.5*(dzmax+dzmin)
      dzm=0.5*(dzmax-dzmin)
      do 32 jz=1,nzhr-1
      dzt(jz)=dzp+dzm*tanh((-width*nnz1/2.)/nnz1)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   32 continue
      else
      nzhr=2*(zhr-zmin)/(dzmax+dzmin)+1
      nnz1=nnz+1-nzhr
      dzp1=0.5*(dzmax+dzmin)
      dzm1=0.5*(dzmax-dzmin)
      dzmax=(2.*(zmax-zhr)-(nnz+1-nzhr)*dzmin)
     1        /float(nnz+1-nzhr)
      dzp=0.5*(dzmax+dzmin)
      dzm=0.5*(dzmax-dzmin)
      ddz=dzp+dzm*tanh((-width*nnz1/2.)/nnz1)
     1   -(dzp1-dzm1*tanh((width*nzhr/2.)/nzhr))
      do 33 jz=1,nzhr-1
      dzt(jz)=ddz+dzp1-dzm1*tanh((width*(jz+1)-width*nzhr/3.)/nzhr)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   33 continue
      endif
      do 34 jz=nzhr,nnz
      dzt(jz)=dzp+dzm*tanh((width*(jz-nzhr)-width*nnz1*2/3.)/nnz1)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   34 continue
      delta=(zmax-zzt(nnz+1))/float(nnz)
      do 35 jz=1,nnz
      dzt(jz)=dzt(jz)+delta
      zzt(jz+1)=zzt(jz)+dzt(jz)
   35 continue
      dzt(nnz+1)=dzt(nnz)
      endif
      if(uniformz.and..not.halfz) goto 60
      if(halfz) then
      if(mz.ne.mmz) then
      do 36 jz=1,mmz
      help(jz)=-zzt(mmz-jz+1)
   36 continue
      do 37 jz=1,mmz
      zzt(jz)=help(jz)
   37 continue
      zzt(mz)=-zzt(nnz)
      do 38 jz=1,mmz
      help(jz)=dzt(mmz-jz+1)
   38 continue
      do 39 jz=1,mmz
      dzt(jz)=help(jz)
   39 continue
      dzt(mz)=dzt(nnz)
      endif
      else
      do 41 jz=1,nnz+1
      help(jz)=-zzt(nnz-jz+2)
      help(jz+nnz)=zzt(jz)
   41 continue
      do 42 jz=1,2*nnz+1
      zzt(jz)=help(jz)
   42 continue
      do 43 jz=1,nnz+1
      help(jz)=dzt(nnz-jz+2)
      help(jz+nnz)=dzt(jz)
   43 continue
      do 44 jz=1,2*nnz+1
      dzt(jz)=help(jz)
   44 continue
      endif
   60 continue
c
      return
      end
c
c
c
      subroutine gridpntnew
      include 'mhd01.for'
      dimension  help(2000),xxh(2000),dxh(2000)
c     if(.not.periodx) then
c------The mesh grid divide into 4 parts: [-30,-10],[-10,0],[0,10],[10,80] for the sun-earth direction
      nnxh=(mx1-1)/2
      nnx=2*nnxh
      mmx=nnx+1
      xxt(1)=xmin0
      dxt(mx1)=0.
      if(.not.halfx) nnx=2*nnxh
      dxi=(xmax-xmin)/float(nnx)
      if(uniformx) then
      nnxh=(mx1-1)/2
      nnx=2*nnxh
      if(.not.halfx) xmin=xmin0
      do 11 jx=1,mmx
      dxt(jx)=dxi
      xxt(jx)=xmin+(jx-1)*dxi
   11 continue
      else
      nnxh=(mx1-1)/8
      nnx1=400
      xxmin=0.
      xxmax=20.
      dxmin=0.02
      width=12.
      dxmax=((xxmax-xxmin)-nnx1*dxmin)/(0.33*nnx1)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      do 12 jx=1,nnx1
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx1*2/3.)/float(nnx1))
      dxt(nnx1-jx+1)=dxh(jx)
   12 continue
      xxmin=0.
      xxmax=10.
      nnx2=300
      dxmin=0.02
      width=8.
      dxmax=((xxmax-xxmin)-nnx2*dxmin)/(0.33*nnx2)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      deltax=dxt(nnx1)-(dxp+dxm*tanh((width-width*nnx2*2/3.)
     1          /float(nnx2)))
      do 13 jx=1,nnx2
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx2*2/3.)/float(nnx2))
      dxt(jx+nnx1)=deltax+dxh(jx)
      dxt(nnx1+2*nnx2-jx+1)=deltax+dxh(jx)
   13 continue
      xxmin=0.
      xxmax=70.
      nnx3=440
      dxmin=0.02
      width=15.
      dxmax=((xxmax-xxmin)-nnx3*dxmin)/(0.33*nnx3)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      deltax=dxt(nnx1+2*nnx2)-(dxp+dxm*tanh((width-width*nnx3*2/3.)
     1          /float(nnx3)))
      do 14 jx=1,nnx3
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx3*2/3.)/float(nnx3))
      dxt(jx+nnx1+2*nnx2)=deltax+dxh(jx)
   14 continue
      xxt(1)=-30.
      do 15 jx=1,mx1-1
      xxt(jx+1)=xxt(jx)+dxt(jx)
   15 continue
      delta=(80.-xxt(mx1))/float(mx1-1)
      do 16 jx=1,mx1-1
      dxt(jx)=dxt(jx)+delta
      xxt(jx+1)=xxt(jx)+dxt(jx)
   16 continue
      endif
      if(nrank.eq.0)then
      open(unit=11,file='gridx.dat',status='unknown',form='formatted')
      write(11,99)(jx,xxt(jx),dxt(jx),jx=1,mx)
   99 format(2(i5,2(1x,f13.5)))
      close(11)
      endif
c
      if(uniformx.and..not.halfx) goto 40
c     if(halfx) then
c     if(mx.ne.mmx) then
c     do 14 jx=1,mmx
c     help(jx)=-xxt(mmx-jx+1)
c  14 continue
c     do 15 jx=1,mmx
c     xxt(jx)=help(jx)
c  15 continue
c     xxt(mx)=-xxt(nnx)
c     do 114 jx=1,mmx
c     help(jx)=dxt(mmx-jx+1)
c 114 continue
c     do 115 jx=1,mmx
c     dxt(jx)=help(jx)
c 115 continue
c     dxt(mx)=dxt(nnx)
c     endif
c     else
c     do 16 jx=1,nnx+1
c     help(jx)=-xxt(nnx-jx+2)
c     help(jx+nnx)=xxt(jx)
c  16 continue
c     do 17 jx=1,2*nnx+1
c     xxt(jx)=help(jx)
c  17 continue
c     do 18 jx=1,nnx+1
c     help(jx)=dxt(nnx-jx+2)
c     help(jx+nnx)=dxt(jx)
c  18 continue
c     do 19 jx=1,2*nnx+1
c     dxt(jx)=help(jx)
c  19 continue
c     endif
c     else
c     dxi=(xmax-xmin)/float(nnx)
c     do 300 jx=1,mx
c     dx(jx)=dxi
c     xxt(jx)=xmin+(jx-1)*dxi
c 300 continue
c     endif
c
   40 nnyh=(my1i-1)/2
      nny=2*nnyh
c     if(.not.periody) then
      mmy=nny+1
      yy(1)=ymin
      dy(my1)=0.
      if(.not.halfy) nny=nnyh
      dyi=(ymax-ymin)/float(nny)
      nnyt=(my1-1)/2
      nnyo=nnyt-nny
      if(uniformy) then
c      if(.not.halfy) ymin=-ymax
      do 20 jy=1,nny+1
      dy(jy)=dyi
      yy(jy)=ymin+(jy-1)*dyi
   20 continue
      else
      width=5.
      if(dymax.eq.dymin) then
      nyhr=(yhr-ymin)/dymin+1
      nny1=nny+1-nyhr
      dymax=(2.*(ymax-yhr)-(nny+1-nyhr)*dymin)
     1        /float(nny+1-nyhr)
      dyp=0.5*(dymax+dymin)
      dym=0.5*(dymax-dymin)
      do 21 jy=1,nyhr-1
      dy(jy)=dyp+dym*tanh((-width*nny1/2.)/nny1)
      yy(jy+1)=yy(jy)+dy(jy)
   21 continue
      else
      nyhr=2*(yhr-ymin)/(dymax+dymin)+1
      nny1=nny+1-nyhr
      dyp1=0.5*(dymax+dymin)
      dym1=0.5*(dymax-dymin)
      dymax=(2.*(ymax-yhr)-(nny+1-nyhr)*dymin)
     1        /float(nny+1-nyhr)
      dyp=0.5*(dymax+dymin)
      dym=0.5*(dymax-dymin)
      ddy=dyp+dym*tanh((-width*nny1/2.)/nny1)
     1   -(dyp1-dym1*tanh((width*nyhr/2.)/nyhr))
      do 121 jy=1,nyhr-1
      dy(jy)=ddy+dyp1-dym1*tanh((width*(jy+1)-width*nyhr/2.)/nyhr)
      yy(jy+1)=yy(jy)+dy(jy)
  121 continue
      endif
      do 22 jy=nyhr,nny
      dy(jy)=dyp+dym*tanh((width*(jy-nyhr)-width*nny1/2.)/nny1)
      yy(jy+1)=yy(jy)+dy(jy)
   22 continue
      delta=(ymax-yy(nny+1))/float(nny)
      do 23 jy=1,nny
      dy(jy)=dy(jy)+delta
      yy(jy+1)=yy(jy)+dy(jy)
   23 continue
      dy(nny+1)=dy(nny)
      endif
cccccccccccccccccccccccccccccc for outer boundary
      dymaxi=dy(nny+1)
      dypo=0.5*(dymaxo+dymaxi)
      dymo=0.5*(dymaxo-dymaxi)
      do 82 jy=1,nnyo
      dy(nny+jy)=dypo+dymo*tanh((5.*jy-5.*nnyo*1/4.)/nnyo)
      yy(nny+1+jy)=yy(nny+jy)+dy(nny+jy)
   82 continue
      deltayo=(ymaxo-yy(nny+1+nnyo))/float(nnyo)
      if(deltayo .le. 0.) print*,'ymaxo should be larger or equal to
     1                            the last yy value, please reduce
     2                            dzmaxo or grid number in z out'
      do 83 jy=1,nnyo
      dy(nny+jy)=dy(nny+jy)+deltayo
      yy(nny+1+jy)=yy(nny+jy)+dy(nny+jy)
   83 continue
      dy(nny+nnyo+1)=dy(nny+nnyo)
cccccccccccccccccccccccccccccc      
c      if(uniformy.and..not.halfy) goto 50
      if(halfy) then
      if(my1.ne.mmy) then
      do 24 jy=1,mmy
      help(jy)=-yy(mmy-jy+1)
   24 continue
      do 25 jy=1,mmy
      yy(jy)=help(jy)
   25 continue
      yy(my1)=-yy(nny)
      do 124 jy=1,mmy
      help(jy)=dy(mmy-jy+1)
  124 continue
      do 125 jy=1,mmy
      dy(jy)=help(jy)
  125 continue
      dy(my1)=dy(nny)
      endif
      else
      do 26 jy=1,nnyt+1
      help(jy)=-yy(nnyt-jy+2)
      help(jy+nnyt)=yy(jy)
   26 continue
      do 27 jy=1,2*nnyt+1
      yy(jy)=help(jy)
   27 continue
      do 28 jy=1,nnyt+1
      help(jy)=dy(nnyt-jy+2)
      help(jy+nnyt)=dy(jy)
   28 continue
      do 29 jy=1,2*nnyt+1
      dy(jy)=help(jy)
   29 continue
      endif
c     else
c     dyi=(ymax-ymin)/float(nny)
c     do 30 jy=1,my
c     dy(jy)=dyi
c     yy(jy)=ymin+(jy-1)*dyi
c  30 continue
c     endif
c
   50 nnzh=(mz1i-1)/2
      nnz=2*nnzh
      mmz=nnz+1
      zzt(1)=zmin
      dzt(mz)=0.
      if(.not.halfz) nnz=nnzh
      dzi=(zmax-zmin)/float(nnz)
      nnzt=(mz1-1)/2
      nnzo=nnzt-nnz
      if(uniformz) then
c      if(.not.halfz) zmin=zmin0
      do 31 jz=1,nnz+1
      dzt(jz)=dzi
      zzt(jz)=zmin+(jz-1)*dzi
   31 continue
      else
      width=20.
      if(dzmax.eq.dzmin) then
      nzhr=(zhr-zmin)/dzmin+1
      nnz1=nnz+1-nzhr
      dzmax=(2.*(zmax-zhr)-(nnz+1-nzhr)*dzmin)
     1        /float(nnz+1-nzhr)
      dzp=0.5*(dzmax+dzmin)
      dzm=0.5*(dzmax-dzmin)
      do 32 jz=1,nzhr-1
      dzt(jz)=dzp+dzm*tanh((-width*nnz1/2.)/nnz1)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   32 continue
      else
      nzhr=2*(zhr-zmin)/(dzmax+dzmin)+1
      nnz1=nnz+1-nzhr
      dzp1=0.5*(dzmax+dzmin)
      dzm1=0.5*(dzmax-dzmin)
      dzmax=(2.*(zmax-zhr)-(nnz+1-nzhr)*dzmin)
     1        /float(nnz+1-nzhr)
      dzp=0.5*(dzmax+dzmin)
      dzm=0.5*(dzmax-dzmin)
      ddz=dzp+dzm*tanh((-width*nnz1/2.)/nnz1)
     1   -(dzp1-dzm1*tanh((width*nzhr/2.)/nzhr))
      do 33 jz=1,nzhr-1
      dzt(jz)=ddz+dzp1-dzm1*tanh((width*(jz+1)-width*nzhr/3.)/nzhr)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   33 continue
      endif
      do 34 jz=nzhr,nnz
      dzt(jz)=dzp+dzm*tanh((width*(jz-nzhr)-width*nnz1*2/3.)/nnz1)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   34 continue
      delta=(zmax-zzt(nnz+1))/float(nnz)
      do 35 jz=1,nnz
      dzt(jz)=dzt(jz)+delta
      zzt(jz+1)=zzt(jz)+dzt(jz)
   35 continue
      dzt(nnz+1)=dzt(nnz)
      endif
cccccccccccccccccccccccccccccc for outer boundary
      dzmaxi=dzt(nnz+1)
      dzpo=0.5*(dzmaxo+dzmaxi)
      dzmo=0.5*(dzmaxo-dzmaxi)
      do 80 jz=1,nnzo
      dzt(nnz+jz)=dzpo+dzmo*tanh((5.*jz-5.*nnzo*1/4.)/nnzo)
      zzt(nnz+1+jz)=zzt(nnz+jz)+dzt(nnz+jz)
   80 continue
      deltazo=(zmaxo-zzt(nnz+1+nnzo))/float(nnzo)
      if(deltazo .le. 0.) print*,'zmaxo should be larger or equal to
     1                            the last zzt value please reduce
     2                            dymaxo or grid number in y out'
      do 81 jz=1,nnzo
      dzt(nnz+jz)=dzt(nnz+jz)+deltazo
      zzt(nnz+1+jz)=zzt(nnz+jz)+dzt(nnz+jz)
   81 continue
      dzt(nnz+nnzo+1)=dzt(nnz+nnzo)
cccccccccccccccccccccccccccccc      
c      if(uniformz.and..not.halfz) goto 60
      if(halfz) then
      if(mz.ne.mmz) then
      do 36 jz=1,mmz
      help(jz)=-zzt(mmz-jz+1)
   36 continue
      do 37 jz=1,mmz
      zzt(jz)=help(jz)
   37 continue
      zzt(mz)=-zzt(nnz)
      do 38 jz=1,mmz
      help(jz)=dzt(mmz-jz+1)
   38 continue
      do 39 jz=1,mmz
      dzt(jz)=help(jz)
   39 continue
      dzt(mz)=dzt(nnz)
      endif
      else
      do 41 jz=1,nnzt+1
      help(jz)=-zzt(nnzt-jz+2)
      help(jz+nnzt)=zzt(jz)
   41 continue
      do 42 jz=1,2*nnzt+1
      zzt(jz)=help(jz)
   42 continue
      do 43 jz=1,nnzt+1
      help(jz)=dzt(nnzt-jz+2)
      help(jz+nnzt)=dzt(jz)
   43 continue
      do 44 jz=1,2*nnzt+1
      dzt(jz)=help(jz)
   44 continue
      endif
   60 continue
c
      return
      end
c
c
c
      subroutine gridpntnew2
      include 'mhd01.for'
      dimension  help(2000),xxh(2000),dxh(2000)
c     if(.not.periodx) then
c------The mesh grid divide into 4 parts: [-30,-10],[-10,0],[0,10],[10,80] for the dayside portion
      nnxh=(mx1-1)/2
      nnx=2*nnxh
      mmx=nnx+1
      xxt(1)=xmin0
      dxt(mx1)=0.
      if(.not.halfx) nnx=2*nnxh
      dxi=(xmax-xmin)/float(nnx)
      if(uniformx) then
      nnxh=(mx1-1)/2
      nnx=2*nnxh
      if(.not.halfx) xmin=xmin0
      do 11 jx=1,mmx
      dxt(jx)=dxi
      xxt(jx)=xmin+(jx-1)*dxi
   11 continue
      else
      nnxh=(mx1-1)/8
      nnx1=400
      xxmin=0.
      xxmax=20.
      dxmin=0.02
      width=12.
      dxmax=((xxmax-xxmin)-nnx1*dxmin)/(0.33*nnx1)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      do 12 jx=1,nnx1
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx1*2/3.)/float(nnx1))
      dxt(nnx1-jx+1)=dxh(jx)
   12 continue
      xxmin=0.
      xxmax=10.
      nnx2=300
      dxmin=0.02
      width=8.
      dxmax=((xxmax-xxmin)-nnx2*dxmin)/(0.33*nnx2)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      deltax=dxt(nnx1)-(dxp+dxm*tanh((width-width*nnx2*2/3.)
     1          /float(nnx2)))
      do 13 jx=1,nnx2
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx2*2/3.)/float(nnx2))
      dxt(jx+nnx1)=deltax+dxh(jx)
      dxt(nnx1+2*nnx2-jx+1)=deltax+dxh(jx)
   13 continue
      xxmin=0.
      xxmax=70.
      nnx3=440
      dxmin=0.02
      width=15.
      dxmax=((xxmax-xxmin)-nnx3*dxmin)/(0.33*nnx3)
      dxp=0.5*(dxmax+dxmin)
      dxm=0.5*(dxmax-dxmin)
      deltax=dxt(nnx1+2*nnx2)-(dxp+dxm*tanh((width-width*nnx3*2/3.)
     1          /float(nnx3)))
      do 14 jx=1,nnx3
      dxh(jx)=dxp+dxm*tanh((width*jx-width*nnx3*2/3.)/float(nnx3))
      dxt(jx+nnx1+2*nnx2)=deltax+dxh(jx)
   14 continue
      xxt(1)=-30.
      do 15 jx=1,mx1-1
      xxt(jx+1)=xxt(jx)+dxt(jx)
   15 continue
      delta=(80.-xxt(mx1))/float(mx1-1)
      do 16 jx=1,mx1-1
      dxt(jx)=dxt(jx)+delta
      xxt(jx+1)=xxt(jx)+dxt(jx)
   16 continue
      endif
      if(nrank.eq.0)then
      open(unit=11,file='gridx.dat',status='unknown',form='formatted')
      write(11,99)(jx,xxt(jx),dxt(jx),jx=1,mx)
   99 format(2(i5,2(1x,f13.5)))
      close(11)
      endif
c
      if(uniformx.and..not.halfx) goto 40
c     if(halfx) then
c     if(mx.ne.mmx) then
c     do 14 jx=1,mmx
c     help(jx)=-xxt(mmx-jx+1)
c  14 continue
c     do 15 jx=1,mmx
c     xxt(jx)=help(jx)
c  15 continue
c     xxt(mx)=-xxt(nnx)
c     do 114 jx=1,mmx
c     help(jx)=dxt(mmx-jx+1)
c 114 continue
c     do 115 jx=1,mmx
c     dxt(jx)=help(jx)
c 115 continue
c     dxt(mx)=dxt(nnx)
c     endif
c     else
c     do 16 jx=1,nnx+1
c     help(jx)=-xxt(nnx-jx+2)
c     help(jx+nnx)=xxt(jx)
c  16 continue
c     do 17 jx=1,2*nnx+1
c     xxt(jx)=help(jx)
c  17 continue
c     do 18 jx=1,nnx+1
c     help(jx)=dxt(nnx-jx+2)
c     help(jx+nnx)=dxt(jx)
c  18 continue
c     do 19 jx=1,2*nnx+1
c     dxt(jx)=help(jx)
c  19 continue
c     endif
c     else
c     dxi=(xmax-xmin)/float(nnx)
c     do 300 jx=1,mx
c     dx(jx)=dxi
c     xxt(jx)=xmin+(jx-1)*dxi
c 300 continue
c     endif
c
   40 nnyh=(my1-1)/2
      nny=2*nnyh
c     if(.not.periody) then
      mmy=nny+1
      yy(1)=ymin
      dy(my1)=0.
      if(.not.halfy) nny=nnyh
      dyi=(ymax-ymin)/float(nny)
      if(uniformy) then
      if(.not.halfy) ymin=-ymax
      do 20 jy=1,mmy
      dy(jy)=dyi
      yy(jy)=ymin+(jy-1)*dyi
   20 continue
      else
      width=10.!5.
      if(dymax.eq.dymin) then
      nyhr=(yhr-ymin)/dymin+1
      nny1=nny+1-nyhr
      dymax=(2.*(ymax-yhr)-(nny+1-nyhr)*dymin)
     1        /float(nny+1-nyhr)
      dyp=0.5*(dymax+dymin)
      dym=0.5*(dymax-dymin)
      do 21 jy=1,nyhr-1
      dy(jy)=dyp+dym*tanh((-width*nny1/2.)/nny1)
      yy(jy+1)=yy(jy)+dy(jy)
   21 continue
      else
      nyhr=2*(yhr-ymin)/(dymax+dymin)+1
      nny1=nny+1-nyhr
      dyp1=0.5*(dymax+dymin)
      dym1=0.5*(dymax-dymin)
      dymax=(2.*(ymax-yhr)-(nny+1-nyhr)*dymin)
     1        /float(nny+1-nyhr)
      dyp=0.5*(dymax+dymin)
      dym=0.5*(dymax-dymin)
      ddy=dyp+dym*tanh((-width*nny1/2.)/nny1)
     1   -(dyp1-dym1*tanh((width*nyhr/2.)/nyhr))
      do 121 jy=1,nyhr-1
c      dy(jy)=ddy+dyp1-dym1*tanh((width*(jy+1)-width*nyhr/2.)/nyhr)
      dy(jy)=ddy+dyp1-dym1*tanh((width*(jy+1)-width*nyhr/3.)/nyhr)
      yy(jy+1)=yy(jy)+dy(jy)
  121 continue
      endif
      do 22 jy=nyhr,nny
c      dy(jy)=dyp+dym*tanh((width*(jy-nyhr)-width*nny1/2.)/nny1)
      dy(jy)=dyp+dym*tanh((width*(jy-nyhr)-width*nny1*2/3.)/nny1)
      yy(jy+1)=yy(jy)+dy(jy)
   22 continue
      delta=(ymax-yy(nny+1))/float(nny)
      do 23 jy=1,nny
      dy(jy)=dy(jy)+delta
      yy(jy+1)=yy(jy)+dy(jy)
   23 continue
      dy(nny+1)=dy(nny)
      endif
      if(uniformy.and..not.halfy) goto 50
      if(halfy) then
      if(my1.ne.mmy) then
      do 24 jy=1,mmy
      help(jy)=-yy(mmy-jy+1)
   24 continue
      do 25 jy=1,mmy
      yy(jy)=help(jy)
   25 continue
      yy(my1)=-yy(nny)
      do 124 jy=1,mmy
      help(jy)=dy(mmy-jy+1)
  124 continue
      do 125 jy=1,mmy
      dy(jy)=help(jy)
  125 continue
      dy(my1)=dy(nny)
      endif
      else
      do 26 jy=1,nny+1
      help(jy)=-yy(nny-jy+2)
      help(jy+nny)=yy(jy)
   26 continue
      do 27 jy=1,2*nny+1
      yy(jy)=help(jy)
   27 continue
      do 28 jy=1,nny+1
      help(jy)=dy(nny-jy+2)
      help(jy+nny)=dy(jy)
   28 continue
      do 29 jy=1,2*nny+1
      dy(jy)=help(jy)
   29 continue
      endif
c     else
c     dyi=(ymax-ymin)/float(nny)
c     do 30 jy=1,my
c     dy(jy)=dyi
c     yy(jy)=ymin+(jy-1)*dyi
c  30 continue
c     endif
c
   50 nnzh=(mz-1)/2
      nnz=2*nnzh
      mmz=nnz+1
      zzt(1)=zmin
      dzt(mz)=0.
      if(.not.halfz) nnz=nnzh
      dzi=(zmax-zmin)/float(nnz)
      if(uniformz) then
      if(.not.halfz) zmin=zmin0
      do 31 jz=1,mmz
      dzt(jz)=dzi
      zzt(jz)=zmin+(jz-1)*dzi
   31 continue
      else
      width=10.!20.
      if(dzmax.eq.dzmin) then
      nzhr=(zhr-zmin)/dzmin+1
      nnz1=nnz+1-nzhr
      dzmax=(2.*(zmax-zhr)-(nnz+1-nzhr)*dzmin)
     1        /float(nnz+1-nzhr)
      dzp=0.5*(dzmax+dzmin)
      dzm=0.5*(dzmax-dzmin)
      do 32 jz=1,nzhr-1
      dzt(jz)=dzp+dzm*tanh((-width*nnz1/2.)/nnz1)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   32 continue
      else
      nzhr=2*(zhr-zmin)/(dzmax+dzmin)+1
      nnz1=nnz+1-nzhr
      dzp1=0.5*(dzmax+dzmin)
      dzm1=0.5*(dzmax-dzmin)
      dzmax=(2.*(zmax-zhr)-(nnz+1-nzhr)*dzmin)
     1        /float(nnz+1-nzhr)
      dzp=0.5*(dzmax+dzmin)
      dzm=0.5*(dzmax-dzmin)
      ddz=dzp+dzm*tanh((-width*nnz1/2.)/nnz1)
     1   -(dzp1-dzm1*tanh((width*nzhr/2.)/nzhr))
      do 33 jz=1,nzhr-1
      dzt(jz)=ddz+dzp1-dzm1*tanh((width*(jz+1)-width*nzhr/3.)/nzhr)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   33 continue
      endif
      do 34 jz=nzhr,nnz
      dzt(jz)=dzp+dzm*tanh((width*(jz-nzhr)-width*nnz1*2/3.)/nnz1)
      zzt(jz+1)=zzt(jz)+dzt(jz)
   34 continue
      delta=(zmax-zzt(nnz+1))/float(nnz)
      do 35 jz=1,nnz
      dzt(jz)=dzt(jz)+delta
      zzt(jz+1)=zzt(jz)+dzt(jz)
   35 continue
      dzt(nnz+1)=dzt(nnz)
      endif
      if(uniformz.and..not.halfz) goto 60
      if(halfz) then
      if(mz.ne.mmz) then
      do 36 jz=1,mmz
      help(jz)=-zzt(mmz-jz+1)
   36 continue
      do 37 jz=1,mmz
      zzt(jz)=help(jz)
   37 continue
      zzt(mz)=-zzt(nnz)
      do 38 jz=1,mmz
      help(jz)=dzt(mmz-jz+1)
   38 continue
      do 39 jz=1,mmz
      dzt(jz)=help(jz)
   39 continue
      dzt(mz)=dzt(nnz)
      endif
      else
      do 41 jz=1,nnz+1
      help(jz)=-zzt(nnz-jz+2)
      help(jz+nnz)=zzt(jz)
   41 continue
      do 42 jz=1,2*nnz+1
      zzt(jz)=help(jz)
   42 continue
      do 43 jz=1,nnz+1
      help(jz)=dzt(nnz-jz+2)
      help(jz+nnz)=dzt(jz)
   43 continue
      do 44 jz=1,2*nnz+1
      dzt(jz)=help(jz)
   44 continue
      endif
   60 continue
c
      return
      end
c
c
c
      subroutine setdt
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      integer jxoo,jyoo,jzoo
      real*8 vxoo,vyoo,vzoo,vpxoo,vpyoo,vpzoo,vaoo,csoo,vfoo
      real*8 x1oo,x2oo,x3oo,x4oo,x5oo,x6oo,x7oo,x8oo
      dtt=100.
      dttold=dtt
      jxoo=3
      jyoo=1
      jzoo=3
      call foreta(time)
      call formu(time)
      do 1 jz=zfirst,z_last
      do 1 jy=1,my-1
      do 1 jx=xfirst,x_last
      if(rr1(jx,jy,jz).ge.1.) then
      vx=x(jx,jy,jz,2)/x(jx,jy,jz,1)
      vy=x(jx,jy,jz,3)/x(jx,jy,jz,1)
      vz=x(jx,jy,jz,4)/x(jx,jy,jz,1)
      va2=(x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2)
     1    /x(jx,jy,jz,1)
      cs2=gamma*x(jx,jy,jz,8)/x(jx,jy,jz,1)
      vpx=abs(vx)+sqrt(abs(va2+cs2))
      vpy=abs(vy)+sqrt(abs(va2+cs2))
      vpz=abs(vz)+sqrt(abs(va2+cs2))
c     dtx=abs(dx(jx))/(vpx/cfl+
c    1   (etan(jx,jy,jz)+fmu0)/dx(jx)/2.)
c      dty=abs(dy(jy))/(vpy/cfl+
c    1   (etan(jx,jy,jz)+fmu0)/dy(jy)/2.)
c      dtz=abs(dz(jz))/(vpz/cfl+
c    1   (etan(jx,jy,jz)+fmu0)/dz(jz)/2.)
      dtx=abs(dx(jx))/(vpx/cfl+
     1   (etan(jx,jy,jz)+fmu(jx,jy,jz))/dx(jx)/2.)
      dty=abs(dy(jy))/(vpy/cfl+
     1   (etan(jx,jy,jz)+fmu(jx,jy,jz))/dy(jy)/2.)
      dtz=abs(dz(jz))/(vpz/cfl+
     1   (etan(jx,jy,jz)+fmu(jx,jy,jz))/dz(jz)/2.)
      dt1=dmin1(dtx,dty)
      dt2=dmin1(dtz,dt1)
      dtt=dmin1(dtt,dt2)
      if(time .ge. 0.0) then
      if(dttold .gt. dtt) then
          dttold=dtt
          jxoo=jx
          jyoo=jy
          jzoo=jz
          vxoo=vx
          vyoo=vy
          vzoo=vz
          vpxoo=vpx
          vpyoo=vpy
          vpzoo=vpz
          vaoo=va2
          csoo=cs2
          vfoo=sqrt(abs(va2+cs2))
          x1oo=x(jx,jy,jz,1)
          x2oo=x(jx,jy,jz,2)
          x3oo=x(jx,jy,jz,3)
          x4oo=x(jx,jy,jz,4)
          x5oo=x(jx,jy,jz,5)
          x6oo=x(jx,jy,jz,6)
          x7oo=x(jx,jy,jz,7)
          x8oo=x(jx,jy,jz,8)
      endif
      endif
      endif
    1 continue
      
      if(dtstatus .eq. 0) then
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,'nrank,dtt,jxoo,jyoo,jzoo='
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,nrank,dtt,jxoo,jyoo,jzoo
      
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,'nrank,dtt,xxoo,yyoo,zzoo='
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,nrank,dtt,xx(jxoo),yy(jyoo),zz(jzoo)
c      if(time .ge. 0.0 .and. time .le. 10002.0) 
c     1  print*,'nrank,dtt,fmu,etan='
c      if(time .ge. 0.0 .and. time .le. 10002.0) 
c     1  print*,nrank,dtt,fmu(jxoo,jyoo,jzoo),etan(jxoo,jyoo,jzoo)
      
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,'nrank,dtt,vxoo,vyoo,vzoo='
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,nrank,dtt,vxoo,vyoo,vzoo
      
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,'nrank,dtt,vaoo,csoo,vfoo='
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,nrank,dtt,vaoo,csoo,vfoo

      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,'nrank,dtt,vpxoo,vpyoo,vpzoo='
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,nrank,dtt,vpxoo,vpyoo,vpzoo
      
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,'nrank,dtt,x1,x2,x3,x4='
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,nrank,dtt,x1oo,x2oo,x3oo,x4oo
      
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,'nrank,dtt,x5,x6,x7,x8='
      if(time .ge. 0.0 .and. time .le. 10002.0) 
     1  print*,nrank,dtt,x5oo,x6oo,x7oo,x8oo
      endif
      
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_allreduce(dtt,dt,1,mpi_double_precision,mpi_min,
     1     mpi_comm_world,ierr)
      return
      end
c
c
c
      subroutine readin
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c     dimension xtt(mx,my,mz,8),prtt(mx,my,mz)
      character*9 contin
      character*3 cn
      do j=0,nsize-1
      contin='mhd'//cn(nst)//cn(j)
      if(nrank.eq.j)then
      open(unit=8,file=contin,status='unknown',form='unformatted')
      read(8)ncase,nstep,time,nst
      read(8)x
      close(8)
      endif
      enddo
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call message4(x,mxpr,my,mzpr,8,myleft,myright,mylow,myupper)
c     call message3(pr,mxpr,my,mzpr,myleft,myright,mylow,myupper)
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      return
      end
c
c
c
      subroutine recrd
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c     dimension xtt(mx,my,mz,8),prtt(mx,my,mz)
c     dimension work4(mxpr,my,mzpr,nsize)
      character*9 output
      character*3 cn
c
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do j=0,nsize-1
      if(nrank.eq.j)then
      output='mhd'//cn(nst)//cn(j)
      open(unit=8,file=output,status='unknown',form='unformatted')
      write(8)ncase,nstep,time,nst
      write(8)x
      close(8)
      endif
      enddo
      return
      end


      subroutine recrdint
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c     dimension xtt(mx,my,mz,8),prtt(mx,my,mz)
c     dimension work4(mxpr,my,mzpr,nsize)
      character*8 output
c
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      
      output='test.dat'
      open(unit=8,file=output,status='unknown',form='formatted')
      write(8,*)1
      write(8,30)((x(jx,101,jz,5),jx=1,mxpr),jz=1,mzpr)
      write(8,*)2
      write(8,30)((x(jx,101,jz,6),jx=1,mxpr),jz=1,mzpr)
      write(8,*)3
      write(8,30)((x(jx,101,jz,7),jx=1,mxpr),jz=1,mzpr)

   30 format(5(1x,e14.5E3))     
      close(8)
      
      return
      end
c
c
c
      subroutine recrdpe
c     to recrd data file by each procs, by ji liang 2019
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
c     dimension xtt(mx,my,mz,8),prtt(mx,my,mz)
c     dimension work4(mxpr,my,mzpr,nsize)
      character*9 output
      character*3 cn
      dimension bl(mxpr,my,mzpr,8)
c
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do j=0,nsize-1,50
      if(nrank.eq.j)then
      do 3 i=1,8
      do 3 jx=xfirst,xlast
      do 3 jy=1,my
      do 3 jz=zfirst,zlast
      if((i.ge.2) .and. (i.le.4)) then
      bl(jx,jy,jz,i)=x(jx,jy,jz,i)/x(jx,jy,jz,1)
      else
      bl(jx,jy,jz,i)=x(jx,jy,jz,i)
      endif
    3 continue
c      output='mhdpe'//cn(nstnew)//cn(j)//'.dat'
      open(unit=8,file='mhdpe'//cn(nstnew)//cn(j)//'.dat'
     1     ,status='unknown',form='formatted')
      write(8,*)3,xlast-xfirst+1,my,zlast-zfirst+1
      write(8,*)time
      write(8,30)(xx(jx),jx=xfirst,xlast),(yy(jy),jy=1,my)
     1          ,(zz(jz),jz=zfirst,zlast)
      write(8,*)1
      write(8,30)(((bl(jx,jy,jz,1),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)2
      write(8,30)(((bl(jx,jy,jz,2),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)3
      write(8,30)(((bl(jx,jy,jz,3),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)4
      write(8,30)(((bl(jx,jy,jz,4),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)5
      write(8,30)(((bl(jx,jy,jz,5),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)6
      write(8,30)(((bl(jx,jy,jz,6),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)7
      write(8,30)(((bl(jx,jy,jz,7),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)8
      write(8,30)(((bl(jx,jy,jz,8),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)9
      write(8,30)(((etan(jx,jy,jz),jx=xfirst,xlast),jy=1,my)
     1                            ,jz=zfirst,zlast)
      write(8,*)10
      write(8,30)(((fmu(jx,jy,jz),jx=xfirst,xlast),jy=1,my)
     1                           ,jz=zfirst,zlast)
   30 format(5(1x,e14.5E3))
      close(8)
      endif
      enddo
      nstnew=nstnew+1
      return
      end
c
c
c
      subroutine recrd2
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      character*8 output
      character*3 cn
      dimension bl(mx,mz,12)
      dimension work200(mxpr,mzpr,nsize),w0tt(mx,mz)
      dimension xl1(mx,mz),zl1(mx,mz)
c     xxend=-10.+0.08*(nst-1)*dtime
c     if(xxend.gt.xmax) xxend=xmax
c     dx0=(xxend-xmin)/float(kx-1)
c     do 1 jx=1,kx
c     xl(jx)=xmin0+float(jx-1)*dx0
c   1 continue
c     zmax0=zmax
c     dz0=(zmax0-zmin0)/float(kz-1)
c     do 2 jz=1,kz
c     zl(jz)=zmin0+float(jz-1)*dz0
c   2 continue
      nyh=my/2+1
      do 20 jx=1,mx
      do 20 jz=1,mz
      xl1(jx,jz)=xxt(jx)
      zl1(jx,jz)=zzt(jz)
   20 continue
c
c     output=cn(nrank)//'.dat'
c     open(unit=8,file=output,status='unknown',form='formatted')
c     do jz=zfirst,zlast
c     write(8,40)(x(jx,nyh,jz,2),jx=1,mxpr)
c     enddo
c  40 format(32(1x,e11.5))
c     close(8)
      do jm=1,8
      do jx=1,mxpr
      do jz=1,mzpr
      work2(jx,jz)=x(jx,nyh,jz,jm)
      enddo
      enddo
      call mpi_gather(work2(1,1),(mxpr)*(mzpr),mpi_double_precision,
     1     work200(1,1,1),(mxpr)*(mzpr),mpi_double_precision,0,
     2     mpi_comm_world,ierr)
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(nrank.eq.0)then
      do i=1,nsize
      if(i.le.(nsize-nprx))then
      if(mod((i),nprx).ne.0)then
      do jz=zfirst,zlast
      do jx=xfirst,xlast
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      else
      do jz=zfirst,zlast
      do jx=xfirst,xlast+mod(mx,nprx)
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      endif
      else
      if(mod((i),nprx).ne.0)then
      do jz=zfirst,zlast+mod(mz,nprz)
      do jx=xfirst,xlast
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      else
      do jz=zfirst,zlast+mod(mz,nprz)
      do jx=xfirst,xlast+mod(mx,nprx)
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      endif
      endif
      enddo
      do 1 jx=1,mx
      do 1 jz=1,mz
      bl(jx,jz,jm)=w0tt(jx,jz)
    1 continue
c     call extrap(w0tt,xxt,zzt,mx,mz,bl(1,1,jm),xl,zl,kx,kz)
      endif
      enddo
c
      do jm=1,3
      do jx=1,mxpr
      do jz=1,mzpr
      work2(jx,jz)=cur(jx,nyh,jz,jm)
      enddo
      enddo
      call mpi_gather(work2(1,1),(mxpr)*(mzpr),mpi_double_precision,
     1     work200(1,1,1),(mxpr)*(mzpr),mpi_double_precision,0,
     2     mpi_comm_world,ierr)
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(nrank.eq.0)then
      do i=1,nsize
      if(i.le.(nsize-nprx))then
      if(mod((i),nprx).ne.0)then
      do jz=zfirst,zlast
      do jx=xfirst,xlast
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      else
      do jz=zfirst,zlast
      do jx=xfirst,xlast+mod(mx,nprx)
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      endif
      else
      if(mod((i),nprx).ne.0)then
      do jz=zfirst,zlast+mod(mz,nprz)
      do jx=xfirst,xlast
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      else
      do jz=zfirst,zlast+mod(mz,nprz)
      do jx=xfirst,xlast+mod(mx,nprx)
      w0tt(mod((i-1),nprx)*(mxpr-4)+jx-2,int((i-1)/nprx)*(mzpr-4)+jz-2)=
     1     work200(jx,jz,i)
      enddo
      enddo
      endif
      endif
      enddo
      do 2 jx=1,mx
      do 2 jz=1,mz
      bl(jx,jz,jm+8)=w0tt(jx,jz)
    2 continue
c     call extrap(w0tt,xxt,zzt,mx,mz,bl(1,1,jm+8),xl,zl,kx,kz)
      endif
      enddo
c
      do 40 i=2,4
      do 40 jx=1,mx
      do 40 jz=1,mz
      bl(jx,jz,i)=bl(jx,jz,i)/bl(jx,jz,1)
   40 continue
c
      if(nrank.eq.0)then
c     output='m2d'//cn(nst)//'.dat'
      open(unit=8,file='./data/m2d'//cn(nst)//'.dat',
     1 status='unknown',form='formatted')
c     write(8,10)'Title="dipole field"'
c     write(8,11)'Variables="X","Z","rho","Vx","Vy","Vz","Bx"
c    1,"By","Bz","P","Jx","Jy","Jz"'
c     write(8,12)'Zone  I=1441  J=801  F=POINT'
c     write(8,30)((xl1(jx,jz),zl1(jx,jz),(bl(jx,jz,i),i=1,11)
c    1,jx=1,mx),jz=1,mz)
c     revised by ji liang in 2018
      write(8,*)2,mx,mz
      write(8,*)time
      write(8,30)xxt,zzt
      write(8,*)1
      write(8,30)((bl(jx,jz,1),jx=1,mx),jz=1,mz)
      write(8,*)2
      write(8,30)((bl(jx,jz,2),jx=1,mx),jz=1,mz)
      write(8,*)3
      write(8,30)((bl(jx,jz,3),jx=1,mx),jz=1,mz)
      write(8,*)4
      write(8,30)((bl(jx,jz,4),jx=1,mx),jz=1,mz)
      write(8,*)5
      write(8,30)((bl(jx,jz,5),jx=1,mx),jz=1,mz)
      write(8,*)6
      write(8,30)((bl(jx,jz,6),jx=1,mx),jz=1,mz)
      write(8,*)7
      write(8,30)((bl(jx,jz,7),jx=1,mx),jz=1,mz)
      write(8,*)8
      write(8,30)((bl(jx,jz,8),jx=1,mx),jz=1,mz)
      write(8,*)9
      write(8,30)((bl(jx,jz,9),jx=1,mx),jz=1,mz)
      write(8,*)10
      write(8,30)((bl(jx,jz,10),jx=1,mx),jz=1,mz)
      write(8,*)11
      write(8,30)((bl(jx,jz,11),jx=1,mx),jz=1,mz)
c     write(8,30)(((bl(jx,jz,ii),jx=1,mx),jz=1,mz),ii=1,11)
c     do 4 ii=1,11
c     write(8,30)bl(:,:,ii)
c   4 continue
   30 format(5(1x,e14.5E3))
c  30 format(1x,e13.5)
   10 format(1x,A20)
   11 format(1x,A72)
   12 format(1x,A28)
      close(8)
      endif
c
      return
      end
c
c
c
      subroutine recrd3
c     parameter(kx=201,ky=61,kz=201)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      character*14 output
      character*3 cn
      dimension bl(mx,my,mz,11)
      dimension xl1(mx,my,mz),yl1(mx,my,mz),zl1(mx,my,mz)
      dimension work4(mxpr,my,mzpr,nsize)
c
c     xmax0=30.
c     dx0=(xmax0-xmin0)/float(kx-1)
c     do 1 jx=1,kx
c     xl(jx)=xmin0+float(jx-1)*dx0
c   1 continue
c     dy0=(ymax-ymin)/float(ky-1)
c     do 2 jy=1,ky
c     yl(jy)=ymin+float(jy-1)*dy0
c   2 continue
c     zmax0=-zmin0
c     if(halfz) zmax0=0.
c     dz0=(zmax0-zmin0)/float(kz-1)
c     do 3 jz=1,kz
c     zl(jz)=zmin0+float(jz-1)*dz0
c   3 continue
c     do 20 jx=1,kx
c     do 20 jy=1,ky
c     do 20 jz=1,kz
c     xl1(jx,jy,jz)=xl(jx)
c     yl1(jx,jy,jz)=yl(jy)
c     zl1(jx,jy,jz)=zl(jz)
c  20 continue
      do 20 jx=1,mx
      do 20 jy=1,my
      do 20 jz=1,mz
      xl1(jx,jy,jz)=xxt(jx)
      yl1(jx,jy,jz)=yy(jy)
      zl1(jx,jy,jz)=zzt(jz)
   20 continue
c
      call current
cccccccccccccccccccccccccccccccccccccccc
c      call diverB
c      do 18 jz=z_first-2,z_last+2
c      do 18 jy=1,my
c      do 18 jx=x_first-2,x_last+2
c      cur(jx,jy,jz,1)=divB(jx,jy,jz)
c      cur(jx,jy,jz,1)=xresd(jx,jy,jz,1)
c      cur(jx,jy,jz,2)=xresd(jx,jy,jz,2)
c      cur(jx,jy,jz,3)=xresd(jx,jy,jz,3)
c   18 continue
cccccccccccccccccccccccccccccccccccccccc
      do jm=1,8
      do jx=1,mxpr
      do jy=1,my
      do jz=1,mzpr
      w(jx,jy,jz)=x(jx,jy,jz,jm)
      enddo
      enddo
      enddo
c
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_gather(w(1,1,1),mxpr*my*mzpr,
     1     mpi_double_precision,work4(1,1,1,1),mxpr*my*mzpr,
     2     mpi_double_precision,0,mpi_comm_world,ierr)
c
      if(nrank.eq.0)then
      do i=1,nsize
        if(i.le.(nsize-nprx))then
        if(mod((i),nprx).ne.0)then
        do jz=zfirst,zlast
        do jy=1,my
        do jx=xfirst,xlast
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        else
        do jz=zfirst,zlast
        do jy=1,my
        do jx=xfirst,xlast+mod(mx,nprx)
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        endif
        else
        if(mod((i),nprx).ne.0)then
        do jz=zfirst,zlast+mod(mz,nprz)
        do jy=1,my
        do jx=xfirst,xlast
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        else
        do jz=zfirst,zlast+mod(mz,nprz)
        do jy=1,my
        do jx=xfirst,xlast+mod(mx,nprx)
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        endif
        endif
      enddo
c     call extrap3(workt(1,1,1),xxt,yy,zzt,mx,my,mz
c    1             ,bl(1,1,1,jm),xl,yl,zl,kx,ky,kz)
      do 1 jx=1,mx
      do 1 jy=1,my
      do 1 jz=1,mz
      bl(jx,jy,jz,jm)=workt(jx,jy,jz)
    1 continue
      endif
      enddo
c
c
      do jm=1,3
      do jx=1,mxpr
      do jy=1,my
      do jz=1,mzpr
      w(jx,jy,jz)=cur(jx,jy,jz,jm)
      enddo
      enddo
      enddo
      call mpi_gather(w(1,1,1),mxpr*my*mzpr,
     1     mpi_double_precision,work4(1,1,1,1),mxpr*my*mzpr,
     2     mpi_double_precision,0,mpi_comm_world,ierr)
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(nrank.eq.0)then
      do i=1,nsize
        if(i.le.(nsize-nprx))then
        if(mod((i),nprx).ne.0)then
        do jz=zfirst,zlast
        do jy=1,my
        do jx=xfirst,xlast
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        else
        do jz=zfirst,zlast
        do jy=1,my
        do jx=xfirst,xlast+mod(mx,nprx)
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        endif
        else
        if(mod((i),nprx).ne.0)then
        do jz=zfirst,zlast+mod(mz,nprz)
        do jy=1,my
        do jx=xfirst,xlast
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        else
        do jz=zfirst,zlast+mod(mz,nprz)
        do jy=1,my
        do jx=xfirst,xlast+mod(mx,nprx)
        workt(mod((i-1),nprx)*(mxpr-4)+jx-2,jy,
     1        int((i-1)/nprx)*(mzpr-4)+jz-2)=work4(jx,jy,jz,i)
        enddo
        enddo
        enddo
        endif
        endif
      enddo
c     call extrap3(workt(1,1,1),xxt,yy,zzt,mx,my,mz
c    1             ,bl(1,1,1,jm+9),xl,yl,zl,kx,ky,kz)
      do 2 jx=1,mx
      do 2 jy=1,my
      do 2 jz=1,mz
      bl(jx,jy,jz,jm+8)=workt(jx,jy,jz)
    2 continue
      endif
      enddo
c
c
      do 3 i=2,4
      do 3 jx=1,mx
      do 3 jy=1,my
      do 3 jz=1,mz
      bl(jx,jy,jz,i)=bl(jx,jy,jz,i)/bl(jx,jy,jz,1)
    3 continue
c
c
      if(nrank.eq.0)then
      output='./data/m3d'//cn(nst)
      open(unit=8,file=output,status='unknown',form='formatted')
c     write(8,10)'Title="dipole field"'
c     write(8,11)'Variables="X","Y","Z","rho","Vx","Vy","Vz","Bx"
c    1,"By","Bz","P","Jx","Jy","Jz"'
c     write(8,12)'Zone  I=201  J=61  K=201  F=POINT'
c     write(8,30)(((xl1(jx,jy,jz),yl1(jx,jy,jz),zl1(jx,jy,jz),
c    1(bl(jx,jy,jz,i),i=1,12),jx=1,mx),jy=1,my),jz=1,mz)
c     revised by ji liang in 2018
      write(8,*)3,mx,my,mz
      write(8,*)time
      write(8,30)xxt,yy,zzt
      write(8,*)1
      write(8,30)(((bl(jx,jy,jz,1),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)2
      write(8,30)(((bl(jx,jy,jz,2),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)3
      write(8,30)(((bl(jx,jy,jz,3),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)4
      write(8,30)(((bl(jx,jy,jz,4),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)5
      write(8,30)(((bl(jx,jy,jz,5),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)6
      write(8,30)(((bl(jx,jy,jz,6),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)7
      write(8,30)(((bl(jx,jy,jz,7),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)8
      write(8,30)(((bl(jx,jy,jz,8),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)9
      write(8,30)(((bl(jx,jy,jz,9),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)10
      write(8,30)(((bl(jx,jy,jz,10),jx=1,mx),jy=1,my),jz=1,mz)
      write(8,*)11
      write(8,30)(((bl(jx,jy,jz,11),jx=1,mx),jy=1,my),jz=1,mz)
c     write(8,30)((((bl(jx,jy,jz,ii),jx=1,mx),jy=1,my),jz=1,mz),ii=1,11)
c     do 4 ii=1,11
c     write(8,30)bl(:,:,:,ii)
c   4 continue
   30 format(11(1x,e10.4))
c   30 format(1x,e13.5)
   10 format(1x,A20)
   11 format(1x,A76)
   12 format(1x,A33)
      close(8)
      endif
c
      return
      end
c
c
c
      character*3 function cn(n)
c
c-----assume that n is no greater than 999
c
c-----separate the digits
c
      n1=n/100
      n2=(n-100*n1)/10
      n3=n-100*n1-10*n2
c
c-----stick together cn using char function
c
      n1=n1+48
      n2=n2+48
      n3=n3+48
      cn(1:1)=char(n1)
      cn(2:2)=char(n2)
      cn(3:3)=char(n3)
c
      return
      end
c
c
c
      subroutine smthxy(t,ddt,km,ms,me,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
c
c three-point zero order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
c
      nsmthx=(nprx-1)*(mxpr-2)
      nsmthy=my-1
      nsmthz=mz-1
c
      if(.not.periodx) then
      if(myright.ne.mpi_proc_null)then
      do 10 k=1,kk
      do 11 m=ms,me
      do 12 jx=x_first,x_last
      theta=3.1415926*(mepx*(mxpr-2)+jx-2)/(nsmthx-3)
      do 12 jz=z_first,z_last
      do 12 jy=2,my-1
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(.5*(1.+cos(theta)))/24.*(
     1         difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   12 continue
      do 13 jx=x_first,x_last
      do 13 jz=z_first,z_last
      do 13 jy=2,my-1
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   13 continue
   11 continue
      call bndry(t,ddt,km)
   10 continue
      else
      if(.not.symmetryx) then
      do 20 k=1,kk
      do 21 m=ms,me
      do 22 jx=x_first,x_last
      theta=3.1415926*(mx-1-(mepx*(mxpr-2)+jx))/(nsmthx-3)
      do 22 jz=z_first,z_last
      do 22 jy=2,my-1
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(.5*(1.+cos(theta)))/24.*(
     1         difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   22 continue
      do 23 jx=x_first,x_last
      do 23 jz=z_first,z_last
      do 23 jy=2,my-1
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   23 continue
   21 continue
      call bndry(t,ddt,km)
   20 continue
      endif
      endif
      endif
c
      if(.not.periody) then
      do 30 k=1,kk
      do 31 m=ms,me
      do 32 jy=2,nsmthy
      theta=3.1415926*(jy-2)/(nsmthy-3)
      do 32 jz=z_first,z_last
      do 32 jx=x_first,x_last
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(.5*(1.+cos(theta)))/24.*(
     1         difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   32 continue
      do 33 jy=2,nsmthy
      do 33 jz=z_first,z_last
      do 33 jx=x_first,x_last
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   33 continue
   31 continue
      call bndry(t,ddt,km)
   30 continue
c
      do 40 k=1,kk
      do 41 m=ms,me
      do 42 jy=my-nsmthy+1,my-1
      theta=3.1415926*(my-1-jy)/(nsmthy-3)
      do 42 jz=z_first,z_last
      do 42 jx=x_first,x_last
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(.5*(1.+cos(theta)))/24.*(
     1         difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   42 continue
      do 43 jy=my-nsmthy+1,my-1
      do 43 jz=z_first,z_last
      do 43 jx=x_first,x_last
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   43 continue
   41 continue
      call bndry(t,ddt,km)
   40 continue
      endif
c
c     if(myupper.ne.mpi_proc_null)then
      do 50 k=1,kk
      do 51 m=ms,me
      do 52 jz=z_first,z_last
      theta=3.1415926*(mepz*(mzpr-2)+jz-2)/(nsmthz-3)
      do 52 jy=2,my-1
      do 52 jx=x_first,x_last
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(.5*(1.+cos(theta)))/24.*(
     1         difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   52 continue
      do 53 jz=z_first,z_last
      do 53 jy=2,my-1
      do 53 jx=x_first,x_last
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   53 continue
   51 continue
      call bndry(t,ddt,km)
   50 continue
c     else
      if(.not.symmetryz) then
      do 60 k=1,kk
      do 61 m=ms,me
      do 62 jz=z_first,z_last
      theta=3.1415926*(mz-1-(mepz*(mzpr-2)+jz))/(nsmthz-3)
      do 62 jy=2,my-1
      do 62 jx=x_first,x_last
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(.5*(1.+cos(theta)))/24.*(
     1         difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   62 continue
      do 63 jz=z_first,z_last
      do 63 jy=2,my-1
      do 63 jx=x_first,x_last
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   63 continue
   61 continue
      call bndry(t,ddt,km)
   60 continue
      endif
c     endif
      return
      end
c
c
c
      subroutine avrg(t,ddt,km,ms,me,caf0,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      real*8 coef_avrg
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      
c coef_weight is how much times larger the edge area to the middle area
      coef_weight=2
      
      do 11 k=1,kk
      do 14 m=ms,me
      do 12 jz=z_first,z_last
      do 12 jy=3,my-2
      do 12 jx=x_first,x_last
      if ((m.ge.5).and.(m.le.7)) then
c      goto 17
c      coef_avrg=0.0
      coef_avrg=1.0/32.0
c      coef_avrg=1.0/32.0+(coef_weight-1)/32.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
c      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(1-caf0)/32.0
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(1-caf0)*coef_avrg
     1         *(1-tanh((xx(jx)+3)/8.))
     1         *(difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
c      if(nrank.eq.0) write(*,*)'avrg jump 17 fail'
c   17 continue
      else
c      coef_avrg=0.0
c      coef_avrg=1.0/8.0
c      coef_avrg=1.0/8.0+(coef_weight-1)/8.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
      
ccccccfor non-uniform coef_avrg=a+(b-a)*0.5*(1+tanh()) varies form a to b
      coef_avrg=1.0/32.0+(1.0/8.0-1.0/32.0)
     1          *0.5*(1+tanh((sqrt(xx(jx)**2+yy(jy)**2+zz(jz)**2)
     2                         -3.)/3.))
      
c      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(1-caf0)/8.0
      w(jx,jy,jz)=xdif(jx,jy,jz,m)+(1-caf0)*coef_avrg
     1         *(difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
      endif
   12 continue
      do 13 jz=z_first,z_last
      do 13 jy=3,my-2
      do 13 jx=x_first,x_last
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   13 continue
   14 continue
      call bndry(t,ddt,km)
      
      goto 27
      if(nrank.eq.0) write(*,*)'avrg jump 27 fail'
      do 24 m=ms,me
      if ((m.ge.5).and.(m.le.7)) then
      goto 24
      else        
      do 22 jz=z_first,z_last
      do 22 jy=3,my-2
      do 22 jx=x_first,x_last
c      coef_avrg=0.0
      coef_avrg=1.0/4.0
c      coef_avrg=1.0/4.0+(coef_weight-1)/4.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
      w(jx,jy,jz)=xdif(jx,jy,jz,m)
c     1         +((1.-tanh((abs(xx(jx)-xmin)-10.)/10.))/8.
c     1         +(1.-tanh((rr0(jx,jy,jz)-6.)/3.))/4.
c     1         +(1.-tanh((rr0(jx,jy,jz)-6.)/3.))*coef_avrg
     1         +(1.-tanh((rr1(jx,jy,jz)-1.3)/0.6))*coef_avrg
     1         *(difc(xdif(jx-1,jy,jz,m),xdif(jx,jy,jz,m)
     2         ,xdif(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdif(jx,jy-1,jz,m),xdif(jx,jy,jz,m)
     4         ,xdif(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdif(jx,jy,jz-1,m),xdif(jx,jy,jz,m)
     6         ,xdif(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   22 continue
      do 23 jz=z_first,z_last
      do 23 jy=3,my-2
      do 23 jx=x_first,x_last
      xdif(jx,jy,jz,m)=w(jx,jy,jz)
   23 continue
      endif
   24 continue
      call bndry(t,ddt,km)
      
   27 continue
      
   11 continue
      return
      end
c
c
c
      subroutine avrgr(t,ddt,km,ms,me,caf0,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      do 11 k=1,kk
      do 12 m=ms,me
      do 12 jz=z_first-1,z_last+1
      do 12 jy=1,my
      do 12 jx=x_first-1,x_last+1
      if(rr1(jx,jy,jz).le.rcr) then
      xdif(jx,jy,jz,m)  = 0.
      else
      xdif(jx,jy,jz,m)=xdif(jx,jy,jz,m)
     1    *tanh((rr1(jx,jy,jz)/1.-1.)**2/0.5)
      endif
   12 continue
   11 continue
      return
      end
c
c
c
      subroutine avrgtt(t,ms,me,caf0,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      real*8 coef_avrg
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      
      coef_weight=2
      
      do 11 k=1,kk
      do 14 m=ms,me
      do 12 jz=z_first,z_last
      do 12 jy=3,my-2
      do 12 jx=x_first,x_last
      if ((m.ge.5).and.(m.le.7)) then
c      goto 17
c      coef_avrg=0.0
      coef_avrg=1.0/32.0
c      coef_avrg=1.0/32.0+(coef_weight-1)/32.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
c      w(jx,jy,jz)=xm(jx,jy,jz,m)+(1-caf0)/32.0
      w(jx,jy,jz)=xm(jx,jy,jz,m)+(1-caf0)*coef_avrg
     1         *(1-tanh((xx(jx)+3)/8.))
     1         *(difc(xm(jx-1,jy,jz,m),xm(jx,jy,jz,m)
     2         ,xm(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xm(jx,jy-1,jz,m),xm(jx,jy,jz,m)
     4         ,xm(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xm(jx,jy,jz-1,m),xm(jx,jy,jz,m)
     6         ,xm(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
c      if(nrank.eq.0) write(*,*)'avrgtt jump 17 fail'
c   17 continue
      else
c      coef_avrg=0.0
c      coef_avrg=1.0/8.0
c      coef_avrg=1.0/8.0+(coef_weight-1)/8.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
      
ccccccfor non-uniform coef_avrg=a+(b-a)*0.5*(1+tanh()) varies form a to b
      coef_avrg=1.0/32.0+(1.0/8.0-1.0/32.0)
     1          *0.5*(1+tanh((sqrt(xx(jx)**2+yy(jy)**2+zz(jz)**2)
     2                         -3.)/3.))
      
c      w(jx,jy,jz)=xm(jx,jy,jz,m)+(1-caf0)/8.0
      w(jx,jy,jz)=xm(jx,jy,jz,m)+(1-caf0)*coef_avrg
     1         *(difc(xm(jx-1,jy,jz,m),xm(jx,jy,jz,m)
     2         ,xm(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xm(jx,jy-1,jz,m),xm(jx,jy,jz,m)
     4         ,xm(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xm(jx,jy,jz-1,m),xm(jx,jy,jz,m)
     6         ,xm(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
      endif
   12 continue
      do 13 jz=z_first,z_last
      do 13 jy=3,my-2
      do 13 jx=x_first,x_last
      xm(jx,jy,jz,m)=w(jx,jy,jz)
   13 continue
   14 continue
      call bndrytt(t,ms,me)
      
      goto 27
      if(nrank.eq.0) write(*,*)'avrgtt jump 27 fail'
      do 24 m=ms,me
      if ((m.ge.5).and.(m.le.7)) then
      goto 24
      else
      do 22 jz=z_first,z_last
      do 22 jy=3,my-2
      do 22 jx=x_first,x_last
c      coef_avrg=0.0
      coef_avrg=1.0/8.0
c      coef_avrg=1.0/8.0+(coef_weight-1)/8.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
      w(jx,jy,jz)=xm(jx,jy,jz,m)
c     1         +((1.-tanh((abs(xx(jx)-xmin)-10.)/10.))/32.
c     1         +(1.-tanh((rr0(jx,jy,jz)-6.)/3.))/8.
c     1         +(1.-tanh((rr0(jx,jy,jz)-6.)/3.))*coef_avrg
     1         +(1.-tanh((rr1(jx,jy,jz)-1.3)/0.6))*coef_avrg
     1         *(difc(xm(jx-1,jy,jz,m),xm(jx,jy,jz,m)
     2         ,xm(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xm(jx,jy-1,jz,m),xm(jx,jy,jz,m)
     4         ,xm(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xm(jx,jy,jz-1,m),xm(jx,jy,jz,m)
     6         ,xm(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   22 continue
      do 23 jz=z_first,z_last
      do 23 jy=3,my-2
      do 23 jx=x_first,x_last
      xm(jx,jy,jz,m)=w(jx,jy,jz)
   23 continue
      endif
   24 continue
      call bndrytt(t,ms,me)
      
   27 continue
      
   11 continue
      return
      end
c
c
c
      subroutine avrgE(t,ms,me,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      real*8 coef_avrg
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      
c coef_weight is how much times larger the edge area to the middle area
      coef_weight=2
      
      do 11 k=1,kk
      do 14 m=ms,me
      do 12 jz=z_first,z_last
      do 12 jy=3,my-2
      do 12 jx=x_first,x_last
c      coef_avrg=0.0
c      coef_avrg=1.0/8.0
c      coef_avrg=1.0/8.0+(coef_weight-1)/8.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
      
ccccccfor non-uniform coef_avrg=a+(b-a)*0.5*(1+tanh()) varies form a to b
      coef_avrg=1.0/32.0+(1.0/8.0-1.0/32.0)
     1          *0.5*(1+tanh((sqrt(xx(jx)**2+yy(jy)**2+zz(jz)**2)
     2                         -3.)/3.))
      
c      w(jx,jy,jz)=efld(jx,jy,jz,m)+(1-caf0)/8.0
      w(jx,jy,jz)=efld(jx,jy,jz,m)+coef_avrg
     1         *(difc(efld(jx-1,jy,jz,m),efld(jx,jy,jz,m)
     2         ,efld(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(efld(jx,jy-1,jz,m),efld(jx,jy,jz,m)
     4         ,efld(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(efld(jx,jy,jz-1,m),efld(jx,jy,jz,m)
     6         ,efld(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   12 continue
      do 13 jz=z_first,z_last
      do 13 jy=3,my-2
      do 13 jx=x_first,x_last
      efld(jx,jy,jz,m)=w(jx,jy,jz)
   13 continue
   14 continue
      call bndryE
c
   11 continue
      return
      end
c
c
c
      subroutine avrgrtt(ms,me,caf0,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      do 11 k=1,kk
      do 12 m=ms,me
      do 12 jz=z_first-1,z_last+1
      do 12 jy=1,my
      do 12 jx=x_first-1,x_last+1
      if(rr1(jx,jy,jz).le.rcr) then
      xm(jx,jy,jz,m)  = 0.
      else
      xm(jx,jy,jz,m)=xm(jx,jy,jz,m)*tanh((rr1(jx,jy,jz)/1.-1.)**2/0.5)
      endif
   12 continue
c     call bndrytt(ms,me)
   11 continue
      return
      end
c
c
c
      subroutine avrgc(ms,me,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      real*8 coef_avrg
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      
      coef_weight=2
      
      do 11 k=1,kk
      do 14 m=ms,me
      do 12 jz=z_first,z_last
      do 12 jy=3,my-2
      do 12 jx=x_first,x_last
c      coef_avrg=0.0
c      coef_avrg=1.0/8.0
c      coef_avrg=1.0/8.0+(coef_weight-1)/8.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
      
ccccccfor non-uniform coef_avrg=a+(b-a)*0.5*(1+tanh()) varies form a to b
      coef_avrg=1.0/32.0+(1.0/8.0-1.0/32.0)
     1          *0.5*(1+tanh((sqrt(xx(jx)**2+yy(jy)**2+zz(jz)**2)
     2                         -3.)/3.))
      
c      w(jx,jy,jz)=cur(jx,jy,jz,m)+1/8.*(
      w(jx,jy,jz)=cur(jx,jy,jz,m)+coef_avrg*(
     1         difc(cur(jx-1,jy,jz,m),cur(jx,jy,jz,m)
     2         ,cur(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(cur(jx,jy-1,jz,m),cur(jx,jy,jz,m)
     4         ,cur(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(cur(jx,jy,jz-1,m),cur(jx,jy,jz,m)
     6         ,cur(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   12 continue
      do 13 jz=z_first,z_last
      do 13 jy=3,my-2
      do 13 jx=x_first,x_last
      cur(jx,jy,jz,m)=w(jx,jy,jz)
   13 continue
   14 continue
      call bndryc
c
   11 continue
      return
      end
c
c
c
      subroutine avrgh(ms,me,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      real*8 coef_avrg
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      
      coef_weight=2
      
      do 11 k=1,kk
      do 14 m=ms,me
      do 12 jz=z_first,z_last
      do 12 jy=3,my-2
      do 12 jx=x_first,x_last
c      coef_avrg=0.0
      coef_avrg=1.0/8.0
c      coef_avrg=1.0/8.0+(coef_weight-1)/8.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
c      w(jx,jy,jz)=xdh(jx,jy,jz,m)+1/8.*(
      w(jx,jy,jz)=xdh(jx,jy,jz,m)+coef_avrg*(
     1         difc(xdh(jx-1,jy,jz,m),xdh(jx,jy,jz,m)
     2         ,xdh(jx+1,jy,jz,m),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(xdh(jx,jy-1,jz,m),xdh(jx,jy,jz,m)
     4         ,xdh(jx,jy+1,jz,m),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(xdh(jx,jy,jz-1,m),xdh(jx,jy,jz,m)
     6         ,xdh(jx,jy,jz+1,m),zz(jz-1),zz(jz),zz(jz+1)))
   12 continue
      do 13 jz=z_first,z_last
      do 13 jy=3,my-2
      do 13 jx=x_first,x_last
      xdh(jx,jy,jz,m)=w(jx,jy,jz)
   13 continue
   14 continue
      call bndryh
c
   11 continue
      return
      end
c
c
c
      subroutine avrgdiv(ms,me,kk)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
      integer status(mpi_status_size)
      real*8 coef_avrg
c
c three-point zero order diffusion (for non-uniform grids)
      difc(fm,f0,fp,xm1,x0,xp1)=
     1   2.*( (fp-f0)*(xp1-x0)-(f0-fm)*(x0-xm1))/(xp1-xm1)
c     difc=h**2*f''(x)=f(x+h)-2f(x)+f(x-h)
c three-point second-order diffusion (for non-uniform grids)
c      difc(fm,f0,fp,xm1,x0,xp1)=
c     1   2.*( (fp-f0)*(x0-xm1)-(f0-fm)*(xp1-x0))/(xp1-xm1)
      
      coef_weight=2
      
      do 11 k=1,kk
      do 14 m=ms,me
      do 12 jz=z_first,z_last
      do 12 jy=3,my-2
      do 12 jx=x_first,x_last
c      coef_avrg=0.0
c      coef_avrg=1.0/8.0
c      coef_avrg=1.0/8.0+(coef_weight-1)/8.0
c     1          *0.5*(1+dtanh((sqrt(yy(jy)**2+zz(jz)**2)-20.)/3.))
c     2          *0.5*(1+dtanh((xx(jx)-20.)/5.))
      
ccccccfor non-uniform coef_avrg=a+(b-a)*0.5*(1+tanh()) varies form a to b
      coef_avrg=1.0/32.0+(1.0/8.0-1.0/32.0)
     1          *0.5*(1+tanh((sqrt(xx(jx)**2+yy(jy)**2+zz(jz)**2)
     2                         -3.)/3.))
      
c      w(jx,jy,jz)=divB(jx,jy,jz)+1/8.*(
      w(jx,jy,jz)=divB(jx,jy,jz)+coef_avrg*(
     1         difc(divB(jx-1,jy,jz),divB(jx,jy,jz)
     2         ,divB(jx+1,jy,jz),xx(jx-1),xx(jx),xx(jx+1))+
     3         difc(divB(jx,jy-1,jz),divB(jx,jy,jz)
     4         ,divB(jx,jy+1,jz),yy(jy-1),yy(jy),yy(jy+1))+
     5         difc(divB(jx,jy,jz-1),divB(jx,jy,jz)
     6         ,divB(jx,jy,jz+1),zz(jz-1),zz(jz),zz(jz+1)))
   12 continue
      do 13 jz=z_first,z_last
      do 13 jy=3,my-2
      do 13 jx=x_first,x_last
      divB(jx,jy,jz)=w(jx,jy,jz)
   13 continue
   14 continue
      call bndrydiv
c
   11 continue
      return
      end
c
c
c
      subroutine funmax(f,fmax,fmin,x,z,mx,mz)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer status(mpi_status_size)
      dimension f(mx,mz),x(mx),z(mz)
      double precision fmax0,fmin0,fmax,fmin,xlmax,zlmax,
     1                 xlmin,jxmin,jzmin
      fmax0=-1000000.
      fmin0=1000000.
      do 10 jz=2,mz-1
      do 10 jx=2,mx-1
      if(f(jx,jz).gt.fmax0) then
      fmax0=f(jx,jz)
      xlmax=x(jx)
      zlmax=z(jz)
      jxmax=jx
      jzmax=jz
      endif
      if(f(jx,jz).lt.fmin0) then
      fmin0=f(jx,jz)
      xlmin=x(jx)
      zlmin=z(jz)
      jxmin=jx
      jzmin=jz
      endif
   10 continue
c     call mpi_barrier(mpi_comm_world,ierr)
      call mpi_allreduce(fmax0,fmax,1,mpi_double_precision,mpi_max,
     1     mpi_comm_world,ierr)
      call mpi_allreduce(fmin0,fmin,1,mpi_double_precision,mpi_min,
     1     mpi_comm_world,ierr)
c     call mpi_barrier(mpi_comm_world,ierr)
      return
      end
c
c
c
      subroutine avct(a,bx,bz,nx,nz,x,z)
      implicit real*8 (a-h,o-z)
c To calculate the vector potential in x-z plane
      dimension a(nx,nz),bx(nx,nz),bz(nx,nz),x(nx),z(nz)
c
      a(1,1)=100.
      do 10 jx=2,nx
      a(jx,1)=a(jx-1,1)-.5*(bz(jx,1)+bz(jx-1,1))*(x(jx)-x(jx-1))
   10 continue
      do 20 jx=1,nx
      do 20 jz=2,nz
      a(jx,jz)=a(jx,jz-1)+.5*(bx(jx,jz)+bx(jx,jz-1))*(z(jz)-z(jz-1))
   20 continue
c
      return
      end
c
c
c
      subroutine message1(x,nx,myleft,myright)
c     Author  --- Lu Xingqiang
c     History --- completed in July, 2009.
      implicit real*8 (a-h,o-z)
c     include 'mhd01.for'
      include 'mpif.h'
      integer nx,myleft,myright
      integer status(mpi_status_size)
      dimension x(1:nx)
c
      call MPI_SENDRECV(x(3),1,MPI_DOUBLE_PRECISION,
     1     myleft,11,x(nx-1),1,MPI_DOUBLE_PRECISION,
     1     myright,11,MPI_COMM_WORLD,status,ierr)
c
      call MPI_SENDRECV(x(4),1,MPI_DOUBLE_PRECISION,
     1     myleft,13,x(nx),1,MPI_DOUBLE_PRECISION,
     1     myright,13,MPI_COMM_WORLD,status,ierr)
c
      call MPI_SENDRECV(x(nx-2),1,MPI_DOUBLE_PRECISION,
     1     myright,12,x(2),1,MPI_DOUBLE_PRECISION,
     1     myleft,12,MPI_COMM_WORLD,status,ierr)

      call MPI_SENDRECV(x(nx-3),1,MPI_DOUBLE_PRECISION,
     1     myright,14,x(1),1,MPI_DOUBLE_PRECISION,
     1     myleft,14,MPI_COMM_WORLD,status,ierr)
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
      return
      end
c
c
c
      subroutine message2(x,mxp,mzp,myleft,myright,mylow,myupper)
c     Author  --- Lu Xingqiang
c     History --- completed in July, 2009.
      implicit real*8 (a-h,o-z)
c     include 'mhd01.for'
      include 'mpif.h'
      integer status(mpi_status_size)
      integer myleft,myright,mylow,myupper
      integer mxp,myp,htype,vtype,ierr
      dimension x(mxp,mzp)
c
      call mpi_type_contiguous(mxp,mpi_double_precision,htype,ierr)
      call mpi_type_commit(htype,ierr)
      call mpi_type_vector(mzp-4,1,(mxp),
     1     mpi_double_precision,vtype,ierr)
      call mpi_type_commit(vtype,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
cccccc锟斤拷锟斤拷一锟斤拷
c-----x锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(x(3,3),1,vtype,
     1     myleft,11,x(mxp-1,3),1,vtype,
     1     myright,11,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(x(4,3),1,vtype,
     1     myleft,15,x(mxp,3),1,vtype,
     1     myright,15,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(x(mxp-2,3),1,vtype,
     1     myright,12,x(2,3),1,vtype,
     1     myleft,12,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(x(mxp-3,3),1,vtype,
     1     myright,16,x(1,3),1,vtype,
     1     myleft,16,MPI_COMM_WORLD,status,ierr)
c-----z锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(x(1,4),1,htype,
     1     mylow,13,x(1,mzp),1,htype,
     1     myupper,13,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(x(1,3),1,htype,
     1     mylow,17,x(1,mzp-1),1,htype,
     1     myupper,17,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(x(1,mzp-3),1,htype,
     1     myupper,14,x(1,1),1,htype,
     1     mylow,14,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(x(1,mzp-2),1,htype,
     1     myupper,18,x(1,2),1,htype,
     1     mylow,18,MPI_COMM_WORLD,status,ierr)
      call mpi_type_free(htype,ierr)
      call mpi_type_free(vtype,ierr)
cccccc锟斤拷锟斤拷锟斤拷锟斤拷
c锟斤拷锟斤拷锟斤拷锟竭斤拷
c     call mpi_send(x(2,1),1,vtype,myleft,10,
c    1     mpi_comm_world,ierr)
c锟斤拷锟斤拷锟揭边斤拷
c     call mpi_send(x(mxp-1,1),1,vtype,myright,20,
c    1     mpi_comm_world,ierr)
c锟斤拷锟斤拷锟铰边斤拷
c     call mpi_send(x(1,mzp-1),1,htype,myupper,30,
c    1     mpi_comm_world,ierr)
c锟斤拷锟斤拷锟较边斤拷
c     call mpi_send(x(1,2),1,htype,mylow,40,
c    1     mpi_comm_world,ierr)
c     return
c锟斤拷锟斤拷锟揭边斤拷
c     call mpi_recv(x(mxp,1),1,vtype,myright,10,
c    1     mpi_comm_world,status,ierr)
c锟斤拷锟斤拷锟斤拷锟竭斤拷
c     call mpi_recv(x(1,1),1,vtype,myleft,20,
c    1     mpi_comm_world,status,ierr)
c     return
c锟斤拷锟斤拷锟较边斤拷
c     call mpi_recv(x(1,1),1,htype,mylow,30,
c    1     mpi_comm_world,status,ierr)
c锟斤拷锟斤拷锟铰边斤拷
c     call mpi_recv(x(1,mzp),1,htype,myupper,40,
c    1     mpi_comm_world,status,ierr)
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     call mpi_type_free((htype,ierror)
c     call mpi_type_free((vtype,ierror)
c
      return
      end
c
c
c
      subroutine message3(x,mxp,myp,mzp,myleft,myright,mylow,myupper)
c     Author  --- Lu Xingqiang
c     History --- completed in July, 2009.
      implicit real*8 (a-h,o-z)
c     include 'mhd01.for'
      include 'mpif.h'
      integer status(mpi_status_size)
      integer myleft,myright,mylow,myupper
      integer mxp,myp,mzp,htype,vtype,ierr
      dimension x(mxp,myp,mzp),work2(mxp,mzp)
c
      call mpi_type_contiguous(mxp,mpi_double_precision,htype,ierr)
      call mpi_type_commit(htype,ierr)
      call mpi_type_vector(mzp-4,1,(mxp),
     1     mpi_double_precision,vtype,ierr)
      call mpi_type_commit(vtype,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do jy=1,myp
      do jz=1,mzp
      do jx=1,mxp
      work2(jx,jz)=x(jx,jy,jz)
      enddo
      enddo
cccccc锟斤拷锟斤拷一锟斤拷
c-----x锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(work2(3,3),1,vtype,
     1     myleft,11,work2(mxp-1,3),1,vtype,
     1     myright,11,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(4,3),1,vtype,
     1     myleft,15,work2(mxp,3),1,vtype,
     1     myright,15,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(mxp-2,3),1,vtype,
     1     myright,12,work2(2,3),1,vtype,
     1     myleft,12,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(mxp-3,3),1,vtype,
     1     myright,16,work2(1,3),1,vtype,
     1     myleft,16,MPI_COMM_WORLD,status,ierr)
c-----z锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(work2(1,4),1,htype,
     1     mylow,13,work2(1,mzp),1,htype,
     1     myupper,13,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,3),1,htype,
     1     mylow,17,work2(1,mzp-1),1,htype,
     1     myupper,17,MPI_COMM_WORLD,status,ierr)
c
      call MPI_SENDRECV(work2(1,mzp-3),1,htype,
     1     myupper,14,work2(1,1),1,htype,
     1     mylow,14,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,mzp-2),1,htype,
     1     myupper,18,work2(1,2),1,htype,
     1     mylow,18,MPI_COMM_WORLD,status,ierr)
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do jz=1,mzp
      do jx=1,mxp
      x(jx,jy,jz)=work2(jx,jz)
      enddo
      enddo
      enddo
      call mpi_type_free(htype,ierr)
      call mpi_type_free(vtype,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
      return
      end
c
c
c
      subroutine message4(x,mxp,myp,mzp,m,myleft,myright,mylow,myupper)
c     Author  --- Lu Xingqiang
c     History --- completed in July, 2009.
      implicit real*8 (a-h,o-z)
c     include 'mhd01.for'
      include 'mpif.h'
      integer status(mpi_status_size)
      integer myleft,myright,mylow,myupper
      integer mxp,myp,mzp,m,htype,vtype,ierr
      dimension x(mxp,myp,mzp,m),work2(mxp,mzp)
c
      call mpi_type_contiguous(mxp,mpi_double_precision,htype,ierr)
      call mpi_type_commit(htype,ierr)
      call mpi_type_vector(mzp-4,1,(mxp),
     1     mpi_double_precision,vtype,ierr)
      call mpi_type_commit(vtype,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do jm=1,m
      do jy=1,myp
      do jz=1,mzp
      do jx=1,mxp
      work2(jx,jz)=x(jx,jy,jz,jm)
      enddo
      enddo
cccccc锟斤拷锟斤拷一锟斤拷
c-----x锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(work2(3,3),1,vtype,
     1     myleft,11,work2(mxp-1,3),1,vtype,
     1     myright,11,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(4,3),1,vtype,
     1     myleft,15,work2(mxp,3),1,vtype,
     1     myright,15,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(mxp-2,3),1,vtype,
     1     myright,12,work2(2,3),1,vtype,
     1     myleft,12,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(mxp-3,3),1,vtype,
     1     myright,16,work2(1,3),1,vtype,
     1     myleft,16,MPI_COMM_WORLD,status,ierr)
c-----z锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(work2(1,4),1,htype,
     1     mylow,13,work2(1,mzp),1,htype,
     1     myupper,13,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,3),1,htype,
     1     mylow,17,work2(1,mzp-1),1,htype,
     1     myupper,17,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,mzp-3),1,htype,
     1     myupper,14,work2(1,1),1,htype,
     1     mylow,14,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,mzp-2),1,htype,
     1     myupper,18,work2(1,2),1,htype,
     1     mylow,18,MPI_COMM_WORLD,status,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do jz=1,mzp
      do jx=1,mxp
      x(jx,jy,jz,jm)=work2(jx,jz)
      enddo
      enddo
      enddo
      enddo
c
      call mpi_type_free(htype,ierr)
      call mpi_type_free(vtype,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
      return
      end
c
c
c
      subroutine message5(x,mxp,myp,mzp,myleft,myright,mylow,myupper)
c     for m=1,x=x(mxp,myp,mzp)
      implicit real*8 (a-h,o-z)
c     include 'mhd01.for'
      include 'mpif.h'
      integer status(mpi_status_size)
      integer myleft,myright,mylow,myupper
      integer mxp,myp,mzp,htype,vtype,ierr
      dimension x(mxp,myp,mzp),work2(mxp,mzp)
c
      call mpi_type_contiguous(mxp,mpi_double_precision,htype,ierr)
      call mpi_type_commit(htype,ierr)
      call mpi_type_vector(mzp-4,1,(mxp),
     1     mpi_double_precision,vtype,ierr)
      call mpi_type_commit(vtype,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do jy=1,myp
      do jz=1,mzp
      do jx=1,mxp
      work2(jx,jz)=x(jx,jy,jz)
      enddo
      enddo
cccccc锟斤拷锟斤拷一锟斤拷
c-----x锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(work2(3,3),1,vtype,
     1     myleft,11,work2(mxp-1,3),1,vtype,
     1     myright,11,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(4,3),1,vtype,
     1     myleft,15,work2(mxp,3),1,vtype,
     1     myright,15,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(mxp-2,3),1,vtype,
     1     myright,12,work2(2,3),1,vtype,
     1     myleft,12,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(mxp-3,3),1,vtype,
     1     myright,16,work2(1,3),1,vtype,
     1     myleft,16,MPI_COMM_WORLD,status,ierr)
c-----z锟斤拷锟斤拷锟斤拷
      call MPI_SENDRECV(work2(1,4),1,htype,
     1     mylow,13,work2(1,mzp),1,htype,
     1     myupper,13,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,3),1,htype,
     1     mylow,17,work2(1,mzp-1),1,htype,
     1     myupper,17,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,mzp-3),1,htype,
     1     myupper,14,work2(1,1),1,htype,
     1     mylow,14,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(work2(1,mzp-2),1,htype,
     1     myupper,18,work2(1,2),1,htype,
     1     mylow,18,MPI_COMM_WORLD,status,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do jz=1,mzp
      do jx=1,mxp
      x(jx,jy,jz)=work2(jx,jz)
      enddo
      enddo
      enddo
c
      call mpi_type_free(htype,ierr)
      call mpi_type_free(vtype,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
      return
      end
c
c
c
      subroutine foreta(t)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
c
      etac=0.d-4
c     cj0=5.0d0
c     alpha0=0.d0
c     crtrn0=0.d0
      do 1 jz=z_first-2,z_last+2
      do 1 jy=1,my
      do 1 jx=x_first-2,x_last+2
c nonlinear resistivity
c     crtrn=cur(jx,jz,3)-cj0
c     crtrn=dmax1(crtrn,crtrn0)
c     ch=alpha0*crtrn*crtrn
c     eta1=dmin1(ch,etac)
c     etan(jx,jy,jz)=eta0
      etan(jx,jy,jz)=eta0
c     1              +0.5d0*eta0*100.*(1+dtanh((dabs(zz(jz))-20.)/3.))
    1 continue
      return
      end
c
c
c
      subroutine formu(t)
      include 'mhd01.for'
      include 'mpif.h'
      include 'mhd02.for'
c
      do 1 jz=z_first-2,z_last+2
      do 1 jy=1,my
      do 1 jx=x_first-2,x_last+2
c     fmu(jx,jy,jz)=fmu0
      fmu(jx,jy,jz)=fmu0
c     1             +0.5d0*fmu0*100.*(1+dtanh((dabs(zz(jz))-20.)/3.))
      frho(jx,jy,jz)=frho0
c     1             +0.5d0*frho0*100.*(1+dtanh((dabs(zz(jz))-20.)/3.))
      fpr(jx,jy,jz)=fpr0
c     1             +0.5d0*fpr0*100.*(1+dtanh((dabs(zz(jz))-20.)/3.))
    1 continue
      return
      end
