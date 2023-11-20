c      parameter (mx1=721,mx=mx1,my1i=61,my1o=40,my1=my1i+my1o,my=my1
c     1           ,mz1i=401,mz1o=40,mz1=mz1i+mz1o,mz=mz1,mz2=2*mz-3)
      parameter (mx1=201,mx=mx1,my1i=31,my1o=20,my1=my1i+my1o,my=my1
     1           ,mz1i=85,mz1o=20,mz1=mz1i+mz1o,mz=mz1,mz2=2*mz-3)
      parameter (nprx=5,nprz=8,nprt=nprx*nprz)
      parameter (mxpr=int(mx/nprx)+4,mzpr=mz/nprz+4)
      implicit real*8 (a-h,o-z)
      logical lrstrt,uniformx,uniformy,uniformz,halfx,halfy,halfz
     1        ,periodz,hall,uyzlimiter,totalE
      logical symmetryx,symmetryy,symmetryz,periody,periodx,viscos,
     1        damp,diffu
      logical mxdat,mxdatfc
      dimension w(mxpr,my,mzpr),cur(mxpr,my,mzpr,3),wf4(mxpr,my,3)
     1          ,wf1(my,mzpr-2,8),wf2(my,mzpr-2,3),xdh(mxpr,my,mzpr,3)
     2          ,wf3(mxpr,my,8),byi(mxpr,my,mzpr),rr0(mxpr,my,mzpr)
     3          ,divB(mxpr,my,mzpr),rr1(mxpr,my,mzpr)
      dimension psi(mxpr,mzpr),bxi(mxpr,my,mzpr),bzi(mxpr,my,mzpr),
     1          rhoi(mxpr,my,mzpr),work2(mxpr,mzpr),vix(mxpr,mzpr)
      dimension bxinfc(my,mzpr),byinfc(my,mzpr),bzinfc(my,mzpr)	 
      dimension xx(mxpr),yy(my),zz(mzpr),cj(mxpr,mzpr),pri(mxpr,my,mzpr)
     1    ,etan(mxpr,my,mzpr),efld(mxpr,my,mzpr,3),xresd(mxpr,my,mzpr,3)
     2    ,fmu(mxpr,my,mzpr),frho(mxpr,my,mzpr),fpr(mxpr,my,mzpr)
      dimension xxt(mx),zzt(mz),dxt(mx),dzt(mz),workt(mx,my,mz)
      dimension bxi2(mx,my,mz),byi2(mx,my,mz),bzi2(mx,my,mz)
      dimension bxpf1(my,mz),bypf1(my,mz),bzpf1(my,mz)
      dimension gx(mxpr,mzpr),gy(mxpr,mzpr),gz(mxpr,mzpr),
     1          feng(mxpr,my,mzpr),fe(mxpr,my),ff(mxpr,mzpr)
      dimension dx(mxpr),dy(my),dz(mzpr),xrest(mzpr,8),nstp(10000)
      dimension ax1(mxpr),bx1(mxpr),cx1(mxpr),dx1(mxpr)
      dimension ax2(mxpr),bx2(mxpr),cx2(mxpr),dx2(mxpr)
      dimension ay1(my),by1(my),cy1(my),dy1(my),cht(mxpr,my,mzpr)
      dimension ay2(my),by2(my),cy2(my),dy2(my)
      dimension az1(mzpr),bz1(mzpr),cz1(mzpr),dz1(mzpr)
      dimension az2(mzpr),bz2(mzpr),cz2(mzpr),dz2(mzpr)
      dimension axm(mxpr),bxm(mxpr),cxm(mxpr)
      dimension axp(mxpr),bxp(mxpr),cxp(mxpr)
      dimension aym(my),bym(my),cym(my)
      dimension ayp(my),byp(my),cyp(my)
      dimension azm(mzpr),bzm(mzpr),czm(mzpr)
      dimension azp(mzpr),bzp(mzpr),czp(mzpr)
      integer nrank,nsize,ierr,myleft,myright,myupper,mylow,
     1        mepx,mepz,x_first,x_last,z_first,z_last,
     2        xfirst,xlast,zfirst,zlast,dtstatus
      common /indexmpi/ nrank,nsize,ierr,myleft,myright,myupper,
     1                  mylow,mepx,mepz,x_first,x_last,z_first,z_last,
     2                  xfirst,xlast,zfirst,zlast,periodx,periodz
      common /logic/ lrstrt,uniformx,uniformy,uniformz,halfx,halfy,hall
     1               ,halfz,symmetryx,symmetryy,symmetryz,periody
     2               ,viscos,damp,diffu,uyzlimiter,totalE,mxdat,mxdatfc
      common /para1/ xmax,xmin,xmin0,ymax,ymin,zmax,zmin,zmin0,xwfront
     1               ,ymaxo,zmaxo
      common /para2/ ymax0,dxmin,dymin,dymax,yhr,dzmin,dzmax,zhr
     1               ,dymaxo,dzmaxo
      common /para3/ rcr	 
      common /cst/ eta0,fmu0,frho0,fpr0,gamma,time,dt,dt0,t0,t00,cfl,
     1             fe0,by0
      common /num/ nsmthx,nsmthy,nsmthz,nsmthxt,nsmthyt,nsmthzt,nst1
      common /run/ nend,nstep,nsp,ncase,nst,nint,nstp,np,nstnew,dtstatus
      common /fourth1/ ax1,bx1,cx1,dx1,ay1,by1,cy1,dy1,az1,bz1,cz1,dz1,
     1                 ax2,bx2,cx2,dx2,ay2,by2,cy2,dy2,az2,bz2,cz2,dz2
      common /bndcoefp/ axp,bxp,cxp,ayp,byp,cyp,azp,bzp,czp
      common /bndcoefm/ axm,bxm,cxm,aym,bym,cym,azm,bzm,czm
      common /coeff1/ cxp1,cxp2,cxm1,cxm2,cyp1,cyp2,cym1,cym2,dtime
      common /coeff2/ czp1,czp2,czm1,czm2,beta0,wids,widp,wid1
      common /pert/ alamda,alpha,caf,pi,wid2,widv,vi0,bxi0
      common /var11/ cur,psi,bxi,bzi,cj,rhoi,pri,etan,xresd,fmu,frho,fpr
      common /var12/ xx1,yy1,zz1,bxi3,byi3,bzi3,bxi2,byi2,bzi2
      common /var13/ bxpf1,bxpf,bypf1,bypf,bzpf1,bzpf
      common /var14/ gx,gy,gz,feng,fe,bxi1,rho0,rho1,efld
      common /var15/ xx,yy,zz,dx,dy,dz,cht,xdh,chall
      common /var16/ xxt,zzt,dxt,dzt
      common /var17/ bxinfc,byinfc,bzinfc
      common /wrk/ w,work2,workt,byi,rr0,xrest,divB,rr1
