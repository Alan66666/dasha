      lrstrt=.false.
c      lrstrt=.true.
c     uniformx=.false.
      uniformx=.true.
      uniformy=.true.
c      uniformy=.false.
      uniformz=.true.
c      uniformz=.false.
      halfx=.false.
      halfy=.false.
      halfz=.false.
      symmetryx=.false.
      symmetryy=.false.
      symmetryz=.false.
      periody=.false.
      periodx=.false.
      periodz=.false.
      viscos=.true.
      damp=.false.
      diffu=.true.
      hall=.false.
      uyzlimiter=.false.
      totalE=.true.!to decide how to calculate -curl E (if including -u1 curl B0 in E)  
      mxdat=.true.  
      mxdatfc=.true.  
      rcr=1.0  
      ncase=1
      nend=5011
      nst=1
      nstnew=1
      nint=1
      dxmin=0.01
      dxmax=0.5
      dymin=0.01!0.2
      dymax=1.0!0.8
      dymaxo=1.0
      dzmin=0.01!0.01
      dzmax=1.0!0.5
      dzmaxo=1.0
      xmin0=-2.1
      xmin=-2.1
      xmax=5.65
      ymin=0.
      ymax=2.5!30.
      ymaxo=2.5
      ymin0=-2.5!-30.
      yhr=0.0
      zmin=0.
      zmax=2.5!30.0
      zmaxo=2.5!50.0
      zmin0=-2.5!-30.0
      zhr=0.0
      pi=3.141592654d0
      nsmthx=2*mx/3-1
      nsmthy=my-1
      nsmthz=mz2/2
      nsmthxt=2*mx/3-1
      nsmthyt=my-1
      nsmthzt=mz2/2
      gamma=5.0/3.0
      eta0=2.d-4
      fmu0=0.!1.d-1
      frho0=0.!1.d-2
      fpr0=0.!1.d-2
c      t0=200./8.
c      t00=200/8
      t0=20.
      t00=20.
      by0=0.0
      fe0=0.
      cfl=2.0
      beta0=0.2
      wids=zmax/20.
      widp=zmax/20.
      wid1=zmax/8.
      wid2=zmax/4.
      widv=2.8
      vi0=0.
      bxi0=1.0
      bxi1=0.5
      rho0=1.0
      rho1=0.
      chall=0.07