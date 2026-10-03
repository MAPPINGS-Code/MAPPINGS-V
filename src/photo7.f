cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'credits.inc'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine photo7
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine photo7 ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******PHOTOIONISATIONMODEL
c     DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION, RADIATION
c     PRESSURE
c
c     space steps derived from optical depth of predicted temperature
c     and ionistation, extrapolated from previous steps
c
c     CHOICE OF EQUILIBRIUM OR FINITE AGE CONDITIONS
c
c     NB. COMPUTATIONS PERFORMED IN SUBROUTINE COMPPH6
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'p7blocks.inc'
c
c           Variables
c
      real*8 blum,ilum,dhlma,dhn,dht,difma,dthre,dti,dtlma
      real*8 epotmi,fin,uinit,qinit,rcin
      real*8 rstsun,starlum
      real*8 scale,dnt,lam
c
      integer*4 j,i,kmax,k,nl,np,idx
      integer*4 elok(mxelem),nel
c
      character carac*36,model*64
      character caract*4,ilgg*4
      character banfil*64
c
      character outsettings*64,outshow*64
c
c           Functions
c
      real*8 densnum,fdilu
      integer*4 mlen
c
      jcon='Y'
      jspot='N'
      dtau0=0.025d0
c
      write (*,10)
c
c     ***INITIAL IONISATION CONDITIONS
c
   10 format(///
     & ' ********************************************************',/,
     & '  Photoionisation model P7 selected:',/,
     & '  Diffuse Field : Radiation Pressure, Exp. API',/,
     & ' ********************************************************'//)
c
      epotmi=iphe
c
c     set up ionisation state
c
      model='protoionisation'
      call popcha (model)
      model='Photo 7'
c
   20 format(a)
c
c     artificial support for low ionisation species
c
      if (expertmode.gt.0) then
c
        carac(1:33)='   '
        j=1
        do i=3,atypes
          j=1+index(carac(1:33),'   ')
          caract='   '
          if (ipote(1,i).lt.epotmi) write (caract,20) elem(i)
          carac(j:j+2)=caract//'   '
        enddo
        i=j+2
   30   write (*,40) carac(1:i)
        ilgg='    '
c
   40    format(//' Allow the elements : ',a/
     & ' to recombine freely? (y/n) : ',$)
c
        read (*,20) ilgg
        call toup (ilgg(1:1), ilgg)
c
        if (ilgg.eq.'Y') goto 50
        if ((ilgg.ne.'N')) goto 30
c
        do i=3,atypes
          arad(2,i)=dabs(arad(2,i))
          if (ipote(1,i).lt.epotmi) arad(2,i)=-arad(2,i)
        enddo
c
   50   continue
c
      endif
c
c     ***CHOOSING PHOTON SOURCE AND GEOMETRICAL PARAMETERS
c
      model='Photo 7'
      call photsou (model)
c
c     Possible to have no photons
c
c      qhlo = dlog(qht+epsilon)
c      rechy = 2.6d-13
c     reclo = dlog(rechy)
c
      call srcsummary (6, 1, soupho)
c
c Default outward only intgration assuming spherical symmetric nebula
c
      jtrans='OUTW'
c
   60 write (*,70)
   70 format(/' Spherical or plane II  geometry (s/p) : ',$)
      read (*,20) jgeo
      call toup (jgeo(1:1), jgeo)
c
      if ((jgeo.ne.'S').and.(jgeo.ne.'P')) goto 60
c
      rstar=1.d0
      astar=rstar*rstar
c
c     ***IF SYMMETRY IS plane parallel :
c
      if (jgeo.eq.'P') then
c
        jtrans='LODW'
c
        write (*,80)
   80  format(/' Radiative Transfer Mode ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    T  :  Two Sided, outward only (symmetrical, default)',/,
     & '    O  :  One Sided, (precursors)',/,
     & ' :: ',$)
        read (*,20) ilgg
        call toup (ilgg(1:1), ilgg)
c
        if ((ilgg.ne.'O').and.(ilgg.ne.'T')) ilgg='T'
        if (ilgg.eq.'T') then
          jtrans='OUTW'
        else
          jtrans='LODW'
        endif
c     ***IF SYMMETRY IS plane parallel :
      endif
c
c     ***IF SYMMETRY IS SPHERICAL :
c
      if (jgeo.eq.'S') then
c
        jtrans='OUTW'
c
c     1st total energy in Inu
c
        call intinu (soupho, qhi, qhei, qheii, qht)
        call inulum (soupho, blum, ilum)
c
        rstar=1.d0
        astar=rstar*rstar
c
        if (blum.gt.0.d0) then
c
   90     write (*,100)
  100 format(/' Define the source size or luminosity ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    R  :  By Radius',/,
     & '    L  :  By Luminosity',/,
     & '    P  :  By Ionizing Photons',/,
     & ' :: ',$)
          read (*,20) ilgg
          call toup (ilgg(1:1), ilgg)
c
          if ((ilgg.ne.'R').and.(ilgg.ne.'L').and.(ilgg.ne.'P')) goto
     &     90
c
          if (ilgg.eq.'R') then
c
  110      format(/,/,' Give photoionisation source radius'
     &            ,'(Rsun=6.957d10 IAU2015)',/
     &            ,' (in solar units (<1.e5) or in cm (>1.e5) : ',$)
            write (*,110)
            read (*,*) rstsun
            if (rstsun.le.0.0d0) goto 90
            if (rstsun.lt.1.d5) then
              rstar=rsun*rstsun
            else
              rstar=rstsun
              rstsun=rstar/rsun
            endif
c
            astar=fpi*rstar*rstar
            blum=astar*blum
            ilum=astar*ilum
c
  120 format(//
     & ' ********************************************************',/,
     & '  Source Total Luminosity              : ', 1pg12.5/
     & '  Source Ionising (13.6eV+) Luminosity : ', 1pg12.5/
     & '  Source Ionising (13.6eV+) Photons    : ', 1pg12.5/
     & ' ********************************************************',/)
            write (*,120) blum,ilum,qht*astar
c
c     end def by radius
c
          endif
c
          if (ilgg.eq.'L') then
c
  130       format(//' Total or Ionising Luminosity (T/I)',$)
            write (*,130)
            read (*,20) ilgg
            call toup (ilgg(1:1), ilgg)
c
            if ((ilgg.ne.'I').and.(ilgg.ne.'T')) ilgg='I'
c
  140       format(//' Give ionising source luminosity ',/,
     & ' (log (<100) or ergs/s (>100)) : ',$)
  150       format(//' Give bolometric source luminosity ',/,
     & ' (log (<100) or ergs/s (>100)) : ',$)
            if (ilgg.eq.'T') then
              write (*,150)
            else
              write (*,140)
            endif
            read (*,*) starlum
            if (starlum.lt.1.d2) starlum=10.d0**starlum
c
c
            if (ilgg.eq.'I') then
              rstsun=dsqrt(starlum/(ilum*fpi))
            else
              rstsun=dsqrt(starlum/(blum*fpi))
            endif
c
c
  160  format(//,
     & ' ********************************************************',/,
     & '  Source Radius : ',1pg12.5,' cm',/,
     & ' ********************************************************',/)
            write (*,160) rstsun
            astar=fpi*rstsun*rstsun
            blum=astar*blum
            ilum=astar*ilum
            write (*,120) blum,ilum,qht*astar
c
            rstar=rstsun
            rstsun=rstar/rsun
            astar=fpi*rstar*rstar
c
c     end def by luminosity
c
          endif
c
          if (ilgg.eq.'P') then
c
  170       format(//' Source ionising photon rate ',/,
     & ' (log (<100) or photons/s (>100)) : ',$)
            write (*,170)
            read (*,*) starlum
            if (starlum.lt.1.d2) starlum=10.d0**starlum
c
c     P/s = qht*fpi*rstsun*rstsun
c
            rstsun=dsqrt(starlum/(qht*fpi))
c
            write (*,160) rstsun
            astar=fpi*rstsun*rstsun
            blum=astar*blum
            ilum=astar*ilum
            write (*,120) blum,ilum,qht*astar
c
            rstar=rstsun
            rstsun=rstar/rsun
            astar=fpi*rstar*rstar
c
c    end source radius with photons
c
          endif
c
c    end blum>0
c
        endif
c
      endif
c
c
  180 format(//,
     & ' ********************************************************',/,
     & '  Setting the Physical Structure  ',/,
     & ' ********************************************************'//)
      write (*,180)
c
c     ***THERMAL STRUCTURE
c
      jthm='S'
      if (expertmode.gt.0) then
c
  190   write (*,200)
  200 format(//,
     & '  Thermal Structure ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    S  : Self-consistent,  (recommended)',/,
     & '    T  : isoThermal, (non-physical, not recommended)',/,
     & ' :: ',$)
        read (*,20) jthm
        call toup (jthm(1:1), jthm)
c
        if ((jthm.ne.'S').and.(jthm.ne.'T').and.(jthm.ne.'F')
     &   .and.(jthm.ne.'I')) goto 190
c
        if (jthm.eq.'T') then
c
  210    format(/,' Set the fixed electron temperature (<=10 as log):'
     &         ,$)
          write (*,210)
          read (*,*) fixtemp
c
          if (fixtemp.le.10.d0) fixtemp=10.d0**fixtemp
c
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***DENSITY BEHAVIOR
c
c default partial pressure of Hydrogen, ignored if not isobaric
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ponk=1.d6
c
  220 write (*,230)
  230 format(/'  Density Structure ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    C  : isoChoric,  (const volume)',/,
     & '    B  : isoBaric, (const pressure)',/,
     & '    F  : Functional form, F(r)',/,
     & ' :: ',$)
      read (*,20) jden
      call toup (jden(1:1), jden)
c
      if ((jden.ne.'C').and.(jden.ne.'B').and.(jden.ne.'F')) goto 220
c
      if (jden.eq.'F') then
        write (*,240)
  240    format(/'  Density Function ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' n(x) = n0*[fx*exp((x-a)/r0) + fp*((x/a)**b)]+c  ',/,
     & ' Give parameters n0 fx fp a b c & r0 (cgs): ',$)
        read (*,*) ofac,xfac,pfac,afac,bfac,cfac,scalen
        dhn=ofac
      endif
c
c     Estimate strom radius for nominal density
c     if we have any photons that is
c
c     if isobaric, get pressure regime
c
      jpres='A'
      if (jden.eq.'B') then
  250   write (*,260)
  260 format(/'  Choose a Pressure Structure ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : isoBaric,  P(r) = P0 + P_rad(r)',/,
     & '    F  : Function,  P(r) = P0*(r/r0)^a + P1 + P_rad(r)',/,
     & ' :: ',$)
        read (*,20) jpres
        call toup (jpres(1:1), jpres)
        if ((jpres.ne.'A').and.(jpres.ne.'F')) goto 250
c
        if (jpres.eq.'A') then
  270    format(/' Give TOTAL gas pressure: P0/k ',
     & ' (n.T = P/k, <10 as log) : ',$)
          write (*,270)
          read (*,*) ponk
          if (ponk.le.10.d0) ponk=10.d0**ponk
c
          tinner=1.d4
c
  280    format(//' Give an estimate of the initial temperature :',//
     & '   : Use 1.0e4 if not known. Initial estimates of density,',/
     & '   : and radius will be uncertain until the model starts.',/
     & '   : Ionisation uncertainty will always make nH estimates',/
     & '   : here approximate, but total pressure will be accurate.',/
     & '   : Initial Q, P/k, T and nH is solved in model,',/
     & '   : so this estimate is not crucial.',//
     & ' [T ~1e4, <10 as log] : ',$)
          write (*,280)
          read (*,*) tinner
          if (tinner.le.10.d0) tinner=10.d0**tinner
c
          dnt=ponk/tinner
          dhn=dnt/2.2d0
c
  290 format(/
     & ' ********************************************************',/,
     & '  Total number density at',1pg10.3,'* ',1pg10.3,/
     & '  Est. Hydrogen number density at',1pg10.3,'* ',1pg10.3/
     & ' ********************************************************',/)
          write (*,290) tinner,dnt,tinner,dhn
c
c     jden = B isobaric
c
        endif
c
        if (jpres.eq.'F') then
          rp0=0.d0
  300    format(/' Give inner gas pressure: P0/k ',
     & ' (n.T = P/k, <10 as log) : ',$)
          write (*,300)
          read (*,*) ponk0
          if (ponk0.le.10.d0) ponk0=10.d0**ponk0
          ponk=ponk0
  310    format(/' Give base gas pressure: P1/k ',
     & ' (n.T = P/k, <10 as log) : ',$)
          write (*,310)
          read (*,*) ponk1
          if (ponk1.le.10.d0) ponk1=10.d0**ponk1
  320    format(/' Give pressure scaling power: p*r^a ',
     & ' (linear, 0.0  to -2.0 ) : ',$)
          write (*,320)
          read (*,*) palpha
c
          tinner=1.d4
c
          write (*,280)
          read (*,*) tinner
          if (tinner.le.10.d0) tinner=10.d0**tinner
c
          dnt=(ponk0+ponk1)/tinner
          dhn=dnt/2.2d0
c
          write (*,290) tinner,dnt,tinner,dhn
c
        endif
c
      endif
c
  330 if (jden.eq.'C') then
  340    format(//' Give hydrogen number density : ',$)
        write (*,340)
        read (*,*) dhn
        tinner=1.0d4
      endif
c
      dht=densnum(dhn)
c
  350 format(//
     & ' ********************************************************',/,
     & '   Mean Ion Density (N):',1pg12.5/
     & '   Mean H  Density (dh):',1pg12.5/
     & ' ********************************************************',/)
c
      write (*,350) dht,dhn
c
  360 format(/' Set filling factor (0<f<=1) : ',$)
      write (*,360)
      read (*,*) fin
      if (((fin.le.0.0d0).or.(fin.gt.1.0d0)).or.(dhn.le.0.0d0)) goto
     &330
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***IF GEOMETRY SPHERICAL :
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ((jgeo.eq.'S')) then
c
c       Spherical Geometry
c
        if (qht.gt.0.d0) then
c
          volstromhb=(qht*astar)/(dhn*dht*alphahb*fin)
          rstromhb=(3.d0*volstromhb/fpi)**0.333333333333d0
          volstromheb=(qheii*astar)/(zion(2)*dhn*dht*alphaheb*fin)
          rstromheb=(3.d0*volstromheb/fpi)**0.333333333333d0
c
          rmax=dmax1(rstromhb,rstromheb)
c
          wdil=fdilu(rstar,(rmax+rstar)*0.5d0)
          qhdnin=qht/dht
          qhdnav=(4.d0*qht/dht)*wdil
c
          unin=qhdnin/cls
          unav=qhdnav/cls
c
          qhdhin=qht/dhn
          qhdhav=qhdnav*(dht/dhn)
c
          uhin=qhdhin/cls
          uhav=qhdhav/cls
c
  370     write (*,380) rstromhb,rstromheb,qht*astar,qheii*astar,qht,
     &     qheii,qhdnin,qhdhin,qhdnav,qhdhav,unin,uhin,unav,uhav
  380 format(//,
     & ' ********************************************************',/,
     & '   Filled Sphere Parameters:',/,
     & ' ********************************************************',/,
     & '  Estimated HII   Stromgren radius:',1pg10.3,' cm.',/,
     & '  Estimated HeIII Stromgren radius:',1pg10.3,' cm.',/,
     & ' ********************************************************',/,
     & '  Photon Luminosity LQH    :',1pg10.3, ' Phot/s',/,
     & '  Photon Luminosity LQHeII :',1pg10.3, ' Phot/s',/,
     & '  Photon Flux       FQH    :',1pg10.3, ' Phot/cm2/s',/,
     & '  Photon Flux       FQHeII :',1pg10.3, ' Phot/cm2/s',/,
     & ' ********************************************************',/,
     & '        QHDN inner   : ',1pg12.5,' cm/s',/,
     & '        QHDH inner   : ',1pg12.5,' cm/s',/,
     & '          <QHDN>     : ',1pg12.5,' cm/s',/,
     & '          <QHDH>     : ',1pg12.5,' cm/s',/,
     & '        U(N) inner   : ',1pg12.5,/,
     & '        U(H) inner   : ',1pg12.5,/,
     & '          <U(N)>     : ',1pg12.5,/,
     & '          <U(H)>     : ',1pg12.5,/,
     & ' ********************************************************',//,
     & ' Set initial radius in terms of distance or Q(N), U(N),',
     & ' Q(H), or U(H) (d/q/n/h/u):',$)
          read (*,20) ilgg
          call toup (ilgg(1:1), ilgg)
c
          if ((ilgg.ne.'Q').and.(ilgg.ne.'H').and.(ilgg.ne.'D')
     &     .and.(ilgg.ne.'U').and.(ilgg.ne.'N')) goto 370
c
          if (ilgg.eq.'U') then
            write (*,390)
  390       format(//,
     & ' Set U(H) at inner radius (<=0 as log) : ',$)
            read (*,*) uinit
c
            if (uinit.le.0) uinit=10.d0**uinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(uhav/uinit)
c
          endif
c
          if (ilgg.eq.'N') then
            write (*,400)
  400       format(//,
     & ' Set U(N) at inner radius (<=0 as log) : ',$)
            read (*,*) uinit
c
            if (uinit.le.0) uinit=10.d0**uinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(unav/uinit)
c
          endif
c
          if (ilgg.eq.'Q') then
            write (*,410)
  410       format(//,
     & ' Set QHDN at inner radius (<= 100 as log) : ',$)
            read (*,*) qinit
c
            if (qinit.le.100) qinit=10.d0**qinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(qhdnav/qinit)
c
          endif
c
          if (ilgg.eq.'H') then
            write (*,420)
  420       format(//,
     & ' Set QHDH at inner radius (<= 100 as log) : ',$)
            read (*,*) qinit
c
            if (qinit.le.100) qinit=10.d0**qinit
c
c     Assume large distance from source and scale by (rmax+rstar)/2
c
            remp=((rmax+rstar)*0.5d0)*dsqrt(qhdhav/qinit)
c
          endif
c
          if (ilgg.eq.'D') then
c
c     give R_0 directly
c
            write (*,430)
  430       format(//' Set the initial radius ',
     & ' (in cm (>1) or in fraction of Stromgren radius) : ',$)
            read (*,*) remp
c
            if (remp.lt.0.0d0) goto 370
            if (remp.le.1) remp=remp*rmax
            if (remp.lt.rstar) remp=rstar
c
          endif
c
          rstromhb=((rstromhb**3.d0)+(remp**3.d0))**0.333333333333d0
          rmax=rstromhb
c
          rstromheb=((rstromheb**3.d0)+(remp**3.d0))**0.333333333333d0
c
          call inulum (soupho, blum, ilum)
c
          wdpl=fdilu(rstar,remp)
          qhdnin=4.d0*qht/dht*wdpl
          unin=qhdnin/cls
          qhdhin=qhdnin*(dht/dhn)
          uhin=qhdhin/cls
          blum=blum*wdpl
          ilum=ilum*wdpl
          wdpl=fdilu(rstar,(rmax+remp)*0.5d0)
          qhdnav=4.d0*qht/dht*wdpl
          unav=qhdnav/cls
          qhdhav=qhdnav*(dht/dhn)
          uhav=qhdhav/cls
  440 format(//,
     & ' ********************************************************',/,
     & '   Partially Filled Sphere Parameters*',/,
     & ' ********************************************************',/,
     & '  Empty inner radius :',1pg12.5,' cm',/,
     & '  Outer HII   radius :',1pg12.5,' cm.',/,
     & '  Outer HeIII radius :',1pg12.5,' cm.',/,
     & ' ********************************************************',/,
     & '        QHDN inner   : ',1pg12.5,' cm/s',/,
     & '        QHDH inner   : ',1pg12.5,' cm/s',/,
     & '          <QHDN>     : ',1pg12.5,' cm/s',/,
     & '          <QHDH>     : ',1pg12.5,' cm/s',/,
     & '        U(N) inner   : ',1pg12.5,/,
     & '        U(H) inner   : ',1pg12.5,/,
     & '          <U(N)>     : ',1pg12.5,/,
     & '          <U(H)>     : ',1pg12.5,/,
     & '   Total intensity   : ',1pg12.5,' erg/s/cm2' /,
     & '   Ionizing intensity: ',1pg12.5,' erg/s/cm2',/,
     & ' ********************************************************'/)
          write (*,440) remp,rstromhb,rstromheb,qhdnin,qhdhin,qhdnav,
     &     qhdhav,unin,uhin,unav,uhav,blum,ilum
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          radin=0.d0
          radout=khuge
c
          write (*,450)
  450    format(//' Volume integration over the whole sphere?',
     & ' (y/n) : ',$)
          read (*,20) ilgg
          call toup (ilgg(1:1), ilgg)
c
          if (ilgg.ne.'N') ilgg='Y'
c
          if (ilgg.eq.'N') then
  460       write (*,470)
  470       format(/' Integration through the line of sight ',
     & 'for a centered ring aperture.',/,
     & 'Give inner and outer radii (cm) :')
            read (*,*) radin,radout
            if ((radin.lt.0.0d0).or.(radin.ge.(0.999d0*radout))) goto
     &       460
          endif
c
c     end with photons (qht>0)
c
        endif
c
c   end Spherical
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***IF GEOMETRY PLANE-PARALLEL :
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jgeo.eq.'P') then
c
        rstar=1.0d0
        astar=1.0d0
        remp=0.0d0
c
c     1st total energy in Inu
c
        call intinu (soupho, qhi, qhei, qheii, qht)
        call inulum (soupho, blum, ilum)
c
        if (qht.gt.0.d0) then
c
c     Determine number of ionising photons
c
          wdpl=0.5d0
          call srcsummary (6, 1, soupho)
c
  480     write (*,490)
  490    format(/,'  Ionizing Flux at inner edge by: ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    B  : Total Flux (erg/s/cm^2)',/,
     & '    I  : Ionising Flux (erg/s/cm^2)',/,
     & '    F  : Ionising Photons FQHI (phot/s/cm^2)',/,
     & '    Q  : Ionisation parameter QHDN ',/,
     & '    H  : Ionisation parameter QHDH ',/,
     & '    N  : Ionisation parameter U(N) ',/,
     & '    U  : Ionisation parameter U(H) ',/,
     & '    X  : No change (Use current flux)',/,
     & ' :: ',$)
          read (*,20) ilgg
          call toup (ilgg(1:1), ilgg)
c
          if (ilgg.eq.'B') then
            write (*,500)
  500       format(/,'Set Bolometric flux at inner edge of cloud:')
            read (*,*) scale
            scale=scale/blum
          elseif (ilgg.eq.'I') then
            write (*,510)
  510       format(/,'Set Ionising flux at inner edge of cloud:')
            read (*,*) scale
            scale=scale/ilum
          elseif (ilgg.eq.'F') then
            write (*,520)
  520       format(/,
     & 'Set Ionising Photon flux at inner edge (<=100 as log):')
            read (*,*) scale
            if (scale.le.100.d0) scale=10**scale
            scale=scale/qht
          elseif (ilgg.eq.'U') then
            write (*,530)
  530     format(/,'Set U(H) at inner edge (<=0 as log):')
            read (*,*) scale
            if (scale.le.0.d0) scale=10**scale
            scale=scale*dhn*cls/qht
          elseif (ilgg.eq.'N') then
            write (*,540)
  540     format(/,'Set U(N) at inner edge (<=0 as log):')
            read (*,*) scale
            if (scale.le.0.d0) scale=10**scale
            scale=scale*dht*cls/qht
          elseif (ilgg.eq.'Q') then
            write (*,550)
  550       format(/,'Set QHDN at inner edge (<=100 as log):')
            read (*,*) scale
            if (scale.le.100.d0) scale=10**scale
            scale=scale*dht/qht
          elseif (ilgg.eq.'H') then
            write (*,560)
  560       format(/,'Set QHDH at inner edge (<=100 as log):')
            read (*,*) scale
            if (scale.le.100.d0) scale=10**scale
            scale=scale*dhn/qht
          elseif (ilgg.eq.'X') then
            scale=1.d0
          else
            goto 480
          endif
c
          do j=1,infph
            soupho(j)=soupho(j)*scale
          enddo
c
          call intinu (soupho, qhi, qhei, qheii, qht)
          call inulum (soupho, blum, ilum)
c
          call srcsummary (6, 1, soupho)
c
          qhdnav=qht/dht
          unav=qhdnav/cls
          qhdhav=qht/dhn
          uhav=qhdhav/cls
  570  format(//,
     & ' ********************************************************',/,
     & ' ** Inner Boundary Ionisation Parameters ****************',/,
     & ' ********************************************************',/,
     & '  Bolometric Flux : ',1pg12.5,' (erg/s/cm^2)',/,
     & '  Ionising Flux   : ',1pg12.5,' (erg/s/cm^2)',/,
     & '  Ionising Phot Flux FQ     : ',1pg12.5,' (phots/s/cm^2)',/,
     & '  Ionisation parameter QHDN : ',1pg12.5,'(cm/s)',/,
     & '  Ionisation parameter QHDH : ',1pg12.5,'(cm/s)',/,
     & '  Ionisation parameter U(N) : ',1pg12.5/
     & '  Ionisation parameter U(H) : ',1pg12.5/
     & ' ********************************************************',/)
          write (*,570) blum,ilum,qht,qhdnav,qhdhav,unav,uhav
c
          volstromhb=qht/(dhn*dht*alphahb*fin)
          rstromhb=volstromhb
          rmax=rstromhb
          volstromheb=qheii/(zion(2)*dhn*dht*alphaheb*fin)
          rstromheb=volstromheb
c
  580     write (*,590)
  590    format(/' Set the geometrical dilution factor (<=0.5) : ',$)
          read (*,*) wdpl
c
          if ((wdpl.le.0.0d0).or.(wdpl.gt.0.5d0)) goto 580
c
          rmax=rmax*wdpl
c
          qhdnin=((2.0d0*qht)*wdpl)/dht
          qhdnav=qhdnin
c
          unin=((2.0d0*qht)*wdpl)/(cls*dht)
          unav=unin
c
          qhdhin=((2.0d0*qht)*wdpl)/dhn
          qhdhav=qhdhin
c
          uhin=((2.0d0*qht)*wdpl)/(cls*dhn)
          uhav=uhin
c
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  600 frlum=0.0d0
      tm00=0.0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***EQUILIBRIUM OR FINITE AGE ASSUMPTION
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  610 format(//,
     & ' ********************************************************',/,
     & '  Ionisation Balance Calculations',/,
     & ' ********************************************************'/)
      write (*,610)
  620 write (*,630)
  630 format(//'  Choose type',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    E  : Equilibrium balance.',/,
     & '    F  : Finite source lifetime.',/,
     & '    P  : Post-Equilibrium decay.',/,
     & '  :: ',$)
      read (*,20) jeq
      call toup (jeq(1:1), jeq)
c
      if ((jeq.ne.'E').and.(jeq.ne.'F').and.(jeq.ne.'P')) goto 620
c
c     forced CIE for given tprofile.
c     useful for hot x-ray bubbles.
c
      if (jeq.eq.'C') then
c
        write (*,640)
  640    format(/'  Fixed Thermal Profile:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',//
     & ' f(x) = t0*[fx*exp((x-ta)/r0) + fp*((x/ta)**tb)]+tc  ',/
     & ' (with r in r0 units , t0<10K as log)',/
     & ' Give parameters t0 ta tb tc r0 : ',$)
        read (*,*) tofac,txfac,tpfac,tafac,tbfac,tcfac,tscalen
      endif
c
c     non equilibrium
c
      if ((jeq.ne.'E').and.(jeq.ne.'C')) then
        rcin=0.9
        if (jgeo.eq.'S') then
          dti=-(dlog(1.0d0-(rcin**3))/(alphahb*dht))
        else
          dti=-(dlog(1.0d0-rcin)/(alphahb*dht))
        endif
        dti=(rmax/3.d10)+dmax1(dti,rmax/3.d10)
        dthre=(1.2/dht)/3d-13
c
c     finite lifetime
c
        if (jeq.eq.'F') then
  650     write (*,660) rmax,dti
  660       format(/'  Source turn on and finite lifetime:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' (Hydrogen Stromgrem radius of the order of',1pg10.3,
     & ' cm',/,' Characteristic time scale for ionisation: ',1pg8.1,
     & ' sec)',/,' Elapsed time since source turned on (sec)',/,
     & ' and give life-time of ionising source (sec) : ',$)
          read (*,*) telap,tlife
          if ((telap.le.0.0d0).or.(tlife.le.0.0d0)) goto 650
          tm00=100.0d0
  670     write (*,680)
  680      format(//' Initial temperature of neutral gas : ',$)
          read (*,*,err=670) tm00
          if (tm00.lt.1.0d0) goto 670
        else
          tlife=0.0d0
  690     write (*,700) rmax,dthre
  700       format(/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Source turn-off and final luminosity',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' (Hydrogen Stromgrem radius of the order of',1pg10.3,'(cm)',/,
     & ' Time scale for recombination: ',1pg8.1,'(sec)',/,
     & ' Give elapsed time since source turned off (sec)',/,
     & ' and fractional final lumnosity of the source (>=0) :: ',$)
          read (*,*) telap,frlum
          if ((telap.le.0.0d0).or.(frlum.lt.0.0d0)) goto 690
        endif
      endif
c
  710 write (*,720)
  720 format(//' Set the step photon absorption fraction: ',$)
      read (*,*) dtau0
      if ((dtau0.le.0.0d0)) goto 710
c
c     standard convergence limits
c
      difma=0.05d0
      dtlma=0.10d0
      dhlma=0.050d0
c
      if (expertmode.gt.0) then
  730    format(/'Alter the convergence criteria ?',
     & ' (not recommended) : ',$)
        write (*,730)
        read (*,20) ilgg
        call toup (ilgg(1:1), ilgg)
        if (ilgg.ne.'Y') ilgg='N'
c
        if (ilgg.eq.'Y') then
          write (*,740)
  740    format(/,' Set the convergence criteria, the maximum changes'
     &   ,'allowed',/
     &   ,' for three parameters: ion population, te and hydrogen'
     &   ,'density.',/
     &   ,' Convergence criteria: difma (0.05),dtlma (0.10)'
     &   ,'dhlma(0.015) : ',$)
c
          read (*,*) difma,dtlma,dhlma
        endif
c
c     end expertmode
c
      endif
c
c     ***PRINT SET UP
c
      if (jeq.eq.'E') then
        write (*,750) dtau0
  750    format(//' Model Summary :',/
     &          ,'    Mode  : Thermal and Ionic Equilibrium ',/
     &          ,'    dTau  :',0pf8.5)
      elseif (jeq.eq.'F') then
        write (*,760) telap,tlife,dtau0
  760 format(//' Model Summary :',/
     & ,'    Mode        :  Non-Equilibrium + Finite Source Life',/
     & ,'    Age         :',1pg10.3,' sec',/
     & ,'    Source Life :',1pg10.3,' sec',/
     & ,'    dTau        :',0pf8.5)
      else
        write (*,770) telap,dtau0
  770 format(//' Model Summary :',/
     & ,'    Mode        :  Equilibrium + Source Switch off',/
     & ,'    Switch off  :',1pg10.3,' sec',/
     & ,'    dTau        :',0pf8.5)
      endif
      if (usekappa) then
  780    format(/
     & ' ********************************************************',/,
     & ' Kappa Electron Distribution Enabled :',/,
     & ' Electron Kappa : ',1pg11.4,/
     & ' ********************************************************',/)
        write (*,780) kappa
      endif
c
  790 format(/
     & ' ********************************************************',/,
     & '   Radiation Field and Parameters*',/,
     & ' ********************************************************',/,
     & t5, ' At estimated T_inner :',1pg10.3,' K',/,
     & t5,' Rsou.',t18,' Remp.',t32,' Rmax',t46,' <DILU>',/,
     & t5, 1pg10.3,t18,1pg10.3,t32,1pg10.3,t46,1pg10.3/
     & t5,' <Hdens>',t18,' <Ndens>',t32,' Fill Factor',/
     & t5, 1pg10.3,t18,1pg10.3,t32,0pf9.6,//
     & ' ********************************************************',/,
     & t5,' FQ tot',t18,' Log FQ',t32,' LQ tot',t46,' Log LQ',/
     & t5, 1pg10.3,t18,0pf8.3,t32,1pg10.3,t46,0pf8.3,/
     & t5,' Q(N) in',t18,' Log Q(N)',t32,' <Q(N)>',t46,' Log<Q(N)>',/
     & t5, 1pg10.3,t18,0pf8.3,t32,1pg10.3,t46,0pf8.3,/
     & t5,' Q(H) in',t18,' Log Q(H)',t32,' <Q(H)>',t46,' Log<Q(H)>',/
     & t5, 1pg10.3,t18,0pf8.3,t32,1pg10.3,t46,0pf8.3,/
     & t5,' U(N) in',t18,' Log U(N)',t32,' <U(N)>',t46,' Log<U(N)>',/
     & t5, 1pg10.3,t18,0pf8.3,t32,1pg10.3,t46,0pf8.3,/
     & t5,' U(H) in',t18,' Log U(H)',t32,' <U(H)>',t46,' Log<U(H)>',/
     & t5, 1pg10.3,t18,0pf8.3,t32,1pg10.3,t46,0pf8.3,/
     & ' ********************************************************',/)
c
      write (*,790) tinner,rstar,remp,rmax,wdpl,dhn,dht,fin,qht,
     &dlog10(qht),astar*qht,dlog10(astar*qht),qhdnin,dlog10(qhdnin),
     &qhdnav,dlog10(qhdnav),qhdhin,dlog10(qhdhin),qhdhav,dlog10(qhdhav),
     &unin,dlog10(unin),unav,dlog10(unav),uhin,dlog10(uhin),uhav,
     &dlog10(uhav)
c
c     ***TYPE OF EXIT FROM THE PROGRAM
c
  800 ielen=1
      jpoen=1
      tend=10.0d0
      fren=0.01d0
      diend=0.0d0
      tauen=0.0d0
c
c
  810 format(//
     & ' ********************************************************',/,
     & '  Boundry Conditions  ',/,
     & ' ********************************************************'/)
      write (*,810)
c
  820 format(/'  Model Ending Boundary Conditions : ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  :   Radiation bounded, HII <',2pf5.2,'%',/,
     & '    B  :   Ionisation bounded, XII < Y',/,
     & '    C  :   Temperature bounded, Tmin.',/,
     & '    D  :   Optical depth limited, Tau',/,
     & '    E  :   Density bounded, Distance from inner edge',/,
     & '    F  :   Column Density limited, atom, ion',/,
     & '    H  :   Total H Column Density limited',/,
     & '       :   ',/,
     & '    R  :   (Reinitialise)',/,
     & '    G  :   (Reset geometry only)',/,
     & ' :: ',$)
      write (*,820) fren
      read (*,20) jend
      call toup (jend(1:1), jend)
c
      if ((jend.ne.'A').and.(jend.ne.'B').and.(jend.ne.'C')
     &.and.(jend.ne.'D').and.(jend.ne.'E').and.(jend.ne.'F')
     &.and.(jend.ne.'R').and.(jend.ne.'G').and.(jend.ne.'H')) goto 800
c
      if (jend.eq.'R') return
c
      if (jend.eq.'G') goto 600
c
      if ((jend.eq.'B').or.(jend.eq.'D')) then
  830   write (*,840)
  840    format(/' Applies to element (Atomic number) : ',$)
        read (*,*) ielen
        if (zmap(ielen).eq.0) goto 830
        ielen=zmap(ielen)
        jpoen=1
        if (arad(2,ielen).le.0.0d0) jpoen=2
      endif
c
      if (jend.eq.'B') then
  850   write (*,860) elem(ielen)
  860    format(/' Give the final ionisation fraction of ',a2,' : ',$)
        read (*,*) fren
        if ((fren.lt.0.0d0).or.(fren.gt.1.0d0)) goto 850
      endif
c
      if (jend.eq.'D') then
  870   write (*,880) elem(ielen)
  880 format(/,
     & ' Give the final optical depth at threshold of ',a2,':',$)
        read (*,*) tauen
        if (tauen.le.0.0d0) goto 870
      endif
c
      if (jend.eq.'C') then
  890   write (*,900)
  900    format(/' Give the final temperature (<10 as log): ',$)
        read (*,*) tend
        if (tend.lt.1.0d0) goto 890
      endif
c
      if (jend.eq.'E') then
  910   write (*,920)
  920    format(/' Give the distance or radius at which the density',
     & ' drops: ',/,
     & ' (in cm (>1E6) or as a fraction of the Stromgren',/,
     & ' radius (<1E6)) : ',$)
        read (*,*) diend
        if (diend.lt.1.d6) diend=diend*rmax
        diend=remp+diend
        if (diend.le.remp) goto 910
      endif
c
      if (jend.eq.'F') then
        write (*,930)
  930    format(/' Give the final column density (<100 as log): ',$)
        read (*,*) colend
        if (colend.lt.1.0d2) colend=10.d0**colend
  940   write (*,950)
  950    format(/' Applies to element (Atomic number) : ',$)
        read (*,*) ielen
        if (zmap(ielen).eq.0) goto 940
        ielen=zmap(ielen)
  960   write (*,970)
  970    format(/' Applies to ion stage : ',$)
        read (*,*) jpoen
        if (jpoen.le.0) goto 960
        if (jpoen.gt.maxion(ielen)) jpoen=maxion(ielen)+1
      endif
c
      if (jend.eq.'H') then
        write (*,980)
  980   format(/' Give the total H column density (<100 as log): ',$)
        read (*,*) colend
        if (colend.lt.1.0d2) colend=10.d0**colend
      endif
c
  990 format(//
     & ' ********************************************************',/,
     & '  Output Requirements  ',/,
     & ' ********************************************************',/)
      write (*,990)
c
      p7pfx='photn'
      nprefix=5
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    New output options menu including reset like in shock5files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c sets unit numbers too
c
c      dynamic unit in wpsou
 1000 jspec='N'
c     dynamic unit in wpsou
      jsou='N'
c
c      bands
      jbnd='N'
c      Lines by wavelength
      jlin='N'
c      rates
      jrat='N'
c      All Ions File
      jall='N'
c      Ionisation of elements
      jiel='N'
c      Multi-level ions
      jiem='N'
c      Ion columns
      jcol='N'
c      first balance , dynamic unit in wbal
      jbal='N'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spectra: F-Lambda plots, ionising and optical default yes
c
 1010 outsettings='A'
c
c     Ion Balances
c
      if (jiel.eq.'Y') then
        call myappend (outsettings, '+B', outshow)
        outsettings=outshow
      endif
c
c    rates file
c
      if (jrat.eq.'Y') then
        call myappend (outsettings, '+C', outshow)
        outsettings=outshow
      endif
c
c    final downstream field
c
      if (jspec.eq.'Y') then
        call myappend (outsettings, '+D', outshow)
        outsettings=outshow
      endif
c
c    Ions/Columns
c
      if (jcol.eq.'Y') then
        call myappend (outsettings, '+E', outshow)
        outsettings=outshow
      endif
c
c    First balance
c
      if (jbal.eq.'Y') then
        call myappend (outsettings, '+F', outshow)
        outsettings=outshow
      endif
c
c     emission bands
c
      if (jbnd.eq.'Y') then
        call myappend (outsettings, '+H', outshow)
        outsettings=outshow
      endif
c
c    monitor lines
c
      if (jlin.eq.'Y') then
        call myappend (outsettings, '+L', outshow)
        outsettings=outshow
      endif
c
c    monitor multi-level ion lines
c
      if (jiem.eq.'Y') then
        call myappend (outsettings, '+M', outshow)
        outsettings=outshow
      endif
c
 1020 format(/,' Output Multi-Option Menu : ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' :: Current: ',a20,'                       ::',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard output - and Reset. ',/,
     & '    B  : Ion balance files.',/,
     & '    C  : Plasma Rates file.',/,
     & '    D  : Final source + Nebula spectra',/,
     & '    E  : Element Ions/Column structure',/,
     & '    F  : First balance',/,
     & '    H  : Emission in x-ray bands.',/,
     & '    L  : Monitor up to ',i3,' lines',/,
     & '    M  : Monitor up to ',i2,' multi-level ions',/,
     & '       :',//,
     & '    R  : Reset',/,
     & '    G  : Everything',/,
     & '    X  : Exit with Settings',/,
     & ' :: ',$)
 1030 write (*,1020) outsettings,mxmonlines,mxmonmlions
      read (*,20) ilgg
      call toup (ilgg(1:1), ilgg)
      write (*,*)
c
      if (ilgg.eq.'Q') ilgg='X'
c
      if ((ilgg.ne.'A').and.(ilgg.ne.'B').and.(ilgg.ne.'C')
     &.and.(ilgg.ne.'D').and.(ilgg.ne.'E').and.(ilgg.ne.'F')
     &.and.(ilgg.ne.'G').and.(ilgg.ne.'H').and.(ilgg.ne.'L')
     &.and.(ilgg.ne.'M').and.(ilgg.ne.'R').and.(ilgg.ne.'X')) goto 1030
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write (*,*) ilgg
      if (ilgg.eq.'A') goto 1000
      if (ilgg.eq.'B') jiel='Y'
      if (ilgg.eq.'C') jrat='Y'
      if (ilgg.eq.'D') jsou='Y'
      if (ilgg.eq.'D') jspec='Y'
      if (ilgg.eq.'E') jcol='Y'
      if (ilgg.eq.'F') jbal='Y'
      if (ilgg.eq.'H') jbnd='Y'
      if (ilgg.eq.'L') jlin='Y'
      if (ilgg.eq.'M') jiem='Y'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ilgg.eq.'R') goto 1000
      if (ilgg.ne.'X') goto 1010
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ilgg.eq.'G') then
c     dynamic unit in wpsou
        jspec='Y'
c     dynamic unit in wpsou
        jsou='Y'
c      bands
        jbnd='Y'
c      Lines by wavelength
        jlin='Y'
c      rates
        jrat='Y'
c      All Ions File
        jall='N'
c      Ionisation of elements
        jiel='Y'
c      Multi-level ions
        jiem='Y'
c      Ion columns
        jcol='Y'
c      first balance , dynamic unit in wbal
        jbal='Y'
      endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Get Output Settings for some types
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      jpfx='psend'
      if ((jsou.eq.'Y').or.(jspec.eq.'Y')) then
c
c     get final field file prefix
c
 1040   format (a64)
 1050   format(//' Prefix for final source files : ',$)
        write (*,1050)
        read (*,1040) jpfx
        nprefix=len(trim(jpfx))
        jpfx=jpfx(1:nprefix)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      njlines=0
      if (jlin.eq.'Y') then
c
 1060  format(//' Select up to ',i3,' lines by A (See spec list):',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
 1070  format(' Number of lines :: ',$)
        write (*,1060) mxmonlines
        write (*,1070)
        read (*,*) nl
        nl=min(max(nl,1),mxmonlines)
        njlines=nl
        do i=1,njlines
          read (*,*) lam
          lam=dabs(lam)
          emlinlist(i)=lam
        enddo
        do i=1,mxmonlines
c      will allow custom bins later
c      0.001A: comfortably covers the +/-0.0005A rounding from reading
c      a wavelength off spec2's 3-decimal-place output, while staying
c      far tighter than any real line-to-line separation encountered in
c      practice - so a wavelength copied from the line list matches
c      exactly, with no risk of picking up an unrelated nearby line.
          emlindeltas(i)=0.001d0
        enddo
        call speclocallineids (emlinlistatom, emlinlistion)
 1080  format(' Found Line: #',i3,' : ', a2,a6,1pg12.5)
        do i=1,njlines
          lam=emlinlist(i)
          write (*,1080) i,elem(emlinlistatom(i)),rom(emlinlistion(i)),
     &     lam
        enddo
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     default monitor all elements
c     used for ion fractions and column densities
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ieln=mxelem
      do idx=1,ieln
        iel(idx)=zmap(idx)
      enddo
c
      if ((jiel.eq.'Y').or.(jcol.eq.'Y')) then
c
 1090   format(//' Monitor ions/columns of max ',i2,' elements:',/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
 1100   format(' Enter number of elements to track : ',$)
 1110   format(/' Elements (Z) : ',$)
        write (*,1090) atypes
        write (*,1100)
        read (*,*) ieln
        ieln=min(max(1,ieln),atypes)
        write (*,1110)
        read (*,*) (iel(i),i=1,ieln)
        do i=1,ieln
          elok(ieln)=0
          do idx=1,atypes
            if (iel(i).eq.mapz(idx)) elok(i)=1
          enddo
        enddo
        nel=ieln
        do i=1,ieln
          if (elok(i).eq.0) nel=idx-1
        enddo
        if ((nel.lt.1).and.(jiel.eq.'Y')) then
          jiel='N'
        endif
        if ((nel.lt.1).and.(jcol.eq.'Y')) then
          jcol='N'
        endif
c
        if ((jiel.eq.'Y').or.(jcol.eq.'Y')) then
c
          ieln=nel
          do i=1,ieln
            iel(i)=zmap(iel(i))
          enddo
 1120      format(/' Monitoring :',30(x,a2),/)
          write (*,1120) (elem(iel(i)),i=1,ieln)
c
        endif
c
        if (jiel.eq.'Y') then
c
 1130      format(/,' Record all ions file (Y/N)? : ',$)
 1140     write (*,1130)
          read (*,20) jall
          call toup (jall(1:1), jall)
          if ((jall.ne.'Y').and.(jall.ne.'N')) goto 1140
c
        endif
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iemn=0
      do idx=1,mxmonmlions
        iem(idx)=idx
      enddo
c
      if (jiem.eq.'Y') then
c
 1150   format(//' Choose up to ',i3,' multi-level ions by id# :',/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::'/)
        write (*,1150) mxmonmlions
 1160   format(//' Number to monitor: ',$)
        write (*,1160)
        read (*,*) iemn
        iemn=min(max(iemn,1),mxmonmlions)
 1170   format(5(2x,i3,': ',a2,a6))
        do i=1,(nfmions/5)+1
          j=(i-1)*5
          kmax=min(j+5,nfmions)-j
          write (*,1170) (j+k,elem(fmatom(j+k)),rom(fmion(j+k)),k=1,
     &     kmax)
        enddo
c
 1180   format(' :: ',$)
        write (*,1180)
        write (*,*)
 1190   format(' Selected: ',i3,': ',a2,a6)
        do i=1,iemn
          read (*,*) idx
          idx=min(max(idx,1),nfmions)
          write (*,1190) idx,elem(fmatom(idx)),rom(fmion(idx))
          iem(i)=idx
        enddo
        write (*,*)
        write (*,*) 'Monitoring Ion Emission:'
        do i=1,iemn
          write (*,'("Ion:",x,i3,": ",a2,a6)') i,elem(fmatom(iem(i))),
     &     rom(fmion(iem(i)))
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      jbfx='p7bal'
      if (jbal.eq.'Y') then
c
c     get ion balance file prefix
c
 1200   format(//' Prefix for first ion balance file : ',$)
        write (*,1200)
        read (*,*) jbfx
        nprefix=len(trim(jbfx))
        jbfx=jbfx(1:nprefix)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     get runname
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 1210 format(//
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Set Name/code for this model: ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' :: ',$)
      write (*,1210)
      read (*,'(a)') runname
      np=mlen(runname)
      runname=runname(1:np)
c
c     remains of old vax batch system (not used)
c
      banfil='INTERACTIVE'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call p7filenames ()
      call createp7files ()
      call p7headers (banfil, dhn, fin)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call compph7 (dhn, fin, banfil, difma, dtlma, dhlma)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine compph7
c! XXXX - add one line purpose here
c! @param [in,out]   real*8       dhn  XXX-meaning
c! @param [in,out]   real*8       fin  XXX-meaning
c! @param [in,out]   real*8    banfil  XXX-meaning
c! @param [in,out]   real*8     difma  XXX-meaning
c! @param [in,out]  Unknown     dtlma  XXX-meaning
c! @param [in,out]   real*8     dhlma  XXX-meaning
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine compph7 (dhn, fin, banfil, difma, dtlma, dhlma)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       PHOTOIONISATION MODEL
c       DIFFUSE FIELD CALCULATED ; 'OUTWARD ONLY' INTEGRATION
c
c       PHOTO7: As for P5 but adds new integration scheme
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'p7blocks.inc'
c
      real*8 poppre(mxion, mxelem), popend(mxion, mxelem)
      real*8 ppre(mxion, mxelem),difma,dtlma,dhlma
      real*8 thist(3),dpw,dxp
      real*8 linfluxes(mxmonlines)
c
      real*8 aadn,agmax,agmu,aper,aperi,apero
      real*8 a,b,c,x,fx,fp,n0,r0,fv,e,f,r,p0,p1
      real*8 de,dedhma,dedhmi,deeq,chcksum
      real*8 dh,dhf,dhl,dhn,dhp,frp
      real*8 en,t_pk,qhratio
      real*8 f5007,f4363,f3727,f6300,f6584,f6720
      real*8 f6716,f6731,f673116
c
      real*8 difg  , difmi
      real*8 difde , difte
      real*8 difpre, difend
c
      real*8 disin,dison,dlos0,dlos1,dlosav
      real*8 dr,drp,gcharge
      real*8 dtau,dtauon,dtaux,dtco,dtco0,dtco1
      real*8 dtoff,dton,dton0,dton1,durmax
      real*8 dv,dvol,eqscal,exl,exla
      real*8 fhi,fhii,fi1,fidhcc,fin,frdw,frx,fthick
      real*8 prescc,prf,protim,rad,radi,recscal,rww
      real*8 t,taux,te1a,tef,telast,tempre
      real*8 tep,tex,texm,tii0,tprop,trea,rdvol,irdvol
      real*8 treach,treach0,treach1,tspr
      real*8 wd,wdil0,wdil1,wei
      real*8 ogas, oonh
      real*8 radpint, pres0, pres1
      real*8 popinttot
      real*8 ffuv,widnu,eng,habing,pahtest(mxinfph)
      real*8 pahsum, pahq0, pahqion, pahqfac
      real*8 q1,q2,q3,q4,q0
      real*8 def
c
      real*8 pop_initial (mxion, mxelem)
      real*8 pop_initial0(mxion, mxelem)
      real*8 pop_initial1(mxion, mxelem)
      real*8 te1_initial
      real*8 dh1_initial
      real*8 de1_initial
c
      integer*4 lut0,m,maxio,mma,n,nidh,niter
      integer*4 i,j,k,ndifatoms,atom, flen, np
      integer*4 idx, nt, itr,init_it
c
      integer*4 fuvmin,fuvmax,pahabsmax
c
      character jjd*4, pollfile*16,tab*4
      character banfil*64
      character imod*4, lmod*4, nmod*4,wmod*4,ispo*4
      character linemod*4, spmod*4, savemod*4
      character fn*128,fspl*128
      character pfx*64,caller*4,sfx*16
      character irfile*64,chargefile*64
      logical iexi
c
c     External Functions
c
      real*8 fdilu,feldens,fradpress,fpressu,frectim,fzgas
c
c     internal functions
c
      real*8 dif,frad,fsight,fring,fpresr
c
      dif(a,b)=dabs(dlog10(a+epsilon)-dlog10(b+epsilon))
      frad(x,n0,fx,fp,a,b,c,r0,fv,e,f)=n0*(fx*dexp((x-a)/r0)+fp*((x/a)**
     &b)+fv*(((r0-x)/e)**(-f)))+c
      fsight(aper,radi)=1.d0-(dcos(dasin(dmin1(1.d0,aper/(radi+(1.d-38*
     &aper)))))*(1.d0-(dmin1(1.d0,aper/(radi+(1.d-38*aper)))**2)))
      fring(aperi,apero,radi)=dmax1(0.d0,fsight(apero,radi)-
     &fsight(aperi,radi))
      fpresr(p0,p1,r0,a,r)=(p0*(r/r0)**a)+p1
c
      jspot='N'
      jcon='Y'
      ispo='SII'
      irfile=' '
      if (graindestructmode.eq.1) then
        grainmode=0
        grainadjust=1
      endif
      tab=char(9)
c
      ndifatoms=2
      q0=qhdhin
      te0=tinner
c
c      nt0 = time()
c      nt1 = 0
c
      mma=mxnsteps-1
      exla=5.d-5
      trea=1.d-5
      dton0=0.0d0
      dton1=0.0d0
      dtoff=0.0d0
      dtaux=0.d0
      drp=0.d0
      dedhmi=0.d0
      dedhma=0.d0
c
      treach0=0.0d0
      treach1=0.0d0
      texm=0.d0
      te0=1.d4
      te1=1.d4
      tep=1.d4
c
      dtco0=0.d0
      dtco1=0.d0
      dlos0=0.d0
      dlos1=0.d0
c
      thist(1)=0.d0
      thist(2)=0.d0
      thist(3)=0.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c open all main P7 files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call appendp7files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Begin Main Calculation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      lut0=0
c
      chargefile=' '
c
c  Set up for PAH existence
c
      fuvmin=1
      fuvmax=1
      pahabsmax=1
      if ((grainmode.eq.1).and.(pahmode.ne.0)) then
        k=1
        do while (photev(k).le.6.d0)
          k=k+1
        enddo
c     FUVMIN = photon bin corresponding to 6 eV energy
        fuvmin=k-1
        do while (photev(k).le.iph)
          k=k+1
        enddo
c     FUVMAX = photon bin corresponding to 13.6 eV energy
        fuvmax=k-1
        do while (photev(k).le.262.012d0)
          k=k+1
        enddo
c     PAHABSMAX = photon bin corresponding to 262.0 eV energy
c       maximal PAH absorption bin
        pahabsmax=k-1
c
c Calculate photon flux in FUV range using Weingartner &
c Draine (2001) ISRF
c
        ffuv=0.d0
        pahsum=0.d0
        do k=fuvmin,fuvmax
          eng=cphotev(k)
          widnu=widbinnu(k)
          pahsum=pahsum+pahnabs(k)*widnu
          if (eng.lt.9.26d0) then
            ffuv=ffuv+2.548d-18*eng**(-1.3322d0)/ev*pahnabs(k)*widnu
          elseif (eng.lt.11.2d0) then
            ffuv=ffuv+1.049d-16*eng**(-3.0d0)/ev*pahnabs(k)*widnu
          else
            ffuv=ffuv+4.126d-13*eng**(-6.4172d0)/ev*pahnabs(k)*widnu
          endif
        enddo
        pahq0=ffuv/pahsum
        write (*,10) pahq0
   10 format(/'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     & 'ISRF photon flux :',1pg12.5,/
     & '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   The routine does some initial setup, guesses at the parameters at
c   the inner boundary of the first spatial step, then for each spatial
c   step it guesses at the parameters at the outer boundary of the step,
c   iterates till these parameters converge, and then iterates over
c   spatial steps till the completion criteria are achieved.
c
c      write(*,*)'***COMPUTATION STARTS HERE'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      fi=fin
      fi1=fin
c
      fidhcc=fin*dhn
c
      dv=1.0d6
      rdvol=(rmax-remp)/300.d0
      irdvol=1.d0/rdvol
c
      dr=0.d0
c
      if (dtau0.gt.0.075d0) dtau0=0.075d0
c
      dtau=dtau0*0.125d0
      difg=1.d0
c
      if (jgeo.eq.'S') then
        vunilog=3.0d0*dlog10(rdvol)
      else
        vunilog=dlog10(rdvol)
      endif
      if ((jeq.eq.'E').or.(jeq.eq.'P')) then
c     Equilibrium structure (E=equilibrium, P=post-equilibrium decay,
        nmod='EQUI'
      else
c     Equilibrium structure (F=finite source)
        nmod='TIM'
      endif
c
c  set convergence criterion parameters
c  agmax = ? ; aadn = ?
c
      if (jden.ne.'B') then
c        Density structure (B=isobaric)
        agmax=0.05d0
        aadn=2.5d0
      else
c     Density structure (C=isochroic, F=function, I=input file)
        agmax=0.80d0
        aadn=5.0d0
      endif
      agmu=agmax
c
      radpint=0.d0
c     usesproto-ionisation
      pres0=fpressu(tinner,dhn,pop0)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*)'***ITERATE, M = STEP NUMBER'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      m=1
c
   20 jspot='N'
      niter=0
      nidh=0
      difmi=1.d30
      wei=(dtau/dtau0)**0.75d0
      dtau=dtau0*wei
      if (dtau.gt.0.075d0) dtau=0.075d0
      wei=3.d0
      agmu=(agmax+(wei*agmu))/(1.d0+wei)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*)'IF M=1 , CHECK CONVERGENCE FOR THE DENSITY'
c      write(*,*)'AND IONIC POP. AT THE INNER BOUNDARY'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (m.eq.1) then
c
        te0=tinner
        rp0=remp
c
        if (graindestructmode.eq.1) te0=3.d5
c
        if (jden.eq.'B') then
c     intialestimate
          if (jpres.eq.'A') dhn=(ponk/te0)/2.2d0
          if (jpres.eq.'F') then
            dhn=(fpresr(ponk0,ponk1,rp0,palpha,remp)/te0)/2.2d0
          endif
        endif
c
c     usesproto-ionisation
        pres0=fpressu(te0,dhn,pop0)
c
        if (jthm.eq.'T') te0=fixtemp
c
        dr=0.d0
        drp=dr
        if (jden.ne.'F') then
c     Density structure (C=isochroic, B=isobaric, I=input file)
c     DH0 = Hydrogen density at inner spatial step boundary
          dh0=dhn
        else
c     Density structure (F=function)
          dh0=frad(remp,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)
        endif
c
c Derive electron density using ionic populations in POP0.
c DE0 = Electron density at the inner spatial step boundary
c DEDHMA = max de/dh [total metallicity of the gas
c          (relative to hydrogen)]
c DEDHMI = min de/dh, allowing for ions prevented from recombining,
c         usually total metallicity of gas for all species with
c               -ve or 0 recombination rates from 2nd ion state???
c   - nope rs
c
        de0=feldens(dh0,pop0)
        dedhmi=0.d0
        dedhma=0.d0
        do j=1,atypes
          dedhma=dedhma+zion(j)
          if (arad(2,j).le.0.d0) dedhmi=dedhmi+zion(j)
        enddo
        call copypop (pop0, pop)
c
c     Derive geometrical dilution factor using photon source radius RSTAR
c     and radial distance from source of the inner step boundary RAD0.
c
        if (jgeo.eq.'S') then
          rad0=remp
          rp0=remp
          dis0=rad0
          wdil0=fdilu(rstar,rad0)
        endif
        if (jgeo.eq.'P') then
          rad0=0.d0
          dis0=remp
          rp0=1.0d0
          wdil0=wdpl
        endif
c
c  Loop back here if hydrogen density change is too great
c
        init_it=0
c
   30   if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh0)
        lmod='DW'
        call localem (te0, de0, dh0)
        call totphot2 (te0, dh0, rad0, dr, dv, wdil0, lmod)
c
c     calculate qhdn into zone
c
        q1=0.0d0
        q2=0.0d0
        q3=0.0d0
        q4=0.0d0
c
c     Integrate local mean intensities to derive number of H/He
c     ionizing photons. Returns q1=QAH1, q2=QAHEI, q3=QAHEII,
c     q4=QATOT.
c
        call intjnu (tphot, q1, q2, q3, q4)
c
        qhdnin=q4/(dh0*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh0)
        uhin=qhdhin/cls
c
        if ((jgeo.eq.'P').and.(jden.eq.'B')) then
          qhratio=q0/qhdhin
          do j=1,infph
            soupho(j)=soupho(j)*qhratio
            tphot(j)=soupho(j)*0.25d0*qhratio
          enddo
          call localem (te0, de0, dh0)
          call totphot2 (te0, dh0, rad0, dr, dv, wdil0, lmod)
          q1=0.0d0
          q2=0.0d0
          q3=0.0d0
          q4=0.0d0
          call intjnu (tphot, q1, q2, q3, q4)
          qhdnin=q4/(dh0*zen)
          unin=qhdnin/cls
          qhdhin=q4/(dh0)
          uhin=qhdhin/cls
        endif
c
        if ((jgeo.eq.'S').and.(jden.eq.'B')) then
          qhratio=qhdhin/q0
          remp=remp*dsqrt(qhratio)
          rp0=remp
          rad0=remp
          dis0=rad0
          wdil0=fdilu(rstar,rad0)
          call localem (te0, de0, dh0)
          call totphot2 (te0, dh0, rad0, dr, dv, wdil0, lmod)
          q1=0.0d0
          q2=0.0d0
          q3=0.0d0
          q4=0.0d0
          call intjnu (tphot, q1, q2, q3, q4)
          qhdnin=q4/(dh0*zen)
          unin=qhdnin/cls
          qhdhin=q4/(dh0)
          uhin=qhdhin/cls
        endif
c
        if (graindestructmode.eq.1) call adjustgrains (te0, qhdhin, 1)
c
        if ((jeq.eq.'E').or.(jeq.eq.'P').and.(jthm.eq.'S')) then
          dton0=1.d33
          call teequi2 (te0, te0, de0, dh0, dton0, nmod)
        else
          distim(1,0)=dis0
          distim(2,0)=dis0/cls
          treach0=distim(2,0)
          dton0=dmax1(0.d0,dmin1(tlife,telap-((2.d0*dis0)/cls)))
          jjd=jden
          jden='C'
          exl=-1.d0
          tex=0.95d0*tm00
          prescc=pres0+fradpress(0.0d0,dhn)
c
          call copypop (pop0, pop)
          call evoltem (tm00, te0, de0, dh0, prescc, dton0, exl, tex,
     &     lut0)
c
          jden=jjd
c
          if (dton0.le.0.d0) then
            write (*,40) dton0
   40         format(/,' Age was incorrectly set, it is shorter than',/,
     & ' needed at the first space step , DTon0 :',1pg10.3,/)
            return
          endif
          dtco0=(gammaEOSU*fpressu(te0,dh0,pop))/(eloss+1.d-36)
        endif
c
        dhp=dh0
        tep=te0
        dlos0=dlos
c
        if (jden.eq.'C') then
          dh0=dhn
        elseif (jden.eq.'F') then
          dh0=frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)
        else
c
c    uses Pfinal=P_init+radPress(dr)
c               =Po+radPint+radPress(dr)
c
c inner edge, rad pres = 0.0, iterate on nh only for
c pres0, used subsequently
c
c     iterate to get initial pressure
c
          pres0=fpressu(te0,dhn,pop)
          if (jpres.eq.'A') dhn=dhn*(ponk*rkb/pres0)
          if (jpres.eq.'F') then
            dhn=dhn*fpresr(ponk0,ponk1,rp0,palpha,dis0)*rkb/pres0
          endif
          dh0=dhn
        endif
c
        dhl=dif(dhp,dh0)
        init_it=init_it+1
        if ((dhl.ge.dhlma).or.((init_it.lt.5).and.(jden.eq.'B'))) goto
     &   30
c
        qhdnin=q4/(dh0*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh0)
        uhin=qhdhin/cls
c
        if (graindestructmode.eq.1) call adjustgrains (te0, qhdhin, 1)
c
        tspr=dton0
c
        tempre=te0
        call copypop (pop, poppre)
        call copypop (pop, ppre)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'END OF INNER BOUNDRY FOR M = 1'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     calculate qhdn into zone dh0
c
        q1=0.0d0
        q2=0.0d0
        q3=0.0d0
        q4=0.0d0
c
c     Integrate local mean intensities to derive number of H/He
c     ionizing photons. Returns q1=QH1, q2=QHEI, q3=QHEII,
c     q4=QTOT.
c
        call intjnu (tphot, q1, q2, q3, q4)
c
        qhdnin=q4/(dh0*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh0)
        uhin=qhdhin/cls
c
      else
c
c     calculate qhdn into zone dh1
c
        q1=0.0d0
        q2=0.0d0
        q3=0.0d0
        q4=0.0d0
c
c     Integrate local mean intensities to derive number of H/He
c     ionizing photons. Returns q1=QAH1, q2=QAHEI, q3=QAHEII,
c     q4=QATOT.
c
        call intjnu (tphot, q1, q2, q3, q4)
c
        qhdnin=q4/(dh1*zen)
        unin=qhdnin/cls
        qhdhin=q4/(dh1)
        uhin=qhdhin/cls
c
      endif
c
c check grains, don't allow to destroy if already formed
c
      if (graindestructmode.eq.1) call adjustgrains (te0, qhdhin, 0)
c
      zgas=fzgas()
      oonh=dlog10(zion(zmap(8))/zsol(zmap(8)))
      ogas=12.d0+dlog10(zion(zmap(8)))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*)'***ESTIMATE OF VALUES AT THE OUTER BOUNDARY'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (poppre, pop)
      call copypop (poppre, popend)
c
      if (m.lt.4) then
c
c     linear extrap from previous 2 steps
c
        te1=te0+(((te0-tep)*te0)/tep)
        te1=dmax1(te0*agmu,dmin1(te0/agmu,te1))
        te1a=te1
      else
        te1a=te1
        te1=thist(3)
      endif
c
      if (jthm.eq.'T') te1=fixtemp
c
      if (te1.lt.0.d0) te1=100.d0
      if (jden.eq.'C') then
        dh1=dhn
      elseif (jden.eq.'F') then
        dh1=frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
c
c     error if pop is very different from equilibrium at te1
c         frp    = fradpress
c         dh1 = (dhn*prescc)/(fpressu(te1,dhn,pop)+frp)
c      presscc=press(prev zone) + new RadPress
c      Press1=Po with new Tf and local pop (Ptilda)
c      Pfinal/Ptilda=nfinal(dh1)/n_init(dhn)
c
        frp=fradpress(dr,dh)
        prescc=pres0+radpint+frp
        if (jpres.eq.'F') then
          prescc=(fpresr(ponk0,ponk1,rp0,palpha,dis0)*rkb)+radpint+frp
        endif
        pres1=fpressu(te1,dhn,pop)
        dh1=dhn*prescc/(pres1)
      endif
      de1=feldens(dh1,pop)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'***DETERMINE SPACE STEP : DR'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      te1_initial=te1
      dh1_initial=dh1
      de1_initial=de1
c
      call copypop (poppre, pop_initial0)
      call copypop (popend, pop_initial1)
      call copypop (pop, pop_initial)
c
   50 call copypop (pop_initial0, poppre)
      call copypop (pop_initial1, popend)
      call copypop (pop_initial, pop)
c
      t=(te0+te1)*0.5d0
      dh=(dh0+dh1)*0.5d0
      de=(de0+de1)*0.5d0
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
      call equion (t, de, dh)
c
      lmod='SO'
      call totphot2 (t, dh, rad0, 0.d0, dv, wdil0, lmod)
      call absdis2 (dh, dtau, dr, rad0, pop)
c
c      dr = dmax1(dr,(1.0d-6*(rmax-remp)))
c
      if (jend.eq.'E') then
        if ((dis0+dr).gt.diend) dr=(diend-dis0)
      endif
      if (jgeo.eq.'S') then
        if (jden.eq.'F') then
          if ((frad(dis0,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &     ffac)).gt.ofac) then
            dpw=dr
            if (pfac.gt.0.0d0) then
              if (bfac.gt.0.0d0) then
                dpw=dis0*((1.1d0**(1/bfac))-1.d0)
              elseif (bfac.lt.0.0d0) then
                dpw=dis0*((0.9d0**(1/bfac))-1.d0)
              endif
            endif
            dxp=dr
            if (xfac.gt.0.d0) then
              dxp=dabs(0.05d0*scalen)
            endif
            dr=1.d0/((1.d0/dpw)+(1.d0/dxp)+(1.d0/dr))
          endif
        endif
        rad1=rad0+dr
        dis1=rad1
        wdil1=fdilu(rstar,rad1)
        if (radout.lt.khuge) then
          dvol=ftpi*((fring(radin,radout,rad1)*((rad1*irdvol)**3))-
     &     (fring(radin,radout,rad0)*((rad0*irdvol)**3)))
        else
          dvol=ftpi*(((rad1*irdvol)**3)-((rad0*irdvol)**3))
        endif
      endif
      if (jgeo.eq.'P') then
        dis1=dis0+dr
        rad1=0.d0
        wdil1=fdilu(0.0d0,0.0d0)
        dvol=dr*irdvol
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'CALCULATION OF IONIC POPULATION AT THE OUTER'
c     of current step
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (jden.eq.'C') then
        dh1=dhn
      elseif (jden.eq.'F') then
        dh1=frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
        frp=fradpress(dr,dh)
        prescc=pres0+radpint+frp
        if (jpres.eq.'F') then
          prescc=(fpresr(ponk0,ponk1,rp0,palpha,dis1)*rkb)+radpint+frp
        endif
        pres1=fpressu(te1,dhn,popend)
        dh1=dhn*prescc/(pres1)
      endif
c
      de1=feldens(dh1,popend)
c
      rww=1.d0
      if (jgeo.eq.'S') then
        rww=(1.d0+((0.75d0*dr)/(rad0+dr)))**2
      else
        rww=1.d0
      endif
c
      t=(te0+(rww*te1))/(1.d0+rww)
      dh=(dh0+(rww*dh1))/(1.d0+rww)
      de=(de0+(rww*de1))/(1.d0+rww)
c
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
c  Determine dust continuum emission for this region
c
      lmod='DW'
      wei=1.d0/(1.d0+rww)
      call averinto (wei, poppre, popend, pop)
      call localem (t, de, dh)
      call totphot2 (t, dh, rad0, dr, dv, wdil1, lmod)
      call zetaeff (dh)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'*EQUILIBRIUM IONISATION AND TEMPERATURE .'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ((jeq.eq.'E').or.((jeq.eq.'P').and.(jthm.eq.'S'))) then
        dton1=1.d33
        call copypop (poppre, pop)
        call teequi2 (te1, tef, de1, dh1, dton1, nmod)
      else
c
c     Finite age, non-equilibrium ionisation. Determine time step.
c     Takes into account ionization front velocity.
c
        distim(1,m)=dis1
        distim(2,m)=distim(2,m-1)+(dr/(viofr+(1.d-35*dr)))
        tprop=distim(2,m)
        dtauon=1.0d0
c
        call absdis2 (dh, dtauon, disin, 0.0d0, pop0)
        fthick=dmin1(disin,dis1-distim(1,0))
        dison=dis1-fthick
        protim=distim(2,0)
        do 60 j=m-1,0,-1
          if (distim(1,j).gt.dison) goto 60
          prf=(dison-distim(1,j))/((distim(1,j+1)-distim(1,j))+epsilon)
          protim=distim(2,j)+(prf*(distim(2,j+1)-distim(2,j)))
          goto 70
   60   continue
   70   continue
c
c
        treach1=dmax1(dis1/cls,protim+(fthick/cls))
        durmax=dmax1(0.d0,tlife+((dis1/cls)-treach1))
        dton1=dmax1(0.d0,dmin1(durmax,(telap-(dis1/cls))-treach1))
        deeq=dmax1(de1/3.d0,exla*dh1)
        eqscal=5.d0*frectim(2.d0*te1,deeq,dh1)
        if (dton1.lt.eqscal) then
c
c      *TIME STEP TOO SHORT , EQUILIBRIUM TEMPERATURE NOT REACHED
c
          jjd=jden
          jden='C'
          exl=-1.d0
          tex=0.95d0*tm00
          call copypop (pop0, pop)
          call evoltem (tm00, tef, de1, dh1, prescc, dton1, exl, tex,
     &     lut0)
          jden=jjd
        else
c
c      *TIME STEP LONG ENOUGH TO ASSUME EQUILIBRIUM TEMPERATURE
c
          dton1=1.d33
          call copypop (pop0, pop)
c     Compute equilibrium temperature and ionization state of the gas
          call teequi2 (te1, tef, de1, dh1, dton1, nmod)
        endif
        dtco1=(gammaEOSU*fpressu(tef,dh1,pop))/(eloss+epsilon)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       write(*,*)'END EQUILIBRIUM IONISATION AND TEMPERATURE .'
c       write(*,*) t,te1,tef
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jthm.eq.'T') tef=fixtemp
c
      if (jden.eq.'C') then
        dhf=dhn
      elseif (jden.eq.'F') then
        dhf=frad(dis1,dhn,xfac,pfac,afac,bfac,cfac,scalen,vfac,efac,
     &   ffac)
      else
        frp=fradpress(dr,dh)
        prescc=pres0+radpint+frp
        if (jpres.eq.'F') then
          prescc=(fpresr(ponk0,ponk1,rp0,palpha,dis1)*rkb)+radpint+frp
        endif
        pres1=fpressu(tef,dhn,pop)
        dhf=dhn*prescc/(pres1)
      endif
      def=feldens(dhf,pop)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   write(*,*)'***COMPARES END VALUES WITH PREVIOUS STEP'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call difhhe (pop, popend, difend)
      call difhhe (pop, poppre, difpre)
c      call difpop(pop, popend, 1e-8, 5, difend)
c      call difpop(pop, poppre, 1e-8, 5, difpre)
      call copypop (pop, popend)
      difend=difend/difma
      difpre=difpre/difma
      difde=dif(def,de1)/dhlma
      difte=dif(tef,te1)/dtlma
c
      difg=dmax1(difpre,difend,difde,difte)
c
      niter=niter+1
      if (difg.le.difmi) then
        texm=te1
        dtaux=dtau
        difmi=difg
      endif
c
      if (niter.gt.1) then
        write (*,80) tef,difte,difde,difpre,difg,difend,dtau,dr
   80  format(t7,0pf8.0,5(f8.3),3(1pg9.2))
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  write(*,*)'AVERAGE QUANTITIES FOR THE SPACE STEP CONSIDERED'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      rww=1.d0
      if (jgeo.eq.'S') then
        rww=(1.d0+((0.75d0*dr)/(rad0+dr)))**2
      else
        rww=1.d0
      endif
      t=(te0+(rww*tef))/(1.d0+rww)
      dh=(dh0+(rww*dhf))/(1.d0+rww)
      de=(de0+(rww*def))/(1.d0+rww)
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      disav=(dis0+(rww*dis1))/(1.d0+rww)
      wei=1.d0/(1.d0+rww)
      call averinto (wei, poppre, popend, pop)
c
      t_pk=fpressu(t,dh,pop)
c
      if (jgeo.eq.'S') then
        rad=disav
        wdil=fdilu(rstar,rad)
      endif
      if (jgeo.eq.'P') then
        wdil=fdilu(0.d0,0.d0)
        rad=0.d0
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  write(*,*)'CHECK CONVERGENCE OF END VALUES WITH PREVIOUS STEP'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      if ((difg.ge.1.d0).and.(difend.gt.0.001d0)) then
        te1=tef
        if (niter.gt.6) then
c     No convergence after 6 iterations
c            write(*, 474)  tef,difte,difde,difpre,difg,difend,dtau,dr
c 474        format(' FINAL',0pf8.0,5(f8.3),3(1pg9.2))
          goto 90
        endif
c
        dtau=dtau*0.3333d0
c
        te1=te1_initial
        dh1=dh1_initial
        de1=de1_initial
c
        goto 50
c     Jump here if no spatial convergence.
   90   continue
      endif
c
      te1=tef
      dlos1=dlos
      telast=t
      frx=1.d0-pop(jpoen,ielen)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'END CHECK CONVERGENCE STEP'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'INTEGRATES COLUMN DENSITY AND END DIFFUSE FIELD'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jtrans.eq.'OUTW') then
        frdw=1.d0
      endif
      if (jtrans.eq.'LOUP') then
        frdw=0.5d0
      endif
      if (jtrans.eq.'LODW') then
        frdw=0.5d0
      endif
      if (jtrans.eq.'DWUP') then
        frdw=0.5d0
      endif
c
      lmod='DW'
      call localem (t, de, dh)
      call totphot2 (t, dh, rad0, dr, dv, wdil1, lmod)
      call zetaeff (dh)
      if ((grainmode.eq.1).and.(irmode.ne.0)) then
        call dusttemp (t, dh, de, dr, m)
      endif
      call newdif2 (t, t, dh, rad0, dr, dv, 0.d0, dv, frdw, jtrans)
c
      imod='COLD'
      call sumdata (t, de, dh, dvol, dr, disav, imod)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'CALCULATION OF AVERAGE QUANTITIES '
c       EQUAL VOLUMES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jgeo.eq.'S') then
        rww=(1.d0+((0.75d0*dr)/(rad0+dr)))**2
      else
        rww=1.d0
      endif
c
      t=(te0+(rww*tef))/(1.d0+rww)
      dh=(dh0+(rww*dhf))/(1.d0+rww)
      de=(de0+(rww*def))/(1.d0+rww)
      en=dh*zen
      if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
      disav=(dis0+(rww*dis1))/(1.d0+rww)
      dlosav=(dlos0+(rww*dlos1))/(1.d0+rww)
      dton=(dton0+(rww*dton1))/(1.d0+rww)
      treach=(treach0+(rww*treach1))/(1.d0+rww)
      dtco=(dtco0+(rww*dtco1))/(1.d0+rww)
      wei=1.d0/(1.d0+rww)
c
      call averinto (wei, poppre, popend, pop)
c
      t_pk=fpressu(t,dh,pop)/rkb
      frp=fradpress(dr,dh)
      radpint=radpint+frp
c
      if (jgeo.eq.'S') then
        rad=disav
        wdil=fdilu(rstar,rad)
      endif
      if (jgeo.eq.'P') then
        wdil=fdilu(0.d0,0.d0)
        rad=0.d0
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    ***CALCULATION OF IONIC POPULATION AND TEMPERATURE WHEN
c       OR IF SOURCE SWITCHED OFF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dtoff=dmax1(0.d0,(telap-((2.d0*disav)/cls))-tlife)
      if ((jeq.ne.'E').and.(dtoff.gt.0.0)) then
c
        if (jeq.eq.'P') recscal=frectim(t,0.8d0*de,dh)
c
        wd=frlum*wdil
        lmod='SO'
        jspot='Y'
        tii0=t
c
        call totphot2 (t, dh, rad, dr, dv, wd, lmod)
        call zetaeff (dh)
        call evoltem (tii0, t, de, dh, prescc, dtoff, exla, tm00, lut0)
c
        if ((fin.ne.1.d0).and.(jden.eq.'B')) fi=dmin1(1.d0,fidhcc/dh)
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*)'COMPUTES SPECTRUM OF THE REGION CONSIDERED'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      imod='REST'
c     Calculate total cooling rate of the plasma
      call cool (t, de, dh)
c     Integrate the line and temperature data and column density
      call sumdata (t, de, dh, dvol, dr, disav, imod)
c
c     record line ratios
c
      if (ox3.ne.0) then
        hoiii(m)=fmbri(8,ox3)
        f5007=hoiii(m)
        f4363=fmbri(10,ox3)
      endif
      if (ox2.ne.0) then
        hoii(m)=(fmbri(1,ox2)+fmbri(2,ox2))
        f3727=hoii(m)
      endif
      if (ox1.ne.0) then
        hoi(m)=fmbri(3,ox1)
        f6300=hoi(m)
      endif
      if (ni2.ne.0) then
        hnii(m)=fmbri(10,ni2)
        f6584=hnii(m)
      endif
      if (su2.ne.0) then
        hsii(m)=(fmbri(1,su2)+fmbri(2,su2))
        f6720=hsii(m)
        f6731=fmbri(1,su2)
        f6716=fmbri(2,su2)
        f673116=fmbri(2,su2)/fmbri(1,su2)
      endif
c
c    Test for PAH existence (if grainmode true)
c
      if ((grainmode.eq.1).and.(pahmode.eq.1)) then
c
c  Calculate habing photodissociation parameter
c
        ffuv=0.d0
        do i=fuvmin,fuvmax-1
          widnu=widbinnu(i)
          eng=cphote(i)
c    photons s-1 cm-2
          ffuv=ffuv+fpi*tphot(i)*widnu/eng
        enddo
        habing=ffuv/(cls*dh)
c
c Calcualte PAH ionization/Local ISRF ratio (Q factor)
c
        pahsum=0.d0
        ffuv=0.d0
        do i=fuvmin,pahabsmax-1
          widnu=widbinnu(i)
          eng=cphote(i)
          pahsum=pahsum+pahnabs(i)*widnu
          ffuv=ffuv+fpi*tphot(i)/eng*pahnabs(i)*widnu
        enddo
        pahqion=ffuv/pahsum
        pahqfac=pahqion/pahq0
c
c      write(*,210) Hab,pahQfac
c 210  format('        HP =',1pg12.5,', PAH ion/ISRF =',1pg12.5,/)
c
c Determine if PAHs exist
c If so, deplete gas of Carbon in PAHs
c
        if (pahactive.eq.0) then
c
          if (pahend.eq.'H') then
            if (habing.lt.pahlimit) then
              pahactive=1
            endif
          endif
          if (pahend.eq.'Q') then
            if (qhdh.lt.pahlimit) then
              pahactive=1
            endif
          endif
          if (pahend.eq.'I') then
            if (pahqfac.lt.pahlimit) then
              pahactive=1
            endif
          endif
c
        endif
c
        if (pahactive.eq.1) then
          atom=zmap(6)
          if (clinpah.eq.1) then
            zion(atom)=zion0(atom)*deltazion(atom)*dion(atom)
          else
            zion(atom)=zion0(atom)*deltazion(atom)*(dion(atom)+(1.d0-
     &       dion(atom))*pahcfrac)
          endif
        endif
c
      endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***OUTPUT QUANTITIES RELATED TO SPACE STEP
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Special output controlled by poll files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      fhi=pop(1,1)
      fhii=pop(2,1)
c
      caller='P7'
      pollfile='balance'
      inquire (file=pollfile,exist=iexi)
      if (((jbal.eq.'Y').and.(m.eq.1)).or.(iexi)) then
        pfx=trim(jbfx)
        np=len(trim(pfx))
        call wbal (caller, pfx, np, pop)
      endif
c
      pollfile='photons'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        pfx='psou'
        np=4
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.d0, tphot)
      endif
c
      pollfile='nebphot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        savemod=lmod
        lmod='NEBL'
        call totphot2 (t, dh, rad, dr, dv, wdil1, lmod)
        pfx='nsou'
        np=4
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.0d0, tphot)
        lmod=savemod
        call totphot2 (t, dh, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='sophot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        savemod=lmod
        lmod='SO'
        call totphot2 (t, dh, rad, dr, dv, wdil1, lmod)
        pfx='ssou'
        np=4
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.0d0, tphot)
        lmod=savemod
        call totphot2 (t, dh, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='loclphot'
      inquire (file=pollfile,exist=iexi)
      if ((iexi)) then
        savemod=lmod
        lmod='LOCL'
        call totphot2 (t, dh, rad, dr, dv, wdil1, lmod)
        pfx='lsou'
        np=4
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.0d0, tphot)
        lmod=savemod
        call totphot2 (t, dh, rad, dr, dv, wdil1, lmod)
      endif
c
      pollfile='speclocal'
      inquire (file=pollfile,exist=iexi)
      if (iexi) then
        pfx='local'
        sfx='lam'
        call newfile (pfx, sfx, fn, flen)
        open (99,file=fn,status='NEW')
        linemod='LAMB'
        spmod='REL'
        call speclocal (99, tloss, eloss, egain, dlos, t, dh, de,
     &   pop(1,1), disav, dr, linemod, spmod)
        close (99)
      endif
c
      if (grainmode.eq.1) then
        pollfile='IRphot'
        inquire (file=pollfile,exist=iexi)
        if ((iexi).and.(irmode.gt.0)) then
          pfx='irsou'
          wmod='REAL'
          do i=1,infph-1
            pahtest(i)=0.5d0*(photev(i+1)+photev(i))*ev*irphot(i)*dr
          enddo
          call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 0.25d0,
     &     pahtest)
        endif
c
        pollfile='pahphot'
        inquire (file=pollfile,exist=iexi)
        if ((iexi).and.(pahactive.eq.1)) then
          pfx='pahs'
          wmod='REAL'
          do i=1,infph
            pahtest(i)=paheng*pahflux(i)
          enddo
          call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 0.25d0,
     &     pahtest)
        endif
c
c   Infrared output (every 1 step)
c
        pollfile='IRlocal'
        inquire (file=pollfile,exist=iexi)
        if ((irmode.ne.0).and.iexi) then
          if (irfile.eq.' ') then
            pfx='IRflux'
            sfx='sou'
            call newfile (pfx, sfx, fn, flen)
            irfile=fn(1:flen)
            open (ir1,file=irfile,status='NEW')
            write (ir1,100) theversion
  100      format('%  Infrared flux per region ',/,
     & '%  MAPPINGS V ',a12,/,
     & '%  given as:',/,
     & '%    energy edge (eV), continuum flux, ',
     & 'IRflux Fnu(erg s-1 cm-2 Hz-1Sr-1)',/,
     & '%')
            close (ir1)
          endif
          do i=1,infph-1
            cnphot(i)=0.d0
          enddo
          call freefree (t, de, dh)
          call freebound (t, de, dh)
          call twophoton (t, de, dh)
          open (ir1,file=irfile,status='OLD',access='APPEND')
          write (ir1,110) m
          write (ir1,*) infph
  110     format(/,' Region ',i4.4)
          do i=1,infph-1
            cnphot(i)=ffph(i)+fbph(i)+p2ph(i)
            write (ir1,120) photev(i),cnphot(i)*dr,irphot(i)*dr
          enddo
  120     format(1pg14.7,' ',1pg14.7,' ',1pg14.7)
          close (ir1)
        endif
c
      endif
c
c     graincharge output
c
      if (grainmode.eq.1) then
        pollfile='graincharge'
        inquire (file=pollfile,exist=iexi)
        if (iexi) then
          if (chargefile.eq.' ') then
            pfx='grpot'
            sfx='ph6'
            call newfile (pfx, sfx, fn, flen)
            chargefile=fn(1:flen)
            open (ir1,file=chargefile,status='NEW')
            write (ir1,130) theversion
  130      format('%  grain charge in each region ',/,
     & '%  MAPPINGS V ',a12,/,
     & '%  given as:',/,
     & '%  m,Av. distance,Temp,dh,de,',/,
     & '%  grain charge(allsizes) (graphite),',/,
     & '%  grain charge(allsizes) silicate',/,
     & '%')
            write (ir1,140) (grainrad(i),i=mindust(1),maxdust(1))
  140      format('gra radii',10(1pg11.4))
            write (ir1,*)
            write (ir1,150) (grainrad(i),i=mindust(2),maxdust(2))
  150      format('sil radii',10(1pg11.4))
            close (ir1)
          endif
          open (ir1,file=chargefile,status='OLD',access='APPEND')
          write (ir1,160) m,disav,t,dh,de
          write (ir1,170) (grainpot(i,1),i=mindust(1),maxdust(1))
          write (ir1,180) (grainpot(i,2),i=mindust(2),maxdust(2))
          close (ir1)
  160      format(/,i3,4(1pg14.5))
  170      format('gra ',30(1pg11.4))
  180      format('sil ',30(1pg11.4))
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Main output
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call appendp7files ()
c
      if ((grainmode.eq.1).and.(pahmode.eq.1)) then
c
  190 format(i4,1pg11.4,3(1x,1pg10.3),13(1pg14.7),6(1x,1pg10.3))
c
        write (luop,190) m,t,dlosav,difg,dtau,disav,(dis0-remp),dr,dh,
     &   de,de+en,fhi,fhii,dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),
     &   hoiii(m),graindgr,zgas,ogas,habing,pahqfac
c
      elseif (graindestructmode.eq.1) then
c
  200 format(i4,1pg11.4,3(1x,1pg10.3),13(1pg14.7),4(1x,1pg10.3))
c
        write (luop,200) m,t,dlosav,difg,dtau,disav,(dis0-remp),dr,dh,
     &   de,de+en,fhi,fhii,dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),
     &   hoiii(m),graindgr,zgas,ogas
c
      else
c
  210 format(i4,1pg11.4,3(1x,1pg10.3),13(1pg14.7))
c
        write (luop,210) m,t,dlosav,difg,dtau,disav,(dis0-remp),dr,dh,
     &   de,de+en,fhi,fhii,dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls)
c
      endif
c
  220 format(i4,20(',',1pg13.6))
      write (lunt,220) m,(dis0-remp),dr,disav,dvol*rdvol**3,t,dh,de,de+
     &en,fhi,fhii,hbeta,f5007,f4363,f3727,f6300,f6584,f6720,f6716,f6731,
     &f673116
c
      if (jeq.eq.'F') write (luop,230) fthick,treach,dtco,dton,dtoff
  230 format(t71,5(1pg9.2))
c
  240 format(t98,5(1pg9.2))
      if (jeq.eq.'P') write (luop,240) recscal,dtoff
c
      if (mod(m,20).eq.1) then
c
        if ((graindestructmode.eq.0).and.(grainmode.eq.1)
     &   .and.(pahmode.eq.1)) then
c
          if (jgeo.eq.'S') then
  250  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103,'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,'<Charge>')
            write (*,250)
          else
  260  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103,'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,'<Charge>')
            write (*,260)
          endif
c
        elseif ((graindestructmode.eq.1).and.(pahmode.eq.1)) then
c
          if (jgeo.eq.'S') then
  270  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,' Zgas',t159,'LogO+12',t167,'<Charge>')
            write (*,270)
          else
  280  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,'Hab.P',t143,'QPAH/ISRF'
     & ,       t151,' Zgas',t159,'LogO+12',t167,'<Charge>')
            write (*,280)
          endif
c
        elseif ((graindestructmode.eq.1).and.(pahmode.eq.0)) then
c
          if (jgeo.eq.'S') then
  290  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,' Zgas',t143,' LogO+12s'
     & ,       t150,'<Charge>')
            write (*,290)
          else
  300  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k'
     & ,       t103, 'Log<QN>', t111,'Log<UH>', t119,'OIII/HB'
     & ,       t127,'Dust/Gas', t135,' Zgas',t143,' LogO+12s'
     & ,       t150,'<Charge>')
            write (*,300)
          endif
c
        else
c
          if (jgeo.eq.'S') then
  310  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<R>',t46,'dr'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k', t103, 'Log<QN>', t111,'Log<UH>'
     & ,       t119,'OIII/HB')
            write (*,310)
          else
  320  format(' Step',t7,'<Te>',t15,'<DLOS>',t24,'dChi',t30,'dTau'
     & ,       t36,'<X>',t46,'dx'
     & ,       t56,'<XHI>',t66,'<nH>',t76,'<ne>',t87,'<nt>'
     & ,       t95,'LogP/k', t103, 'Log<QN>', t111,'Log<UH>'
     & ,       t119,'OIII/HB')
            write (*,320)
          endif
        endif
      endif
c
  330 format(i4,x,0pf8.0,(1pg9.2),0pf6.3,0pf6.3
     & ,6(1pg10.3),5(0pf8.3),2(0pf7.3,x),6(0pf8.3))
c
      gcharge=0.5d0*(avgpot(1)+avgpot(2))
c
      if ((graindestructmode.eq.0).and.(grainmode.eq.1)
     &.and.(pahmode.eq.1)) then
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),f5007/hbeta,
     &   graindgr,habing,pahqfac,gcharge
c
      elseif ((graindestructmode.eq.1).and.(pahmode.eq.1)) then
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),
     &   dlog10(qhdn/cls),f5007/hbeta,graindgr,habing,pahqfac,zgas,ogas,
     &   gcharge
c
      elseif ((graindestructmode.eq.1).and.(pahmode.eq.0)) then
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),
     &   dlog10(qhdn/cls),f5007/hbeta,graindgr,zgas,ogas,gcharge
c
      else
c
        write (*,330) m,t,dlosav,difg,dtau,disav,dr,fhi,dh,de,de+en,
     &   dlog10(t_pk),dlog10(qhdn),dlog10(qhdh/cls),f5007/hbeta
c
      endif
c
      if (jiel.eq.'Y') then
  340  format(1x,i4,',', 10(1pg12.5,', '),31(1pg12.5,', '))
        do idx=1,ieln
          write (luions(idx),340) m,disav,(dis0-remp),dr,t,de,dh,de+en,
     &     dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),(pop(j,iel(idx)),
     &     j=1,maxion(iel(idx)))
        enddo
      endif
c
      if (jcol.eq.'Y') then
        do idx=1,ieln
          write (lucols(idx),340) m,disav,(dis0-remp),dr,t,de,dh,de+en,
     &     dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),(popint(j,iel(idx)
     &     ),j=1,maxion(iel(idx)))
          write (lucols(idx),*) ''
        enddo
      endif
c
      if (jiem.eq.'Y') then
  350  format(1x,i4,',',11(1pg12.5,', '),31(1pg12.5,', '))
        do idx=1,iemn
          i=iem(idx)
          nt=nfmtrans(i)
          write (luemiss(idx),350) m,disav,(dis0-remp),dr,dvol*rdvol**3,
     &     t,de,dh,de+en,dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),
     &     hbeta,(fmbri(itr,iem(idx)),itr=1,nt)
        enddo
      endif
c
      if (jlin.eq.'Y') then
  360  format(1x,i4,',', 11(1pg12.5,', '),31(1pg12.5,', '))
        call speclocallines (linfluxes)
        write (lulin,360) m,disav,(dis0-remp),dr,dvol*rdvol**3,t,de,dh,
     &   de+en,dlog10(qhdh),dlog10(qhdh/cls),dlog10(qhdn),hbeta,
     &   (linfluxes(itr),itr=1,njlines)
        write (lulin,*) ''
      endif
c
      if (jall.eq.'Y') then
        write (lusl,370)
        write (lusl,380) m,t,de,dh,fi,dvol*rdvol,disav,(dis0-remp),dr
  370    format(//' Step   Te Ave.(K)   ',
     & '  ne(cm^-3)   ',
     & '  nH(cm^-3)   ',
     & '  fill. Fact. ',
     & '  dVol.  ',
     & ' Dist.Ave.(cm)',
     & ' Din-Remp.(cm)',
     & '    dR (cm)   ')
  380    format(1x,i4,8(1pg14.7)/)
      endif
c
      wmod='PROP'
      call wmodel (lupf, t, de, dh, disav, wmod)
      wmod='LOSS'
      call wmodel (lups, t, de, dh, disav, wmod)
c
      if (jall.eq.'Y') then
        call wionabal (lusl, pop)
        call wionabal (lusl, popint)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***END OUTPUT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call closep7files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***TEST ENDING CONDITIONS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (((jend.eq.'A').or.(jend.eq.'B')).and.(frx.le.fren)) goto 410
      if ((jend.eq.'C').and.(telast.le.tend)) goto 410
      if ((jend.eq.'E').and.(dis1.ge.diend)) goto 410
c
      do 390 n=1,ionum
c      find first (lowest energy) cross section that matches species
        if ((atpho(n).eq.ielen).and.(ionpho(n).eq.jpoen)) goto 400
  390 continue
c sigpho contains new threshold values if verner in force.
  400 taux=sigpho(n)*popint(jpoen,ielen)
      if ((jend.eq.'D').and.(taux.ge.tauen)) goto 410
      if (jend.eq.'F') then
        if (jpoen.gt.maxion(ielen)) then
          popinttot=0.d0
          do i=1,maxion(ielen)
            popinttot=popinttot+popint(i,ielen)
          enddo
          if (popinttot.ge.colend) goto 410
        else
c         write (*,*) jpoen,ielen,popint(jpoen,ielen),colend
          if (popint(jpoen,ielen).ge.colend) goto 410
        endif
      endif
      if (jend.eq.'H') then
        popinttot=0.d0
        popinttot=popinttot+popint(1,1)+popint(2,1)
        if (popinttot.ge.colend) goto 410
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Check Visual extinction condition or end when HII<1%
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if (jend.eq.'H') then
c       a_v=dustav*dustint
c       if (a_v.gt.a_vend) goto 410
c     endif
c
      if (((((de1/dh1)-dedhmi)/dedhma).le.exla).and.(grainmode.eq.0))
     &goto 410
c
      pollfile='terminate'
      inquire (file=pollfile,exist=iexi)
      if (iexi) goto 410
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***RESET INNER BOUNDARY QUANTITIES FOR NEXT SPACE STEP
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (popend, poppre)
c
      tep=te0
      te0=tef
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Record electron temperatures in THIST for quadratic extrapolation
c     when m > 3.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      thist(1)=thist(2)
      thist(2)=thist(3)
      thist(3)=tef
c
      de0=def
      dh0=dhf
      dton0=dton1
      treach0=treach1
      dtco0=dtco1
      fi=fi1
      rad0=rad1
      dis0=dis1
      wdil0=wdil1
      drp=dr
      dlos0=dlos1
c
      if (difg.lt.0.1d0) dtau=0.1d0*dtau/difg
      if (difg.gt.0.75d0) dtau=dtau/(1.d0+(difg-0.75d0)*5.d0)**2.d0
      dtau=dmax1(1.0d-5,dtau)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*)'LOOP BACK AND INCREMENT M, DO NEXT STEP.....'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      m=m+1
c
c     stop things getting silly....
c
      if (m.gt.mma) goto 410
c
c     loop back
c
      goto 20
c
  410 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write(*,*)'MODEL ENDED ; OUTPUT RESULTS **************'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call appendp7files ()
c
      write (luop,420)
      write (lunt,420)
  420 format(//' Model ended',t20,'Tfinal',t28,'DISfin',t38,'FHXF',t48,
     &'TAUXF',t58,'Ending')
      write (luop,430) t,dis1,frx,taux,jend
  430 format(' ###########' ,t18,0pf8.0,3(1pg10.3),4x,a4,//)
      if ((jeq.ne.'E').and.(frlum.gt.0.0)) write (luop,440) frlum
  440 format(/' Final fractional luminosity of source after ',
     &'turn off :',1pg10.3)
      if ((radout.lt.khuge).and.(jgeo.eq.'S')) write (luop,450) radin,
     &radout
  450 format(/' NB. ***************Integration through the line',
     &'of sight for a ring aperture of radii :',2(1pg10.3))
      write (lupf,460) tempre,tspr
  460 format(//t4,'Preionisation conditions for step#0 ',
     &'for all elements   (TEpr :',0pf8.0,' Time :',1pg9.2,' )  :'/)
      write (lupf,470) (elem(i),i=1,atypes)
  470 format(' ',t8,16(4x,a2,4x)/)
      maxio=0
      do j=1,atypes
        if (maxio.le.maxion(j)) maxio=maxion(j)
      enddo
      do j=1,maxio
        chcksum=0.0d0
        do i=1,atypes
          chcksum=chcksum+ppre(j,i)
        enddo
        if (chcksum.gt.0.0d0) then
          write (lupf,480) rom(j),(ppre(j,i),i=1,atypes)
  480    format(' ',a6,t8,16(1pg10.3))
        endif
      enddo
c
      savemod=lmod
      if (jspec.eq.'Y') then
c
c and as nu-nufnu
c
        caller='P7'
        pfx=trim(jpfx)
        np=len(trim(pfx))
c
c     Down stream nebula only photon field and source
c
        wmod='LFLM'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.0d0, tphot)
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.0d0, tphot)
        wmod='NFNU'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.0d0, tphot)
c
      endif
c
      if (jsou.eq.'Y') then
c
c      write nebula spectrum in lam - flam (ergs/s/cm2/A)
c
        caller='P7'
        pfx='flam_'
        np=len(trim(pfx))
c
c      write nebula spectrum into lam - flam (ergs/s/cm2/A)
c
        do j=1,infph
         tp1(j)=tphot(j)
        enddo
c
        savemod=lmod
        lmod='NEBL'
        call totphot2 (te1, dh1, dis0, dr, 0.0d0, wdil0, lmod)
        do j=1,infph
         tp2(j)=tphot(j)
        enddo
        lmod=savemod
c
        savemod=lmod
        lmod='SO'
        call totphot2 (te1, dh1, dis0, dr, 0.0d0, wdil0, lmod)
        do j=1,infph
         tp3(j)=tphot(j)
        enddo
        lmod=savemod
c
        savemod=lmod
        lmod='NEBC'
        call totphot2 (te1, dh1, dis0, dr, 0.0d0, wdil0, lmod)
        do j=1,infph
         tp4(j)=tphot(j)
        enddo
        lmod=savemod
c
        call wplam4 (caller, pfx, np, t, de, dh, dr, 1.d0)
c
        lmod=savemod
        call totphot2 (te1, dh1, dis0, dr, 0.0d0, wdil0, lmod)
c
      endif
c
c    ***COMPUTE AVERAGE SPECTRUM AND WRITES IT IN FILE PHN
c
      call avrdata
c
      call wrsppop (luop)
c
      linemod='LAMB'
      spmod='REL'
c
      call spec2 (lusp, linemod, spmod)
c
      call closep7files ()
c
      write (*,490) banfil,filnam
  490 format(//a12,' Output created &&&&&& File : ',a/)
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine p7headers
c! XXXX - add one line purpose here
c! @param [in,out] character    banfil  XXX-meaning
c! @param [in,out] character       dhn  XXX-meaning
c! @param [in,out]   real*8       fin  XXX-meaning
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine p7headers (banfil, dhn, fin)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine to get the file header buisness out where
c     it can be worked on, and/or modified for other headers.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'p7blocks.inc'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           Variables
c
      real*8 dhn,fin
      real*8 dht,zi(mxelem)
c
      integer*4 i, j, at, io
      integer*4 idx, nt, nl, itr
      integer*4  nr, nf, nb
c
      character abundtitle*64
      character banfil*64
      character el*4
c
      real*8 lamvac,lamair
      real*8 densnum, fnair
      integer*4 mlen
c
      dht=densnum(dhn)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call appendp7files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***WRITES ON FILE INITIAL PARAMETERS
c
      nr=mlen(runname)
      nb=mlen(banfil)
      nf=mlen(filnam)
   10 format(//,
     & ' PHOTO 7: Photoionisation Model',/,
     & ' ============================================',/,
     & ' (Diffuse Field , Radiation Pressure)',//
     & ' Calculated by MAPPINGS V ',a12,/,
     & ' Run   :, ',a/
     & ' Input :, ',a/
     & ' Output:, ',a/)
      write (luop,10) theversion,runname,banfil,filnam
      write (lupf,10) theversion,runname,banfil,filnam
      write (lups,10) theversion,runname,banfil,filnam
      write (lusp,10) theversion,runname,banfil,filnam
      write (lunt,10) theversion,runname,banfil,filnam
      if (jall.eq.'Y') then
        write (lusl,10) theversion,runname(1:nr),banfil(1:nb),
     &   filnam(1:nf)
      endif
      if (jlin.eq.'Y') then
        write (lulin,10) theversion,runname(1:nr),banfil(1:nb),
     &   filnam(1:nf)
      endif
      if (jiel.eq.'Y') then
        do idx=1,ieln
          write (luions(idx),10) theversion,runname(1:nr),banfil(1:nb),
     &     filnam(1:nf)
        enddo
      endif
      if (jcol.eq.'Y') then
        do idx=1,ieln
          write (lucols(idx),10) theversion,runname(1:nr),banfil(1:nb),
     &     filnam(1:nf)
        enddo
      endif
      if (jiem.eq.'Y') then
        do idx=1,iemn
          write (luemiss(idx),10) theversion,runname(1:nr),banfil(1:nb),
     &     filnam(1:nf)
        enddo
      endif
      write (lupf,20)
   20 format(' Plasma properties for ',
     &'each model step,',/,
     &' and the pre-ionisation conditions used.'/ )
   30 format(' Nebula structure summary,',/,
     & ' and the strong line emissivities.',/,
     & ' Weight by shell volumes to get relative luminosities.'/)
      write (lunt,30)
   40 format(' Ionisation fraction for ',
     &'all atomic elements as a function of distance,',/,
     &' and the abundances used. See photn for model details.'/ )
      if (jall.eq.'Y') then
        write (lusl,40)
      endif
   50 format(' Line emissivities (erg/cm^3/s/sr) for'
     & ' ',i3,' lines as a function of distance or radius.',/,
     & ' Weight by shell volumes x4pi to get luminosities.'/)
      if (jlin.eq.'Y') then
        write (lulin,50) njlines
      endif
   60   format(' Line emissivities (erg/cm^3/s/sr) for'
     & ' multi-level atom models as a function of distance.',/,
     & ' Weight by shell volumes x4pi to get luminosities.'/)
      if (jiem.eq.'Y') then
        do idx=1,iemn
          write (luemiss(idx),60)
        enddo
      endif
      if (jiel.eq.'Y') then
        do idx=1,ieln
          el=elem(iel(idx))
          write (luions(idx),'(" Ionic Fraction Structure : X")')
          write (luions(idx),'(/" Element : ",a2/)') el
        enddo
      endif
      if (jcol.eq.'Y') then
        do idx=1,ieln
          el=elem(iel(idx))
          write (lucols(idx),'(" Ionic Column Density Structure : N")')
          write (lucols(idx),'(/" Element : ",a2/)') el
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     abundances file header
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      abundtitle=' Initial Abundances :'
      do i=1,atypes
        zi(i)=zion0(i)*deltazion(i)
      enddo
c
      call dispabundances (luop, zi, abundtitle)
c
      abundtitle=' Gas Phase Abundances :'
      call dispabundances (luop, zion, abundtitle)
c
      call wabund (lupf)
      call wabund (lups)
      call wabund (lusp)
      call wabund (lunt)
c
      if (jall.eq.'Y') then
        call wabund (lusl)
      endif
      if (jlin.eq.'Y') then
        call wabund (lulin)
      endif
      if (jiem.eq.'Y') then
        do idx=1,ieln
          call wabund (luemiss(idx))
        enddo
      endif
      if (jcol.eq.'Y') then
        do idx=1,ieln
          call wabund (lucols(idx))
        enddo
      endif
      if (jiel.eq.'Y') then
        do idx=1,ieln
          call wabund (luions(idx))
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   70     format(
     & '#[1],      <R> [2],   DeltaR [3],       dR [4],      Vol [5]',
     & ',      <T> [6],     <ne> [7],     <nH> [8], <nTotal> [9]',
     & ', <logQH> [10], <logUH> [11], <logQN> [12]',$)
   80     format(
     & '#[1],      <X> [2],   DeltaX [3],       dX [4],      Vol [5]',
     & ',      <T> [6],     <ne> [7],     <nH> [8], <nTotal> [9]',
     & ', <logQH> [10], <logUH> [11], <logQN> [12]',$)
   90     format(
     & '#[1],      <R> [2],   DeltaR [3],       dR [4]',
     & ',      <T> [5],     <ne> [6],     <nH> [7], <nTotal> [8]',
     & ',  <logQH> [9], <logUH> [10], <logQN> [11]',$)
  100     format(
     & '#[1],      <X> [2],   DeltaX [3],       dX [4]',
     & ',      <T> [5],     <ne> [6],     <nH> [7], <nTotal> [8]',
     & ',  <logQH> [9], <logUH> [10], <logQN> [11]',$)
  110     format(', <HB> [13]',$)
c   95     format(',    <HB> [12]',$)
  120     format(',',1pg14.6,'[',i3,']',$)
  130     format(',',1pg14.6,'[',i3,']',$)
  140     format(', X_',a2,a6,'[',i3,']',$)
  150     format(', N_',a2,a6,'[',i3,']',$)
c
  160     format(/,'ion:',12(',             '),$)
  170     format(', ',a2,a6,' ',$)
      if (jlin.eq.'Y') then
        write (lulin,160)
        write (lulin,170) (elem(emlinlistatom(itr)),
     &   rom(emlinlistion(itr)),itr=1,njlines)
        write (lulin,*)
        if (jgeo.eq.'S') then
          write (lulin,70)
        else
          write (lulin,80)
        endif
        write (lulin,110)
        write (lulin,130) (emlinlist(itr),itr+13,itr=1,njlines)
        write (lulin,*)
      endif
      if (jiem.eq.'Y') then
        do idx=1,iemn
          i=iem(idx)
          at=fmatom(i)
          io=fmion(i)
          nt=nfmtrans(i)
          nl=fmnl(i)
          write (luemiss(idx),60)
          write (luemiss(idx),'(/" Emissivity for Element Ion:",a2,a6)')
     &      elem(at),rom(io)
  180    format(' Atomic Level Data for Ion: ',/,
     & '  Levels:',i4)
          write (luemiss(idx),180) ,nl
  190    format ( '  Level,  Term,      g,   E (cm^-1)  ')
  200    format ( '#  n  , i, j,    lambda_vac,    lambda_std ')
  210    format ( '#   ',i3,',', a6,',     ',i2.2,',',1pg14.6)
  220    format ( '# ',i2.2,'  ',2(',',i2.2),2(',', 1pg14.6))
          write (luemiss(idx),190)
          do j=1,nl
            write (luemiss(idx),210) j,nfmterm(j,i),wim(j,i),eim(j,1,i)/
     &       (plk*cls)
          enddo
  230    format('   Lines:',i3)
          write (luemiss(idx),230) ,nt
          write (luemiss(idx),200)
          do j=1,nt
            lamvac=1.0d8*fmlam(j,i)
            lamair=lamvac/fnair(lamvac)
            write (luemiss(idx),220) j,nfmlower(j,i),nfmupper(j,i),
     &       lamvac,lamair
          enddo
          write (luemiss(idx),*)
          if (jgeo.eq.'S') then
            write (luemiss(idx),70)
          else
            write (luemiss(idx),80)
          endif
          write (luemiss(idx),110)
          write (luemiss(idx),120) (1.0d8*fmlam(itr,i)/fnair(1.0d8*
     &     fmlam(itr,i)),itr+13,itr=1,nt)
          write (luemiss(idx),*)
        enddo
      endif
c
      if (jcol.eq.'Y') then
        do idx=1,ieln
          el=elem(iel(idx))
          if (jgeo.eq.'S') then
            write (lucols(idx),90)
          else
            write (lucols(idx),100)
          endif
          write (lucols(idx),150) (el,rom(j),j+11,j=1,maxion(iel(idx)))
          write (lucols(idx),*)
        enddo
      endif
c
      if (jiel.eq.'Y') then
        do idx=1,ieln
          el=elem(iel(idx))
          write (luions(idx),'(" Ion Fractions : ")')
          write (luions(idx),'(" Element : ",a2/)') elem(iel(idx))
          if (jgeo.eq.'S') then
            write (luions(idx),90)
          else
            write (luions(idx),100)
          endif
          write (luions(idx),140) (el,rom(j),j+11,j=1,maxion(iel(idx)))
          write (luions(idx),*)
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jeq.eq.'E') then
        write (luop,240) dtau0
  240 format( ' Model Summary :',//
     &,t5,' Mode  : Thermal and Ionic Equilibrium ',/
     &,t5,' dTau  :',0pf8.5,/)
      elseif (jeq.eq.'F') then
        write (luop,250) telap,tlife,dtau0
  250 format( ' Model Summary :',//
     &,t5,' Mode        :  Non-Equilibrium + Finite Source Life ',/
     &,t5,' Age         :',1pg10.3,' sec',/
     &,t5,' Source Life :',1pg10.3,' sec',/
     &,t5,' dTau        :',0pf8.5,/)
      else
        write (luop,260) telap,dtau0
  260 format( ' Model Summary :',//
     &,t5,' Mode        :  Equilibrium + Source Switch off ',/
     &,t5,' Switch off  :',1pg10.3,' sec',/
     &,t5,' dTau        :',0pf8.5,/)
      endif
      if (usekappa) then
  270 format( ' Kappa Electron Distribution Enabled :',//
     & ,t5,' Electron Kappa : ',1pg11.4,/)
        write (luop,270) kappa
      endif
c
  280 format(/
     & ' ********************************************************',/,
     & ' ** Radiation Field and Parameters **********************',/,
     & ' ********************************************************',/,
     & t5, ' At estimated T_inner :',1pg10.3,' K',/,
     & t5,' Rsou.',t20,' Remp.',t35,' Rmax',t48,' <DILU>',/,
     & t5, 1pg10.3,t20,1pg10.3,t35,1pg10.3,t48,1pg10.3/
     & t5,' <Hdens>',t20,' <Ndens>',t35,' Fill Factor',/
     & t5, 1pg10.3,t20,1pg10.3,t35,0pf9.6,//
     & ' ********************************************************',/,
     & t5,' FQ tot',t20,' Log FQ',t35,' LQ tot',t48,' Log LQ',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t48,0pf8.3,/
     & t5,' Q(N) in',t20,' Log Q(N)',t35,' <Q(N)>',t48,' Log<Q(N)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t48,0pf8.3,/
     & t5,' Q(H) in',t20,' Log Q(H)',t35,' <Q(H)>',t48,' Log<Q(H)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t48,0pf8.3,/
     & t5,' U(N) in',t20,' Log U(N)',t35,' <U(N)>',t48,' Log<U(N)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t48,0pf8.3,/
     & t5,' U(H) in',t20,' Log U(H)',t35,' <U(H)>',t48,' Log<U(H)>',/
     & t5, 1pg10.3,t20,0pf8.3,t35,1pg10.3,t48,0pf8.3,/
     & ' ********************************************************'/)
c
      write (luop,280) tinner,rstar,remp,rmax,wdpl,dhn,dht,fin,qht,
     &dlog10(qht),astar*qht,dlog10(astar*qht),qhdnin,dlog10(qhdnin),
     &qhdnav,dlog10(qhdnav),qhdhin,dlog10(qhdhin),qhdhav,dlog10(qhdhav),
     &unin,dlog10(unin),unav,dlog10(unav),uhin,dlog10(uhin),uhav,
     &dlog10(uhav)
c
  290 format(t6,'MOD',t10,'Temp.',t22,'Alpha',t31,'Turn-on',
     & t39,'Cut-off',/,
     & t6,a2,x,1pg10.3,x,0pf8.3,x,0pf8.3,x,0pf8.3,//,
     & t6,'Zstar',t16,'FQHI',t27,'FQHEI',t38,'FQHEII',/,
     & t6,0pf8.4,1x,3(1pg10.3,x),/)
      write (luop,290) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      if (grainmode.eq.1) then
        write (luop,330) galpha,amin(1),amax(1),amin(2),amax(2),
     &   graindens(1),graindens(2),bgrain,yinf
c
  310  format(  t3,' Projected dust area/H atom (graphite,silicate):',/
     &         ,t5,'  ',1pg12.4,1x,' (cm^2)',1x,1pg12.4,1x,' (cm^2)' )
        write (luop,310) siggrain(1),siggrain(2)
c
  320  format( t3,' Composition Mass Ratios:  ',/
     &        ,t5,'    Gas/H     :',1pg12.4/
     &        ,t5,'    Dust/Gas  :',1pg12.4/
     &        ,t5,'    PAH/Gas   :',1pg12.4/)
c
        write (luop,320) grainghr,graindgr,grainpgr
c
      endif
c
  330 format( ' Dust Parameters :',//
     & ,t5,' Alpha :',1pg10.3,/
     & ,t5,' Graphite min:       max:',/
     & ,t5,3x,1pg10.3,x,1pg10.3,/
     & ,t5,' Silicate min:       max:',/
     & ,t5,3x,1pg10.3,x,1pg10.3,/
     & ,t5,' Graph dens  Sil dens',/
     & ,t5,3x,1pg10.3,x,1pg10.3,/
     & ,t5,' Bgrain      Yinf ',/
     & ,t5,3x,1pg10.3,x,1pg10.3)
c
  340 format(   /,' Model Parameters:',//
     & ,t6,'Jden',t12,'Jgeo',t18,'Jend',t23,'Ielen',t29,
     & 'Jpoen',t38,'Fren',/
     & ,t6,3(a4,2x),2(i2,4x),0pf6.4,//
     & ,t6,'Tend',t12,'DIend',t22,'TAUen',t31,'Jeq',t38,'Teini',/
     & ,t5,0pf6.0,2(1pg10.3),x,a4,x,0pf7.1,//,
     & '************************************',
     & '************************************',/)
      write (luop,340) jden,jgeo,jend,ielen,jpoen,fren,tend,diend,tauen,
     &jeq,tm00
c
  350 format(/' Description of each space step:')
      write (luop,350)
c
      if ((grainmode.eq.1).and.(pahmode.eq.1)) then
        if (jgeo.eq.'S') then
  360 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     &  t40,'dTau',t50,'<R>',t64,'Delta R',t78,'dr',t93,'nH',
     & t107,'ne',t120,'nt',t134,'XHI',t148,'XHII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)',t232,'Hab.P.',t246,'QPAH/ISRF')
          write (luop,360)
        else
  370 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     &  t40,'dTau',t50,'<X>',t64,'Delta X',t78,'dx',t93,'nH',
     & t107,'ne',t120,'nt',t134,'XHI',t148,'XHII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)',t232,'Hab.P.',t246,'QPAH/ISRF')
          write (luop,370)
        endif
      else
        if (jgeo.eq.'S') then
  380 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     &  t40,'dTau',t50,'<R>',t64,'Delta R',t78,'dr',t93,'nH',
     & t107,'ne',t120,'nt',t134,'XHI',t148,'XHII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)')
          write (luop,380)
        else
  390 format(/,t3,'#',t8,'<Te>',t18,'<DLOS>',t29,'dChi',
     &  t40,'dTau',t50,'<X>',t64,'Delta X',t78,'dx',t93,'nH',
     & t107,'ne',t120,'nt',t134,'XHI',t148,'XHII',t162,'Log<P/k>'
     & t176,'Log<Q(H)>',t190,'Log<U(H)>',t204,'Log<Q(N)>',
     & t218,'Log<U(N)')
          write (luop,390)
        endif
      endif
c
      if (jeq.eq.'F') write (luop,400)
  400 format(t73,'Fthick',t82,'Treach',t91,'DTCO',t100,'Dton',t109,
     &'DToff')
      if (jeq.eq.'P') write (luop,410)
  410 format(t100,'RECscal',t109,'DToff')
      write (*,420) filnam
  420 format(/' Output in file : ',a/)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phsem file header
c
c     contains structure summary and emissivity of strong lines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  430  format(
     & '#[1],      <R> [2],   DeltaR [3]',
     & ',       dR [4],      Vol [5]',$)
  440  format(
     & '#[1],      <X> [2],   DeltaX [3]',
     & ',       dX [4],      Vol [5]',$)
  450  format(
     & ',      <T> [6],     <ne> [7]',
     & ',     <nH> [8], <nTotal> [9]',
     & ',    X_HI [10],   X_HII [11]',
     & ', e_HBeta [12],  e_5007 [13]',
     & ',  e_4363 [13],  e_3727+[14]',
     & ',  e_6300 [15],  e_6584 [16]',
     & ',  e_6725+[17],  e_6731 [18]',
     & ',  e_6716 [19], 6731/16 [20]')
      if (jgeo.eq.'S') then
        write (lunt,430)
      else
        write (lunt,440)
      endif
      write (lunt,450)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phapn file header
c
c     new format: contains plasma properties in line, for
c     easy plotting later.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  460 format(' Dist.',t14,'  Te',t28,'  de',t42,'  dh',t56,
     &'  en',t70,'  XHI',t84,'  Density',t98,'  Pressure',t112,
     &'  Flow',t126,'  Ram Press.',t140,'  Sound Spd',t154,
     &'  Slab depth',t168,'  Alfen Spd.',t182,
     &'  Mag. Field',t196,'  Grain Pot.')
      write (lupf,460)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     phlss file header
c
c     contains plasma cooling/heating in line, for
c     easy plotting later.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call wmodel (lups, 0.d0, 0.0d0, 0.d0, 0.0d0, 'LOSH')
c
      call closep7files ()
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine p7filenames
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine p7filenames ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     set initial file unit numbers and default mode flags
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'p7blocks.inc'
c
      integer*4 i, flen
      integer*4 nt, nl, at, io, idx, ie
c
      character fn*128
      character pfx*64,sfx*16
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c units
c
c     photn filna default
      luop=20
c     phapn filnb default
      lupf=21
c     phlss filnc default
      lups=22
c     spec filnd default
      lusp=23
c     phsem filnt default
      lunt=24
c     bands filpb jbnd
      lupb=25
c     lines single file filin jlin
      lulin=26
c     rates file filrt jrat
      lurt=27
c     allions  filio jall
      lusl=28
c     luions, filn jiel
      do i=1,atypes
        luions(i)=32+i
      enddo
c     luemiss filnem jiem
      do i=1,atypes
        luemiss(i)=32+atypes+i
      enddo
c     jcol, lucols, filncol
      do i=1,mxmonmlions
        lucols(i)=32+atypes+atypes+i
      enddo
      ir1=lucols(mxmonmlions)+1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 20
      fn=' '
      pfx='photn'
      sfx='ph7'
      call newfile (pfx, sfx, fn, flen)
      filna=fn(1:flen)
      filnam=filna
c 21
      fn=' '
      pfx='phapn'
      sfx='ph7'
      call newfile (pfx, sfx, fn, flen)
      filnb=fn(1:flen)
c 22
      fn=' '
      pfx='phlss'
      sfx='ph7'
      call newfile (pfx, sfx, fn, flen)
      filnc=fn(1:flen)
c 23
      fn=' '
      pfx='spec'
      sfx='csv'
      call newfile (pfx, sfx, fn, flen)
      filnd=fn(1:flen)
c 24
      fn=' '
      pfx='phsem'
      sfx='ph7'
      call newfile (pfx, sfx, fn, flen)
      filnt=fn(1:flen)
c 25
      fn=' '
      pfx='bands'
      sfx='csv'
      call newfile (pfx, sfx, fn, flen)
      filpb=fn(1:flen)
c 26
      fn=' '
      pfx='lines'
      sfx='csv'
      call newfile (pfx, sfx, fn, flen)
      filin=fn(1:flen)
c 27
      fn=' '
      pfx='rates'
      sfx='csv'
      call newfile (pfx, sfx, fn, flen)
      filrt=fn(1:flen)
c 28
      fn=' '
      pfx='ions_'
      sfx='ph7'
      call newfile (pfx, sfx, fn, flen)
      filio=fn(1:flen)
c
      sfx='csv'
      do i=1,ieln
        ie=iel(i)
        fn=' '
        pfx=elem(ie)
        if (elem_len(ie).eq.1) then
          pfx=elem(ie)
          pfx=pfx(1:1)//'_ion'
        else
          pfx=elem(ie)
          pfx=pfx(1:2)//'_ion'
        endif
        call newfile (pfx, sfx, fn, flen)
        filn(i)=fn(1:flen)
c
        if (elem_len(ie).eq.1) then
          pfx=elem(ie)
          pfx=pfx(1:1)//'_col'
        else
          pfx=elem(ie)
          pfx=pfx(1:2)//'_col'
        endif
        call newfile (pfx, sfx, fn, flen)
        filncol(i)=fn(1:flen)
      enddo
c     emmissivity structure
      do i=1,iemn
        idx=iem(i)
c
        at=fmatom(idx)
        io=fmion(idx)
        nt=nfmtrans(idx)
        nl=fmnl(idx)
        fn=' '
c
        pfx=elem(at)//rom(io)
        pfx=trim(pfx)//'_'//trim(jpfx)
c
        sfx='csv'
        call newfile (pfx, sfx, fn, flen)
        filnem(i)=fn(1:flen)
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine createp7files
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine createp7files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'p7blocks.inc'
c
      integer*4 i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Create All files NEW
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      photn
      open (luop,file=filna,status='NEW')
c      phapn
      open (lupf,file=filnb,status='NEW')
c      phlss
      open (lups,file=filnc,status='NEW')
c      spec
      open (lusp,file=filnd,status='NEW')
c      phsem
      open (lunt,file=filnt,status='NEW')
      if (jbnd.eq.'Y') then
c     bands
        open (lupb,file=filpb,status='NEW')
      endif
      if (jlin.eq.'Y') then
c     lines
        open (lulin,file=filin,status='NEW')
      endif
      if (jrat.eq.'Y') then
c     lines
        open (lurt,file=filrt,status='NEW')
      endif
      if (jall.eq.'Y') then
c     ions_
        open (lusl,file=filio,status='NEW')
      endif
      if (jiel.eq.'Y') then
        do i=1,ieln
          open (luions(i),file=filn(i),status='NEW')
        enddo
      endif
      if (jiem.eq.'Y') then
        do i=1,iemn
          open (luemiss(i),file=filnem(i),status='NEW')
        enddo
      endif
      if (jcol.eq.'Y') then
        do i=1,ieln
          open (lucols(i),file=filncol(i),status='NEW')
        enddo
      endif
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine appendp7files
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine appendp7files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'p7blocks.inc'
c
      integer*4 i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open All files APPEND
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      logical unitopen
c
c     photn
      inquire (unit=luop,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (luop,file=filna,status='OLD',access='APPEND')
      endif
c     phapn
      inquire (unit=lupf,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lupf,file=filnb,status='OLD',access='APPEND')
      endif
c     phlss
      inquire (unit=lups,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lups,file=filnc,status='OLD',access='APPEND')
      endif
c     spec
      inquire (unit=lusp,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lusp,file=filnd,status='OLD',access='APPEND')
      endif
c     phsem
      inquire (unit=lunt,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lunt,file=filnt,status='OLD',access='APPEND')
      endif
      if (jbnd.eq.'Y') then
c     bands
      inquire (unit=lupb,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lupb,file=filpb,status='OLD',access='APPEND')
      endif
      endif
      if (jlin.eq.'Y') then
c     lines
      inquire (unit=lulin,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lulin,file=filin,status='OLD',access='APPEND')
      endif
      endif
c     rates
      if (jrat.eq.'Y') then
      inquire (unit=lurt,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lurt,file=filrt,status='OLD',access='APPEND')
      endif
      endif
      if (jall.eq.'Y') then
c     old allions now ions_
      inquire (unit=lusl,opened=unitopen)
      if (unitopen.eqv..false.) then
        open (lusl,file=filio,status='OLD',access='APPEND')
      endif
      endif
c     individual elements
      if (jiel.eq.'Y') then
        do i=1,ieln
          inquire (unit=luions(i),opened=unitopen)
          if (unitopen.eqv..false.) then
            open (luions(i),file=filn(i),status='OLD',access='APPEND')
          endif
        enddo
      endif
c     emmissivity structure
      if (jiem.eq.'Y') then
        do i=1,iemn
          inquire (unit=luemiss(i),opened=unitopen)
          if (unitopen.eqv..false.) then
          open (luemiss(i),file=filnem(i),status='OLD',access='APPEND')
          endif
        enddo
      endif
c     column densities
      if (jcol.eq.'Y') then
        do i=1,ieln
        inquire (unit=luemiss(i),opened=unitopen)
        if (unitopen.eqv..false.) then
          open (lucols(i),file=filncol(i),status='OLD',access='APPEND')
        endif
        enddo
      endif
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine closep7files
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine closep7files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 'p7blocks.inc'
c
      integer*4 i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Close All P7 files If Open ignores closed units
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      logical unitopen
c
      inquire (unit=luop,opened=unitopen)
      if (unitopen) close (luop)
      inquire (unit=lupf,opened=unitopen)
      if (unitopen) close (lupf)
      inquire (unit=lups,opened=unitopen)
      if (unitopen) close (lups)
      inquire (unit=lusp,opened=unitopen)
      if (unitopen) close (lusp)
      inquire (unit=lunt,opened=unitopen)
      if (unitopen) close (lunt)
      inquire (unit=lupb,opened=unitopen)
      if (unitopen) close (lupb)
      inquire (unit=lulin,opened=unitopen)
      if (unitopen) close (lulin)
      inquire (unit=lurt,opened=unitopen)
      if (unitopen) close (lurt)
      inquire (unit=lusl,opened=unitopen)
      if (unitopen) close (lusl)
c
      do i=1,ieln
        inquire (unit=luions(i),opened=unitopen)
        if (unitopen) close (luions(i))
      enddo
c
      do i=1,iemn
        inquire (unit=luemiss(i),opened=unitopen)
        if (unitopen) close (luemiss(i))
      enddo
c
      do i=1,ieln
        inquire (unit=lucols(i),opened=unitopen)
        if (unitopen) close (lucols(i))
      enddo
c
      inquire (unit=ir1,opened=unitopen)
      if (unitopen) close (ir1)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
