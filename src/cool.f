cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'credits.inc'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine cool
c! XXXX - add one line purpose here
c! @param [in,out]   real*8         t  XXX-meaning
c! @param [in,out]   real*8        de  XXX-meaning
c! @param [in,out]   real*8        dh  XXX-meaning
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine cool (t, de, dh)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c*******COMPUTES TOTAL COOLING RATE OF PLASMA , TEMP : T
c     ELECTRON DENSITY : DE , HYDROGEN DENSITY : DH
c     RETURNS TLOSS : TOTAL COOLING RATE (ERG.CM-3.S-1)
c     CALL SUBR.  HYDRO,RESON,INTER,FORBID,FREFRE,COLOSS,
c                 PHEAT,NETGAIN,ALLRATES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
c
c total losses and gain
c
      real*8 tll, tgg
c
c effective  losses and gains including closed box reombination approx
c and allowing charge exchange to be +ve or -ve cooling
c
      real*8 ell, egg, ett
c
      character jjmod*4
      real*8 t, de, dh
      real*8 csum
      integer*4 i, j
c
      tll=0.d0
      tgg=0.d0
c
c clear element/ion cooling totals
c
      do j=1,atypes
        coolz(j)=0.d0
        heatz(j)=0.d0
        do i=1,maxion(j)
          coolzion(i,j)=0.d0
          heatzion(i,j)=0.d0
        enddo
      enddo
      oiii5007loss=0.d0
c
c    ***COMPUTES NEW RATES IF TEMP. OR PHOTON FIELD HAVE CHANGED
c
      if ((t.le.mintemp).or.(dh.le.0.d0)) then
        write (*,*) 'Cool out of range',t,de,dh
        stop
      endif
c
      jjmod='ALL'
      call allrates (t, jjmod)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c cooling
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***New hydrogen + helium II line cooling
c    ***Also calls 2 photon and He I calcs
c
      hheloss=0.d0
      hloss=0.d0
      if (hhecollmode.eq.0) call hhecoll (t, de, dh)
      call hydro (t, de, dh)
c includes He coll losses
      tll=hloss
c
c Heavy element recombination, negligible cooling
c
      call recom_cii (t, de, dh)
      call recom_nii (t, de, dh)
      call recom_oi (t, de, dh)
      call recom_oii (t, de, dh)
      call recom_neii (t, de, dh)
c
c    LEGACY RESONANCE,INTERCOMBINATION AND FORBIDDEN LINES
c
      call fine3 (t, de, dh)
      tll=tll+f3loss
      call inter (t, de, dh)
      tll=tll+fslos
c      call inter2 (t, de, dh)
c
c multi-levels must be called before the reson routines
c
      fmloss=0.d0
      call multilevel (t, de, dh)
      tll=tll+fmloss
      feloss=0.d0
      call multiiron (t, de, dh)
      tll=tll+feloss
c
c     rloss=0.0d0
c     call reson (t, de, dh)
c     call reson2 (t, de, dh)
c
      xr3loss=0.d0
      xrlloss=0.d0
      call reson3 (t, de, dh)
      call resonl (t, de, dh)
      tll=tll+xr3loss+xrlloss
c
c      call helif (t, de, dh)
c
c
c     tll=tll+rloss+fslos+fmloss+xr3loss+xrlloss
c     tll=tll+f3loss+feloss
cc
c
c    ***FREE-FREE COOLING
c
      fflos=0.d0
      call frefre (t, de, dh)
      tll=tll+fflos
c
c
c    ***COLLISIONAL IONISATION LOSSES
c
      cmplos=0.d0
c     colos
      call coloss (de, dh)
      tll=tll+colos
c
c      tll=hloss+rloss+fslos+fmloss+xrloss+xr3loss+xrlloss+xiloss
c      tll=tll+fflos+colos
c      tll=tll+f3loss+feloss+gcool
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c heating
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***RECOMB. COOLING (PLUS ESTIMATE OF ON THE SPOT HEATING IF APPL)
c    can be heating or cooling see below
c
      rngain=0.d0
      call netgain (t, de, dh)
c     rngain can be + or -
c
c    ***PHOTOIONISATION HEATING
c
      pgain=0.d0
      call pheat (de, dh)
      tgg=tgg+pgain
c
c grain / pah heating and cooling
c
      chgain=0.d0
      call hgrains (t, de, dh)
      tll=tll+gcool
      tgg=tgg+gheat
      paheat=0.d0
      call hpahs (t, de, dh)
      tgg=tgg+paheat
c
c    ***CHARGE EXCHANGE HEATING
c
      chgain=0.d0
      call cheat (dh)
c     chgain can be + or -
c
c    ***COSMIC RAY HEATING
c
c      call cosmic (t, de, dh)
      call cosmic (dh)
      tgg=tgg+cosgain
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    ***(SEPARATELY) EFFECTIVE LOSS AND GAIN  :  ELOSS,EGAIN
c
      ell=tll
c
      egg=tgg
c
      if (rngain.lt.0.0d0) then
        ell=ell-rngain
      else
        egg=egg+rngain
      endif
c
      if (chgain.lt.0.0d0) then
        ell=ell-chgain
      else
        egg=egg+chgain
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ett=ell+egg
      if (ett.gt.0.d0) then
        dlos=(ell-egg)/ett
      else
        dlos=1.0d0
      endif
c
      csum=0.d0
      do j=1,atypes
        do i=1,maxion(j)
          if (dabs(coolzion(i,j)).lt.epsilon) coolzion(i,j)=0.d0
        enddo
        if (dabs(coolz(j)).lt.epsilon) coolz(j)=0.d0
        csum=csum+coolz(j)
      enddo
c
      tloss=tll
      tgain=tgg
      egain=egg
      eloss=ell
c
      if (alphacoolmode.eq.1) then
        tloss=de*de*alphac0*((1.0d-6*t)**alphaclaw)
        tgain=0.d0
        egain=0.d0
        eloss=tloss
      endif
c
      if (dabs(tloss).lt.epsilon) tloss=0.d0
      if (dabs(tgain).lt.epsilon) tgain=0.d0
      if (dabs(egain).lt.epsilon) egain=0.d0
      if (dabs(eloss).lt.epsilon) eloss=0.d0
c
c       write(*,*) 'Cool:',tloss, tgain, egain, eloss, csum
c
      return
      end
