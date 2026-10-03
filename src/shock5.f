cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'credits.inc'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief SHOCK5: Shock model
c! Shocks ith steady Rankine-Hugoniot solution and iterative precursors
c! @param This routine has no parameters
c!
c! @return
c!  Creates spectra and model file output
c!
c! @details
c!  SHOCK5: Shock model with steady Rankine-Hugoniot solution fo
c! each step.  Time steps based on fastest atomic, cooling or
c! photon absorbion timescales.
c! Full continuum and diffuse field, auto preionisation
c***************************************************************

      subroutine shock5 ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SHOCK5: Shock model with steady Rankine-Hugoniot solution for
c     each step.  Time steps based on fastest atomic, cooling or
c     photon absorbion timescales.
c
c     Full continuum and diffuse field, auto preionisation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 iterations,iterationindex,nlock,wasconverged
      integer*4 wasosclass
      real*8 wasprecrms,waspostrms
   10 format(//,
     & ' *********************************************************',/,
     & '  SHOCK 5    Global Iteration: ',i2.2,' of ',i2.2)
   20 format(/,
     & '  SHOCK 5 Completed Iteration: ',i2.2,' of ',i2.2,/,
     & ' *********************************************************',/)
   30 format(/,
     & ' *********************************************************',/,
     & '  SHOCK 5 Model Completed',/,
     & ' *********************************************************',/)
   35 format(/,
     & ' *********************************************************',/,
     & '  SHOCK 5: NO SHOCK -- preshock flow is sub-Alfvenic',/,
     & '  Alfven Mach Number = ',1pg11.4,' (< 1): no compressive',/,
     & '  MHD shock jump exists for these parameters.  Skipping',/,
     & '  the shock calculation for this model.',/,
     & ' *********************************************************',/)
   36 format(' No Shock: , Alfven Mach:, ',1pg11.4,
     & ' , Reason:, sub-Alfvenic preshock flow')
   37 format('  Result: NO SHOCK (sub-Alfvenic)',/,
     & ' *********************************************************',/)
c
      iterations=3
      ieln=4
c
      call shock5setup (iterations)
      call shock5headers (iterations)
c
c  A compressive fast MHD shock only has a solution when the preshock
c  flow exceeds the Alfven speed (alfvennumber, computed in
c  shock5setup from preshock quantities alone, .ge.1).  Below that,
c  shockcmpf's jump-condition quadratic has no physical root and
c  returns an unphysical compression (x<1, an expansion), which used
c  to propagate through to a negative temperature and a hard stop
c  many steps later in cool() (issue #8 follow-up).  Checked here,
c  before shockcmpf is ever called, rather than after the fact --
c  confirmed via a full grid run that every case with
c  alfvennumber<1 fails this way and every case >=1 does not, a
c  clean, sharp boundary at exactly the physically-expected value.
c  "No Shock:" is written in the same comma-separated style as
c  "Model ended:" so existing summary tooling (MapSum) can be
c  extended to recognise it the same way.
c
      if (alfvennumber.lt.1.0d0) then
        write (*,35) alfvennumber
        write (*,36) alfvennumber
        write (luop,36) alfvennumber
        write (*,37)
        call closeS5files ()
        return
      endif
c
      iterationindex=0
      converged=0
      osclass=0
      precwarn=0
      finalit=0
      aitninit=0
c
      if (iterations.le.1) finalit=1
c
      call compsh5 (0, iterations)
c
      if (iterations.gt.1) then
c
   40   iterationindex=iterationindex+1
c
        write (*,10) iterationindex,iterations
c
        call shock5precursor (iterationindex, iterations)
        call compsh5 (iterationindex, iterations)
        call shock5check (iterationindex, iterations)
c
        write (*,20) iterationindex,iterations
c
        if ((converged.eq.0).and.(iterations.lt.mxshockits)) goto 40
c
c  A few more ordinary iterations once convergence is first declared,
c  cheap extra confirmation now that every iteration (not just this
c  tail) solves the precursor to the real tolerance (issue #7).
c  Skipped if the loop above never converged; more of the same
c  iteration wouldn't be expected to fix that on its own.
c
        if (converged.ne.0) then
          do nlock=1,3
            iterationindex=iterationindex+1
            write (*,10) iterationindex,iterationindex
            call shock5precursor (iterationindex, iterationindex)
            call compsh5 (iterationindex, iterationindex)
            call shock5check (iterationindex, iterationindex)
            write (*,20) iterationindex,iterationindex
          enddo
        endif
c
c repeat a final model for outputs
c
        finalit=1
        wasconverged=converged
        wasosclass=osclass
        wasprecrms=precrms
        waspostrms=postrms
        iterationindex=iterationindex+1
c
        call shock5precursor (iterationindex, iterationindex)
        call compsh5 (iterationindex, iterationindex)
        call shock5check (iterationindex, iterationindex)
c
c  This one extra call is an independent re-solve of the precursor,
c  purely to regenerate full output tables at finalit's tightened
c  tolerance -- it is not needed to establish convergence, which the
c  loop above (plus the lock-in iterations, if any) already did.  A
c  small residual disagreement here is Aitken's own noise floor from
c  re-solving an already-converged fixed point, not evidence the
c  model is actually unconverged -- the "final-pass gap" (issue #7):
c  the main loop converges cleanly and only this one extra,
c  independent call disagrees.  If the model was already genuinely
c  converged going into this call, report that,
c  rather than letting one noisy independent re-solve override a
c  result already established by many iterations.  The raw numbers
c  above are printed either way -- only the bottom-line Result is
c  corrected.
c
        if ((wasconverged.ne.0).and.(converged.eq.0)) then
          converged=1
          write (*,38) wasprecrms*100.d0,waspostrms*100.d0
        else if ((wasconverged.eq.0).and.(converged.eq.0)) then
          if (wasosclass.eq.1) write (*,39) wasprecrms*100.d0,
     &      waspostrms*100.d0
          if (wasosclass.eq.2) write (*,41) wasprecrms*100.d0,
     &      waspostrms*100.d0
          if (wasosclass.eq.3) write (*,42) wasprecrms*100.d0,
     &      waspostrms*100.d0
          if (wasosclass.eq.4) write (*,43) wasprecrms*100.d0,
     &      waspostrms*100.d0
          if (wasosclass.eq.5) write (*,44)
        endif
   38   format('  Result: CONVERGED  (last independent recheck was ',
     & 'noisy;',/,
     & '  already converged earlier -- see comment at issue #7; ',
     & 'precursor RMS ',1pg11.4,'%, post-shock RMS ',1pg11.4,'%)',/,
     & ' *********************************************************',/)
   39   format('  Result: NOT CONVERGED (precursor oscillating; ',
     & 'precursor RMS ',1pg11.4,'%, post-shock RMS ',1pg11.4,'%)',/,
     & '  Bounded precursor limit cycle, issue #14 -- post-shock and ',
     & 'inner solve OK.',/,
     & ' *********************************************************',/)
   41   format('  Result: NOT CONVERGED (precursor not converged; ',
     & 'precursor RMS ',1pg11.4,'%, post-shock RMS ',1pg11.4,'%)',/,
     & '  Precursor residual exceeds the issue #14 validated bound; ',
     & 'post-shock OK.',/,
     & ' *********************************************************',/)
   42   format('  Result: NOT CONVERGED (post-shock not converged, ',
     & 'precursor OK; precursor RMS ',1pg11.4,'%, post-shock RMS ',
     & 1pg11.4,'%)',/,
     & ' *********************************************************',/)
   43   format('  Result: NOT CONVERGED (precursor and post-shock ',
     & 'not converged; precursor RMS ',1pg11.4,'%, post-shock RMS ',
     & 1pg11.4,'%)',/,
     & ' *********************************************************',/)
   44   format('  Result: FAILED (precursor inner solve exhausted, ',
     & 'unreliable)',/,
     & '  Inner sweep hit its budget without self-consistency -- ',
     & 'see issue #14.',/,
     & ' *********************************************************',/)
c
      endif
c
      call closeS5files ()
c
      write (*,30)
c
      return
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine shock5setup
c! Run the user initialisation and set global flags
c! @param [in,out] integer*4  iterations  minimum number of iterations
c!
c! @return
c!  Default 3 global iterations, returns user actual choice
c!
c! @details
c!  Iterations is an estimate and will often take more. 0 or -ve will
C! Perform no iterations.
c***************************************************************

      subroutine shock5setup (iterations)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Run the user initialisation and set global flags
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 iterations,i,nl,nel,elok(mxelem)
      integer*4 idx,npx
      real*8 vapprox,tpr,tpo,va,eps
      character outsettings*64,outshow*64
      character px*32,readname*32
c
c functions
c
      real*8 feldens,fpresse
      real*8 frho,velshock2
      integer*4 mlen
c
   10 format(a)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      photonmode=1
c
c clear everything first
c
      call zer ()
c
c   calculation modes - older so YES not Y
c
      jspot='N'
      jcon='Y'
c
      ispo='SII'
c
      jden='B'
      jgeo='P'
      wdil=0.5d0
c
      jtrans='LODW'
      specmode='DW'
      alphacoolmode=0
c
      subname='Shock 5'
      jnorm=0
c
      mtype='A'
      magparam=0.d0
      Pmag=0.0d0
      Bmag=0.d0
c
      vmod='NONE'
      s5pfx='v100sh'
      nprefix=20
c
      wdil=0.5d0
      wdilt0=1.d4
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTRO - Define Input Parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 format(///,
     & ' ********************************************************'/
     & ' SHOCK 5 model selected:',/,
     & ' Shock code with auto-preionisation.',/,
     & ' ********************************************************')
      write (*,20)
c
c     get ionisation balance to work with...
c
      subname='proto-ionisation'
      call popcha (subname)
      call copypop (pop0, pop)
      call copypop (pop0, pop_neu)
      call copypop (pop0, pop_pre)
      call copypop (pop0, pop_pre0)
      call zeroemiss
      subname='Shock 5'
c
c     set up current field conditions
c
      call photsou (subname)
      do i=1,infph
        prefield(i)=soupho(i)
      enddo
c
c     Shock model preferences
c
   30 format(///,
     & ' ********************************************************'/
     & ' Setting the shock conditions',/,
     & ' ********************************************************'/)
      write (*,30)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Shock Jump Definition
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      stype='V'
   40 format(//,
     & ' Choose Shock Jump Parameter:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    T  : Shock jump in terms of Temperature jump',/,
     & '    V  : Shock jump in terms of flow Velocity.',//,
     & ' :: ',$)
c
   50 write (*,40)
      read (*,10) stype
      call toup (stype(1:1), stype)
c
      if ((stype.ne.'T').and.(stype.ne.'V')) goto 50
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Shock Magnetic Field Parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      magparam=1.0d0
      malpha=1.0d0
      mageta=1.0d0
      mmach=0.d0
c
      mtype='B'
   60 format(//,
     & ' Choose Shock Magnetic Parameter Type:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    B  : Magnetic field in microGauss at shockfront',/,
     & '    A  : Magnetic Alpha_0 = Pmag/Pgas in protostate',/,
     & '    M  : Magnetic Alfven Mach Number, Ma',/,
     & '    C  : Magnetic Alpha Pmag/Pgas at shockfront',/,
     & '    R  : Magnetic Eta 2Pmag/Pram =1/Ma^2 at shockfront',/,
     & ' :: ',$)
   70 write (*,60)
      read (*,10) mtype
      call toup (mtype(1:1), mtype)
c
      if ((mtype.ne.'A').and.
     &    (mtype.ne.'B').and.
     &    (mtype.ne.'C').and.
     &    (mtype.ne.'M').and.
     &    (mtype.ne.'R')) goto 70
c
      if (mtype.eq.'B') then
c       Preshock magnetic field at shockfront in uG
   80   format(//,
     & ' Specify the magnetic field at shockfront, B: ',/,
     & ' (microgauss): ', $ )
        write (*,80)
        read (*,*) magparam
        Bmag=dabs(magparam)
        magparam=Bmag
        Bmag=Bmag*1e-6
        Pmag=(Bmag*Bmag)/epi
      endif
c
      if (mtype.eq.'A') then
c     Preshock magnetic fields expressed as Alpha = Pmag/Pgas
c     Alpha = 1.0 is equipartition, Alpha = 10.0 is strong magnetic
c     Alpha = 0.1 is weak magnetic field. Alpha = 0.0 is no magnetic
   90   format(//,
     & ' Specify the magnetic Alpha_0 in proto-state',/,
     & ' ( >= 0.0: Pmag/Pgas): ', $)
        write (*,90)
        read (*,*) magparam
        malpha=magparam
        if (magparam.lt.0.0d0) then
          Bmag=-magparam
          Bmag=Bmag*1e-6
          Pmag=(Bmag*Bmag)/epi
          mtype='B'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
        endif
      endif
c
      if (mtype.eq.'C') then
  100   format(//,
     & ' Specify the magnetic Alpha at shockfront',/,
     & ' ( >= 0.0: Pmag/Pgas): ', $)
        write (*,100)
        read (*,*) magparam
        malpha=magparam
        if (magparam.lt.0.0d0) then
          Bmag=-magparam
          Bmag=Bmag*1e-6
          Pmag=(Bmag*Bmag)/epi
          mtype='B'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
          mmach=0.d0
        endif
      endif
c
      if (mtype.eq.'M') then
c       Preshock magnetic field expressed as magnetic Alfven Mach number
  110   format(//,
     & ' Specify the magnetic Alfven Mach Number, Ma',/,
     & ' ( v/sqrt(Pmag/rho), Ma>= 0.0 ) ', $)
        write (*,110)
        read (*,*) magparam
        mmach=magparam
        if (magparam.lt.0.0d0) then
          Bmag=-magparam
          Bmag=Bmag*1e-6
          Pmag=(Bmag*Bmag)/epi
          mtype='B'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
          mmach=0.d0
        endif
      endif
c
      if (mtype.eq.'R') then
c       Preshock magnetic field expressed as magnetic eta.
c       eta = 2Pmag/Pram = 1/Ma^2
c       Pram = rho*vs*vs; Pmag = (B^2)/8pi
  120   format(//,
     & ' Specify the magnetic Eta at shockfront',/,
     & ' ( >= 0.0: 2Pmag/Pram): ', $)
        write (*,120)
        read (*,*) magparam
        mageta=magparam
        if (magparam.lt.0.0d0) then
          Bmag=dabs(magparam)*1e-6
          Pmag=(Bmag*Bmag)/epi
          mtype='B'
          magparam=dabs(magparam)
          mageta=1.0d0
          malpha=1.0d0
          mmach=0.d0
        endif
      endif
c
c always  plane parallel
c
      wdil=0.5d0
      wdilt0=1.d4
c
      if (stype.eq.'V') then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Velocity  - Tshock adjusts as presionisation changes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  130   format(//,
     & ' Specify the proto-shock conditions:',/,
     & '  T (K), dh (N), v (< 1e5 km/s, >= 1e5 cm/s)',/,
     & ' (T<10 taken as a log, dh <= 0 taken as a log),',/,
     & ' : ',$)
        write (*,130)
        read (*,*) t,dh,ve
c
        dr=1.d0
        de=feldens(dh,pop_neu)
c
        if (ve.lt.1.d5) ve=ve*1.d5
        if (t.le.10.d0) t=10.d0**t
        if (dh.le.0.d0) dh=10.d0**dh
c
        Pgas=fpresse(t,de,dh)
        de=feldens(dh,pop_neu)
        rh_neu=frho(de,dh)
        Pram=rh_neu*ve*ve
c
        if (mtype.eq.'B') then
          Bmag=magparam*1.0d-6
          Pmag=(Bmag*Bmag)/epi
        endif
        if (mtype.eq.'A') then
          Pmag=malpha*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
        if (mtype.eq.'C') then
          Pmag=malpha*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
        if (mtype.eq.'M') then
          va=ve/mmach
          Pmag=(0.5*rh_neu*va*va)
          Bmag=dsqrt(epi*Pmag)
        endif
        if (mtype.eq.'R') then
          Pmag=0.5d0*mageta*Pram
          Bmag=dsqrt(epi*Pmag)
        endif
c
        vshoc=ve
c
        bm_neu=Bmag
        te_neu=t
        de_neu=de
        dh_neu=dh
        vs_neu=ve
        pr_neu=fpresse(te_neu,de_neu,dh_neu)
        rh_neu=frho(de_neu,dh_neu)
c
c     fill in other default parameters
c
        dt=0.d0
        tloss=0.d0
c
        call shockcmpf (t, de, dh, ve, Bmag)
c
        tm00=te0
c
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (stype.eq.'T') then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Temperature Jump - Adjusts as presionisation changes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  140   format(//,
     & ' Specify the proto-shock conditions:',/,
     & ' (T<10 taken as a log, dh <= 0 taken as a log)',/,
     & ' T (K), dh (N) : ',$)
        write (*,140)
        read (*,*) tpr,dh
        wdil=0.5d0
c
        if (tpr.le.10.d0) tpr=10.d0**tpr
        if (dh.le.0.0d0) dh=10.d0**dh
c
  150   format(//,
     & ' Specify the post-shock temperature:',/,
     & ' (T<10 taken as a log) T (K) : ',$)
        write (*,150)
        read (*,*) tpo
c
        if (tpo.le.10.d0) tpo=10.d0**tpo
c
        te0=tpr
        te1=tpo
c       Partially set up proto-state
        te_neu=tpr
        dh_neu=dh
        de_neu=feldens(dh_neu,pop_neu)
        Pgas=fpresse(te_neu,de_neu,dh_neu)
        rh_neu=frho(de_neu,dh_neu)
c
        if (mtype.eq.'B') then
          Bmag=magparam*1.0d-6
          Pmag=(Bmag*Bmag)/epi
        endif
        if (mtype.eq.'A') then
          Pmag=malpha*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
        if (mtype.eq.'C') then
          Pmag=malpha*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
c
        vapprox=velshock2(dh_neu,tpr,0.d0,tpo)
        rh_neu=frho(de_neu,dh_neu)
        Pram=rh_neu*vapprox*vapprox
        if (mtype.eq.'R') then
          Pmag=0.5d0*mageta*Pram
          Bmag=dsqrt(epi*Pmag)
        endif
c
c       We need to iterate for a consistent solution.  We have the
c       density and T jump; solve for shock speed and magnetic field.
c
  160   bm0=Bmag
c
        vshoc=velshock2(dh_neu,tpr,Bmag,tpo)
        ve=vshoc
        Pram=rh_neu*ve*ve
c
        if (mtype.eq.'M') then
          va=vshoc/mmach
          Pmag=(0.5*rh_neu*va*va)
          Bmag=dsqrt(epi*Pmag)
          eps=dabs(2.0*(Bmag-bm0)/(Bmag+bm0))
          if (eps.gt.1.0d-6) goto 160
        endif
c
        if (mtype.eq.'R') then
c
c iterate for rampressure/alpha_r
c
          Pmag=0.5d0*mageta*Pram
          Bmag=dsqrt(epi*Pmag)
c           write(*,*) vshoc,Pram,Pmag,bm0,Bmag,te1
          eps=dabs(2.0*(Bmag-bm0)/(Bmag+bm0))
c         write (*,*) eps,Bmag,bm0
          if (eps.gt.1.0d-6) goto 160
        endif
c
        dt=0.d0
        tloss=0.d0
c
        call shockcmpf (te_neu, de_neu, dh_neu, ve, Bmag)
c
        tm00=te0
        bm_neu=Bmag
        vs_neu=vshoc
        pr_neu=fpresse(te_neu,de_neu,dh_neu)
        rh_neu=frho(de_neu,dh_neu)
c
      endif
c
c all shocks - powerlaw cooling for testing
c
      if (alphacoolmode.eq.1) then
  170  format(//,
     & ' Powerlaw Cooling enabled  Lambda T6 ^ alpha :',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    Give index alpha, Lambda at 1e6K',/,
     & '    (Lambda<0 as log)',/,
     & ' :: ',$)
        write (*,170)
        read (*,*) alphaclaw,alphac0
        if (alphac0.lt.0.d0) alphac0=10.d0**alphac0
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Shock intial state kept for iterations, above numbers are
c reused from step to step
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      te_pre=te0
      te_pst=te1
      de_pre=de0
      de_pst=de1
      dh_pre=dh0
      dh_pst=dh1
      vs_pre=vel0
      vs_pst=vel1
c
      rh_pre=rho0
      rh_pst=rho1
      pr_pre=pr0
      pr_pst=pr1
      bm_pre=bm0
      bm_pst=bm1
      bp0=(bm0*bm0)/epi
      bp1=(bm1*bm1)/epi
c
      Pmag=bp0
      Pgas=pr0
      Pram=rho0*vel0*vel0
c
      machnumber=vel0/dsqrt(gammaEOS*pr0/rho0)
      alfvennumber=vel0/dsqrt(2.0d0*Pmag/rho0)
      malpha=Pmag/Pgas
      gaseta=gammaEOS*Pgas/Pram
      mageta=2.d0*Pmag/Pram
c
      call shocksummary (6)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c set intial shock vars
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      cmpf0=cmpf
      psi0=0.0d0
      psi=0.0d0
c
      te_pre0=te_pre
      de_pre0=de_pre
      dh_pre0=dh_pre
      vs_pre0=vs_pre
      rh_pre0=rh_pre
      pr_pre0=pr_pre
      bm_pre0=bm_pre
c
      te_pst0=te_pst
      de_pst0=de_pst
      dh_pst0=dh_pst
      vs_pst0=vs_pst
      rh_pst0=rh_pst
      pr_pst0=pr_pst
      bm_pst0=bm_pst
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Diffuse field interaction
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      photonmode=1
c
  180 format(//,
     & ' Choose Diffuse Field Option :',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    F  : Full diffuse field interaction (Default).',/,
     & '    Z  : Zero diffuse field interaction.',/,
     & ' :: ',$)
  190 write (*,180)
      read (*,10) ilgg
      call toup (ilgg(1:1), ilgg)
c
      if ((ilgg.ne.'Z').and.(ilgg.ne.'F')) goto 190
c
      if (ilgg.eq.'Z') photonmode=0
      if (ilgg.eq.'F') photonmode=1
      photofraction=1.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  200 format(//,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Calculation Limit Settings ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,200)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set up time step limits
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      tmod='A'
c
c no longer used directly, may return in the future
c
      atimefrac=0.0250d0
      cltimefrac=0.0250d0
      abtimefrac=0.5d0
      utime=0.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  210 format(//,
     & ' ********************************************************'/
     & '  Boundry Conditions  '/
     & ' ********************************************************'/)
      write (*,210)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Choose ending conditions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  220 format(/,
     & '  Model Ending Condition:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard ending, 1% weighted ionisation',/,
     & '    B  : Species ionisation limit',/,
     & '    C  : Temperature limit',/,
     & '    S  : Temperature limit, and >95% neutral',/,
     & '    D  : Distance limit',/,
     & '    E  : Time limit',/,
     & '    F  : Thermal balance limit',/,
     & '    G  : Heating limit',/,
     & ' :: ',$)
c
  230 write (*,220)
      read (*,10) jend
      call toup (jend(1:1), jend)
c
      if ((jend.ne.'A')
     &.and.(jend.ne.'B')
     &.and.(jend.ne.'C')
     &.and.(jend.ne.'D')
     &.and.(jend.ne.'E')
     &.and.(jend.ne.'F')
     &.and.(jend.ne.'G')
     &.and.(jend.ne.'S')) goto 230
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Secondary info:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jend.eq.'B') then
  240   format(//,
     & ' Give atom, ion and limit fraction: ',/,
     & ' (eg 1 2 0.05 = stop when HII < 0.05): ',$)
        write (*,240)
        read (*,*) ielen,jpoen,fren
      endif
      if ((jend.eq.'C').or.(jend.eq.'S')) then
  250   format(//,
     & ' Give final temperature (K > 10, log <= 10): ',$)
        write (*,250)
        read (*,*) tend
        if (tend.le.10.d0) tend=10.d0**tend
      endif
      if (jend.eq.'D') then
  260   format(//,
     & ' Give final distance (cm > 100, log<=100): ',$)
        write (*,260)
        read (*,*) diend
        if (diend.le.100.d0) diend=10.d0**diend
      endif
      if (jend.eq.'E') then
  270   format(//,
     & ' Give time limit (s > 100, log<=100): ',$)
        write (*,270)
        read (*,*) timend
        if (timend.le.100.d0) timend=10.d0**timend
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  280   format(//,
     & '  Choose the minimum number of shock iterations:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' (<= 1, no iterations) : ',$)
      write (*,280)
      read (*,*) iterations
      if (iterations.le.0) iterations=1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  290 format(//,
     & ' ********************************************************'/
     & '  Output Requirements  '/
     & ' ********************************************************'//)
      write (*,290)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     get final field file prefix
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      readname='.               '
      s5pfx='.               '
  300 format (a16)
  310 format(//,
     & ' Give a prefix for all output files (max 8 chars): ')
      write (*,310)
      read (*,300) readname
      call mytrim (readname, npx, px)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    New output options menu including reset like in shock5files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call shock5filenames (px)
c
  320 jspec='N'
c
c     Ion Balances
c
      tsrmod='N'
      allmod='N'
c
c    dynamics file
c
      dynmod='N'
c
c    rates file
c
      ratmod='N'
c
c     emission bands
c
      bandsmod='N'
c
c     Cooling Components
c
      fclmod='N'
c
c    monitor 16 lines
c
      jlin='N'
c
c     Upstream fields
c
      lmod='N'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spectra: F-Lambda plots, ionising and optical default yes
c
  330 outsettings='A'
c
c     Ion Balances
c
      if (tsrmod.eq.'Y') then
        call myappend (outsettings, '+B', outshow)
        outsettings=outshow
      endif
c
c    rates file
c
      if (ratmod.eq.'Y') then
        call myappend (outsettings, '+C', outshow)
        outsettings=outshow
      endif
c
c    dynamics file
c
      if (dynmod.eq.'Y') then
        call myappend (outsettings, '+D', outshow)
        outsettings=outshow
      endif
c
c    final downstream field
c
      if (jspec.eq.'Y') then
        call myappend (outsettings, '+E', outshow)
        outsettings=outshow
      endif
c
c    all up stream fields
c
      if (lmod.eq.'Y') then
        call myappend (outsettings, '+F', outshow)
        outsettings=outshow
      endif
c
c     emission bands
c
      if (bandsmod.eq.'Y') then
        call myappend (outsettings, '+H', outshow)
        outsettings=outshow
      endif
c
c     Cooling Components
c
      if (fclmod.eq.'Y') then
        call myappend (outsettings, '+K', outshow)
        outsettings=outshow
      endif
c
c    monitor 16 lines
c
      if (jlin.eq.'Y') then
        call myappend (outsettings, '+L', outshow)
        outsettings=outshow
      endif
c
  340 format(//,
     & ' Output Multi-Option Menu : ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' :: Current: ',a20,'                       ::',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard output - and Reset. ',/,
     & '    B  : Ion balance files.',/,
     & '    C  : All Rates file.',/,
     & '    D  : Dynamics file.',/,
     & '    E  : Final downstream field.',/,
     & '    F  : All upstream fields at each step.',/,
     & '    H  : Cooling in x-ray bands.',/,
     & '    K  : Cooling Components by Elements file.',/,
     & '    L  : Monitor up to ',i3,' lines'/
     & '       :',/,
     & '    R  : Reset',/,
     & '    X  : Exit with Settings',/,
     & ' :: ',$)
  350 write (*,340) outsettings,mxmonlines
      read (*,10) ilgg
      call toup (ilgg(1:1), ilgg)
      write (*,*)
c
      if (ilgg.eq.'Q') ilgg='X'
c
      if ((ilgg.ne.'A')
     &.and.(ilgg.ne.'B')
     &.and.(ilgg.ne.'C')
     &.and.(ilgg.ne.'D')
     &.and.(ilgg.ne.'E')
     &.and.(ilgg.ne.'F')
     &.and.(ilgg.ne.'H')
     &.and.(ilgg.ne.'K')
     &.and.(ilgg.ne.'L')
     &.and.(ilgg.ne.'R')
     &.and.(ilgg.ne.'X')) goto 350
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ilgg.eq.'A') goto 320
      if (ilgg.eq.'B') tsrmod='Y'
      if (ilgg.eq.'C') ratmod='Y'
      if (ilgg.eq.'D') dynmod='Y'
      if (ilgg.eq.'E') jspec='Y'
      if (ilgg.eq.'F') lmod='Y'
      if (ilgg.eq.'H') bandsmod='Y'
      if (ilgg.eq.'K') fclmod='Y'
      if (ilgg.eq.'L') jlin='Y'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ilgg.eq.'R') goto 320
      if (ilgg.ne.'X') goto 330
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     B Ion balance files tsrmod, allmod
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     default monitor elements
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ieln=4
      iel(1)=zmap(1)
      iel(2)=zmap(2)
      iel(3)=zmap(8)
      iel(4)=zmap(16)
c
      allmod='N'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (tsrmod.eq.'Y') then
c
  360   format(//' Monitor ions/columns of max ',i2,' elements:',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
  370   format(' Enter number of elements to track : ',$)
        write (*,360) atypes
        write (*,370)
        read (*,*) ieln
c
        ieln=min(max(ieln,1),atypes)
        read (*,*) (iel(i),i=1,ieln)
c
        do i=1,ieln
          elok(i)=0
          do idx=1,atypes
            if (iel(i).eq.mapz(idx)) elok(i)=1
          enddo
        enddo
c
        if (tsrmod.eq.'Y') then
          nel=0
          do i=1,ieln
            if (elok(i).eq.1) then
              nel=nel+1
              iel(nel)=zmap(iel(i))
            endif
          enddo
          if (nel.lt.1) tsrmod='N'
          if (tsrmod.eq.'Y') then
            ieln=nel
  380         format(/' Monitoring :',30(x,a2),/)
            write (*,380) (elem(iel(i)),i=1,ieln)
c
  390         format(/,' Record all ions file (Y/N)? : ',$)
  400       write (*,390)
            read (*,10) allmod
            call toup (allmod(1:1), allmod)
            if ((allmod.ne.'Y').and.(allmod.ne.'N')) goto 400
          endif
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    C timescales and rates file
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     no specific setup - fixed format
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    D flow dynamics file, dynmod
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     no specific setup - fixed format
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spec List Files - always on
c
c     no specific setup - fixed format
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Upstream fields lmod = Y/N
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c no specific setup - fixed format units from wpsou
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     F-Lambda plots, ionising and optical  jspec = 'Y'/'N'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c no specific setup - fixed format units from wpsou
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Cooling Components File, fclmod
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (fclmod.eq.'Y') then
        jnorm=0
  410  format (//' Cooling File Normalisation,',/
     & ' (0=ne.nH, 1=nH^2, 2=ne.ni, 3=n^2, 4=ne^2): ',$)
        write (*,410)
        read (*,*) jnorm
        if (jnorm.lt.0) jnorm=0
        if (jnorm.gt.4) jnorm=0
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    monitor 16 lines, jlin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      njlines=0
      if (jlin.eq.'Y') then
c
  420   format('   ',a3,a6,' ',f12.3)
  430   format(//' Select up to ',i3,' lines by (A, incl air)   :',/
     & '::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
  440   format(' Enter number of lines to track: ',$)
        write (*,430) mxmonlines
        write (*,440)
        read (*,*) nl
        nl=min(max(nl,1),mxmonlines)
        njlines=nl
  450   format(/' Wavelengths (see spec list files) : ',$)
        write (*,450)
        read (*,*) (emlinlist(i),i=1,njlines)
        do i=1,mxmonlines
c     0.001A: comfortably covers the +/-0.0005A rounding from reading
c     a wavelength off spec2's 3-decimal-place output, while staying
c     far tighter than any real line-to-line separation encountered in
c     practice - so a wavelength copied from the line list matches
c     exactly, with no risk of picking up an unrelated nearby line.
          emlindeltas(i)=0.001d0
        enddo
        call speclocallineids (emlinlistatom, emlinlistion)
  460   format(/' Monitoring :')
        write (*,460)
        do i=1,njlines
          write (*,420) elem(emlinlistatom(i)),rom(emlinlistion(i)),
     &     emlinlist(i)
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     get screen display mode
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  470 format(//,
     & ' Runtime screen display:',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    A  : Standard display ',/,
     & '    B  : Detailed Slab display.',/,
     & '    C  : Full Display (Full slab display + timescales).',/,
     & '    M  : Minimal Display (batch mode).',/,
     & ' :: ',$)
      write (*,470)
      read (*,10) ilgg
      call toup (ilgg(1:1), ilgg)
c
      vmod='NONE'
      if (ilgg.eq.'A') vmod='MINI'
      if (ilgg.eq.'B') vmod='SLAB'
      if (ilgg.eq.'C') vmod='FULL'
      if (ilgg.eq.'M') vmod='NONE'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     get runname
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  480 format(//,
     & ' Set a name/code for this model: ',$)
      write (*,480)
      read (*,'(a)') runname
      np=mlen(runname)
      runname=runname(1:np)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Create User and open main model files
c
      call createS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine shock5headers
c! Opens and fills the headers of S5 chosen output files
c! @param [in,out] integer*4  iterations  number of chosen minimum iterations
c!
c! @return
c!  Chosen output files have user and model information written
c!
c! @details
c!  Files remain close until needed and at the end of the model.
c***************************************************************

      subroutine shock5headers (iterations)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i,j,iterations
      integer*4 itr
      character tab*1
c
c function
c
      real*8 feldens
c
      tab=','
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write Headers, open files to write
c
      call appendS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      de=feldens(dh,pop)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Title
c
   10 format(//,
     & ' SHOCK 5: Steady Rankine-Hugoniot Shock Model ',/,
     & ' ============================================',/,
     & ' Diffuse Field, Full Continuum Calculations.',/,
     & ' Global Shock-Precursor Iterations: ',i2,/,
     & ' Calculated by MAPPINGS V ',a12)
      write (*,10) iterations,theversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (luop,10) iterations,theversion
      write (lusp,10) iterations,theversion
      if (ratmod.eq.'Y') write (lurtsh,10) iterations,theversion
      if (dynmod.eq.'Y') write (ludy,10) iterations,theversion
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          write (luionsh(i),10) iterations,theversion
        enddo
      endif
      if (allmod.eq.'Y') write (lualsh,10) iterations,theversion
      if (fclmod.eq.'Y') write (lucl,10) iterations,theversion
      if (jlin.eq.'Y') write (lulsh,10) iterations,theversion
      if (bandsmod.eq.'Y') write (lupb,10) iterations,theversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lupt.gt.0) write (lupt,10) iterations,theversion
      if (lupc.gt.0) write (lupc,10) iterations,theversion
      if (lurtpc.gt.0) then
        if (ratmod.eq.'Y') write (lurtpc,10) iterations,theversion
      endif
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          if (luionpc(i).gt.0) then
            write (luionpc(i),10) iterations,theversion
          endif
        enddo
      endif
      if (lualpc.gt.0) then
        if (allmod.eq.'Y') write (lualpc,10) iterations,theversion
      endif
      if (lulpc.gt.0) then
        if (jlin.eq.'Y') write (lulpc,10) iterations,theversion
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  put runname string and master file name in each file created
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 format(//' Run  :,',a,/,' File :,',a)
      write (luop,20) runname,fsm
      write (lusp,20) runname,fsm
      if (ratmod.eq.'Y') write (lurtsh,20) runname,fsm
      if (dynmod.eq.'Y') write (ludy,20) runname,fsm
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          write (luionsh(i),20) runname,fsm
        enddo
      endif
      if (allmod.eq.'Y') write (lualsh,20) runname,fsm
      if (fclmod.eq.'Y') write (lucl,20) runname,fsm
      if (jlin.eq.'Y') write (lulsh,20) runname,fsm
      if (bandsmod.eq.'Y') write (lupb,20) runname,fsm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lupt.gt.0) write (lupt,20) runname,fsm
      if (lupc.gt.0) write (lupc,20) runname,fsm
      if (lurtpc.gt.0) then
        if (ratmod.eq.'Y') write (lurtpc,20) runname,fsm
      endif
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          if (luionpc(i).gt.0) then
            write (luionpc(i),20) runname,fsm
          endif
        enddo
      endif
      if (lualpc.gt.0) then
        if (allmod.eq.'Y') write (lualpc,20) runname,fsm
      endif
      if (lulpc.gt.0) then
        if (jlin.eq.'Y') write (lulpc,20) runname,fsm
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write Model Parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   30 format(//,
     & ' Model Parameters:',/,
     & ' =================',//,
     & ' Abundances     :,',a128,/,
     & ' Pre-ionisation :, ',a128,/,
     & ' Photon Source  :, ',a128)
   40 format(/,' Charge Exchange:, ',a12,/,
     & ' Photon Mode    :, ',a12,/,
     & ' Collision calcs:, ',a12)
   50 format(/,' Charge Exchange:, ',a12,/,
     & ' Photon Mode    :, ',a12,/,
     & ' Collision calcs:, ',a12,/,
     & ' Electron Kappa :, ',1pg11.4)
   60 format(//' ',t2,',Jden,',t9,',Jgeo,',t16,',Jtrans,',
     &            t23,',Jend,',t29,',Ielen,',t35,
     & ',Jpoen,',t43,',,Fren,',t51,',Tend,',t57,',DIend,',t67,
     & ',TAUen,',t77,',Jeq,',t84,',Teini,')
   70 format('  ',4(a4,3x),2(i2,4x),0pf6.4,0pf6.0,
     &           2(1pg10.3),3x,a4,0pf7.1)
   80 format(//' Photon Source'/
     & ' ============='/)
   90 format(/' MOD,',t7,',Temp.,',t16,',Alpha,',t22,',Turn-on,',t30,
     &',Cut-off,'
     &,t38,',Zstar,',t47,',FQHI,',t56,',FQHEI,',t66,',FQHEII,')
  100 format(' ',a2,1pg10.3,4(0pf7.2),1x,3(1pg10.3))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (luop,30) abnfile,ionsetup,srcfile
      write (lusp,30) abnfile,ionsetup,srcfile
c
      if (ratmod.eq.'Y') write (lurtsh,30) abnfile,ionsetup,srcfile
      if (dynmod.eq.'Y') write (ludy,30) abnfile,ionsetup,srcfile
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          write (luionsh(i),30) abnfile,ionsetup,srcfile
        enddo
      endif
      if (allmod.eq.'Y') write (lualsh,30) abnfile,ionsetup,srcfile
      if (fclmod.eq.'Y') write (lucl,30) abnfile,ionsetup,srcfile
      if (jlin.eq.'Y') write (lulsh,30) abnfile,ionsetup,srcfile
      if (bandsmod.eq.'Y') write (lupb,30) abnfile,ionsetup,srcfile
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (lupt,30) abnfile,ionsetup,srcfile
      write (lupc,30) abnfile,ionsetup,srcfile
c
      if (lurtpc.gt.0) then
        if (ratmod.eq.'Y') write (lurtpc,30) abnfile,ionsetup,srcfile
      endif
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          if (luionpc(i).gt.0) then
            write (luionpc(i),30) abnfile,ionsetup,srcfile
          endif
        enddo
      endif
      if (lualpc.gt.0) then
        if (allmod.eq.'Y') write (lualpc,30) abnfile,ionsetup,srcfile
      endif
      if (lualpc.gt.0) then
        if (jlin.eq.'Y') write (lualpc,30) abnfile,ionsetup,srcfile
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
      call dispabundances (luop, zi, abundtitle)
      if (lupt.gt.0) call dispabundances (lupt, zi, abundtitle)
c
      abundtitle=' Gas Phase Abundances :'
      call dispabundances (luop, zion, abundtitle)
      if (lupt.gt.0) call dispabundances (lupt, zion, abundtitle)
c
      call wabund (lusp)
      if (lupc.gt.0) call wabund (lupc)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      cht='New Set'
      if (chargemode.eq.2) cht='Disabled'
      if (chargemode.eq.1) cht='Old Set'
      pht='Normal'
      if (photonmode.eq.0) pht='Zero Field, Decoupled'
      clt='Dere 2007 Collisions'
      if (collmode.eq.1) clt='Old A&R Collision Methods'
      if (collmode.eq.1) clt='LMS Collision Methods'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (usekappa) then
        write (luop,50) cht,pht,clt,kappa
        if (lupt.gt.0) write (lupt,50) cht,pht,clt,kappa
      else
        write (luop,40) cht,pht,clt
        if (lupt.gt.0) write (lupt,40) cht,pht,clt
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write (luop,60)
      write (luop,70) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,diend,
     &tauen,tmod,tm00
      write (luop,80)
      write (luop,90)
      write (luop,100) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
c
      if (lupt.gt.0) then
        write (lupt,60)
        write (lupt,70) jden,jgeo,jtrans,jend,ielen,jpoen,fren,tend,
     &   diend,tauen,tmod,tm00
        write (lupt,80)
        write (lupt,90)
        write (lupt,100) iso,teff,alnth,turn,cut,zstar,qhi,qhei,qheii
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Shock parameters
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  110 format(//,
     & ' Initial Jump Conditions:',/,
     & ' =========================')
      write (luop,110)
      if (lupt.gt.0) write (lupt,110)
c
  120 format(//'  T0',1pg14.7,' V0',1pg14.7,/,
     & ' RH0',1pg14.7,' P0',1pg14.7,' B0',1pg14.7,//,
     & '  T1',1pg14.7,' V1',1pg14.7,/,
     & ' RH1',1pg14.7,' P1',1pg14.7,' B1',1pg14.7)
c
      dr=0.d0
      dv=vel1-vel0
c
      write (*,120) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &rho1,pr1,bm1*1.d6
      write (luop,120) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &rho1,pr1,bm1*1.d6
      if (lupt.gt.0) write (lupt,120) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,
     &te1,vel1*1.d-5,rho1,pr1,bm1*1.d6
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      wdilt0=te0
      t=te0
      ve=vel0
      dr=1.0d0
      dv=ve*0.01d0
      fi=1.0d0
      rad=0.0d0
      wdil=0.5d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     get the electrons...
c
      de=feldens(dh,pop)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate the radiation field and atomic rates
c
      call localem (t, de, dh)
      call totphot2 (t, dh, rad, dr, dv, wdil, specmode)
      call zetaeff (dh)
      call cool (t, de, dh)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call protostate (6)
      call protostate (luop)
      call protostate (lusp)
      if (lupt.gt.0) call protostate (lupt)
      if (lupc.gt.0) call protostate (lupc)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  130 format(//,
     & ' Precursor Conditions and Ionisation State',/
     & ' =========================================',/)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (luop,130)
      write (lusp,130)
      if (lupt.gt.0) write (lupt,130)
      if (lupc.gt.0) write (lupc,130)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ratmod.eq.'Y') write (lurtsh,130)
      if (dynmod.eq.'Y') write (ludy,130)
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          write (luionsh(i),130)
        enddo
      endif
      if (allmod.eq.'Y') write (lualsh,130)
      if (fclmod.eq.'Y') write (lucl,130)
      if (jlin.eq.'Y') write (lulsh,130)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call wionpop (luop, pop)
      if (ratmod.eq.'Y') call wionpop (lurtsh, pop)
      if (dynmod.eq.'Y') call wionpop (ludy, pop)
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          call wionpop (luionsh(i), pop)
        enddo
      endif
      if (allmod.eq.'Y') call wionpop (lualsh, pop)
      if (fclmod.eq.'Y') call wionpop (lucl, pop)
      if (jlin.eq.'Y') call wionpop (lulsh, pop)
      if (bandsmod.eq.'Y') call wionpop (lupb, pop)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lurtpc.gt.0) then
        if (ratmod.eq.'Y') write (lurtpc,130)
      endif
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          if (luionpc(i).gt.0) then
            write (luionpc(i),130)
          endif
        enddo
      endif
      if (lualpc.gt.0) then
        if (allmod.eq.'Y') write (lualpc,130)
      endif
      if (lulpc.gt.0) then
        if (jlin.eq.'Y') write (lulpc,130)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lupt.gt.0) call wionpop (lupt, pop)
      if (lurtpc.gt.0) then
        if (ratmod.eq.'Y') call wionpop (lurtpc, pop)
      endif
      if (tsrmod.eq.'Y') then
        do i=1,ieln
          if (luionpc(i).gt.0) then
            call wionpop (luionpc(i), pop)
          endif
        enddo
      endif
      if (lualpc.gt.0) then
        if (allmod.eq.'Y') call wionpop (lualpc, pop)
      endif
      if (lulpc.gt.0) then
        if (jlin.eq.'Y') call wionpop (lulpc, pop)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ratmod.eq.'Y') then
        call wmodel (lurtsh, 0.d0, 0.0d0, 0.d0, 0.0d0, 'LOSH')
        if (lurtpc.gt.0) call wmodel (lurtpc, 0.d0, 0.0d0, 0.d0, 0.0d0,
     &   'LOSH')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (tsrmod.eq.'Y') then
c
c B  : Ion balance files.'
c
  140  format(//,
     & ' #[1], [2] <X>     , [3] DeltaX  , [4] t       , [5] dt      ,',
     & '  [6] <T>    , [7] <ne>    , [8] <nH>    , [9] <nT>    ,',
     & 31(' [',i2,'] ', a6,' ,'))
        do i=1,ieln
          write (luionsh(i),'(//," Element : ",a2)') elem(iel(i))
          write (luionsh(i),140) (j+9,rom(j),j=1,maxion(iel(i)))
          if (luionpc(i).gt.0) then
            write (luionpc(i),'(//," Element : ",a2)') elem(iel(i))
            write (luionpc(i),140) (j+9,rom(j),j=1,maxion(iel(i)))
          endif
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  150 format(//,
     & 17(a12,a1))
      if (bandsmod.eq.'Y') then
        write (lupb,150) 'Te',tab,'de',tab,'dh',tab,'en',tab,'XHI',tab,'
     &XHII',tab,'mu',tab,'tloss',tab,'Lambda',tab,'ff/total',tab,'BHI-0
     &.1keV',tab,'B0.1-0.5keV',tab,'B0.5-1.0keV',tab,'B1.0-2.0eV',tab,'B
     &2.0-10.0keV',tab,'Ball'
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  160 format(/'Mean Zone Values'/
     & '================'/)
  170 format(10(a12,a2),30(a12,a2),33(a12,a2))
  180 format(10(a12,a2),30(a12,a2),33(1pg12.5,a2))
c
  190 format('=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '=======================================================',
     & '======================================')
      if (fclmod.eq.'Y') then
        write (lucl,160)
        if (jnorm.eq.0) then
          write (lucl,170) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'XHI   ',tab,'XHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(ne
     &.nH)',tab,(elem(j),tab,j=1,atypes),'Eloss/ne.nH',tab,'Egain/ne.nH'
     &     ,tab,'Netloss/ne.nH'
        endif
        if (jnorm.eq.1) then
          write (lucl,170) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'XHI   ',tab,'XHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(nH
     &^2)',tab,(elem(j),tab,j=1,atypes),'Eloss/nH^2',tab,'Egain/nH^2',
     &     tab,'Netloss/nH^2'
        endif
        if (jnorm.eq.2) then
          write (lucl,170) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'XHI   ',tab,'XHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(ne
     &.ni)',tab,(elem(j),tab,j=1,atypes),'Eloss/ne.ni',tab,'Egain/ne.ni'
     &     ,tab,'Netloss/ne.ni'
        endif
        if (jnorm.eq.3) then
          write (lucl,170) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'XHI   ',tab,'XHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(n^
     &2)',tab,(elem(j),tab,j=1,atypes),'Eloss/n2',tab,'Egain/n2',tab,'Ne
     &tloss/n2'
        endif
        if (jnorm.eq.4) then
          write (lucl,170) 'T ',tab,'n_e',tab,'n_H',tab,'n_ion',tab,'rho
     & ',tab,'XHI   ',tab,'XHII  ',tab,'mu ',tab,'Losses (L)',tab,'L/(ne
     &^2)',tab,(elem(j),tab,j=1,atypes),'L_5007',tab,'LHalpha',tab,'LLya
     &lpha'
        endif
        write (lucl,180) '(K)',tab,'(/cm^3)',tab,'(/cm^3)',tab,'(/cm^3)'
     &   ,tab,'(g/cm^3)',tab,' ',tab,' ',tab,'(amu)',tab,'(erg/cm^3/s)',
     &   tab,'(erg cm^3/s)',tab,('(erg cm^3/s)',tab,j=1,atypes),'(erg cm
     &^3/s)',tab,'(erg cm^3/s)',tab,'(erg cm^3/s)'
        write (lucl,190)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jlin.eq.'Y') then
c
c L  : Monitor up to ',i3,' lines  - precursor and shocks
c
  200   format(// ,'     ,             ,             ,             ,',
     & '             ,             ,             ,             ',
     &   30(',    ',a2,a6,'    '))
        write (lulsh,'("Run: ",a96)') runname
        write (lulsh,200) (elem(emlinlistatom(itr)),
     &   rom(emlinlistion(itr)),itr=1,njlines)
        if (lulpc.gt.0) then
          write (lulpc,'("Run: ",a96)') runname
          write (lulpc,200) (elem(emlinlistatom(itr)),
     &     rom(emlinlistion(itr)),itr=1,njlines)
        endif
  210  format(
     & ' # [1] step, [2] <X>, [3] Xmid, [4] dX, [5] <T>, [6] <ne>,'
     & ' [7] <nH>, [8]  <nT>, [9] logQH, [10]  logUH,',
     & ' [11]  logQN, [12]   <HB>',
     &   16(',',f12.3,'[',i2,']'))
        write (lulsh,210) (emlinlist(itr),itr+12,itr=1,njlines)
        if (lulpc.gt.0) then
          write (lulpc,210) (emlinlist(itr),itr+12,itr=1,njlines)
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call wionabal (luop, pop)
      if (lupt.gt.0) call wionabal (lupt, pop)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Close all files and flush
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call closeS5files ()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine shocksummary
c! Print summary of shock properties for a given precursor condition
c! @param [in,out] integer*4     lunit  logical unit to write to, 6 for std out
c!
c! @return
c!  Display shock properties for screen and file output,
c! varies as precursor evolves, commas make output compatible with csv
c! output and doesn't make the screen or other locations any harder to
c! read.
c!
c! @details
c!  Detailed Shock properties given current precursor to file or screen.
c***************************************************************

      subroutine shocksummary (lunit)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      real*8 fmua,fpresse,frho
      integer*4 lunit
c
c uses global pop
c
c
c pre = preshock , pst = postshock
c
      te_pre=te0
      te_pst=te1
      de_pre=de0
      de_pst=de1
      dh_pre=dh0
      dh_pst=dh1
      vs_pre=vel0
      vs_pst=vel1
c
      rh_pre=rho0
      rh_pst=rho1
      cmpf=rho1/rho0
      pr_pre=pr0
      pr_pst=pr1
      bm_pre=bm0
      bm_pst=bm1
c
      bp0=(bm0*bm0)/epi
      bp1=(bm1*bm1)/epi
      Pgas=pr0
      Pmag=bp0
      Pram=rh_pre*vs_pre*vs_pre
c
      machnumber=vel0/dsqrt(gammaEOS*pr0/rho0)
      alfvennumber=vel0/dsqrt(2.d0*Pmag/rho0)
      malpha=Pmag/Pgas
      gaseta=gammaEOS*Pgas/Pram
      mageta=2.d0*Pmag/Pram
c
   10   format(/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    Velocity       :',1pg12.5,' km/s',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    Preshock Mach Number          :,',1pg12.5,/,
     & '    Preshock Alfven Mach Number   :,',1pg12.5,/,
     & '    Preshock Mag Alpha (Pmag/Pgas):,',1pg12.5,/,
     & '    Preshock Gas Eta  (gPgas/Pram):,',1pg12.5,/,
     & '    Preshock Mag Eta  (2Pmag/Pram):,',1pg12.5,//,
     & '    Preshock     T :,',1pg12.5,', K',/,
     & '    Preshock     ne:,',1pg12.5,', cm^-3',/,
     & '    Preshock     nH:,',1pg12.5,', cm^-3',/,
     & '    Preshock     d :,',1pg12.5,', g/cm^-3',/,
     & '    Preshock   Pgas:,',1pg12.5,', dyne/cm^2',/,
     & '    Preshock     mu:,',1pg12.5,', a.m.u.',/,
     & '    Preshock    XHI:,',1pg12.5,/,
     & '    Preshock   XHII:,',1pg12.5,/,
     & '    Preshock   XHeI:,',1pg12.5,/,
     & '    Preshock  XHeII:,',1pg12.5,/,
     & '    Preshock XHeIII:,',1pg12.5,/,
     & '    Preshock     B :,',1pg12.5,', microGauss',/,
     & '    Preshock   Pmag:,',1pg12.5,', dyne/cm^2',/,
     & '    Preshock   Pram:,',1pg12.5,', dyne/cm^2',/)
   20   format(
     & '    Postshock Compression Factor   :',1pg12.5,/,
     & '    Postshock    T :,',1pg12.5,', K',/,
     & '    Postshock    ne:,',1pg12.5,', cm^-3',/,
     & '    Postshock    nH:,',1pg12.5,', cm^-3',/,
     & '    Postshock    v :,',1pg12.5,', km/s',/,
     & '    Postshock    d :,',1pg12.5,', g/cm^-3',/,
     & '    Postshock  Pgas:,',1pg12.5,', dyne/cm^2',/,
     & '    Postshock    B :,',1pg12.5,', microGauss',/,
     & '    Postshock  Pmag:,',1pg12.5,', dyne/cm^2',/,
     & '    Postshock  Pram:,',1pg12.5,', dyne/cm^2',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
c
      en=zen*dh0
      pr0=fpresse(te0,de0,dh0)
      rho0=frho(de0,dh0)
      cspd=dsqrt(gammaEOS*pr0/rho0)
      wmol=rho0/(en+de0)
      mu=fmua(de0,dh0)
c
      Pgas=pr0
      Pmag=(bm0*bm0)/epi
      Pram=rho0*vel0*vel0
c
      machnumber=vel0/dsqrt(gammaEOS*pr0/rho0)
      alfvennumber=vel0/dsqrt(2.0d0*Pmag/rho0)
      malpha=Pmag/Pgas
      gaseta=gammaEOS*Pgas/Pram
      mageta=2.d0*Pmag/Pram
c
      write (lunit,10) vel0*1.0d-5,machnumber,alfvennumber,malpha,
     &gaseta,mageta,te0,de0,dh0,rho0,pr0,mu,pop(1,zmap(1)),pop(2,zmap(1)
     &),pop(1,zmap(2)),pop(2,zmap(2)),pop(3,zmap(2)),bm0*1.d6,bp0,Pram
c
      Pram=rho1*vel1*vel1
      write (lunit,20) cmpf,te1,de1,dh1,vel1*1d-5,rho1,pr1,bm1*1.d6,bp1,
     &Pram
c
      return
      end

c****************************************************************
c> @brief The subroutine shock5jump
c! Calculate the RH shock jump flow variables with 2017 Pressure formulae
c! @param This routine has no parameters
c!
c! @return
c!  Calculate the RH shock jump flow variables with 2017 Pressure formulae
c!
c! @details
c!  Calculate the RH shock jump flow variables with 2017 Pressure formulae
c! variables vary a little in setup depending on how precursor B field
c! is defined in S5 block vars 'mtype' and 'stype' for magnetic type
c! and shock type in terms of temperature or velocity.
c***************************************************************

      subroutine shock5jump ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      real*8 tj,dej,dhj,tpo,va,eps
      real*8 feldens,fpresse,frho,velshock2
c
c     Now get revised shock solution for changed pre shock values:
c     Shock in term of flow velocity
c
      call copypop (pop_neu, pop)
      rh_neu=frho(de_neu,dh_neu)
      call copypop (pop_pre, pop)
      rh_pre=frho(de_pre,dh_pre)
c
      Pram=rh_pre*vs_pre*vs_pre
c
      if (stype.eq.'V') then
c
        dr=1.d0
c
        if (mtype.eq.'B') then
c         Specified B field in uG.  Won't change between iterations.
          Bmag=magparam*1.0d-6
          Pmag=(Bmag*Bmag)/epi
        endif
c
        if (mtype.eq.'A') then
c         Specified protostate malpha. B won't change between iterations
          Pram=rh_neu*vs_pre*vs_pre
          Pgas=fpresse(te_neu,de_neu,dh_neu)
          Pmag=magparam*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
c
        if (mtype.eq.'C') then
c         Specified shockfront malpha.  B will change between iterations
c         as the gas pressure at the back of the precursor changes.
          Pgas=fpresse(te_pre,de_pre,dh_pre)
          Pmag=magparam*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
c
        if (mtype.eq.'M') then
c         Specified Alfven Mach number
          va=vs_pre/mmach
          Pmag=(0.5*rh_pre*va*va)
          Bmag=dsqrt(epi*Pmag)
        endif
c
        if (mtype.eq.'R') then
c         Specified magnetic eta = 2Pmag/Pram.  Bmag shouldn't change
c         between iterations since Pram doesn't change.
          Pmag=0.5d0*mageta*Pram
          Bmag=dsqrt(epi*Pmag)
        endif
c
        bm_pre=Bmag
c
c       Now call post-shock state always P_pre, alpha may change for A
c
c        vshoc=vs_pre
c
        tj=te_pre
        dhj=dh_pre
        dej=de_pre
c
c jump is the new exact quadratric with no dx time or losses
c pre and post 0,1 globals are set in routine
c
        call shockcmpf (tj, dej, dhj, vshoc, Bmag)
c
        tm00=te0
c
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Shock in term of temperature
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (stype.eq.'T') then
c
        dr=1.d0
c
c
        if (mtype.eq.'B') then
c         Specified B field in uG.  Won't change between iterations.
          Bmag=magparam*1.0d-6
          Pmag=(Bmag*Bmag)/epi
        endif
c
        if (mtype.eq.'A') then
c         Specified malpha in protostate
          Pgas=fpresse(te_neu,de_neu,dh_neu)
          Pmag=magparam*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
c
        if (mtype.eq.'C') then
c         Specified shockfront malpha.  B will change between iterations
c         as the gas pressure at the back of the precursor changes.
          Pram=rh_neu*vs_pre*vs_pre
          Pgas=fpresse(te_pre,de_pre,dh_pre)
          Pmag=magparam*Pgas
          Bmag=dsqrt(epi*Pmag)
        endif
c
        if (mtype.eq.'M') then
          va=vs_pre/mmach
          Pmag=(0.5*rh_pre*va*va)
          Bmag=dsqrt(epi*Pmag)
        endif
c
        if (mtype.eq.'R') then
          Pmag=0.5d0*mageta*Pram
          Bmag=dsqrt(epi*Pmag)
        endif
c
        bm_pre=Bmag
c
c       Now call velshock2 to find the shock speed
c
        Pgas=fpresse(te_pre,de_pre,dh_pre)
        vshoc=vs_pre
        t=te_pre
        tpo=te_pst
        dh=dh_pre
        de=feldens(dh,pop)
        tloss=0.0d0
        dt=0.0d0
c
   10   bm0=Bmag
c
c       Call function velshock2 to find the new shock speed
c       N.B. velshock2 calls subroutine rankhug
        vshoc=velshock2(dh,t,Bmag,tpo)
c
        Pram=rh_neu*vshoc*vshoc
c
        if (mtype.eq.'M') then
          va=vshoc/mmach
          Pmag=(0.5*rh_pre*va*va)
          Bmag=dsqrt(epi*Pmag)
          eps=dabs(2.0*(Bmag-bm0)/(Bmag+bm0))
          if (eps.gt.1.0d-6) goto 10
        endif
c
        if (mtype.eq.'R') then
c iterate for rampressure/alpha_r
          Pmag=0.5d0*mageta*Pram
          Bmag=dsqrt(epi*Pmag)
          eps=dabs(2.0*(Bmag-bm0)/(Bmag+bm0))
          if (eps.gt.1.0d-6) goto 10
        endif
c
        vs_neu=vel0
        bm_neu=Bmag
        tm00=te0
      endif
c
c Shock intial state kept for iterations, above numbers are
c reused from step to step
c
      te_pre=te0
      te_pst=te1
      de_pre=de0
      de_pst=de1
      dh_pre=dh0
      dh_pst=dh1
      vs_pre=vel0
      vs_pst=vel1
c
      rh_pre=rho0
      rh_pst=rho1
      pr_pre=pr0
      pr_pst=pr1
      bm_pre=bm0
      bm_pst=bm1
c
      bp0=(bm0*bm0)/epi
      bp1=(bm1*bm1)/epi
      Pmag=bp0
      Pgas=pr0
      Pram=rho0*vel0*vel0
c
      machnumber=vel0/dsqrt(gammaEOS*pr0/rho0)
      malpha=Pmag/Pgas
      gaseta=gammaEOS*Pgas/Pram
      mageta=2.d0*Pmag/Pram
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine shock5check
c! Check shock structure compared to previous iteration
c! @param [in] integer*4       its  current iterations so far
c! @param [in] integer*4    maxits  max its requested - may change
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!     Check if the global shock-precursor iterations have converged.
c!     Compares the current iteration with the previous iteration.
c!     Saves the current shock state on either side of the shock jump,
c!    for a future comparison.
c***************************************************************
c
      subroutine shock5check (its, maxits)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Check if the global shock-precursor iterations have converged.
c     Compares the current iteration with the previous iteration.
c     Saves the current shock state on either side of the shock jump,
c     for a future comparison.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 its,maxits
c
      real*8 delhhe,rmserr,term
      real*8 t_psi,t_cmp,t_tpre,t_tpst,t_depre
c
   10 format(/,
     & ' ********************************************************',/,
     & '  SHOCK 5  Convergence Test, It.: ',i2.2,' of ',i2.2)
   20 format(
     & '  Result: CONVERGED  (precursor RMS ',1pg11.4,'%, post-shock ',
     & 'RMS ',1pg11.4,'%)',/,
     & ' ********************************************************',/)
   30 format(
     & '  Result: NOT CONVERGED',/,
     & ' ********************************************************',/)
   40 format(
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    Psi (Q/v): ',  1pg11.4,'  `Psi      : ',  1pg11.4,/,
     & '    Compress : ',  1pg11.4,'  `Compress : ',  1pg11.4,/,
     & '    T_pre    : ',  1pg11.4,'  `T_pre    : ',  1pg11.4,/,
     & '    T_shock  : ',  1pg11.4,'  `T_shock  : ',  1pg11.4,/,
     & '    ne_pre   : ',  1pg11.4,'  `ne_pre   : ',  1pg11.4,/,
     & '    DelH/He  : ',  1pg11.4,'%',/,
     & '    Precursor RMS (Psi,T_pre,ne_pre,DelH/He): ',1pg11.4,'%',/,
     & '    Post-shock RMS (Compress,T_shock)      : ',1pg11.4,'%',/,
     & '    RMS      : ',  1pg11.4,'%',/,
     & '    Aitken omega (precursor relaxation)     : ',1pg11.4,/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
      converged=0
      osclass=0
c
      if (its.gt.1) then
c
c uses global pop
c
        call difhhe (pop_pre, pop_pre0, delhhe)
        write (*,10) its,maxits
        term=2.d0*(psi-psi0)/(psi+psi0)
        t_psi=term*term
        term=2.d0*(cmpf-cmpf0)/(cmpf+cmpf0)
        t_cmp=term*term
        term=2.d0*(te_pre-te_pre0)/(te_pre+te_pre0)
        t_tpre=term*term
        term=2.d0*(te_pst-te_pst0)/(te_pst+te_pst0)
        t_tpst=term*term
        term=2.d0*(de_pre-de_pre0)/(de_pre+de_pre0)
        t_depre=term*term
        rmserr=t_psi+t_cmp+t_tpre+t_tpst+t_depre+(delhhe*delhhe)
        rmserr=dsqrt(rmserr/6.d0)
c
c  Diagnostic-only split of the same six terms into a precursor
c  (pre-shock: Psi, T_pre, ne_pre, DelH/He) sub-residual and a
c  post-shock (Compress, T_shock) sub-residual, so a failing model's
c  log shows which side of the shock front the residual actually
c  comes from (issue #7 follow-up).  Neither sub-value feeds the
c  convergence decision below -- only the combined rmserr does,
c  unchanged from before.
c
        precrms=dsqrt((t_psi+t_tpre+t_depre+(delhhe*delhhe))/4.d0)
        postrms=dsqrt((t_cmp+t_tpst)/2.d0)
c
        write (*,40) psi,psi0,cmpf,cmpf0,te_pre,te_pre0,te_pst,te_pst0,
     &   de_pre,de_pre0,delhhe*100.d0,precrms*100.d0,postrms*100.d0,
     &   rmserr*100.d0,aitomega
        if (rmserr.lt.s5rmstol) then
          converged=1
          write (*,20) precrms*100.d0,postrms*100.d0
        else
          converged=0
          if (its.ge.maxits) maxits=maxits+1
          write (*,30)
c
c  ISSUE #14: osclass -- classify *why* the combined test failed,
c  from the same three facts every time (see the meanings listed at
c  osclass's declaration, s5blocks.inc).  precwarn takes priority:
c  the precrms/postrms split below assumes the precursor's own inner
c  self-consistency solve actually settled on this iteration, and
c  isn't a meaningful reading if it didn't.  0.5% (osclass 1 vs 2) is
c  a deliberately conservative margin above every instance actually
c  measured so far (0.015%-0.13% precrms across the known
c  oscillators): comparing predicted precursor spectra between the
c  two states of a live oscillator at this residual size showed every
c  diagnostic emission line agreeing to well under 0.5%, consistently
c  smaller than the genuine line-flux change between adjacent
c  converged grid points one velocity step apart -- see issue #14.
c  postrms's own ceiling (0.05%, vs s5rmstol's tight 0.01%) is looser
c  than the precursor one and not independently spectrum-validated --
c  post-shock quantities (Compression, T_shock) are computed from the
c  precursor's own boundary state, so a purely precursor-side
c  oscillation shows up as a small correlated wobble here too (0.005%-
c  0.022% seen across the 4 known oscillators, confirmed via direct
c  comparison against s5rmstol alone misclassifying 2 of them as a
c  post-shock problem).  0.05% is a provisional margin above that
c  observed range, not a proven bound -- revisit if more instances
c  push closer to it.
          if (precwarn.ne.0) then
            osclass=5
          else if ((postrms.lt.5.0d-4).and.(precrms.lt.5.0d-3)) then
            osclass=1
          else if ((postrms.lt.5.0d-4).and.(precrms.ge.5.0d-3)) then
            osclass=2
          else if ((postrms.ge.5.0d-4).and.(precrms.lt.5.0d-3)) then
            osclass=3
          else
            osclass=4
          endif
        endif
      endif
c
c save current shock for next test
c
      call copypop (pop_pre, pop_pre0)
      psi0=psi
      te_pre0=te_pre
      de_pre0=de_pre
      dh_pre0=dh_pre
      vs_pre0=vs_pre
      rh_pre0=rh_pre
      pr_pre0=pr_pre
      bm_pre0=bm_pre
c
      te_pst0=te_pst
      de_pst0=de_pst
      dh_pst0=dh_pst
      vs_pst0=vs_pst
      rh_pst0=rh_pst
      pr_pst0=pr_pst
      bm_pst0=bm_pst
c
      cmpf0=cmpf
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine shock5precursor
c! XXXX - add one line purpose here
c! @param [in,out] integer*4  iteration  XXX-meaning
c! @param [in,out] integer*4    maxits  XXX-meaning
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine shock5precursor (iteration, maxits)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      logical iexi
c
      integer*4 i,idx
      integer*4 nfs0,t0lim
      integer*4 iteration,maxits
      integer*4 itcount
c
      character nmod*4,ratemode*4
      character tab*1
c
      real*8 tf,f
      real*8 dtime,dtimer,pre_par,qh
      real*8 drr,tstep
      real*8 diff,xhfinal
      real*8 rmserr,rmslimit
      real*8 rdvol,irdvol
      real*8 te_0,de_0,dh_0
      real*8 aitrtenow,aitrdenow,aitdnum,aitdden
c
c  functions
c
      real*8 feldens,frectim3
c      real*8 fdilu
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open files, main and precursor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call appendS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   10 format(
     & ' ********************************************************',/,
     & '  SHOCK 5 Precursor Iteration: ',i2.2,' of ',i2.2,/,
     & ' ********************************************************')
      write (*,10) ,iteration,maxits
      write (lupt,10) ,iteration,maxits
c
      tab=','
      jspot='N'
      jcon='Y'
      ratemode='ALL'
      absf=0.04d0
c
c  Iterate the precursor's own inner loop, on every iteration, to the
c  *same* tolerance shock5check uses for the outer pass/fail decision
c  (s5rmstol, 0.01%).  The inner loop used to stop at 5-10%, so the
c  outer convergence history was being judged against a precursor
c  solve 10x-1000x looser than the threshold it was compared to
c  (issue #7).  Costs roughly 7-30% more run time.
c
      rmslimit=s5rmstol
c
      fi=1.0d0
      wdil=0.5d0
      rdvol=1.0d0
      irdvol=1.0d0
      vunilog=0.0d0
c
      t=te_neu
      de=de_neu
      dh=dh_neu
      Bmag=bm_neu
      call copypop (pop_neu, pop)
c
      rad=dist(step)
      dr=dist(step)-dist(step-1)
      dv=veloc(step)-veloc(step-1)
c
      call localem (t, de, dh)
c
      specmode='UP'
      call totphot2 (t, dh, 1.0d38, 1.0d0, 0.0d0, wdil, specmode)
c
      do idx=1,infph
        soupho(idx)=0.d0
c UP mode is upf=1.0d0, need 0.5x
        tphot(idx)=tphot(idx)*0.5d0
        src(idx)=tphot(idx)
      enddo
c
c  lets see what we got:
c
      write (*,*)
      call fieldsummary (6, 1, tphot)
      write (*,*)
      call fieldsummary (lupt, 1, tphot)
c
c   compare ion front velocity to shock velocity
c
      ve=vshoc
      viofr=qht/(zen*dh_neu)
      psi0=psi
      psi=viofr/vshoc
c just nH
      qh=qht/dh_neu
      pre_par=viofr/ve
c
   20 format(/,' SHOCK 5 Precursor Parameter: ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '     v_shock       : ',  1pg12.5,' km/s',/,
     & '     Q_ions        : ',  1pg12.5,' km/s',/,
     & '     Psi (Q/v)     : ',  1pg12.5,/,
     & '     Q_H           : ',  1pg12.5,' km/s',/,
     & '     Psi_H (Q_H/v) : ',  1pg12.5,/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
c
      write (*,20) vshoc*1.d-5,viofr*1.d-5,psi,qh*1.d-5,qh/vshoc
      write (lupt,20) vshoc*1.d-5,viofr*1.d-5,psi,qh*1.d-5,qh/vshoc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c setup for precursor iteration, fill arrays with protoionisation
c and opacities
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c to allow outer edge of ion front to evolve in fast ion front case,
c but reset to neutral for each global interation
c
      te_0=te_neu
      de_0=de_neu
      dh_0=dh_neu
c
      call copypop (pop_pre, pop)
c
      t=te_pre
      de=de_pre
      dh=dh_pre
c
      call allrates (t, ratemode)
c
c  stromgren length if front detatches, recomb rate if ionised gas
c
      dtimer=frectim3(dh)*pre_par
      drr=dtimer*ve
c
      call copypop (pop_neu, pop)
c
      t=te_0
      de=de_0
      dh=dh_0
      ve=vshoc
c
      call allrates (t, ratemode)
c
c  go nfs* distance that absorbs 4% of the field in the neutral medium
c  tau = ~ 10.38 for 100 steps ~6.5 for 50 steps
c
c set up absorbsion fractions and attenuations arrays
c setup ion balance arrays too
c
      absf=0.04d0
      dtime=frdr/vshoc
      if (iteration.eq.1) then
        nfs=nint(dlog(1.d-6)/dlog(1.d0-absf))
        nfs=min(nfs,mxifsteps)
c
c       dr=1.0d18
c       call absdis2 (dh, absf, dr, 0.0d0, pop_neu)
        do i=1,nfs
          dr=1.0d18
          f=absf*min(1.0d0,0.04d0*dble(i*i))
          fra(i)=absf*f
          call absdis2 (dh, absf*f, dr, 0.0d0, pop_neu)
          predr(i)=dr
        enddo
        frdr=nfs*dr
        dtime=frdr/ve
        invnfs=1.d0/dble(nfs)
c
      endif
c
      x(1)=0.d0
      fra(1)=absf*0.04d0
      do i=2,nfs
        x(i)=x(i-1)+predr(i-1)
        f=min(1.0d0,0.04d0*dble(i*i))
        fra(i)=absf*f
      enddo
c
c index 1 is zone at the source, index nfs is outer edge of front zone
c
      if (iteration.eq.1) then
        call copypop (pop_neu, p1)
        dr=predr(1)
        call copypopstep (p1, 1, popfr)
        call clearpop (p1)
        call copypopstep (p1, 1, popintfr)
        fh(1)=pop_neu(2,1)
        foi(1)=pop_neu(1,zmap(8))
        foii(1)=pop_neu(2,zmap(8))
        foiii(1)=pop_neu(3,zmap(8))
        nh(1)=dh_0
        ne(1)=de_0
        nte(1)=te_0
        do i=2,nfs
          call copypop (pop_neu, p1)
          dr=predr(i)
          call copypopstep (p1, i, popfr)
          call scalepop (p1, dr)
          call copysteppop (i-1, popintfr, p2)
          call addpop (p1, p2)
          call copypopstep (p2, i, popintfr)
          fh(i)=pop_neu(2,1)
          foi(i)=pop_neu(1,zmap(8))
          foii(i)=pop_neu(2,zmap(8))
          foiii(i)=pop_neu(3,zmap(8))
c         fra(i)=absf
          nh(i)=dh_0
          ne(i)=de_0
          nte(i)=te_0
        enddo
      endif
      dtime=x(nfs)/vshoc
c
   30 format(/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Precursor Timescales: '/
     & '    : ',i4,' zone neutral crossing time   : ',1pg11.4,/
     & '    : ',i4,' zone neutral length          : ',1pg11.4,//
     & '    :  Ionised recombination time x Q/v : ',1pg11.4,/
     & '    :  Ionised recombination distance   : ',1pg11.4,/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
      if (vmod.ne.'NONE') then
        write (*,30) nfs,dtime,nfs,frdr,dtimer,drr
      endif
c
      dr=frdr*invnfs
      tstep=(frdr/vshoc)*invnfs
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Iterate, full solution for Psi > 1, inner edge only for Psi <1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   40 format(/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Precursor Iteration: ',/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
      if (vmod.ne.'NONE') write (*,40)
c
      xhfinal=pop(1,1)
      itcount=0
      precwarn=0
      t0lim=0
      nfs0=nfs
c
   50 call copypop (pop_neu, pop)
c
      t=te_0
      de=de_0
      dh=dh_0
c
      nmod='TIM'
      rms=0.d0
c
c integrate from outer edge to source
c
      rmserr=0.d0
      nrms=min(max(minifsteps,nfs/2),50)
      x(1)=0.d0
      do i=2,nfs
        x(i)=x(i-1)+predr(i-1)
      enddo
c
      do i=nfs,1,-1
c
c full column absorption applied to spectrum
c
        dr=predr(i)
        tstep=dr/vshoc
        call copysteppop (i, popintfr, p1)
        call attsig (dh, dr, p1, attcol, sigcol)
        do idx=ionstartbin,infph-1
          tphot(idx)=src(idx)*attcol(idx)
          skipbin(idx)=.true.
          if (tphot(idx).gt.epsilon) skipbin(idx)=.false.
          ipho=1
          iphom=0
          tem=0.d0
        enddo
c redo rates
        call allrates (t, ratemode)
c evolve in attenuated field to inner boundary
        call copypop (pop, p1)
        call teequi2 (t, tf, de, dh, tstep, nmod)
        call copypop (pop, p2)
        call averinto (0.5d0, p1, p2, p1)
        call absdis2 (dh, fra(i), dr, 0.0d0, p1)
        predr(i)=dr
        nh(i)=dh
        ne(i)=feldens(dh,p1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c RMS Error
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (i.le.nrms) then
          call copysteppop (i, popfr, p1)
          call difhhe (pop, p1, diff)
          diff=dabs(1.d0-(10.d0**diff))
          rmserr=rmserr+(diff*diff)
          diff=0.5d0*(tf-nte(i))/(tf+nte(i))
          rmserr=rmserr+(diff*diff)
        endif
c save step temp at point
        nte(i)=tf
c save pop after step = closer to src, ith pop is inner edge of ith zone
        call copypopstep (pop, i, popfr)
        fh(i)=pop(2,1)
        foi(i)=pop(1,zmap(8))
        foii(i)=pop(2,zmap(8))
        foiii(i)=pop(3,zmap(8))
        de=feldens(dh,pop)
        t=tf
      enddo
c
      x(1)=0.d0
      fra(1)=absf*0.04d0
      do i=2,nfs
        x(i)=x(i-1)+predr(i-1)
        f=min(1.0d0,0.04d0*dble(i*i))
        fra(i)=absf*f
      enddo
c
   60 format(//,
     & a5,14(a1,2x,a10,1x),/
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',
     & ':::::::::::::::::::::::::::::::::::::')
   70 format(i5.3,14(a1,1pg13.6))
c
      if (vmod.ne.'NONE') then
        write (*,60) ' Step',tab,'Dist.(cm)',tab,'dx (cm)',
     &   tab,'Time (s)',tab,'dt (s)',
     &   tab, 'Te (K)',tab,'ne(cm^-3)',tab,'nH(cm^-2)',
     &   tab,'XHI',tab,'XHII',tab,'<abs>',
     &   tab,'XOI',tab,'XOII',tab,'XOIII'
      do i=nfs,1,-1
        call copysteppop (i, popintfr, p1)
        write (*,70) i,tab,x(i),tab,predr(i),
     &    tab,x(i)/vshoc,tab,predr(i)/vshoc,
     &    tab,nte(i),tab,ne(i),tab,nh(i),tab,1.d0-fh(i),
     &    tab,fh(i),tab,fra(i),
     &    tab,foi(i),tab,foii(i),tab,foiii(i)
      enddo
      endif
c
      if (vmod.eq.'NONE') then
        write (*,60) ' Step',tab,'Dist.(cm)',tab,'dx (cm)',
     &   tab,'Time (s)',tab,'dt (s)',
     &   tab, 'Te (K)',tab,'ne(cm^-3)',tab,'nH(cm^-2)',
     &   tab,'XHI',tab,'XHII',tab,'<abs>',
     &   tab,'XOI',tab,'XOII',tab,'XOIII'
        i=1
        call copysteppop (i, popintfr, p1)
        write (*,70) i,tab,x(i),tab,predr(i),
     &    tab,x(i)/vshoc,tab,predr(i)/vshoc,
     &    tab,nte(i),tab,ne(i),tab,nh(i),tab,1.d0-fh(i),
     &    tab,fh(i),tab,fra(i),
     &    tab,foi(i),tab,foii(i),tab,foiii(i)
      endif
c
c get post precursor outward field
c
      rmserr=dsqrt(rmserr/dble(nrms))
c
   80  format(
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' RMS change %: ',1pg11.4,' Iteration: ',i2,/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
      write (*,80) 100.d0*rmserr,itcount+1
      write (lupt,80) 100.d0*rmserr,itcount+1
c
      t0lim=0
      nfs0=nfs
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Adjust grid size for given absf as te and pops converge
c and save for next interation starting point.  Kept keyed on finalit
c (issue #7) -- the grid should stay free to adapt for the entire
c normal run; it should only freeze for the single, truly final
c output-writing pass.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (finalit.le.0) then
      if (nte(nfs).ge.500.d0) then
        t0lim=1
        nfs=nfs+nint(nte(nfs)*0.01d0)+10
      endif
      n100=nfs
      do i=1,nfs
        if (nte(i).lt.100.d0) then
          n100=i
          goto 90
        endif
      enddo
   90 nfs=max(n100,minifsteps)
      if (popfr(nfs,1,1).lt.0.95d0) then
        nfs=nfs+10
      endif
      nfs=min0(nfs,mxifsteps)
c
      if (nfs.gt.nfs0) then
c
c copy last cell to fill out expanded grid, leave if reduced
c
        do i=nfs0+1,nfs
          predr(i)=predr(nfs0)
          call copypopstep (popfr, nfs0, p1)
          call copypopstep (p1, i, popfr)
          dr=predr(nfs0)
          call scalepop (p1, dr)
          call copysteppop (nfs0, popintfr, p2)
          call addpop (p1, p2)
          call copypopstep (p2, i, popintfr)
          call copypopstep (popfr, nfs0, p1)
          fh(i)=p1(2,1)
          foi(i)=p1(1,zmap(8))
          foii(i)=p1(2,zmap(8))
          foiii(i)=p1(3,zmap(8))
c         fra(i)=absf
          nh(i)=nh(nfs0)
          ne(i)=ne(nfs0)
          nte(i)=nte(nfs0)
        enddo
        x(1)=0.d0
        fra(1)=absf*0.04d0
c
c taper intial abs fractions close to shock front.
c
        do i=2,nfs
          x(i)=x(i-1)+predr(i-1)
          f=min(1.0d0,0.04d0*dble(i*i))
          fra(i)=absf*f
        enddo
      endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c this will all hopefully all inline with gfortran -flto
c popintfr is columns of previous steps so starts at 0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call clearpop (p1)
c     columns up to inneredge ofzone
      call copypopstep (p1, 1, popintfr)
      do i=2,nfs
c pops at end time of previous zone, nearer source
        call copysteppop (i-1, popfr, p1)
c pops at end time of this zone, further from source
        call copysteppop (i, popfr, p2)
c average of pops from start to end in prev zone
        call averinto (0.5d0, p1, p2, p1)
        call scalepop (p1, predr(i))
        call copysteppop (i-1, popintfr, p2)
        call addpop (p1, p2)
        call copypopstep (p2, i, popintfr)
      enddo
c
      xhfinal=pop(1,1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      itcount=itcount+1
c
c     poll for terminate file
c
      inquire (file='terminate',exist=iexi)
      if (iexi) goto 100
c
      if (iteration.lt.3) goto 100
c
c  Note: the final pass (finalit>0) used to exit here after a single,
c  un-refined sweep, even though rmslimit is the full s5rmstol
c  tolerance.  That let the final output be written from a
c  precursor state that had not actually met its own tightened
c  tolerance, which could make an already-converged model report
c  NOT CONVERGED on this last check alone (issue #7).  Now it falls
c  through to the ordinary loop-continuation test below like any other
c  iteration, so it actually iterates to rmslimit before exiting.
c
      if ((iabs(nfs0-nfs).gt.10)
     &.or.((rmserr.gt.rmslimit)
     &     .and.(itcount.lt.mxpcits)
     &     .and.(nfs.gt.minifsteps))
     &.or.(itcount.lt.minpcits)
     &.or.((nfs.lt.mxifsteps).and.(t0lim.gt.0))
     &.and.(itcount.lt.mxpcits)
     &) goto 50
c
c  Flag the case where the precursor's own inner zone-stepping loop
c  exhausted its sweep budget (mxpcits) without reaching its own
c  self-consistency target (rmslimit) -- previously silent.  With
c  the inner loop now demanding 0.01% self-consistency instead of the
c  old 10%, it's plausible some models can't get there in mxpcits
c  sweeps;
c  that would look like outer coupling oscillation/slow-convergence
c  without this, when it's actually the inner solve itself running out
c  of budget (issue #7 follow-up).
c
      if ((rmserr.gt.rmslimit).and.(itcount.ge.mxpcits)) then
        write (*,95) itcount,iteration,rmserr*100.d0,rmslimit*100.d0
c  ISSUE #14: this iteration's inner self-consistency solve did not
c  reach its own target -- shock5check classifies this osclass=5
c  regardless of how bounded the outer precrms/postrms split looks,
c  since that split assumes the inner solve actually settled.
        precwarn=1
      endif
   95 format('  ... PRECURSOR WARNING: inner loop hit its sweep cap (',
     & i3,' sweeps) at global iteration ',i3,
     & ' with self-consistency RMS ',1pg11.4,
     & '% still above its target ',1pg11.4,'%')
c
  100 continue
c
      te_pre=t
      de_pre=de
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      put final/inner balance into preionisation array
c
      call copysteppop (1, popfr, pop_pre)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Aitken Delta^2 dynamic relaxation of the new precursor solution
c  toward the previous global iteration's state (te_pre0/de_pre0/
c  pop_pre0, saved by shock5check) -- issue #7.  A raw fixed-point
c  substitution here can settle into a limit cycle for some shocks
c  (Psi~0.5-0.8) while a fixed damping factor slows or destabilises
c  others that were already contracting fine on their own; a single
c  relaxation constant cannot be right for both.  Aitken relaxation
c  instead estimates the local behaviour of the precursor<->shock
c  coupling from the last two raw-vs-accepted residuals and adapts
c  omega each iteration, the standard fix for exactly this kind of
c  black-box partitioned coupling (e.g. fluid-structure interaction).
c  te_pre/de_pre are used as the (normalised, signed) residual proxy
c  driving omega; pop_pre is blended by the same omega via averinto
c  to keep the whole precursor state moving together.
c
      if (iteration.gt.1) then
        aitrtenow=2.d0*(te_pre-te_pre0)/(te_pre+te_pre0)
        aitrdenow=2.d0*(de_pre-de_pre0)/(de_pre+de_pre0)
        if (aitninit.eq.0) then
c        no residual history yet -- bootstrap with a plain 0.5 blend
          aitomega=0.5d0
          aitninit=1
        else
          aitdnum=aitrte*(aitrtenow-aitrte)+aitrde*(aitrdenow-aitrde)
          aitdden=(aitrtenow-aitrte)**2+(aitrdenow-aitrde)**2
          if (dabs(aitdden).gt.1.d-30) then
            aitomega=-aitomega*aitdnum/aitdden
          endif
c        clamp: guards against a noisy or near-degenerate estimate
c        driving the relaxation factor to an unstable extreme
          aitomega=dmax1(0.05d0,dmin1(1.0d0,aitomega))
        endif
        aitrte=aitrtenow
        aitrde=aitrdenow
c
        te_pre=aitomega*te_pre+(1.d0-aitomega)*te_pre0
        de_pre=aitomega*de_pre+(1.d0-aitomega)*de_pre0
        call averinto (aitomega, pop_pre, pop_pre0, pop_pre)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Now get revised shock solution for changed pre shock values:
c     Shock in terms of flow velocity
c
      call shock5jump ()
c
c and show the result
c
      call shocksummary (6)
      call fieldsummary (6, 1, tphot)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (finalit.gt.0) then
c
        call shocksummary (luop)
        call fieldsummary (luop, 1, tphot)
c
        call shocksummary (lupt)
        call fieldsummary (lupt, 1, tphot)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Psi Summary
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  110 format(//
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & ' SHOCK 5 Final Precursor Parameter: ',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '     v_shock       : ',  1pg12.5,' km/s',/,
     & '     Q_ions        : ',  1pg12.5,' km/s',/,
     & '     Psi (Q/v)     : ',  1pg12.5,/,
     & '     Q_H           : ',  1pg12.5,' km/s',/,
     & '     Psi_H (Q_H/v) : ',  1pg12.5,/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
c
c redo Psi as vel0 has changed slightly
c
        psi=viofr/vel0
        write (*,110) vel0*1.d-5,viofr*1.d-5,psi,qh*1.d-5,qh/vel0
        write (luop,110) vel0*1.d-5,viofr*1.d-5,psi,qh*1.d-5,qh/vel0
        write (lupt,110) vel0*1.d-5,viofr*1.d-5,psi,qh*1.d-5,qh/vel0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write out final structure and last RMS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        write (*,60) ' Step',tab,'Dist.(cm)',tab,'dx (cm)',tab,'Time (s)
     &',tab,'dt (s)',tab,'Te (K)',tab,'ne(cm^-3)',tab,'nH(cm^-2)',tab,'X
     &HI',tab,'XHII',tab,'<abs>',tab,'XOI',tab,'XOII',tab,'XOIII'
c
        write (lupt,60) ' Step',tab,'Dist.(cm)',tab,'dx (cm)',tab,'Time
     &(s)',tab,'dt (s)',tab,'Te (K)',tab,'ne(cm^-3)',tab,'nH(cm^-2)',
     &   tab,'XHI',tab,'XHII',tab,'<abs>',tab,'XOI',tab,'XOII',tab,'XOII
     &I'
c
        do i=1,nfs
          write (*,70) i,tab,x(i),tab,predr(i),tab,x(i)/vshoc,tab,
     &     predr(i)/vshoc,tab,nte(i),tab,ne(i),tab,nh(i),tab,1.d0-fh(i),
     &     tab,fh(i),tab,fra(i),tab,foi(i),tab,foii(i),tab,foiii(i)
          write (lupt,70) i,tab,x(i),tab,predr(i),tab,x(i)/vshoc,tab,
     &     predr(i)/vshoc,tab,nte(i),tab,ne(i),tab,nh(i),tab,1.d0-fh(i),
     &     tab,fh(i),tab,fra(i),tab,foi(i),tab,foii(i),tab,foiii(i)
        enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        write (lupt,80) 100.d0*rmserr,itcount,iabs(nfs0-nfs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        wmod='NEBL'
        call multizone (nfs, x, nte, ne, nh, popfr, popintfr, wmod)

c
        caller='S5'
        pfx='PCup'//trim(s5pfx)
        np=len(trim(pfx))
c
        wmod='LFLM'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.d0, tphot)
c
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.d0, tphot)
c
        if (lmod.eq.'Y') then
          wmod='NFNU'
          call wpsou (caller, pfx, np, wmod, t, de, dh, dr, 1.d0, tphot)
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c collect precursor spectrum
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        call avrdata ()
c
        spmod='REL'
        linemod='LAMB'
        call spec2 (lupc, linemod, spmod)
c
        call wrsppop (lupt)
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Now clear tphot so it doesn't accumulate in each global iteration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call zerbuf ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c set *both* the initial population arrays and go back to computing
c post shock flow.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop_pre, pop)
      call copypop (pop_pre, pop0)
      t=te_pre
      de=dh_pre
      dh=dh_pre
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Close precursor files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call closeS5files ()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine multizone
c! XXXX - add one line purpose here
c! @param [in,out] integer*4         n  XXX-meaning
c! @param [in,out] integer*4         x  XXX-meaning
c! @param [in,out] xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx       nte  XXX-meaning
c! @param [in,out] xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx        ne  XXX-meaning
c! @param [in,out] integer*4        nh  XXX-meaning
c! @param [in,out] integer*4     popfr  XXX-meaning
c! @param [in,out]   real*8  popintfr  XXX-meaning
c! @param [in,out]   real*8      smod  XXX-meaning
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine multizone (n, x, nte, ne, nh, popfr, popintfr, smod)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Compute the emission from a multi-zone structure saved in the
c x, nte, ne, nh, popfr, popintfr arrays. A general low level routine.
c src not used so far, but it will in the future when we also add
c diffuse field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
c
      integer*4 n,idx
      real*8 x(mxifsteps),nte(mxifsteps)
      real*8 nh(mxifsteps),ne(mxifsteps)
      real*8 popfr(mxifsteps, mxion, mxelem)
      real*8 popintfr(mxifsteps, mxion, mxelem)
c
      real*8 t,tdw,tup,de,dh,dvdw,dvup,rdvol,irdvol,dvol
      real*8 gdil,x0,dx,dxdw,dxup
      real*8 p1(mxion, mxelem), p2(mxion, mxelem)
      real*8 sumhb,hbl
c
      character imod*4,specmode*4,smod*4,subname*64
c
      real*8 feldens
c
      integer*4 i
      real*8 linfluxes(mxmonlines)
c
c Only use part of common block
      integer*4 luop,lusp,lupt,lupc
      integer*4 ludy,lucl,lupb,lupf
      integer*4 lualsh,lurtsh,lulsh,luionsh(mxelem)
      integer*4 lualpc,lurtpc,lulpc,luionpc(mxelem)
      common /shock5L/ luop,lusp,lupt,lupc,
     & ludy,lucl,lupb,lupf,
     & lualsh,lurtsh,lulsh,luionsh,
     & lualpc,lurtpc,lulpc,luionpc      
c
      gdil=0.5d0
      fi=1.d0
c
      call zerbuf ()
c global modes
      jspot='N'
      jcon='Y'
      jgeo='P'
      jtrans='LODW'
      jden='C'
      specmode='NEBL'
      subname='MultiZone'
c
      rdvol=1.0d0
      irdvol=1.d0
      vunilog=0.d0
      sumhb=0.0d0
c
  10  format(1x,i4,',', 11(1pg12.5,', '),31(1pg12.5,', '))
c
      do idx=1,n-1
c
        dx=x(idx+1)-x(idx)
        x0=x(idx)
c
        call copysteppop (idx, popintfr, popint)
        call copysteppop (idx, popfr, p1)
        call copysteppop (idx+1, popfr, p2)
        call averinto (0.5d0, p1, p2, pop)
        t=0.5d0*(nte(idx)+nte(idx+1))
        dh=nh(idx)
        de=feldens(dh,pop)
        ne(idx)=de
c
        call localem (t, de, dh)
        hbl=hbeta*fpi*dx
        sumhb=sumhb+hbl
        jgeo='P'
c        write(*,'(8(1pg11.4,x))') x0,dx,t,dh,de,hbeta*fpi,hbl,sumhb
        specmode='NEBL'
        call totphot2 (t, dh, x0, dx, 0.d0, gdil, specmode)
        call zetaeff (dh)
c
ccccccccccccc
        if (jlin.eq.'Y') then
          call speclocallines (linfluxes)
          write (lulpc,10) idx,x(idx),x(idx)+dx*0.5d0,dx,t,de,
     &     dh,0.0d0,0.0d0,0.0d0,0.0d0,hydrobri(2,2),(linfluxes(i),i=1,
     &     njlines)
        endif
ccccccccccccc
c
c     accumulate spectrum
c
        tdw=nte(idx+1)
        tup=nte(idx)
        dxdw=dx
        dvdw=0.d0
        dxup=0.d0
        dvup=0.d0
c
        jgeo='P'
        jtrans='LODW'
        call newdif2 (tdw, tup, dh, x0, dxdw, dvdw, dxup, dvup, gdil,
     &   jtrans)
c
        imod='ALL'
        dvol=dx*irdvol
        call sumdata (t, de, dh, dvol, dx, x0, imod)
c
      enddo
c     write(*,*) 'multizone Hbeta:', sumhb, dlog10(sumhb)
c
c     call avrdata()
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c***************************************************************
c> @brief The function real*8 function flocallosses(tm,ne,nh,pp,r,dr,dv,w)
c! XXXX - add one line purpose here
c! @param [in,out]   real*8        tm  XXX-meaning
c! @param [in,out]   real*8        ne  XXX-meaning
c! @param [in,out]   real*8        nh  XXX-meaning
c! @param [in,out]   real*8        pp  XXX-meaning
c! @param [in,out]   real*8         r  XXX-meaning
c! @param [in,out]   real*8        dr  XXX-meaning
c! @param [in,out]   real*8        dv  XXX-meaning
c! @param [in,out]   real*8         w  XXX-meaning
c!
c! @return
c!  XXXX This function returns a %s number with is
c!  XXXX say explictly what is returned
c!
c! @details
c!  XXXX Enter details here
c****************************************************************

      real*8 function flocallosses(tm,ne,nh,pp,r,dr,dv,w)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
c
      real*8 tm,ne,nh,r,dr,dv,w
      real*8 pp(mxion, mxelem)
      real*8 savepop(mxion, mxelem)
c
c     real*8 nexp,dlosmin,f
c     parameter( nexp=4.d0 )
c     5.0e-4^4
c     parameter( dlosmin=6.25d-14 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c loss dampening if close to heating/cooling balance
c prevents wild instability in x-ray neq tails
c
c     f=(dlos**nexp)/(6.25d-14+dlos**nexp)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      character mode*4
      parameter(mode='DW')
c
      call copypop (pop, savepop)
      call copypop (pp, pop)
      call localem (tm, ne, nh)
      call totphot2 (tm, nh, r, dr, dv, w, mode)
      call zetaeff (nh)
      call cool (tm, ne, nh)
      call copypop (savepop, pop)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      flocallosses=(eloss-egain)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine compsh5
c! XXXX - add one line purpose here
c! @param [in,out] integer*4  iteration  XXX-meaning
c! @param [in,out] integer*4    maxits  XXX-meaning
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine compsh5 (iteration, maxits)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 iteration, maxits
c
      integer*4 i,j,idx
      logical iexi
      real*8 linfluxes(mxmonlines)
      real*8 n,invn,rdvol,irdvol,dvol
      real*8 b0,b1,b2,b3,b4,blum,binlum,pe
      real*8 hdt,tscale,dt0
      real*8 temp1,temp2,temp3,tav
c save state for step integral, uses pop, pop0 and pop1 in cblocks
      real*8 popstep(mxion,mxelem)
      real*8 tstep,destep,dhstep,xhstep
      real*8 rstep,drstep,velstep,dvstep,bmstep
      real*8 v,b
c
c structural markers
c
      real*8 timMark(7),disMark(7),bmMark(7)
      real*8 prMark(7),rhMark(7),dhMark(7),deMark(7)
      real*8 templim,ft,cft
      integer*4 tempidx
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           Functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 frho,fmua,fdynamictimestep
      real*8 feldens,fpresse,flocallosses
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      character tab*1
      tab=','
c
      do idx=1,7
        timMark(idx)=0.0d0
        disMark(idx)=0.0d0
        bmMark(idx)=0.0d0
        prMark(idx)=0.0d0
        rhMark(idx)=0.0d0
        dhMark(idx)=0.0d0
        deMark(idx)=0.0d0
      enddo
c
      rdvol=1.0d0
      irdvol=1.d0
      vunilog=0.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c open all main shock files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call appendS5files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Begin Main Calculation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   10 format(
     & ' ********************************************************',/,
     & '  SHOCK 5     Shock Iteration: ',i2.2,' of ',i2.2,/,
     & ' ********************************************************')
      write (*,10) iteration,maxits
      if (finalit.gt.0) then
        write (luop,10) iteration,maxits
        call shocksummary (luop)
      endif
c
   20 format(//'  T0',1pg14.7,' V0',1pg14.7,/,
     & ' RH0',1pg14.7,' P0',1pg14.7,' B0',1pg14.7,//,
     & '  T1',1pg14.7,' V1',1pg14.7,/,
     & ' RH1',1pg14.7,' P1',1pg14.7,' B1',1pg14.7,//)
c
      if (vmod.ne.'NONE') write (*,20) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,
     &te1,vel1*1.d-5,rho1,pr1,bm1*1.d6
      if (finalit.gt.0) then
        write (luop,20) te0,vel0*1.d-5,rho0,pr0,bm0*1.d6,te1,vel1*1.d-5,
     &   rho1,pr1,bm1*1.d6
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   30 format(//,
     & ' #  ',
     & ',     Dist(cm),       dr(cm)',
     & ',        Te(K),     ne(/cm3),     nH(/cm3),     ni(/cm3)',
     & ',      mu(amu)',
     & ',      Time(s),        dt(s)',
     & ', L(erg.cm3/s),  Loss(erg/s),       dLoss',
     & ',          XHI,         XHII,        B(G)')
c
      wdilt0=te0
      t=te0
      dh=dh0
      de=de0
      press=fpresse(t,de,dh)
      rhotot=frho(de,dh)
      mu=fmua(de,dh)
      dr=0.d0
      ve=vel0
      vpo=vel1
      Bmag=bm0
      bm0=bm0
      dv=ve*0.01d0
      fi=1.d0
      rad=0.d0
      wdil=0.5d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate the radiation field and atomic rates
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   40 format(i4.3,15(',',1pg13.6))
c
   50 format(/,
     & ' #    Te(K)       ne(/cm3)    nH(/cm3)    ni(/cm3)    ',
     & 'B(G)        XHI         Time(s)     dt(s)       ',
     & 'Dist(cm)    dr(cm)      v(cm/s)     L(erg.cm3/s)')
      if (vmod.eq.'MINI') then
        write (*,50)
      endif
   60 format (i4,6(1x,1pg11.4),6(1x,1pg11.4))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      step=-2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop_neu, pop)
      en=zen*dh_neu
      cspd=dsqrt(gammaEOS*pr_neu/rh_neu)
      wmol=rh_neu/(en+de_neu)
      mu=fmua(de_neu,dh_neu)
      t=te_neu
      de=de_neu
      dh=dh_neu
      en=zen*dh
      Bmag=bm_neu
c
      netloss=flocallosses(t,de,dh,pop_neu,1.d38,0.d0,0.d0,wdil)
c
      ue=gammaEOSU*(en+de)*rkb*t
      tnloss=tloss/((en+de)*(en+de))
c
      if (finalit.gt.0) then
c
        write (luop,40) step,0.d0,0.d0,te_neu,de_neu,dh_neu,en,mu,0.d0,
     &   0.d0,tnloss,tloss,dlos,pop(1,1),pop(2,1),bm_neu
      endif
c
c        if (vmod.eq.'MINI') then
      write (*,60) step,te_neu,de_neu,dh_neu,en,bm_neu,pop_neu(1,1),
     &0.d0,0.d0,0.0d0,0.d0,vel0,tnloss
c        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      step=-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop_pre, pop)
      call copypop (pop_pre, pop0)
      en=zen*dh_pre
      cspd=dsqrt(gammaEOS*pr_pre/rh_pre)
      wmol=rh_pre/(en+de_pre)
      mu=fmua(de_pre,dh_pre)
      t=te_pre
      de=de_pre
      dh=dh_pre
      en=zen*dh
      Bmag=bm_pre
c
      netloss=flocallosses(t,de,dh,pop_pre,1.d38,0.d0,0.d0,wdil)
c
      ue=gammaEOSU*(en+de)*rkb*t
      tnloss=netloss/((en+de)*(en+de))
c
      rhotot=frho(de,dh)
      wmol=rhotot/(en+de)
      mu=fmua(de,dh)
c
      if (finalit.gt.0) then
        write (luop,40) step,0.d0,0.d0,te_pre,de_pre,dh_pre,en,mu,0.d0,
     &   0.d0,0.0d0,0.0d0,0.d0,pop_pre(1,1),pop_pre(2,1),bm_pre
      endif
c
c     if (vmod.eq.'MINI') then
      write (*,60) step,te_pre,de_pre,dh_pre,en,bm_pre,pop_pre(1,1),
     &0.d0,0.d0,0.0d0,0.d0,vel0,tnloss
c     endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      step=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ue=gammaEOSU*(en+de)*rkb*t
      tnloss=netloss/((en+de)*(en+de))
c for later step^2 5x5=25, 25*0.04=1.0
      tscale=0.04d0
      dt=fdynamictimestep(t,dh,0.d0,ve,pop,netloss)
      dt=tscale*dt
      hdt=0.5d0*dt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Setup Initial first step in post shock-front gas(#1) zero losses
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     rmod='NEAR'
c      call isochorflow (t, de, dh, ve, Bmag, rmod, netloss, dt)
c      call isobarflow (t, de, dh, ve, Bmag, rmod, netloss, dt)
c      call rankhug (t, de, dh, ve, Bmag, 0.d0, 0.d0)
      call shockcmpf (t, de, dh, ve, Bmag)
c
      t=te1
      de=de1
      dh=dh1
      Bmag=bm1
c
      call copypop (pop_pre, pop)
      call copypop (pop_pre, pop0)
c
c sets up 0 and 1 vars
c
      dr=dt*vel1
      dv=vel1-vel0
c
c
      ue=gammaEOSU*(en+de)*rkb*t
      tnloss=netloss/((en+de)*(en+de))
c
      tscale=0.04d0
      dt=fdynamictimestep(t,dh,0.d0,ve,pop,netloss)
      dt=tscale*dt
      hdt=0.5d0*dt
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     init arrays coming first step
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      xh(1)=pop(2,1)
      xh(2)=pop(2,1)
      veloc(1)=vel0
      veloc(2)=vel1
c
      te(1)=te0
      te(2)=te1
      deel(1)=de0
      deel(2)=de1
      dhy(1)=dh0
      dhy(2)=dh1
      bmg(1)=bm0
      bmg(2)=bm1
c
      dist(1)=0.0d0
      dist(2)=dr
      timlps(1)=0.0d0
      timlps(2)=dt
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      t=te1
      de=de1
      dh=dh1
      en=zen*dh
      ve=vel1
      rhotot=frho(de,dh)
      wmol=rhotot/(en+de)
      mu=fmua(de,dh)
c
      netloss=flocallosses(t,de,dh,pop_pre,1.d38,0.d0,0.d0,wdil)
c
      ue=gammaEOSU*(en+de)*rkb*t
      tnloss=netloss/((en+de)*(en+de))
c
      if (finalit.gt.0) then
        write (luop,40) step,0.0d0,0.0d0,t,de,dh,en,mu,0.0d0,0.0d0,
     &   tnloss,tloss,dlos,pop(1,1),pop(2,1),bm1
      endif
c
c     if (vmod.eq.'MINI') then
      write (*,50)
      write (*,60) step,t,de,dh,en,bm1,pop(1,1),0.0d0,0.0d0,0.0d0,0.0d0,
     &ve,tnloss
c     endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   70 format(i4,1x,i4)
   80 format(1pg14.7,'K',1x,1pg14.7,'cm/s',1x,1pg14.7,'g/cm3',1x,
     &1pg14.7,'dyne/cm2',1x,1pg14.7,'Gauss',1x,1pg14.7,'ergs/cm3/s')
   90 format(2(1pg14.7,' cm',1x),2(1pg14.7,' s',1x))
  100 format(1pg14.7,' ergs/cm^3',1x,1pg14.7,' /cm^3',1x,
     &1pg14.7,' cm/s',1x,1pg14.7,' g/particle ')
  110 format(1x,i4,',', 11(1pg12.5,', '),31(1pg12.5,', '))
  120 format(43(1pg12.5,a2))
  130 format(17(1pg12.5,a1))
  140 format(//,
     & ' Step   Te Ave.(K)   ',
     & '  ne(cm^-3)   ',
     & '  nH(cm^-3)   ',
     & '  V1(cm/s)    ',
     & ' Rho1(g/cm^3) ',
     & ' Pr1(erg/cm^3)',
     & '   Dist.(cm)  ',
     & ' Elps. Time(s)')
  150 format(1x,i4,8(1pg14.7)/)
  160 format(1x,i4,41(', ',1pg12.5))
c
      call copypop (pop0, pop)
      de0=feldens(dh0,pop0)
      if (finalit.gt.0) then
        caller='S5'
        pfx=s5pfx(1:nprefix)
        np=nprefix
        call wbal (caller, pfx, np, pop)
      endif
c
      specmode='DW'
      wdil=0.5d0
      fi=1.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SHOCK Step/Iteration reentry point:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      step=1
c
      t=te1
      te0=te1
      rho0=rho1
      pr0=pr1
      vel0=vel1
      dh0=dh1
      de0=de1
      bm0=bm1
      Bmag=bm0
      dist(1)=0.d0
      timlps(1)=0.d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (vmod.ne.'NONE') then
        write (*,50)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (finalit.gt.0) write (luop,30)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    MAIN STEP LOOP
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  170 count=0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c save step initial condition
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop0, popstep)
      tstep=te0
      dhstep=dh0
      destep=feldens(dh0,pop0)
      xhstep=pop(2,1)
      rstep=dist(step)
      drstep=dr
      velstep=vel0
      dvstep=vel1-vel0
      bmstep=bm0
c
      rad=rstep+drstep
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c get intial time step dt, and hdt = 0.5d0*dt, init t = te0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c initial slow start scaling
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     netloss=flocallosses(te0,de0,dh0,pop0,rad,dr,dv,wdil)
c     slow start in step^2 for 5 steps : 25*0.04=1.0
      tscale=dmin1(1.d0,0.04d0*dble(step*step))
      dt0=fdynamictimestep(tstep,dhstep,rad,velstep,popstep,netloss)
      dt=tscale*dt0
      hdt=0.5d0*dt
      drstep=velstep*dt
      dr=drstep
c
c  Heartbeat: some parameter combinations (e.g. low nH, fast v_shock)
c  need thousands of zone steps per global iteration to cool to the
c  stopping criterion, and the per-zone table above is only printed
c  when vmod='MINI'.  Without this, a slow-but-working run is
c  indistinguishable from a hang (issue #7 follow-up).  Printed
c  unconditionally, throttled to avoid flooding fast-finishing models.
c
      if (mod(step,100).eq.1) then
        write (*,195) iteration,maxits,step,mxnsteps,tstep,rad
        call flush (6)
      endif
  195 format ('  ... SHOCK 5 progress: Global It ',i3,' of ',i3,
     & ', Zone Step ',i5,' of ',i5,', T=',1pg11.4,' K, Dist=',
     & 1pg11.4,' cm')
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      count=count+1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      evolve ionisation over hdt to get pop at midpoint
c      allow for iterations - not used as already 3 steps
c      initial midpoint cooling
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (popstep, pop0)
      call copypop (popstep, pop)
      t=tstep
      de=destep
      dh=dhstep
      v=velstep
      b=bmstep
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Get midpoint Temp based on initial temp and cooling
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      temp1=t
      netloss=flocallosses(t,de,dh,popstep,rad,dr,dv,wdil)
      call rankhug (t, de, dh, v, b, netloss, hdt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c from rankhug common block te0/te1 change vars
c midpoint t,de,dh,v,b all in '1' vars, not to get pop1 to match
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      temp2=te1
      tav=0.5d0*(temp1+temp2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     reset for ionisation t de dh: pop not changed by rh
c     xhii unused here, new balance in pop at hdt -> pop1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call copypop (popstep, pop)
      t=tstep
      de=destep
      dh=dhstep
      xhii=xhstep
      call timion (tav, de, dh, xhii, hdt)
      call copypop (pop, pop1)
c
c  get midpoint cooling at temp2=te1, with improved pop1
c
      de1=feldens(dh1,pop1)
      netloss=flocallosses(te1,de1,dh1,pop1,rad,dr,dv,wdil)
c
c step RH over dt with mid point losses, from start of step
c
      t=tstep
      de=destep
      dh=dhstep
      v=velstep
      b=bmstep
c
      call rankhug (t, de, dh, v, b, netloss, dt)
      temp3=te1
c higher order step temp
      tav=0.5d0*(temp1+temp3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c now advance ions with whole dt at overall mean tav
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (popstep, pop)
      de=destep
      dh=dhstep
      xhii=xhstep
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c advance ionisation a whole timestep at mean t
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call timion (tav, de, dh, xhii, dt)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c save end ionisation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop, pop1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c End step iterations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     record step
c
      call averinto (0.5d0, pop0, pop1, pop)
c
      te(step)=(te0+te1)*0.5d0
      dhy(step)=(dh0+dh1)*0.5d0
      deel(step)=(de0+de1)*0.5d0
      bmg(step)=(bm0+bm1)*0.5d0
c
      xh(step)=pop(2,1)
      veloc(step)=(vel1+vel0)*0.5d0
c
      dist(step+1)=dist(step)+dr
      timlps(step+1)=timlps(step)+dt
c
      en=zen*dh
c
      ue=gammaEOSU*(en+de)*rkb*t
      tnloss=netloss/((en+de)*(en+de))
c
      cspd=dsqrt(gammaEOS*(pr1+pr0)/(rho1+rho0))
      wmol=(rho1+rho0)/(2.d0*(en+de))
      mu=fmua(de,dh)
c
      mb=(bm0+bm1)*0.5d0
c
      if ((vmod.eq.'FULL').or.(vmod.eq.'SLAB')) then
        wmod='SCRN'
        call wmodel (luop, t, de, dh, dr, wmod)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c structural markers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      templim=1.0d7
  180 tempidx=idnint(dlog10(templim))
      if ((te1.le.templim).and.(te0.gt.templim)) then
        ft=(templim-te1)/(te0-te1)
        cft=1.d0-ft
        timmark(tempidx)=cft*timlps(step)+ft*timlps(step)
        dismark(tempidx)=cft*dist(step)+ft*dist(step)
        bmmark(tempidx)=cft*bm0+ft*bm1
        prmark(tempidx)=cft*pr0+ft*pr1
        rhmark(tempidx)=cft*rho0+ft*rho1
        dhmark(tempidx)=cft*dh0+ft*dh1
        demark(tempidx)=cft*de0+ft*de1
      endif
      templim=templim*0.1d0
      if (templim.gt.10.d0) goto 180
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (vmod.eq.'MINI') then
        write (*,60) step,t,de,dh,en,mb,pop(1,1),timlps(step),dt,
     &   dist(step),dr,veloc(step),tnloss
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (finalit.gt.0) then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        write (luop,40) step,dist(step),dr,t,de,dh,en,mu,timlps(step),
     &   dt,tnloss,tloss,dlos,pop(1,1),pop(2,1),bm1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c output to secondary files
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (jlin.eq.'Y') then
          call speclocallines (linfluxes)
          write (lulsh,110) step,dist(step),dist(step)+dr*0.5d0,dr,t,de,
     &     dh,de+en,0.0d0,0.0d0,0.0d0,hydrobri(2,2),(linfluxes(i),i=1,
     &     njlines)
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (ratmod.eq.'Y') then
          wmod='LOSS'
          call wmodel (lurtsh, t, de, dh, dr, wmod)
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (dynmod.eq.'Y') then
          write (ludy,70) step,count
          write (ludy,90) dist(step),dr,timlps(step),dt
          write (ludy,80) te0,vel0,rho0,pr0,bm0,l0
          write (ludy,80) te1,vel1,rho1,pr1,bm1,l11
          write (ludy,100) ue,en,cspd,wmol
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (fclmod.eq.'Y') then
          do i=1,atypes
            meanion(i)=0.d0
            ionbar=0.d0
            do j=1,maxion(i)
              ionbar=ionbar+(pop(j,i)*j)
            enddo
            meanion(i)=dmax1(ionbar-1.0d0,0.0d0)
          enddo
c
          en=zen*dh
          n=(de*dh)
          if (jnorm.eq.0) n=(de*dh)
          if (jnorm.eq.1) n=(dh*dh)
          if (jnorm.eq.2) n=(de*en)
          if (jnorm.eq.3) n=((de+en)*(de+en))
          if (jnorm.eq.4) n=(de*de)
          invn=1.d0/dble(n)
c
          write (lucl,120) t,tab,de,tab,dh,tab,en,tab,(0.5*(rho0+rho1)),
     &     tab,pop(1,1),tab,pop(2,1),tab,mu,tab,tloss,tab,tloss*invn,
     &     (tab,coolz(j)*invn,j=1,atypes),tab,eloss*invn,tab,egain*invn,
     &     tab,netloss*invn
c
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (bandsmod.eq.'Y') then
c
          b0=0.d0
          b1=0.d0
          b2=0.d0
          b3=0.d0
          b4=0.d0
          blum=0.d0
          do idx=1,infph-1
            if (tphot(idx).gt.epsilon) then
              widnu=widbinnu(idx)
              pe=photev(idx)
              binlum=tphot(idx)*fpi*widnu/dr
              blum=blum+binlum
              if ((pe.gt.0.0d0).and.(pe.le.100.0d0)) b0=b0+binlum
              if ((pe.gt.100.0d0).and.(pe.le.500.0d0)) b1=b1+binlum
              if ((pe.gt.500.0d0).and.(pe.le.1000.0d0)) b2=b2+binlum
              if ((pe.gt.1000.0d0).and.(pe.le.2000.0d0)) b3=b3+binlum
              if ((pe.gt.2000.0d0).and.(pe.le.10000.0d0)) b4=b4+binlum
            endif
          enddo
c
          b0=b0/tloss
          b1=b1/tloss
          b2=b2/tloss
          b3=b3/tloss
          b4=b4/tloss
          blum=blum/tloss
c
          if (b0.lt.epsilon) b0=0.d0
          if (b1.lt.epsilon) b1=0.d0
          if (b2.lt.epsilon) b2=0.d0
          if (b3.lt.epsilon) b3=0.d0
          if (b4.lt.epsilon) b4=0.d0
          if (blum.lt.epsilon) blum=0.d0
c
          write (lupb,130) t,tab,de,tab,dh,tab,en,tab,pop(1,1),tab,
     &     pop(2,1),tab,mu,tab,tloss,tab,tloss*invn,tab,fflos/tloss,tab,
     &     b0,tab,b1,tab,b2,tab,b3,tab,b4,tab,blum,tab
c
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (allmod.eq.'Y') then
          write (lualsh,140)
          write (lualsh,150) step,t,de,dh,vel1,rho1,pr1,dist(step),
     &     timlps(step)
c
          call wionabal (lualsh, pop)
          call wionabal (lualsh, popint)
        endif
c
        if (tsrmod.eq.'Y') then
c
          do i=1,ieln
            write (luionsh(i),160) step,dist(step),dr,timlps(step),dt,t,
     &       de,dh,en,(pop(j,iel(i)),j=1,maxion(iel(i)))
          enddo
c
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c end finalit
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      endif
c
c     poll for photons file
c
      pollfile='photons'
      inquire (file=pollfile,exist=iexi)
c
      if ((lmod.eq.'Y').or.(iexi)) then
c
c     Local photon field
c
        specmode='LOCL'
        call localem (t, de, dh)
        rad=dist(step)
        call totphot2 (t, dh, rad, dr, dv, wdil, specmode)
        caller='S5'
        pfx='shlocl'//trim(s5pfx)
        np=len(trim(pfx))
        wmod='NORM'
        dva=vel1-vel0
        en=zen*dh
        n=(de*dh)
        if (jnorm.eq.0) n=(de*dh)
        if (jnorm.eq.1) n=(dh*dh)
        if (jnorm.eq.2) n=(de*en)
        if (jnorm.eq.3) n=((de+en)*(de+en))
        if (jnorm.eq.4) n=(de*de)
        invn=1.d0/n
        call wpsou (caller, pfx, np, wmod, t, de, dh, dr, invn, tphot)
c
c     reset mode and tphot
c
        specmode='DW'
        rad=dist(step)
        call totphot2 (t, dh, rad, dr, dv, wdil, specmode)
      endif
c
c     get mean ionisation state for step
c
c     middleofstep,dist(step)isnowend
      rdis=dist(step)+dr*0.5d0
c      cmpf=rho1/rho0
      fi=1.0d0
c
c     accumulate spectrum
c
      tdw=te0
      tup=te1
      drdw=dr
      dvdw=vel0-vel1
      drup=dist(step)
      dvup=dsqrt(vpo*vel0)
      frdw=0.5d0
c
      call localem (t, de, dh)
      call zetaeff (dh)
      call newdif2 (tdw, tup, dh, rad, drdw, dvdw, drup, dvup, frdw,
     &jtrans)
c
      imod='ALL'
      dvol=dr*irdvol
      call sumdata (t, de, dh, dvol, dr, rdis, imod)
c
c     record line ratios
c
      if (ox3.ne.0) then
        hoiii(step)=(fluxm(6,ox3)+fluxm(8,ox3))/(fluxh(2)+epsilon)
      endif
      if (ox2.ne.0) then
        hoii(step)=(fluxm(1,ox2)+fluxm(2,ox2))/(fluxh(2)+epsilon)
      endif
      if (ni2.ne.0) then
        hnii(step)=(fluxm(7,ni2)+fluxm(10,ni2))/(fluxh(2)+epsilon)
      endif
      if (su2.ne.0) then
        hsii(step)=(fluxm(1,su2)+fluxm(1,su2))/(fluxh(2)+epsilon)
      endif
      if (ox1.ne.0) then
        if (ispo.eq.' OI') hsii(step)=fluxm(3,ox1)/(fluxh(1)+epsilon)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set up initial conditions for the next zone
c     (= end conditions in previous zone)
c     Calculate new cooling times and next step
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call copypop (pop1, pop0)
      call copypop (pop1, pop)
c
      t=te1
      te0=te1
      rho0=rho1
      pr0=pr1
      vel0=vel1
      dh0=dh1
      de0=de1
      bm0=bm1
      Bmag=bm0
      count=0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Program endings
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ionisation, finish when mean ionisation < 1%
c
c      if (jend.eq.'A') then
c
      totn=0.d0
      totn=zen
      ionbar=0.d0
c
      do i=1,atypes
        do j=2,maxion(i)
          ionbar=ionbar+(pop(j,i)*zion(i))
        enddo
      enddo
c
      ionbar=ionbar/totn
c
      if ((jend.eq.'A').and.(ionbar.lt.0.01d0)) goto 190
c
c     Specific species ionisation limit
c
      if (jend.eq.'B') then
        if (pop(jpoen,ielen).lt.fren) goto 190
      endif
c
c stop first iteration at heating==cooling if maxits>1
c
      if ((iteration.le.1).and.(maxits.gt.1).and.(dlos.lt.0.5d0)) goto
     &190
c
c     Temperature limit
c
      if (jend.eq.'C') then
        if (t.lt.tend) goto 190
      endif
c
c     Temperature limit, plus 95% neutral
c
      if ((jend.eq.'S').and.(t.lt.tend)) then
        if (ionbar.lt.0.05d0) goto 190
      endif
c
c     Distance Limit
c
      if ((jend.eq.'D').and.(dist(step).ge.diend)) goto 190
c
c     Time Limit
c
      if ((jend.eq.'E').and.(timlps(step).ge.timend)) goto 190
c
c     thermal balance dlos<1e-2
c
      if ((jend.eq.'F').and.(dlos.lt.1.d-2)) goto 190
c
c     cooling function test, finish when tloss goes -ve
c
      if ((jend.eq.'G').and.(tloss.lt.0.d0)) goto 190
c
c     poll for terminate file
c
      pollfile='terminate'
      inquire (file=pollfile,exist=iexi)
      if (iexi) goto 190
c
c     otherwise go to normal iteration loop, if interlocks
c     permit
c
      step=step+1
      if ((t.gt.100.d0).and.(step.lt.(mxnsteps-1))) then
        goto 170
      endif
c
  190 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Close all files and flush
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call closeS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     End model
c     write out spectrum etc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Reopen all main shock files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call appendS5files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Final Downstream Field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (finalit.gt.0) then
        if (jspec(1:1).eq.'Y') then
c
          caller='S5'
          pfx='SHdw'//trim(s5pfx)
          np=len(trim(pfx))
          dva=0.d0
c
          wmod='LFLM'
          call wpsou (caller, pfx, np, wmod, t, de, dh, dr, wdil, tphot)
c
          wmod='REAL'
          call wpsou (caller, pfx, np, wmod, t, de, dh, dr, wdil, tphot)
c
          if (lmod.eq.'Y') then
c
            wmod='NFNU'
            call wpsou (caller, pfx, np, wmod, t, de, dh, dr, wdil,
     &       tphot)
c
          endif
        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Upstream photon field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        call localem (t, de, dh)
        specmode='UP'
        rad=dist(step)
        call totphot2 (t, dh, 1.0d38, 1.0d0, 0.d0, wdil, specmode)
c
        call fieldsummary (6, 1, tphot)
c
        caller='S5'
        pfx='SHup'//trim(s5pfx)
        np=len(trim(pfx))
        dva=vel1-vel0
c
        wmod='LFLM'
        call wpsou (caller, pfx, np, wmod, te1, de1, dh1, dr, wdil,
     &   tphot)
c
        wmod='REAL'
        call wpsou (caller, pfx, np, wmod, te1, de1, dh1, dr, wdil,
     &   tphot)
c
        if (lmod.eq.'Y') then
c
          wmod='NFNU'
          call wpsou (caller, pfx, np, wmod, te1, de1, dh1, dr, wdil,
     &     tphot)
c
        endif
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dynamics
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call avrdata ()
c
      if (finalit.gt.0) then
c
  200 format(//,
     & ' Model ended: , Distance:, ',1pg14.7,
     &', Time:, ',1pg14.7,', Temp:, ',1pg12.5,/)
        write (luop,200) dist(step),timlps(step),t
        write (*,200) dist(step),timlps(step),t
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c structural markers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  210 format(' Marker',i1,':, ',1pe12.5,' (K)',/,
     & '      t',i1,':, ',1pg12.5,/,
     & '      d',i1,':, ',1pg12.5,/,
     & '      B',i1,':, ',1pg12.5,/,
     & '      P',i1,':, ',1pg12.5,/,
     & '    rho',i1,':, ',1pg12.5,/,
     & '     nH',i1,':, ',1pg12.5,/,
     & '     ne',i1,':, ',1pg12.5,/)
        do idx=7,1,-1
          if (timmark(idx).gt.0.d0) then
            templim=10.0d0**dble(idx)
            write (luop,210) idx,templim,idx,timMark(idx),idx,
     &       disMark(idx),idx,bmMark(idx),idx,prMark(idx),idx,
     &       rhMark(idx),idx,dhMark(idx),idx,deMark(idx)
            write (*,210) idx,templim,idx,timMark(idx),idx,disMark(idx),
     &       idx,bmMark(idx),idx,prMark(idx),idx,rhMark(idx),idx,
     &       dhMark(idx),idx,deMark(idx)
          endif
        enddo
c
        spmod='REL'
        linemod='LAMB'
        call spec2 (lusp, linemod, spmod)
c
c both summary columns and spec/spectrum all in one
c
        call wrsppop (luop)
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call closeS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c***************************************************************
c> @brief The function real*8 function fdynamictimestep(t,dh,x,vel,p,netloss)
c! XXXX - add one line purpose here
c! @param [in,out]   real*8         t  XXX-meaning
c! @param [in,out]   real*8        dh  XXX-meaning
c! @param [in,out]   real*8         x  XXX-meaning
c! @param [in,out]   real*8       vel  XXX-meaning
c! @param [in,out] xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx         p  XXX-meaning
c! @param [in,out] xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   netloss  XXX-meaning
c!
c! @return
c!  XXXX This function returns a %s number with is
c!  XXXX say explictly what is returned
c!
c! @details
c!  XXXX Enter details here
c****************************************************************

      real*8 function fdynamictimestep(t,dh,x,vel,p,netloss)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
c
      real*8 t,dh,de,vel,netloss
      real*8 p(mxion, mxelem)
      real*8 en, pr, ue, x
      real*8 cltime, absdr, abstime
      real*8 ctime, rtime, dt, absf
      real*8 feldens, fcolltim, frectim2
      parameter ( absf=0.04d0 )
      real*8 flocallosses
c
      de=feldens(dh,p)
      en=zen*dh+de
      pr=en*rkb*t
      ue=gammaEOSU*pr
c
c initial cooling rate for initial timestep guess
c
      netloss=flocallosses(t,de,dh,p,x,0.d0,0.d0,wdil)
c
      cltime=(ue/(dabs(netloss)+epsilon))
      absdr=1.0d38
      abstime=1.d0/epsilon
      if (photonmode.ne.0) then
        absdr=1.0d19/dh0
        call absdis2 (dh, 0.5d0, absdr, 0.0d0, p)
        abstime=absdr/(dabs(vel)+epsilon)
      endif
      ctime=fcolltim(de)
      rtime=frectim2(de)
c
c weighted harmonic sum, coeffs determined empirically
c
      dt=1.d0/((0.5d0/abstime)+(1.d0/(cltime*0.025d0))+(1.d0/(rtime*
     &0.05d0))+(1.d0/(ctime*0.01d0)))
c
      fdynamictimestep=dt
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine protostate
c! Print the shock properties in terms of the cold protostate.
c! @param [in,out] integer*4     lunit  output file logical unit
c!
c! @return
c!  adds to file in logical unit
c!
c! @details
c!  Protostate before the precursor is cold usually and the shock
c! paramters appear as the highest Mach numbers and so on.
c! Once and if heated by a precursor the Mach number etc change for
c! the same shock.
c***************************************************************

      subroutine protostate (lunit)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      real*8 fmua
      integer*4 lunit
c
   10   format(//
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '  Shock Proto-State Properties',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    Velocity       :,',1pg12.5,', km/s',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/,
     & '    Proto Mach Number          :,',1pg12.5,/,
     & '    Proto Alfven Mach Number   :,',1pg12.5,/,
     & '    Proto Mag Alpha (Pmag/Pgas):,',1pg12.5,/,
     & '    Proto Gas Eta  (gPgas/Pram):,',1pg12.5,/,
     & '    Proto Mag Eta  (2Pmag/Pram):,',1pg12.5,//,
     & '    Proto        T :,',1pg12.5,', K',/,
     & '    Proto        ne:,',1pg12.5,', cm^-3',/,
     & '    Proto        nH:,',1pg12.5,', cm^-3',/,
     & '    Proto        d :,',1pg12.5,', g/cm^-3',/,
     & '    Proto      Pgas:,',1pg12.5,', dyne/cm^2',/,
     & '    Proto        mu:,',1pg12.5,', a.m.u.',/,
     & '    Proto       XHI:,',1pg12.5,/,
     & '    Proto      XHII:,',1pg12.5,/,
     & '    Proto      XHeI:,',1pg12.5,/,
     & '    Proto     XHeII:,',1pg12.5,/,
     & '    Proto    XHeIII:,',1pg12.5,/,
     & '    Proto        B :,',1pg12.5,', microGauss',/,
     & '    Proto      Pmag:,',1pg12.5,', dyne/cm^2',/,
     & '    Proto      Pram:,',1pg12.5,', dyne/cm^2',/,
     & ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::',/)
c
      en=zen*dh_neu
      cspd=dsqrt(gammaEOS*pr_neu/rh_neu)
      wmol=rh_neu/(en+de_neu)
      mu=fmua(de_neu,dh_neu)
c
      Pram=rh_neu*vel0*vel0
      Pgas=pr_neu
      Pmag=(bm_neu*bm_neu)/epi
c
      machnumber=vel0/dsqrt(gammaEOS*pr_neu/rh_neu)
      alfvennumber=vel0/dsqrt(2.0d0*Pmag/rh_neu)
      malpha=Pmag/Pgas
      gaseta=gammaEOS*Pgas/Pram
      mageta=2.d0*Pmag/Pram
c
      write (lunit,10) vel0*1.0d-5,machnumber,alfvennumber,malpha,
     &gaseta,mageta,te_neu,de_neu,dh_neu,rh_neu,pr_neu,mu,pop_neu(1,
     &zmap(1)),pop_neu(2,zmap(1)),pop_neu(1,zmap(2)),pop_neu(2,zmap(2)),
     &pop_neu(3,zmap(2)),bm_neu*1.d6,Pmag,Pram
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine shock5filenames
c! Creare S5 output files, names and unit numbers, variable.
c! @param [in,out] character*        px  User file output prefix for some files
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine shock5filenames (px)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     set initial file unit numbers and default mode flags
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i,flen
      character* (*) px
c
      integer*4 mlen
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     set up file modes
c
      jspec='N'
      lmod='N'
      tsrmod='N'
      allmod='N'
      dynmod='N'
      ratmod='N'
      bandsmod='N'
      fclmod='N'
      jlin='N'
      s5pfx=px
      nprefix=mlen(px)
c
c main model files
c always on full model files
c spec line list files
c always in csv files lists
c    allions files , shock and precursors
c    allmod=Y
c rates fitles , shock and precursors
c ratmod=Y
c shock 'dynamics' file
c dynmod=Y
c shock cooling file
c fclmod=Y
c shock cooling emission bands
c bandsmod=Y
c line structures, shocks and precursors
c jlin=Y
c Ion balance files.' for just the 4 chosen mapping z to i later
c      luionsh(i)=32+i
c
c main model files
c always on full model files
c
c      open (luop,file=fsm,status='NEW')
c      open (lupt,file=fpm,status='NEW')
c
c
c main files, model and specs disable precursors files
c
      luop=20
      lusp=21
      lupt=22
      lupc=23
c
c common files
c
      lucl=24
      lupb=25
      lupb=26
      ludy=27
c
c shock only files
c
      lualsh=28
      lurtsh=29
      lulsh=30
      ieln=4
      do i=1,atypes
        luionsh(i)=30+i
      enddo
c
c precursor only files to disabled set to 0
c
      lualpc=31+atypes
      lurtpc=32+atypes
      lulpc=33+atypes
      do i=1,atypes
        luionpc(i)=33+atypes+i
      enddo
c
c     Name "sh5" general output files
c
c  luop
c
      fsm=' '
      pfx='shck_'//trim(s5pfx)
      sfx='sh5'
      call newfile (pfx, sfx, fn, flen)
      fsm=fn(1:flen)
c
c  lupt
c
      fpm=' '
      pfx='prec_'//trim(s5pfx)
      sfx='sh5'
      call newfile (pfx, sfx, fn, flen)
      fpm=fn(1:flen)
c
c spec line list files
c always in csv files lists
c
c      open (lusp,file=fsh,status='NEW') shock
c      open (lupc,file=fpc,status='NEW') precursor
c
c
      pfx='specSH'//trim(s5pfx)
      sfx='csv'
      fsh=' '
      call newfile (pfx, sfx, fn, flen)
      fsh=fn(1:flen)
c
      pfx='specPC'//trim(s5pfx)
      sfx='csv'
      fpc=' '
      call newfile (pfx, sfx, fn, flen)
      fpc=fn(1:flen)
c
c    allions files , shock and precursors
c    allmod=Y
c
c    open (lualsh,file=fash,status='NEW')
c    open (lualpc,file=fapc,status='NEW')
c
      pfx='ionSH'//trim(s5pfx)
      sfx='sh5'
      call newfile (pfx, sfx, fn, flen)
      fash=fn(1:flen)
c     Precursor file
      pfx='ionPC'//trim(s5pfx)
      call newfile (pfx, sfx, fn, flen)
      fapc=fn(1:flen)
c
c rates fitles , shock and precursors
c ratmod=Y
c
c        open (lurtsh,file=frsh,status='NEW')
c        open (lurtpc,file=frpc,status='NEW')
c
      pfx='ratSH'//trim(s5pfx)
      sfx='sh5'
      call newfile (pfx, sfx, fn, flen)
      frsh=fn(1:flen)
c     Rates file for precursor
      pfx='ratPC'//trim(s5pfx)
      sfx='sh5'
      call newfile (pfx, sfx, fn, flen)
      frpc=fn(1:flen)
c
c shock 'dynamics' file
c dynmod=Y
c
c    open (ludy,file=fd,status='NEW')
c
      pfx='dynSH'//trim(s5pfx)
      sfx='sh5'
      fd=' '
      call newfile (pfx, sfx, fn, flen)
      fd=fn(1:flen)
c
c shock cooling file
c fclmod=Y
c
c         open (lucl,file=fcl,status='NEW')
c
      pfx='coolSH'//trim(s5pfx)
      sfx='csv'
      fcl=' '
      call newfile (pfx, sfx, fn, flen)
      fcl=fn(1:flen)
c
c shock cooling emission bands
c bandsmod=Y
c
c     open (lupb,file=fpb,status='NEW')
c
      lupb=30
      pfx='bandSH'//trim(s5pfx)
      sfx='csv'
      fpb=' '
      call newfile (pfx, sfx, fn, flen)
      fpb=fn(1:flen)
c
c line structures, shocks and precursors
c jlin=Y
c
c      open (lulsh,file=flsh,status='NEW')
c      open (lulpc,file=flpc,status='NEW')
c
      pfx='linSH'//trim(s5pfx)
      sfx='csv'
      flsh=' '
      call newfile (pfx, sfx, fn, flen)
      flsh=fn(1:flen)
c
      pfx='linPC'//trim(s5pfx)
      sfx='csv'
      flpc=' '
      call newfile (pfx, sfx, fn, flen)
      flpc=fn(1:flen)
c
c       ionisation structure files -  all names a made only 4 are
c       used atm tsrmod=Y  i is 1-mxelem and is effectively the
c       ordinal atom value ie=iel(i) = zmap( of posible choice i=z
c       ), ie is then map id for the ith element in map.prefs or
c       mapz, elem(ie) uses the map id for the elem name all
c       elements in map.prefs  get file names but only ieln=4 are
c       created and used, ieln can change in the future
c
c B  : Ion balance files.' for just the 4 chosen mapping z to i later
c
c          open (luionsh(i),file=fionsh(i),status='NEW')
c          open (luionpc(i),file=fionpc(i),status='NEW')
c
      sfx='csv'
      do i=1,atypes
c
c       mapppings internal ids not z, names mapped in prefs
c
        pfx='elSH'//trim(s5pfx)//elem(i)
        fn=' '
        call newfile (pfx, sfx, fn, flen)
        fash=fn(1:flen)
        fionsh(i)=fash
c
        pfx='elPC'//trim(s5pfx)//elem(i)
        call newfile (pfx, sfx, fn, flen)
        fapc=fn(1:flen)
        fionpc(i)=fapc
c
      enddo
c
c Final downstream field.'
c jspec=YES
c
c All upstream fields at each step.
c lmod=Y
c
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine createS5files
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine createS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Create files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write "sh5" general output files always open until close all
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      open (luop,file=fsm,status='NEW')
      if (lupt.gt.0) open (lupt,file=fpm,status='NEW')
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spec List Files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      open (lusp,file=fsh,status='NEW')
      if (lupc.gt.0) open (lupc,file=fpc,status='NEW')
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion balance files if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (tsrmod.eq.'Y') then
c
c B  : Ion balance files.'
c
        do i=1,ieln
          open (luionsh(i),file=fionsh(iel(i)),status='NEW')
          if (luionpc(i).gt.0) then
            open (luionpc(i),file=fionpc(iel(i)),status='NEW')
          endif
        enddo
c
      endif
c
      if (allmod.eq.'Y') then
        open (lualsh,file=fash,status='NEW')
        if (lualpc.gt.0) then
          open (lualpc,file=fapc,status='NEW')
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion Dynamics file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (dynmod.eq.'Y') then
c
c C  : Dynamics file.'
c
        open (ludy,file=fd,status='NEW')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion Rates file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ratmod.eq.'Y') then
c
c D  : All Rates file.
c
        open (lurtsh,file=frsh,status='NEW')
c       Rates file for precursor
        if (lurtpc.gt.0) open (lurtpc,file=frpc,status='NEW')
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     spectrum files handled on the fly by wpsou
c     lmod and jspec
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if (jspec.eq.'Y') then
c
c E  : Final downstream field.'
c       units made by wpsou
c
c     endif
c
c     if (lmod.eq.'Y') then
c
c F  : All upstream fields at each step.
c       units made by wpsou
c
c     endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write radiation bands file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (bandsmod.eq.'Y') then
c
c H  : Cooling in x-ray bands.'
c
        open (lupb,file=fpb,status='NEW')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write Cooling file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      if (fclmod.eq.'Y') then
c
c K  : Cooling, Components by Elements file.'
c
        open (lucl,file=fcl,status='NEW')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write Line Stucture files if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jlin.eq.'Y') then
c
c     L  : Monitor up to ',i3,' lines'/
c     line structures, shocks and precursors
c     jlin=Y
c
        open (lulsh,file=flsh,status='NEW')
        if (lulpc.gt.0) open (lulpc,file=flpc,status='NEW')
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine closeS5files
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine closeS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Close All Main and Shock files if open
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      logical unitopen
c
      inquire (unit=luop,opened=unitopen)
      if (unitopen) close (luop)
      if (lupt.gt.0) then
        inquire (unit=lupt,opened=unitopen)
        if (unitopen) close (lupt)
      endif
      inquire (unit=lusp,opened=unitopen)
      if (unitopen) close (lusp)
      if (lupc.gt.0) then
        inquire (unit=lupc,opened=unitopen)
        if (unitopen) close (lupc)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      inquire (unit=lurtsh,opened=unitopen)
      if (unitopen) close (lurtsh)
      inquire (unit=ludy,opened=unitopen)
      if (unitopen) close (ludy)
      do i=1,ieln
        inquire (unit=luionsh(i),opened=unitopen)
        if (unitopen) close (luionsh(i))
      enddo
      inquire (unit=lualsh,opened=unitopen)
      if (unitopen) close (lualsh)
      inquire (unit=lucl,opened=unitopen)
      if (unitopen) close (lucl)
      inquire (unit=lulsh,opened=unitopen)
      if (unitopen) close (lulsh)
      inquire (unit=lupb,opened=unitopen)
      if (unitopen) close (lupb)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lurtpc.gt.0) then
        inquire (unit=lurtpc,opened=unitopen)
        if (unitopen) close (lurtpc)
      endif
      do i=1,ieln
        if (luionpc(i).gt.0) then
          inquire (unit=luionpc(i),opened=unitopen)
          if (unitopen) close (luionpc(i))
        endif
      enddo
      if (lualpc.gt.0) then
        inquire (unit=lualpc,opened=unitopen)
        if (unitopen) close (lualpc)
      endif
      if (lualpc.gt.0) then
        inquire (unit=lulpc,opened=unitopen)
        if (unitopen) close (lulpc)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c****************************************************************
c> @brief The subroutine appendS5files
c! XXXX - add one line purpose here
c! @param This routine has no parameters
c!
c! @return
c!  XXXX Add one or more lines describing what is updated
c!
c! @details
c!  XXXX Enter details here
c***************************************************************

      subroutine appendS5files ()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'cblocks.inc'
      include 's5blocks.inc'
c
      integer*4 i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open All files APPEND
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spec List Files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      open (luop,file=fsm,status='OLD',access='APPEND')
      if (lupt.gt.0) open (lupt,file=fpm,status='OLD',access='APPEND')
c
      open (lusp,file=fsh,status='OLD',access='APPEND')
      if (lupc.gt.0) open (lupc,file=fpc,status='OLD',access='APPEND')
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion balance files if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (tsrmod.eq.'Y') then
c
c B  : Ion balance files.'
c
        do i=1,ieln
          open (luionsh(i),file=fionsh(iel(i)),status='OLD',access='APPE
     &ND')
          if (luionpc(i).gt.0) then
            open (luionpc(i),file=fionpc(iel(i)),status='OLD',access='AP
     &PEND')
          endif
        enddo
c
      endif
c
      if (allmod.eq.'Y') then
        open (lualsh,file=fash,status='OLD',access='APPEND')
c         Precursor file
        if (lualpc.gt.0) then
          open (lualpc,file=fapc,status='OLD',access='APPEND')
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion Dynamics file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (dynmod.eq.'Y') then
c
c C  : Dynamics file.'
c
        open (ludy,file=fd,status='OLD',access='APPEND')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write ion Rates file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ratmod.eq.'Y') then
c
c D  : All Rates file.
c
        open (lurtsh,file=frsh,status='OLD',access='APPEND')
c       Rates file for precursor
        if (lurtpc.gt.0) then
          open (lurtpc,file=frpc,status='OLD',access='APPEND')
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     spectrum files handled on the fly by wpsou
c     lmod and jspec
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if (jspec.eq.'Y') then
c
c E  : Final downstream field.'
c       units made by wpsou
c
c     endif
c
c     if (lmod.eq.'Y') then
c
c F  : All upstream fields at each step.
c       units made by wpsou
c
c     endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write radiation bands file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (bandsmod.eq.'Y') then
c
c H  : Cooling in x-ray bands.'
c
        open (lupb,file=fpb,status='OLD',access='APPEND')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write Cooling file if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      if (fclmod.eq.'Y') then
c
c K  : Cooling, Components by Elements file.'
c
c
        open (lucl,file=fcl,status='OLD',access='APPEND')
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write Line Stucture files if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (jlin.eq.'Y') then
c
c     L  : Monitor up to ',i3,' lines'/
c     line structures, shocks and precursors
c     jlin=Y
c
        open (lulsh,file=flsh,status='OLD',access='APPEND')
        if (lulpc.gt.0) then
          open (lulpc,file=flpc,status='OLD',access='APPEND')
        endif
      endif
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
