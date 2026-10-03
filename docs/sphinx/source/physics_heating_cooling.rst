.. _physics_heating_cooling:

##################################################################
Physical Processes: Heating, Cooling, and Thermal Balance
##################################################################

This page describes how MAPPINGS determines the electron temperature
of a zone — which heating and cooling channels compete, how they are
assembled into a single energy balance, and how that balance is
actually solved (as an equilibrium, or evolved forward in time). It is
about the physics only; the routines here are called from the same
zone-stepping loops documented in :doc:`code_photo` and :doc:`code_s5`.
For the ionisation processes referenced below, see
:doc:`physics_ionization`; for the line emissivities, see
:doc:`physics_lines`; for the dust/PAH terms, see :doc:`physics_dust`.

.. contents:: Contents
   :local:
   :depth: 2

-------------------------------------------------
Overview
-------------------------------------------------

``cool.f`` is the master per-zone assembler: given a trial temperature,
electron density, and hydrogen density, it calls out to every cooling
and heating routine in turn, accumulates two running totals — ``tll``
(total losses) and ``tgg`` (total gains) — and reduces them to a single
normalised imbalance,

.. code-block:: none

   dlos = (ell - egg) / (ell + egg)

where ``ell``/``egg`` are "effective" versions of the totals that
additionally route the recombination net-gain and charge-exchange
terms to whichever side (loss or gain) they actually act on for that
zone, since both can be heating *or* cooling depending on conditions.
``dlos = 0`` means heating balances cooling; the temperature solvers
described further down exist purely to drive ``dlos`` to zero (for
equilibrium) or to integrate the imbalance forward in time (for
non-equilibrium gas).

-------------------------------------------------
Cooling channels
-------------------------------------------------

``cool.f`` calls, in order:

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Channel
     - Routine
     - Notes
   * - H/He recombination + collisional line losses
     - ``hydro`` (``hloss``)
     - Same routine that computes H/He recombination-line brightness
       in :doc:`physics_lines` — line emission and cooling are the
       same energy leaving the gas, computed once.
   * - Heavy-element recombination lines
     - ``recom_cii``/``recom_nii``/``recom_oi``/``recom_oii``/
       ``recom_neii`` (``heavyrec.f``)
     - Noted in the source as a "negligible cooling" contribution but
       included regardless; also the emissivity source for the
       corresponding heavy-element recombination line fluxes.
   * - Fine-structure (3-level atom) cooling
     - ``fine3`` (``f3loss``)
     - Legacy 3-level treatment for simple fine-structure transitions.
   * - Intercombination line cooling
     - ``inter`` (``fslos``)
     -
   * - Forbidden/semi-forbidden (CEL) cooling
     - ``multilevel`` (``fmloss``)
     - The main collisionally-excited-line coolant; same routine
       documented in :doc:`physics_lines` as the optically-thin line
       emissivity source.
   * - Iron-ion cooling
     - ``multiiron`` (``feloss``)
     - Separate from ``multilevel`` owing to the size of the iron
       model ion (up to ~900 levels).
   * - Trapped resonance-line cooling
     - ``reson3``/``resonl`` (``xr3loss``/``xrlloss``)
     - The cooling contribution of the same radiatively-trapped
       resonance lines whose escape-probability transfer is documented
       in :doc:`physics_lines` (``transferline.f``).
   * - Free-free cooling
     - ``frefre`` (``fflos``), in ``freefree.f``
     - Distinct entry point in the same file as ``freefree``, which
       computes the free-free continuum emissivity for the diffuse
       field.
   * - Collisional-ionisation cooling
     - ``coloss`` (``colos``)
     - Each collisional ionisation event removes kinetic energy equal
       to the ion's ionisation potential from the free-electron pool:
       ``cl = col(ion,atom) * ipote(ion,atom) * n_ion * n_e * n_H``,
       using the same rates computed for :doc:`physics_ionization`.
   * - Grain collisional cooling
     - ``hgrains`` (``gcool``)
     - See :doc:`physics_dust`.

-------------------------------------------------
Heating channels
-------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Channel
     - Routine
     - Notes
   * - Photoionisation heating
     - ``pheat`` (``pgain``)
     - Interpolates the tabulated per-ion heating rates ``heaph``
       computed in ``phion.f`` (:doc:`physics_ionization`) as a
       function of the local electron fraction, rather than
       re-integrating the photon field directly.
   * - Grain + PAH photoelectric heating
     - ``hgrains``/``hpahs`` (``gheat``, ``paheat``)
     - See :doc:`physics_dust`.
   * - Charge-exchange heating
     - ``cheat`` (``chgain``)
     - Uses the energy defect of each charge-exchange reaction
       (``chxrde``, in units of k\ :sub:`B`·10\ :sup:`4` K); can be
       net heating or net cooling depending on the reaction, using the
       same rates set up for :doc:`physics_ionization`
       (``allrates``/``charex``).
   * - Cosmic-ray heating
     - ``cosmic`` (``cosgain``)
     - The same routine documented as a secondary-ionisation source in
       :doc:`physics_ionization` — cosmic rays both ionise and heat.
   * - Recombination net-gain (on-the-spot heating vs. recombination
       cooling)
     - ``netgain`` + ``spotap`` (``rngain``)
     - See below.

Recombination net-gain and the on-the-spot approximation
============================================================

``netgain.f`` computes the balance between the cooling from
free electrons recombining (losing their kinetic energy) and, when the
on-the-spot approximation is enabled (``jspot = 'Y'``), an offsetting
heating term from ``spotap`` representing photons produced by
recombinations to excited states being immediately reabsorbed rather
than transported — the same on-the-spot assumption discussed in
:doc:`physics_ionization`'s recombination section. The result
(``rngain``) follows **Seaton (1959, MNRAS 119, 81)** and
**Hummer & Seaton (1964, MNRAS 127, 230)**, and can be positive or
negative — ``cool.f`` routes it to whichever side of the energy budget
it actually belongs on for that zone (see Overview).

-------------------------------------------------
Dormant heating/cooling channels
-------------------------------------------------

Two fully-implemented channels are never invoked anywhere in the
current source (verified by an exhaustive grep for their call sites).
As with the Bowen-fluorescence and precursor-line-monitor gaps
documented elsewhere (:doc:`physics_lines`, :doc:`physics_output_spectra`),
these would need source changes to become active; no run configuration
currently activates them.

Compton heating/cooling
==========================

``compton.f`` is complete and self-contained: it implements the
standard non-relativistic Compton heating/cooling balance following
**Krolik, McKee & Tarter (1981)**, using the Thomson cross-section
(σ\ :sub:`T` = 6.6524×10\ :sup:`-25` cm²) and the mean photon energy of
the local field ``tphot``:

.. code-block:: none

   cmplos = (sigma_T * F / (m_e c^2)) * n_e * (4kT - <h*nu>)

(Compton cooling scales with 4kT, Compton heating with the mean photon
energy; ``cmplos`` is the net of the two.) It even has a built-in
``expertmode`` diagnostic print. The routine ``compton(t, de)`` is
simply never called — not from ``cool.f``, not from anywhere else —
with no gating flag or input mode to check first. It reads as code
that was disconnected from the call graph at some point, rather than
a feature deliberately left unfinished.

``timtqui`` — a distinct regime, not a redundant copy of ``teequi2``
========================================================================

Unlike Compton above, ``timtqui.f`` is not simply disconnected
code — it targets a genuinely different physical question from the
two active temperature solvers. All three use the same
arctanh(dlos)=B+A·ln(T) iterative scheme, but differ in what ionisation
state they hold the trial temperature against:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Solver
     - Question being solved
   * - ``teequi2`` (active)
     - What T and ionisation state are mutually consistent *at
       equilibrium* (calls ``equion`` at each trial T)?
   * - ``evoltem`` (active)
     - How do T and ionisation evolve *together, forward in real
       time* (calls ``timion`` once, over the actual elapsed time)?
   * - ``timtqui`` (unused)
     - If ionisation is whatever a fixed elapsed time ``tstep`` of
       time-dependent evolution produces, what T is thermally
       self-consistent with *that* — resetting the population to its
       initial state and re-calling ``timion`` at every trial T?

That third question is a legitimate quasi-equilibrium regime (useful
when ionisation lags temperature and a self-consistent instantaneous
temperature is wanted for a partially-relaxed, not fully time-evolved,
ionisation snapshot) — it just has no caller anywhere in the current
source.

-------------------------------------------------
A simplified alternative: parametrised cooling
-------------------------------------------------

``cool.f`` also supports bypassing the entire physical heating/cooling
calculation above via ``alphacoolmode``. When enabled, the total
cooling rate is instead a pure power law,

.. code-block:: none

   tloss = de^2 * alphac0 * (T / 1e6 K)^alphaclaw
   tgain = 0 ; egain = 0 ; eloss = tloss

with no heating term at all — useful for controlled numerical
experiments (e.g. testing the hydrodynamics in isolation from the
detailed atomic cooling physics) rather than a physically realistic
model.

-------------------------------------------------
How temperature is actually solved for
-------------------------------------------------

Three different routines turn the ``cool``-computed imbalance into an
actual temperature, depending on whether the gas is assumed to be in
equilibrium or evolving in time:

``teequi2`` — the modern equilibrium/finite-age solver
===========================================================

Used by P6/P7 for equilibrium and finite-age zones, and by both shock
models for the initial equilibrium state (``shock4.f``, ``shock5.f`` —
see :doc:`physics_shocks`). It alternates trial calls to ``equion``
(ionisation balance at the trial temperature) and ``cool`` (net
imbalance ``dlos`` at that temperature and ionisation state), and
fits

.. code-block:: none

   arctanh(dlos) = B + A * ln(T)

to successive trial points, solving for the temperature where this
predicts ``dlos = 0`` — an accelerated secant-like iteration in
log-temperature space rather than a plain bisection.

``teequi`` — the older equilibrium solver
=============================================

A simpler predecessor to ``teequi2``, still actively used by
``phocrv.f`` (a standalone cooling-curve utility), by ``preion.f`` (the
shock precursor's initial-condition calculation — see
:doc:`physics_shocks`), and by the single-zone models (``sinsla.f``,
``slab.f``).

``evoltem`` — forward time-dependent evolution
==================================================

Rather than solving for equilibrium, ``evoltem`` integrates
temperature and ionic populations forward together over an actual
elapsed time (calling ``timion`` for the ionisation side), for gas that
has not had time to reach equilibrium — used by P6/P7's finite-age and
post-equilibrium modes. It exits early if the electron fraction or
temperature drop below supplied thresholds (e.g. once the gas has
effectively recombined to neutral).

-------------------------------------------------
Data flow summary
-------------------------------------------------

.. code-block:: none

   Trial T, n_e, n_H at a zone
           |
           |  cool.f calls, per channel:
           |    hydro, heavyrec, fine3, inter, multilevel, multiiron,
           |    reson3/resonl, frefre, coloss, hgrains        -> tll
           |    netgain(+spotap), pheat, hgrains/hpahs, cheat,
           |    cosmic                                        -> tgg
           v
   dlos = (ell - egg) / (ell + egg)
           |
           |  dlos = 0 sought by:                dlos evolved forward by:
           |    teequi2 (P6/P7, both shock          evoltem (P6/P7 finite-age /
           |    models' initial state) or            post-equilibrium; integrates
           |    teequi (single-zone models,          T + ionic populations together
           |    phocrv, shock precursor init)         via timion over elapsed time)
           v                                     v
   Converged equilibrium T, ionisation state      Time-evolved T, ionisation state

-------------------------------------------------
Key routines
-------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Routine
     - Role
   * - ``cool``
     - Assembles every heating and cooling channel into ``tll``/``tgg``
       and the normalised imbalance ``dlos``.
   * - ``netgain``
     - Recombination cooling vs. on-the-spot heating balance
       (``rngain``); calls ``spotap``.
   * - ``pheat``
     - Turns ``phion.f``'s tabulated ``heaph`` rates into the actual
       photoionisation heating term.
   * - ``cheat``
     - Charge-exchange heating/cooling from reaction energy defects.
   * - ``coloss``
     - Collisional-ionisation cooling.
   * - ``teequi2``
     - Modern equilibrium/finite-age temperature solver; iterative fit
       in log-T / arctanh(``dlos``) space.
   * - ``teequi``
     - Older equilibrium solver, still used by single-zone models,
       ``phocrv.f``, and the shock precursor's initial conditions.
   * - ``evoltem``
     - Forward time-dependent evolution of temperature and ionic
       populations together, via ``timion``.
   * - ``compton``
     - Non-relativistic net Compton heating (Krolik, McKee & Tarter
       1981) — implemented, currently unused (see Developer Note).

.. todo :: Document ``timion``'s ODE-stepping mechanics itself
   (``sdifeq.f``, ``spotap.f``, ``ionab.f``, ``poputil.f``) — this page
   covers what temperature solver calls it and why, not how the
   non-equilibrium ionisation integration itself works.
