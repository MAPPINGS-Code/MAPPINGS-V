.. _code_s5:

##################################
Code Operation: S5 Shock Model
##################################

This page describes how the S5 shock model operates internally.  For the
user-facing inputs and outputs see :doc:`models` and :doc:`outputs`, or
:doc:`walkthrough_s5` for a full worked example.  For the broader code
structure see :doc:`code_overview`.  For how the radiative processes used
here compare to the photoionisation models, see :doc:`physics_shocks`.

The S5 model computes a steady-state, plane-parallel, magnetised radiative
shock wave and its photo-ionised precursor.  The implementation lives in
``src/shock5.f``, which is a single file of ~4600 lines containing the
main dispatcher and all supporting subroutines.

.. contents:: Contents
   :local:
   :depth: 2

-------------------------------------------------
Top-level structure
-------------------------------------------------

The entry point is ``subroutine shock5()``, ~170 lines.  It orchestrates
four phases:

1. **Setup** — ``shock5setup`` collects all user input and computes the
   initial Rankine-Hugoniot jump conditions.
2. **Sub-Alfvénic check** — before any shock jump is attempted, the
   preshock Alfvén Mach number (computed in ``shock5setup`` from preshock
   quantities alone) is checked.  If it is below 1, no compressive fast
   MHD shock solution exists for these parameters at all, and the model
   exits cleanly with a "No Shock" message rather than attempting a jump
   calculation that has no physical root — see :ref:`s5_no_shock` below.
3. **Iteration loop** — alternates between integrating the post-shock
   cooling zone and recomputing the precursor, repeating until convergence.
   See :ref:`s5_convergence` for how the loop actually converges — the
   iteration count alone is not the whole story.
4. **Final output pass** — once convergence is declared, a few more
   ordinary iterations run to let the state settle further, then one
   last pass with the ``finalit`` flag set writes all output files at
   full detail.

The default minimum number of global iterations is 3.  The loop
continues until the convergence flag is set or ``mxshockits`` (20) is
reached.

-------------------------------------------------
Phase 1 — Setup (``shock5setup``)
-------------------------------------------------

``shock5setup`` is the longest routine (~1000 lines).  It handles all
user interaction and prepares the initial conditions.

1. Calls :doc:`popcha` — the user sets the initial ionisation balance
   of the upstream (proto-shock) gas.  The result is copied into four
   population arrays: ``pop_neu`` (neutral/proto-shock), ``pop_pre``
   (precursor), and their saved copies ``pop0`` and ``pop_pre0``.
2. Calls :doc:`photsou` — the user defines any ambient radiation field
   (can be set to zero for a pure shock).  The field is stored as
   ``prefield`` and added to the shock's own radiation in later steps.
3. Prompts for the shock parameters in order:

   - Shock jump type: velocity ``V`` (most common) or post-shock
     temperature ``T``
   - Magnetic field parameterisation: ``B``, ``A``, ``M``, ``C``, or
     ``R`` (see :doc:`models`)
   - Pre-shock velocity and hydrogen density (or pre- and post-shock
     temperature for ``T`` mode)
   - Diffuse field mode: full or zeroed
   - Stopping condition and limit value
   - Minimum number of iterations
   - Output file prefix and which optional output tables to write

4. Calls ``shockcmpf`` — applies the **Rankine-Hugoniot jump
   conditions** to compute the immediate post-shock state (T\ :sub:`1`,
   n\ :sub:`1`, v\ :sub:`1`, B\ :sub:`1`) from the pre-shock state.
   This becomes the starting point for the cooling zone integration.
5. Calls ``shocksummary`` to print the pre- and post-shock state to
   the terminal and output files, including the preshock Alfvén Mach
   number and the ratio of magnetic to gas pressure.

.. _s5_no_shock:

-------------------------------------------------
No shock: sub-Alfvénic preshock flow
-------------------------------------------------

A compressive fast MHD shock only has a solution when the preshock flow
speed exceeds the Alfvén speed.  ``shock5()`` checks this immediately
after ``shock5setup`` returns, using the Alfvén Mach number computed
there from preshock quantities alone — before ``shockcmpf`` (the jump
solver) is ever called.

If the Alfvén Mach number is below 1, ``shockcmpf``'s jump-condition
quadratic has no physical root: it returns an unphysical compression
factor (an *expansion*, not a compression) that would otherwise
propagate through to a negative temperature and a hard stop deep inside
the cooling-zone integration, many steps later, in ``cool()``.  Rather
than let that happen, MAPPINGS detects the condition upfront and exits
the model cleanly:

.. code-block:: text

   ********************************************************
    SHOCK 5: NO SHOCK -- preshock flow is sub-Alfvenic
    Alfven Mach Number =  0.5341     (< 1): no compressive
    MHD shock jump exists for these parameters.  Skipping
    the shock calculation for this model.
   ********************************************************

    No Shock: , Alfven Mach:,  0.5341     , Reason:, sub-Alfvenic preshock flow

    Result: NO SHOCK (sub-Alfvenic)

The ``No Shock:`` line is written in the same comma-separated style as
``Model ended:`` (see :doc:`outputs`) so it can be parsed the same way,
and the ``Result:`` line means ``grep "Result:"`` behaves consistently
across every outcome — see :ref:`s5_convergence_outcomes` for the full
set, eight in total.  No ``.sh5``/``.csv`` structure files are produced
for a no-shock model, since no shock structure exists to write.

-------------------------------------------------
Phase 2a — Cooling zone integration (``compsh5``)
-------------------------------------------------

``compsh5`` integrates the post-shock cooling zone zone by zone,
marching downstream (away from the shock front).  Before the main loop
it writes three header rows to the output file recording the proto-shock
state (step −2), the precursor state (step −1), and the immediate
post-shock state (step 0).

The main step loop then advances from step 1 onward.  Each step:

Step 1 — Save initial state
============================

The temperature, density, ionisation populations (``pop0``), distance,
velocity, and magnetic field at the start of the step are saved.  These
will be used to restart the step if needed.

Step 2 — Choose the timestep
==============================

``fdynamictimestep`` returns the smallest of the atomic recombination
timescale, the cooling timescale, and the photon absorption timescale.
This is then multiplied by a scale factor (0.04 × step²) that ramps up
slowly for the first five steps, preventing instabilities immediately
behind the shock front, before reaching the full adaptive step size.

Step 3 — Half-step predictor
==============================

A half-step (``hdt = dt/2``) is taken to obtain mid-step estimates:

- ``timion`` advances the ionisation balance by ``hdt`` at the initial
  temperature.
- ``rankhug`` integrates the MHD flow equations (Rankine-Hugoniot in
  differential form) by ``hdt`` to get a mid-step temperature.

Step 4 — Full-step corrector
==============================

The mid-step net cooling rate is evaluated with ``flocallosses`` using
the improved mid-step ionisation state.  Then:

- ``rankhug`` integrates the full flow step ``dt`` from the start of
  the step using the mid-point loss rate (a second-order
  predictor-corrector scheme).
- ``timion`` advances the ionisation balance over the full ``dt`` at
  the mean temperature of the step.

Step 5 — Record the step
=========================

Start and end values are averaged to give zone-centred quantities:
temperature ``te``, electron density ``deel``, hydrogen density ``dhy``,
magnetic field ``bmg``, flow velocity ``veloc``.  Cumulative distance
``dist`` and elapsed time ``timlps`` are updated:

- ``dist`` (cm) is the downstream distance of the current zone from the
  shock front (``dist(1) = 0`` at the front), accumulated by adding
  each zone's width ``dr`` (``shock5.f:3673``).  This is what the
  ``D`` ending condition compares against.
- ``timlps`` (s) is the elapsed **flow time** of the shocked gas
  parcel since it crossed the front — a Lagrangian clock carried with
  the gas, derived from the same per-step ``dt`` used to integrate the
  ionisation and flow equations (``dr = dt * vel1``,
  ``shock5.f:3379,3674``) — not wall-clock runtime of the MAPPINGS
  process itself.  This is what the ``E`` ending condition compares
  against.

Structural markers are recorded at the distances where the temperature
passes through 10\ :sup:`7`, 10\ :sup:`6`, 10\ :sup:`5`, …, 10\ :sup:`2` K.

If ``finalit`` is set, the step is written to the main structure file
(``.sh5``) and to any optional output tables (rates, dynamics, cooling,
ion fractions, line monitor) that were selected at setup.

Step 6 — Check stopping condition
===================================

This condition (``jend``) and its limit value are chosen by the user
during setup (see :doc:`models`, S5 step 7, for the input prompts as
the user sees them).  Each step, the chosen criterion is tested against
the current zone's state:

- ``A`` — 1% weighted ionisation fraction (default)
- ``B`` — specific ion fraction below a threshold
- ``C`` or ``S`` — temperature below a limit (``S`` also requires >95%
  neutral)
- ``D`` — distance limit ``diend`` (cm downstream of the shock front) reached
- ``E`` — elapsed flow-time limit ``timend`` (s of Lagrangian post-shock
  flow time, not wall-clock runtime) reached
- ``F`` — thermal balance reached
- ``G`` — heating limit

When the condition is met the loop exits and the emission spectrum,
line fluxes, and other integrated quantities are assembled and written.

Beyond the user-chosen ``jend`` letter, three further conditions can also
end the zone loop early, checked every step regardless of which stopping
criterion was selected:

- **Hard zone cap** — the loop only continues while ``step`` is below
  ``mxnsteps`` (4096, set in ``const.inc``) and the temperature is above
  100 K.  Reaching either limit silently ends the model even if the
  chosen ``jend`` condition was never satisfied; a model that hits this
  cap has not reached its intended stopping condition and should be
  treated as incomplete.
- **``terminate`` poll file** — if a file literally named ``terminate``
  exists in the run directory when a step completes, the model stops
  immediately, independent of ``jend``.  This offers a way to abort a
  long-running shock model in place without killing the process.
- **First-iteration early exit** (S5 only) — on the *first* precursor
  iteration, if more than one iteration was requested (``maxits`` > 1),
  the zone loop stops as soon as ``dlos`` (the thermal-balance residual)
  drops below 0.5, well short of the requested ``jend`` condition.  This
  gives ``shock5precursor`` a cheap first-pass structure to seed the
  precursor calculation; the full ``jend`` condition is honoured on
  subsequent iterations once the precursor loop is iterating for real.

A consequence for ``D`` and ``E`` specifically: because the hard zone
cap and 100 K temperature floor are checked in the same statement that
otherwise allows the loop to continue (``shock4.f:1957``,
``shock5.f:4037``), setting ``diend``/``timend`` to a value the model
will never physically reach does not produce an infinite or unbounded
run.  The model instead terminates on whichever of the two overriding
conditions is met first — typically the temperature floor, once the
post-shock gas has radiatively cooled to 100 K, or (for a shock whose
post-shock gas is held near some higher equilibrium and never cools
that far) the 4096-zone cap — with ``jend`` itself never having been
satisfied.  Since this whole zone loop is re-run on every global
precursor iteration (up to ``mxshockits``, 20), an unreachable
``diend``/``timend`` multiplies this wasted cost by however many
iterations the outer loop takes.

Key physics routines called each step
======================================

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Routine
     - Role
   * - ``rankhug``
     - Solves the Rankine-Hugoniot MHD conservation equations in
       differential form; advances T, ρ, v, B given the net cooling rate.
       When magnetic pressure strongly dominates gas pressure, the
       pressure-balance step can be numerically ill-conditioned enough
       that a tiny, physically-correct compression drives the result
       non-physical; if that happens, ``rankhug`` retries with a
       reduced effective cooling term (leaving the timestep itself
       untouched) until the result is physical, printing a notice when
       it does so.
   * - ``timion``
     - Time-dependent multi-species ionisation balance solver; evolves
       all ion fractions over a timestep at a given temperature and density.
   * - ``flocallosses``
     - Net cooling rate = total radiative losses − photoionisation
       heating − all other heating terms.
   * - ``localem`` / ``totphot2``
     - Compute the local emissivity and the radiation field at the
       current position; ``totphot2`` propagates photons in the
       downstream (``DW``) direction.
   * - ``averinto``
     - Averages start and end ionisation populations to give the
       zone-centred mean.


-------------------------------------------------------
Phase 2b — Precursor computation (``shock5precursor``)
-------------------------------------------------------

After each cooling-zone integration, ``shock5precursor`` computes the
**photo-ionised precursor** — the upstream gas that is pre-ionised by
UV and X-ray radiation escaping forward from the shock.

1. Takes the accumulated downstream radiation field from ``compsh5``
   and calls ``totphot2`` in upstream (``UP``) mode to propagate it
   into the neutral pre-shock gas.  A factor of 0.5 is applied to the
   source field to account for the half-space geometry: only half the
   photons leaving the shock face travel in the upstream direction.
2. Computes the **precursor parameter** ψ = Q(ions) / v\ :sub:`shock`
   — the ratio of the ionising photon flux to the shock velocity.
   ψ ≫ 1 indicates a fully radiation-dominated precursor; ψ < 1
   means only partial pre-ionisation.
3. Integrates the precursor zone, evolving the ionisation state of the
   upstream gas step by step using the same ``timion`` + ``rankhug``
   machinery as the cooling zone.
4. Updates ``pop_pre`` (the precursor ionisation balance) and the
   precursor thermodynamic state (``te_pre``, ``de_pre``, ``dh_pre``,
   ``bm_pre``), which become the new pre-shock boundary condition for
   the next ``compsh5`` call.

.. _s5_convergence:

How the precursor↔shock loop actually converges
=================================================

The precursor calculation above and the cooling-zone integration
(``compsh5``) depend on each other: the precursor's radiation field
comes from the shock's own downstream emission, and the shock's
preshock boundary condition comes from the precursor's ionisation
state.  Global iteration alternates between the two until they agree.
Two things make that agreement reliable rather than a raw
fixed-point substitution, which for some shocks can settle into a
non-decaying oscillation instead of converging:

**Consistent precision throughout.**  The precursor's own inner
zone-stepping loop always solves to the same tolerance
(``s5rmstol``, 0.01%) that the outer convergence check demands, from
the very first global iteration — not just on the final pass.  A
loose inner tolerance held only until the very end can let the outer
loop believe it has converged when the underlying precursor solve was
never actually self-consistent to that precision.

**Aitken Δ² dynamic relaxation.**  Rather than accepting each new
precursor solve outright (which can oscillate) or blending it with a
fixed, guessed fraction of the previous iteration (which helps some
shocks and measurably hurts others whose coupling was already
well-behaved), the relaxation weight applied each iteration is
computed from how the last two iterations actually moved — the
standard technique for this kind of black-box partitioned coupling,
also used for e.g. fluid-structure interaction.  It is printed each
iteration as ``Aitken omega`` in the convergence-test block below.

Once the loop first reports convergence, a few more ordinary
iterations run to let the state settle further, then a final pass
(``finalit=1``) writes all output tables.  That final pass is still
one *more*, independent re-solve of the precursor, and can show a
small residual against the immediately-preceding converged state —
this is the relaxation scheme's own noise floor from re-solving an
already-converged point, not evidence the model is actually
unconverged.  If the loop had already genuinely converged before this
one extra call, ``shock5check`` reports the true prior result rather
than letting one noisy independent re-solve flip the final line to
``NOT CONVERGED``.  The raw numbers from that last check are still
printed either way — only the bottom-line ``Result:`` is corrected. A
model whose main loop never actually converges is not affected by
this and still correctly reports one of the outcomes described in
:ref:`s5_convergence_outcomes`.

**Persistent precursor oscillation.**  For a small minority of models
— empirically always precursor-dominated, and typically low-to-mid
shock velocity — the precursor↔shock coupling settles into a genuine,
verified period-2 limit cycle rather than converging: the precursor's
boundary state (Ψ, T\ :sub:`pre`, n\ :sub:`e,pre`) alternates between
two nearby values indefinitely instead of approaching a single fixed
point. This is a real eigenvalue-crossing instability in the discrete
precursor↔shock map itself, not an artifact of the Aitken relaxation
scheme — bypassing Aitken entirely (raw fixed-point substitution)
shows the same period-2 cycle, just larger in amplitude. It traces to
a mesh↔ionisation feedback loop concentrated in the precursor's
He-ionisation-transition zones: the per-zone mesh width array
(``predr``) persists across outer iterations (a ``COMMON`` array, not
reset each pass), so each iteration's mesh adapts to the *previous*
iteration's ionisation state, closing a feedback loop with no inertia.

A direct comparison of the two alternating states' predicted precursor
emission-line spectra (using the line-summary tooling in
``py_progs/mappings``) shows the oscillation doesn't matter for the
actual output: every diagnostic line agreed to well under 0.5% between
the two states, consistently smaller than the genuine spectral
difference between adjacent, cleanly-converged grid points one
velocity step apart.

Given that, MAPPINGS does not attempt to force this class of model to
converge. Instead ``shock5check`` classifies *why* the combined test
failed, and a distinct, honest ``Result:`` line is reported — see
:ref:`s5_convergence_outcomes` below.

.. note::

   **Radiation transport in the precursor: absorption only.**
   Within the precursor zone the radiation field at each upstream step
   is computed as:

   .. code-block:: none

      tphot(ν) = src(ν) × attcol(ν)

   where ``src`` is the field at the shock face and ``attcol`` is the
   Beer-Lambert photoelectric attenuation through the integrated column
   ``popintfr`` (ion populations × dr, accumulated from the shock face
   outward).  **No geometrical dilution factor is applied.**  The
   ``fdilu`` function, which computes the proper dilution factor
   W = ½(1 − √(1 − (r\ :sub:`source`/r)\ :sup:`2`)) and is used by
   the photoionisation models for spherical geometry, is present in the
   code but explicitly commented out in ``shock5precursor``.

   The plane-parallel shock geometry makes this self-consistent for an
   infinite planar shock sweeping through uniform gas: in that limit a
   uniform emitting sheet produces a constant-intensity field at all
   upstream distances and there is no beam spreading.  However,
   geometrical dilution can be significant in practice.  When the
   emitting region (the cooling zone) subtends a small solid angle as
   seen from a distant precursor zone — for instance in fast shocks
   with compact cooling zones, or when modelling the bow shock around
   an isolated object — the 1/r\ :sup:`2` fall-off of the radiation
   field is neglected.  This leads to an overestimate of the ionising
   flux reaching the outer precursor, and therefore an overestimate of
   the degree of pre-ionisation at large upstream distances.

-------------------------------------------------
Phase 2c — Convergence check (``shock5check``)
-------------------------------------------------

After each ``shock5precursor`` + ``compsh5`` pair, ``shock5check``
computes the RMS fractional change, relative to the previous
iteration, of six quantities: the ionisation parameter Ψ (Q/v), the
compression factor, the precursor temperature, the post-shock
temperature, the precursor electron density, and the H/He ion
fraction change.  If the combined RMS is below ``s5rmstol`` (0.01%)
the flag ``converged`` is set, which exits the iteration loop.
Otherwise another iteration begins, up to ``mxshockits`` (20) global
iterations.

Two further diagnostics are computed alongside the six raw quantities:
a split of the same six terms into a **precursor** (pre-shock: Ψ,
T\ :sub:`pre`, n\ :sub:`e,pre`, ΔH/He) sub-RMS (``precrms``) and a
**post-shock** (compression, T\ :sub:`shock`) sub-RMS (``postrms``),
showing which side of the shock front a lingering residual actually
comes from (in practice, almost always the precursor); and the Aitken
relaxation weight ``ω`` used that iteration (see :ref:`s5_convergence`
above). A separate check tracks whether the precursor's own inner
zone-stepping loop exhausts its sweep budget without reaching its own
self-consistency target — a different, nested convergence test from
the outer coupling, which (see below) takes priority over it: if the
inner solve itself never settled, the outer ``precrms``/``postrms``
split isn't a meaningful reading regardless of what it shows.

.. _s5_convergence_outcomes:

Outcome classification
=================================================

When the main loop finishes — either because it converged, or because
it exhausted ``mxshockits`` — MAPPINGS reports one of eight mutually
exclusive, jointly exhaustive outcomes on the bottom-line ``Result:``
line. Every outcome also carries the actual precursor/post-shock RMS
percentages, so the classification is self-verifiable without digging
through the log.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - ``Result:`` line
     - Meaning
   * - ``NO SHOCK (sub-Alfvenic)``
     - No compressive shock solution exists — see :ref:`s5_no_shock`.
   * - ``CONVERGED``
     - Ordinary clean convergence.
   * - ``CONVERGED (last independent recheck was noisy; already
       converged earlier)``
     - Converged, but the mandatory final independent re-solve showed
       a small residual on its own — Aitken's noise floor, not a real
       disagreement (see :ref:`s5_convergence` above).
   * - ``NOT CONVERGED (precursor oscillating)``
     - The characterised, bounded period-2 limit cycle described in
       :ref:`s5_convergence` — ``postrms`` is tight, ``precrms`` is
       within the spectrum-validated 0.5% ceiling, and the precursor's
       inner solve settled normally. Safe to use as-is; not an
       ordinary convergence, but not a problem either.
   * - ``NOT CONVERGED (precursor not converged)``
     - Post-shock side and inner solve are fine, but ``precrms``
       exceeds the validated ceiling — an uncharacterised precursor
       problem.
   * - ``NOT CONVERGED (post-shock not converged, precursor OK)``
     - The precursor coupling is fine, but the post-shock side itself
       has not converged — never yet observed in practice, but now
       distinguishable from a precursor problem rather than folded
       into the same bucket.
   * - ``NOT CONVERGED (precursor and post-shock not converged)``
     - Neither side converged.
   * - ``FAILED (precursor inner solve exhausted, unreliable)``
     - The precursor's own inner zone-stepping loop hit its sweep
       budget without reaching self-consistency. Takes priority over
       the five outcomes above, since the outer split isn't a
       meaningful reading if the inner solve itself never settled.

For example, a model in the characterised oscillating regime reports:

.. code-block:: text

    Result: NOT CONVERGED (precursor oscillating; precursor RMS  3.2138E-02%, post-shock RMS  1.3283E-02%)
    Bounded precursor limit cycle -- post-shock and inner solve OK.

The post-shock ceiling used above (0.05%) is deliberately looser than
the tight ``s5rmstol`` (0.01%) used for the combined test: post-shock
quantities are computed from the precursor's own boundary state, so a
purely precursor-side oscillation shows up as a small correlated
wobble on the post-shock side too, and gating on the tight threshold
directly misclassified known precursor-only cases as post-shock
problems during validation. 0.05% is a provisional margin above the
0.005%–0.022% actually observed across the known oscillating models,
not an independently validated bound the way the 0.5% precursor
ceiling is.

When the loop exits, a few more ordinary iterations run (see
:ref:`s5_convergence`), then one final pass runs with ``finalit=1``
so that all output tables are written at full detail.

-------------------------------------------------
Data flow summary
-------------------------------------------------

.. code-block:: none

   proto-shock gas (pop_neu, T, nH, v, B)
           |
           | shock5setup: preshock Alfven Mach number
           v
   Alfven Mach < 1?  --yes-->  No Shock: exit cleanly (no jump exists)
           |
           no
           v
           | shockcmpf (Rankine-Hugoniot jump)
           v
   post-shock state (T1, n1, v1, B1)
           |
           | compsh5: step-by-step cooling zone integration
           |   each step:
           |     rankhug  (MHD flow)
           |     timion   (ionisation balance)
           |     flocallosses (net cooling)
           v
   shock structure + emission spectrum
           |
           | radiation field escaping upstream
           v
   shock5precursor: ionises upstream gas
           |
           v
   updated pop_pre → new pre-shock state
   (Aitken-relaxed toward the previous iteration's state)
           |
           +----( iterate until converged, s5rmstol throughout )----+
           |
           v
   a few more ordinary iterations to settle further
           |
           v
   final pass (finalit=1) → write all output files
