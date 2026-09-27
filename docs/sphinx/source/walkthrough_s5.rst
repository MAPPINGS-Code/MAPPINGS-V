.. _walkthrough_s5:

##################################
Walkthrough: An S5 Shock Model
##################################

This page runs one concrete S5 shock model end to end — the actual
prompts answered, the actual files it produces, and what is inside
them — so the whole path from "run a shock model" to "here is my line
flux as a function of position" is in one place. :doc:`models`
(prompt order), :doc:`code_s5` (what the solver does internally), and
:doc:`outputs` (full file/column reference) remain the authoritative
pages for those topics; this page does not restate them, it walks one
worked example and points back to them for the general case.

.. contents:: Contents
   :local:
   :depth: 1

-------------------------------------------------
The scenario
-------------------------------------------------

A 200 km/s, unmagnetised (B = 0) radiative shock into solar-abundance
gas at n\ :sub:`H` = 1 cm\ :sup:`−3`, with MAPPINGS' auto-iterated
photo-ionised precursor and the full internal diffuse radiation field.
B = 0 is the more typical choice for a first shock model — it isolates
the radiative/hydrodynamic physics without also having to pick a
magnetic field strength or geometry, and is a more natural default than
parameterising the field via α (magnetic-to-gas pressure ratio) and
setting α = 1 (equipartition), which is a comparatively specific,
strongly magnetised case. The proto-shock gas is started at 1000 K
(rather than a lower guess), the model runs to a 1%-weighted-ionisation
ending condition rather than a fixed temperature floor, and only the
standard output files are requested — a minimal, representative script
rather than one exercising every optional output table.

-------------------------------------------------
The script
-------------------------------------------------

Piped to ``map52`` non-interactively (see :doc:`inputs`, "Running from
a script"). This is a real, tested script — download
:download:`s5_shock200.mv <../../../lab/examples/s5_shock200.mv>` and
its abundance file,
:download:`s5_shock_abund.txt <../../../lab/examples/s5_shock_abund.txt>`,
and run it as shown below. It also ships with every install, in
``~/mappings520/lab/examples/``.

.. literalinclude:: ../../../lab/examples/s5_shock200.mv
   :language: none

Each of these is one instance of the general prompt sequence documented
in :doc:`popcha` (ionisation balance), :doc:`photsou` (radiation
field — zeroed here with ``X`` immediately, since the shock supplies its
own), and :doc:`models` (S5 jump type/field/velocity/density/stopping
prompts) — this page just shows concrete answers.

A few of these are worth spelling out:

* ``B`` at "Magnetic field parameterisation" asks for the field strength
  directly, in μG; ``0`` gives a purely hydrodynamic, unmagnetised
  shock. See :doc:`models` for the other parameterisations (``A`` for
  α = P\ :sub:`mag`/P\ :sub:`gas`, ``M`` for Alfvén Mach number, etc.)
  if you want a magnetised case instead.
* The velocity/density line (``1000.0 1.0 200.0``) answers "Specify the
  proto-shock conditions: T (K), dh (N), v (km/s)". T ≤ 10 is taken as a
  log, so ``1000.0`` (> 10) is read literally as 1000 K; ``2.0`` would
  give the same thing via 10²=100 K. This is the *proto-shock* gas
  far upstream, distinct from the *pre-shock* state immediately ahead of
  the jump (~9600 K in this run), which the precursor solver computes
  self-consistently regardless of the proto-shock guess — see
  "Does the starting temperature matter?" below.
* ``A`` at "Model Ending Condition" is "Standard ending, 1% weighted
  ionisation" — the default stopping rule, needing no follow-up value
  (unlike ``C``/``S``, which ask for a temperature). See :doc:`code_s5`
  for the other lettered stopping criteria.
* In the Output Multi-Option Menu, ``A`` means "Standard output — and
  Reset" (clears any options set so far), so ``A`` then ``X`` produces
  the plain standard output set with none of the optional tables (ion
  balance files, rates, cooling, etc.) — see :doc:`outputs` for what
  each of the other letters adds.
* The output-file prefix (``v200s_``) is read with a fixed-width
  read, not a free-form one — unlike almost every other prompt in this
  script, MAPPINGS does **not** ignore trailing text on this line, so it
  must contain the prefix alone. A trailing underscore keeps generated
  filenames readable, since MAPPINGS inserts no separator of its own —
  see :doc:`outputs`, "Notes on File Naming".

Run it (from wherever you saved the two downloaded files, or from
``~/mappings520/lab/examples/`` if using the shipped copy)::

    map52 < s5_shock200.mv

-------------------------------------------------
Does the starting temperature matter?
-------------------------------------------------

Starting the proto-shock gas at 1000 K instead of a lower guess (e.g.
100 K) does converge: this run's pre-shock temperature settled to
9566.6 K, post-shock temperature 568 937 K, compression 3.9370. The
precursor solver recomputes the
pre-shock equilibrium state from radiative balance each iteration, so
the proto-shock guess only affects the very first iteration's starting
point, not the converged answer.

It does, however, affect *how many* global iterations are needed: this
run required more than the requested minimum of 3, reaching
"CONVERGED" only after MAPPINGS auto-extended the loop to 5, then ran
one further pass (6 of 6) to write output at the tightened tolerance —
the same behaviour documented in :doc:`code_s5`, "Phase 2c —
Convergence check". A closer initial guess can save a couple of
iterations; a distant one does not prevent convergence, it just costs
a bit more runtime.

-------------------------------------------------
Checking for convergence and other failures
-------------------------------------------------

**A shock model that finishes without a crash has not necessarily
converged.** MAPPINGS never treats non-convergence as a fatal error: if
a run isn't converged within the requested number of iterations, it
silently extends the iteration cap by one and tries again, up to a hard
ceiling of 20 global iterations (see :doc:`outputs`, "Convergence
diagnostics and standard output"). Whatever happens, it then runs a few
more ordinary iterations plus one final pass and writes full output
regardless — same exit code, same file set, no warning anywhere in the
``.sh5``/``.csv`` files themselves. Concretely, always check the
following before trusting a shock result:

1. **Capture stdout.** None of this is recorded in the output files, so
   redirect it when you run the model::

       ./map52 < model.mv | tee run.log

2. **Check the last "SHOCK 5 Convergence Test" block**, not just that
   the program exited cleanly::

       grep "Result:" run.log | tail -1

   It should say ``Result: CONVERGED``, along with the actual
   precursor/post-shock RMS percentages. This line is reliable for a
   model whose main loop genuinely converges — the mandatory final
   independent re-solve no longer overrides an already-established
   result with a spurious disagreement (see :doc:`code_s5`, "How the
   precursor↔shock loop actually converges"). If the main loop never
   converges, the line instead classifies *why* — see :doc:`code_s5`,
   :ref:`s5_convergence_outcomes`, for the full set of eight possible
   outcomes. Not all of them mean the same thing: ``NOT CONVERGED
   (precursor oscillating)`` is a characterised, bounded limit cycle
   validated (via a direct predicted-spectrum comparison) to leave the
   physical output unaffected, and is safe to use as-is; the other
   ``NOT CONVERGED``/``FAILED`` variants are uncharacterised and the
   physical state MAPPINGS wrote out should be treated with more
   suspicion. A result of ``NO SHOCK`` means no compressive shock
   solution exists for these parameters at all (sub-Alfvénic preshock
   flow) — this is not a failure, just a model with no shock to
   compute.

3. **Look at how many iterations it actually took**
   (``grep "Convergence Test" run.log``). A run that converges
   comfortably under the 20-iteration ceiling is on solid ground. A run
   that only reaches "CONVERGED" at or near iteration 20, or that
   prints a long, non-decreasing run of ``Aitken omega`` values that
   never settle down (see :doc:`outputs`), is at or near the edge of
   what the solver can currently handle — worth treating with more
   suspicion and, if practical, rerunning with a less demanding
   parameter combination to see if it settles faster and to the same
   answer. Watching the per-iteration RMS trend is also informative: a
   smoothly, monotonically shrinking RMS is a genuinely converging
   model that just needs more headroom; an RMS that plateaus or bounces
   between values without shrinking further is a real oscillation that
   more iterations will not fix — but see :doc:`code_s5`,
   :ref:`s5_convergence_outcomes`: a bounded oscillation on the
   precursor side alone is now identified and reported as such, and
   does not need a rerun.

4. **Confirm you're reading the actual final output file.** When
   MAPPINGS extends the iteration count beyond what you requested, it
   can start a fresh, higher-numbered output file part-way through,
   leaving the original low-numbered ``.sh5``/``.csv`` files short
   (header only, no "Model ended" line) and the real output in
   ``..._0002.*`` instead. This isn't guaranteed to happen on every run
   that extends — this page's own example run also extended from 3 to
   6 iterations, without splitting the file — but it's worth checking
   for. The reliable check, independent of the convergence question and
   useful as a general "did this run actually finish" test too, is::

       grep -l "Model ended" shck_*.sh5

   A ``.sh5`` file without that line — whether because of the sequence-
   number split above or a genuine mid-run failure — is not the model's
   final state and should not be used.

-------------------------------------------------
What got created
-------------------------------------------------

With no optional output letters selected (just ``A`` then ``X``), this
run produced only the always-written files:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - File
     - Contents
   * - ``shck_v200s_0001.sh5``
     - Post-shock cooling-zone structure (see below).
   * - ``prec_v200s_0001.sh5``
     - Precursor-zone structure, same format.
   * - ``specSHv200s_0001.csv`` / ``specPCv200s_0001.csv``
     - Integrated line list for the shock / precursor (see below).
   * - ``SHupv200s_0001.lam`` / ``.sou``; ``PCupv200s_0001.lam`` / ``.sou``
     - Upstream radiation field snapshots. Written unconditionally —
       not gated by any output-menu option (see :doc:`outputs`).
   * - ``v200s_0001.bln``
     - Final ionisation balance table (not a spectrum — see :doc:`outputs`).

No ``elSH``/``elPC`` (needs ``B``), ``bandSH`` (needs ``H``), ``coolSH``
(needs ``K``), ``linSH``/``linPC`` (needs ``L``), or ``SHdw*`` (needs
``E``) — this is the minimal S5 output set. Every pattern and its full
column layout is catalogued in :doc:`outputs` ("Shock Models — S5").

-------------------------------------------------
A peek inside
-------------------------------------------------

``shck_v200s_0001.sh5`` — after the header, abundance table, and the
Rankine-Hugoniot jump summary (preshock T = 9566.6 K, postshock T =
568 937 K, compression factor 3.937 — close to the strong-shock limit
of 4 for a purely hydrodynamic, unmagnetised shock, since there is no
magnetic pressure resisting compression), one row per downstream step::

    607, 7.568605E+17, 7.556891E+11, 465.576, 93.0247, 6410.99, 6963.32,
    1.24183, 3.916991E+11, 2.422357E+08, 1.362000E-27, 5.168346E-19,
    6.997252E-02, 0.989372, 1.062804E-02, 0.00000

Columns are step, distance (cm), dr (cm), T (K), n\ :sub:`e`,
n\ :sub:`H`, n\ :sub:`i`, μ, elapsed time, dt, normalised net cooling,
total cooling, d\ :sub:`los`, X(H⁰), X(H⁺), B (G) — the full column
reference is in :doc:`outputs`. The B column is 0.00000 throughout, as
expected for this unmagnetised case. By step 607 (of 622) the gas has
cooled from the 568 937 K post-shock plateau to ~466 K and is almost
fully recombined (X(H⁺) ≈ 0.011), because the ``A`` ending condition
(1% weighted ionisation) lets the gas keep cooling well past a fixed
500 K floor. The full run reached 369.8 K before stopping.

``specSHv200s_0001.csv`` — the integrated shock line list. The high
post-shock temperature reaches ionisation stages a slow shock never
does — the list opens with O VII X-ray lines rather than He II/Mg II::

     Lambda(A),  E (eV) , Flux (HB=1.0), Species  , Kind, Accuracy (1-5)
     =======================================================================
          17.396,  7.12717E+02,  1.37964E-07, O  VII   , CC  ,    3
          21.602,  5.73948E+02,  6.57202E-05, O  VII   , CC  ,    3

-------------------------------------------------
Getting line strength as a function of position
-------------------------------------------------

The above files give you the *integrated* structure and the
*integrated* line list, but not individual line strengths at each
position. For that, add output option ``L`` ("Monitor lines"). Starting
from the same minimal script, only the output-menu section changes —
because ``B`` and ``K`` are not selected, none of their follow-up
prompts (element tracking, cooling normalisation) appear, so the
sequence is shorter than it would be with those options active.
Download :download:`s5_shock200_linemonitor.mv
<../../../lab/examples/s5_shock200_linemonitor.mv>` (uses the same
abundance file as above) — it differs from the base script only in the
output-prefix/menu block:

.. literalinclude:: ../../../lab/examples/s5_shock200_linemonitor.mv
   :language: none
   :lines: 16-21

The 11 monitored lines here are a fairly typical optical diagnostic
set: the [O II] 3726/3729 and [S II] 6716/6731 density-sensitive
doublets, [O III] 5007 and [S III] 9069/9531 (ionisation/excitation),
[O I] 6300 and [N II] 6583 (low-ionisation), and Hβ/Hα for the Balmer
decrement — a mix of recombination and collisionally-excited lines
spanning a wide range of critical densities and excitation
temperatures, rather than just the three strongest lines in the
spectrum.

This produces ``linSHv200sl_0001.csv`` (post-shock zone) and
``linPCv200sl_0001.csv`` (precursor zone), one row per step with the
local flux of each selected line::

     # [1] step, ... ,    3726.032[13],    3728.815[14],    4861.333[15],    5006.843[16],    6300.304[17],    6562.819[18],    6583.454[19],    6716.440[20],    6730.816[21],    9068.620[22],    9530.619[23],
        1,  0.0000    , ...,  1.24381E-24,  1.78622E-24,  1.75670E-23,  4.59217E-24,  2.91859E-27,  5.37271E-23,  6.07674E-26,  1.04168E-26,  7.59840E-27,  3.24782E-26,  8.15837E-26,

Verified end to end for this scenario: 286 rows in
``linPCv200sl_0001.csv``, 622 in ``linSHv200sl_0001.csv``, all 11
species correctly identified (confirmed via the setup-time echo: e.g.
``O  II  3726.032``, ``S  III  9530.619``, not defaulted or
misidentified). Every column carries real, position-varying flux
somewhere in the run; how much of the run depends on the physics, not
on the matching — Hβ, Hα, and [O I] are nonzero in all 622 shock rows
(present across the whole temperature range), while lines needing a
specific, narrower ionisation stage — [O III], the [S II]/[S III]
doublets, [N II] — are nonzero in roughly 520–605 of the 622 rows,
zero only where the gas is too hot or too cool for that particular ion
to exist, which is correct physics rather than a matching failure.
Both files share the same selected-line list and header, so a line's
strength can be traced continuously from the precursor through the
shock front and into the cooling zone. See :doc:`physics_output_spectra`
for the underlying mechanism, including a note on when this
position-resolved precursor output was added, and the note below on
how precisely a requested wavelength needs to be specified.

.. note::

   The wavelengths above are given to 3–4 decimal places, matching
   MAPPINGS' own tabulated values — at minimum the same precision
   :doc:`outputs`' line-list files report wavelengths to. Matching
   requires this level of precision: a value that misses by more than
   0.001 Å stops the run immediately with an error listing every real
   line within 1 Å, so a wrong or imprecise wavelength can be
   corrected immediately rather than silently producing a zero-flux
   column. See
   :doc:`physics_output_spectra`, "How a requested wavelength gets
   matched", for why this precision is required and what the error
   looks like.

-------------------------------------------------
Adapting this to your own model
-------------------------------------------------

* **Proto-shock conditions** (``1000.0 1.0 200.0``) — temperature (K,
  or log if ≤ 10), hydrogen density (cm⁻³, or log if ≤ 0), and velocity
  (km/s if < 1e5, else cm/s), in that order.
* **Magnetic parameterisation** (``B`` then ``0``) — the field strength
  in μG; switch to ``A`` for α = P\ :sub:`mag`/P\ :sub:`gas`, or
  ``M``/``C``/``R`` for the other parameterisations (see :doc:`models`).
* **Ending condition** (``A``) — switch to ``C``/``S`` for a temperature
  floor (with a follow-up value), ``D``/``E`` for a distance/time limit,
  or ``B`` for a specific ion-fraction threshold; see :doc:`code_s5` for
  what each does internally.
* **Output options** — add ``B`` for per-element ion balance files,
  ``H`` for X-ray/UV cooling bands, ``K`` for cooling-by-element, ``D``
  for rates/timescales, ``C`` for the dynamics file (note the C/D
  menu-label swap documented in :doc:`outputs`), or ``L`` as shown above
  for position-resolved line monitoring — each adds its own follow-up
  prompts before the output menu can be exited with ``X``.
