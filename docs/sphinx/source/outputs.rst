.. _outputs:

###############
Output Files
###############

MAPPINGS produces a variety of output files depending on the model type run.
All file names are constructed from a **prefix** and a zero-padded four-digit
**sequence number** (e.g. ``photn0001.ph6``, ``photn0002.ph6``, …), so that
multiple models in a single run produce distinct files.  For shock models the
user is also prompted for an additional string that is embedded in the prefix
(e.g. ``shck_v100s_0001.sh5``).

The file extensions group files by model type and content as described in the
sections below.  For what the line-flux and continuum-spectrum numbers in
these files physically represent (a volume-integrated total vs. a
radiatively-transferred emergent field), see :doc:`physics_output_spectra`.

.. contents:: Contents
   :local:
   :depth: 2

---------------------------------
Photoionization Models
---------------------------------

Photoionization models are run by the ``photo6`` (mode 6) and ``photo7``
(mode 7) solvers and produce files with the extensions ``.ph6`` and ``.ph7``
respectively, together with a set of ``.csv`` files and optional photon source
(``.sou``) and spectrum (``.lam``) files.  For a worked example tracing one
concrete run's choices to the files below, see :doc:`walkthrough_p6`.

Main Output Files
=================

These files are always written for each photoionization model.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Filename pattern
     - Contents
   * - ``photn<N>.ph6`` / ``photn<N>.ph7``
     - Primary output.  Contains the full nebular structure zone by zone,
       emission line fluxes, ionic fractions, temperatures, densities, and
       an integrated line list.
   * - ``phapn<N>.ph6`` / ``phapn<N>.ph7``
     - Aperture-summed photon output, integrating the emission over a
       specified projected aperture.
   * - ``phlss<N>.ph6`` / ``phlss<N>.ph7``
     - Line and spectrum summary; a condensed form of the line list and
       broadband fluxes suitable for quick inspection.
   * - ``phsem<N>.ph6`` / ``phsem<N>.ph7``
     - Spectral energy map; the full spectral energy distribution resolved
       across the nebular zones.
   * - ``ions_<N>.ph6`` / ``ions_<N>.ph7``
     - Ion fraction profiles as a function of depth through the nebula for
       all tracked species.

Optional Tabular Output Files
==============================

These ``.csv`` files are written when the corresponding output option is
selected in the run menu.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Filename pattern
     - Contents
   * - ``spec<N>.csv``
     - The full discrete line list (via ``spec2``), one row per computed
       emission line: wavelength (Å, air), photon energy (eV), flux
       relative to H\ :math:`\beta` = 1, species (atom + ion stage),
       line-type code, and an index — six columns, not two. Same format
       as ``specSH<pfx>``/``specPC<pfx>`` (S5) and ``lines<N>.csv``
       below; this is a volume-integrated total with no wavelength
       binning/resolution limit, unlike ``.lam``/``.sou``/``.nfn``
       (see below).
   * - ``lines<N>.csv``
     - Line flux table listing every computed emission line, its
       identification, and its flux relative to H\ :math:`\beta`.
   * - ``bands<N>.csv``
     - Integrated fluxes in standard broadband photometric filters.
   * - ``rates<N>.csv``
     - Reaction rate coefficients and heating/cooling rates by process.
   * - ``flam_<N>.csv``
     - F\ :math:`\lambda` spectrum (via ``wplam4``) in five columns:
       wavelength (Å, air), total model flux, source-only flux,
       nebula-only flux, and nebula continuum-only flux, all in
       erg s\ :sup:`−1` cm\ :sup:`−2` Å\ :sup:`−1`. Despite "flam" in
       the name this is on the same fixed, log-uniform wavelength grid
       as ``.lam`` (see the caveat below) and its "total"/"nebula"
       columns include emission lines, not just continuum.
   * - ``<X>_ion<N>.csv``
     - Per-element ion fraction profiles, where ``<X>`` is a one- or
       two-character element abbreviation.
   * - ``<X>_col<N>.csv``
     - Per-element column density profiles.

Radiation Field Files
======================

These files record the internal radiation field at selected points and can be
read back as input photon sources for downstream models, allowing models to be
chained.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Filename pattern
     - Contents
   * - ``psou<N>.sou``
     - Primary (stellar or AGN) ionising photon source vector in
       J\ :math:`_\nu` units (1/4π sr).
   * - ``nsou<N>.sou``
     - Diffuse nebular photon source vector accumulated over the model.
   * - ``ssou<N>.sou``
     - Secondary photon source vector (e.g. two-photon continuum).
   * - ``lsou<N>.sou``
     - Local radiation field source vector at the illuminated face.
   * - ``irsou<N>.sou``
     - Infrared radiation field source vector.
   * - ``IRflux<N>.sou``
     - Total IR flux source vector including dust emission.
   * - ``pahs<N>.sou``
     - PAH (polycyclic aromatic hydrocarbon) emission source vector.
   * - ``local<N>.lam``
     - Local radiation field as wavelength (Å) versus
       F\ :math:`\lambda` (erg s\ :sup:`−1` cm\ :sup:`−2` Å\ :sup:`−1`).
   * - ``grpot<N>.ph6``
     - Grain potential output — the electrostatic grain charge as a
       function of depth.


---------------------------------
Shock Models — S5
---------------------------------

The S5 shock solver (``shock5.f``) computes steady-state, magnetised, radiative shock waves
and their photo-ionised precursors.  When starting a shock run the user is asked for a short
prefix string (up to 8 characters) that is embedded in all output filenames, for example
``v100s``, giving files such as ``shck_v100s_0001.sh5``.

Most output files share a common plain-text header block containing the MAPPINGS version, run
name, abundance file, ionisation setup, radiation source, shock parameters, and jump conditions.
After the header, the format varies by file type as described below.  For a worked example
tracing one concrete run's choices to the files below, see :doc:`walkthrough_s5`.

Always-created files
======================

Four files are written for every S5 model regardless of which options are selected.

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Filename pattern
     - Contents
   * - ``shck_<pfx>_<N>.sh5``
     - Main shock cooling-zone structure (see below).
   * - ``prec_<pfx>_<N>.sh5``
     - Precursor zone structure in the same format (see below).
   * - ``specSH<pfx>_<N>.csv``
     - Integrated emission-line spectrum from the shock: wavelength (Å) and
       flux relative to H\ :math:`\beta`, one line per row.
   * - ``specPC<pfx>_<N>.csv``
     - Integrated emission-line spectrum from the precursor, same format.

Structure of ``shck_*.sh5``
-----------------------------

After the header block the file contains one comma-separated data row per spatial step through
the post-shock cooling zone.  Negative step indices cover the initial state:

* **step −2** — far-upstream neutral gas (pre-shock boundary condition)
* **step −1** — pre-shock ionised gas immediately upstream of the shock jump
* **step 0** — immediately post-shock gas (high-temperature plateau after the jump)
* **steps 1, 2, …** — the cooling zone, marching downstream from the shock face

The 16 columns in each row are:

.. list-table::
   :header-rows: 1
   :widths: 8 24 68

   * - Col
     - Quantity
     - Notes
   * - 1
     - Step number
     - Integer; negative values label the pre/jump states listed above
   * - 2
     - Distance (cm)
     - Cumulative distance downstream from the shock face
   * - 3
     - Step width dr (cm)
     -
   * - 4
     - Temperature T (K)
     -
   * - 5
     - n\ :sub:`e` (cm\ :sup:`−3`)
     - Electron density
   * - 6
     - n\ :sub:`H` (cm\ :sup:`−3`)
     - Total hydrogen nuclei
   * - 7
     - n\ :sub:`i` (cm\ :sup:`−3`)
     - Total positive ions
   * - 8
     - μ (amu)
     - Mean molecular weight
   * - 9
     - ``timlps`` (s)
     - Elapsed time since the shock front — the **shock age** at this step
   * - 10
     - dt (s)
     - Time-step width
   * - 11
     - Net cooling / n² (erg cm\ :sup:`3` s\ :sup:`−1`)
     - Normalised net cooling rate; positive = net cooling
   * - 12
     - Total cooling (erg cm\ :sup:`−3` s\ :sup:`−1`)
     -
   * - 13
     - d\ :sub:`los` (cm)
     - Cooling column
   * - 14
     - X(H\ :sup:`0`)
     - Neutral hydrogen fraction n(HI)/n\ :sub:`H`
   * - 15
     - X(H\ :sup:`+`)
     - Ionised hydrogen fraction n(HII)/n\ :sub:`H`
   * - 16
     - B (G)
     - Magnetic field strength

After the per-step rows the file contains three further sections:

* **"Model ended" summary** — the final distance, shock age (``timlps``), and temperature at
  the point where the model stopped.  This line is also printed to the terminal.
* **Temperature structural markers** — for each temperature decade from 10\ :sup:`7` K down
  to 10\ :sup:`2` K that the cooling gas passed through, the file records: the threshold
  temperature, the elapsed time, distance from the shock face, magnetic field, pressure,
  mass density, n\ :sub:`H`, and n\ :sub:`e`.  These are also printed to the terminal.
* **Integrated emission line list** — the complete set of line fluxes relative to H\ :math:`\beta`.

Structure of ``prec_*.sh5``
-----------------------------

The format is identical to the shock structure file, but the rows record the precursor zone.
During each convergence iteration the precursor solver also writes a diagnostic grid table
(step, distance, dx, time, Te, n\ :sub:`e`, n\ :sub:`H`, XHI, XHII, XOI, XOII, XOIII)
to both this file and to the terminal.

Optional output files
=======================

Additional files are enabled through the **Output Multi-Option Menu** presented during setup.
Enter one letter at a time and exit with **X**.

.. list-table::
   :header-rows: 1
   :widths: 8 18 38 36

   * - Key
     - Flag set
     - Shock file(s) created
     - Precursor counterpart
   * - ``B``
     - ``tsrmod``
     - ``elSH<pfx><El>_<N>.csv`` — per tracked element (up to 4)
     - ``elPC<pfx><El>_<N>.csv``
   * -
     - ``allmod`` (sub-option of B)
     - ``ionSH<pfx>_<N>.sh5`` — all ion fractions, all elements
     - ``ionPC<pfx>_<N>.sh5``
   * - ``C``
     - ``dynmod``
     - ``dynSH<pfx>_<N>.sh5`` — dynamics file (see note below)
     - none
   * - ``D``
     - ``ratmod``
     - ``ratSH<pfx>_<N>.sh5`` — rates and timescales
     - ``ratPC<pfx>_<N>.sh5``
   * - ``H``
     - ``bandsmod``
     - ``bandSH<pfx>_<N>.csv`` — X-ray/UV band fluxes
     - none
   * - ``K``
     - ``fclmod``
     - ``coolSH<pfx>_<N>.csv`` — cooling by element
     - none
   * - ``L``
     - ``jlin``
     - ``linSH<pfx>_<N>.csv`` — monitored line emissivities
     - ``linPC<pfx>_<N>.csv``
   * - ``E``
     - ``jspec``
     - ``SHdw<pfx>_<N>.lam`` / ``.sou`` — final downstream field (via ``wpsou``)
     - —
   * - ``F``
     - ``lmod``
     - Adds ``.nfn`` (ν vs. νF\ :sub:`ν`) alongside every ``.lam``/``.sou``
       pair below, including the always-written ones
     - —

.. note::

   **The menu labels for C and D are swapped relative to their actual effect.**
   The menu displays "C = All Rates file" and "D = Dynamics file", but the code
   (``shock5.f`` lines 904–905) sets ``dynmod`` on C and ``ratmod`` on D.  In practice:

   * pressing **C** produces ``dynSH*.sh5`` (the dynamics file)
   * pressing **D** produces ``ratSH*.sh5`` / ``ratPC*.sh5`` (the rates/timescales file)

Selecting **B** (ion balance files) triggers two follow-up prompts:

1. The number of elements to monitor (1–4), then their atomic numbers Z — e.g. ``4`` then
   ``1 2 8 16`` for H, He, O, S.
2. "Record all ions file (Y/N)?" — answering Y additionally creates the wide
   ``ionSH`` / ``ionPC`` files for all elements.

.. note::

   **``SHup<pfx>_<N>.lam``/``.sou`` and ``PCup<pfx>_<N>.lam``/``.sou`` are
   written unconditionally** — the upstream radiation field snapshot is not
   gated by ``E``, ``F``, or any other output-menu option; it is written on
   every run regardless of which optional output tables are selected. Only
   the *downstream* field (``SHdw*``) requires ``E``, and only the extra
   ``.nfn`` detail file for any of the three requires ``F``.

   For what these ``.lam``/``.sou``/``.nfn`` files actually contain
   (continuum and lines together, on a fixed-resolution grid) and their
   limits for observational comparison, see the caveat under
   "Continuum and Spectrum Files" below.

Format of the optional files
==============================

``elSH<pfx><El>`` / ``elPC<pfx><El>``  (CSV, option B)
---------------------------------------------------------

One file per tracked element, named with the element symbol appended directly
after the prefix (e.g. ``elSHv100s_O0001.csv`` for prefix ``v100s_``).  No
separator is inserted before the element symbol or before the sequence
number — any underscore visible in the filename comes from what was typed
at the "prefix for all output files" prompt, which is why the shock prefix
is conventionally given with a trailing underscore (e.g. ``v100s_`` rather
than ``v100s``).  Comma-separated, one row per step.  Columns:

   step, distance (cm), dr (cm), timlps (s), dt (s), T (K),
   n\ :sub:`e` (cm\ :sup:`−3`), n\ :sub:`H` (cm\ :sup:`−3`),
   n\ :sub:`i` (cm\ :sup:`−3`),
   then one column per ion stage of the element (neutral fraction first, then singly
   ionised, doubly ionised, …).

``ionSH<pfx>`` / ``ionPC<pfx>``  (``.sh5``, sub-option of B)
--------------------------------------------------------------

Written only when "Record all ions file" is answered Y.  Contains every ion fraction for
every element in the model at every step.  These are the widest S5 output files.  Per-step
rows include: step, T, n\ :sub:`e`, n\ :sub:`H`, velocity, density, pressure, distance,
``timlps``, then all ion fractions.

``dynSH<pfx>``  (``.sh5``, option C)
--------------------------------------

Not a flat table — each step is a block of five text lines:

1. Step number and solver iteration count
2. Distance (cm), step width (cm), timlps (s), dt (s)
3. Pre-shock state: T (K), velocity (cm/s), ρ (g/cm³), pressure (dyne/cm²), B (G),
   cooling luminosity (erg/cm³/s)
4. Post-shock state: same quantities after the jump
5. Internal energy density (erg/cm³), total number density (cm\ :sup:`−3`), sound speed
   (cm/s), mean molecular weight (g/particle)

``coolSH<pfx>``  (CSV, option K)
----------------------------------

Cooling by element at each step.  Two header lines name each column and give units.
Data columns:

   T (K), n\ :sub:`e`, n\ :sub:`H`, n\ :sub:`i`, ρ (g/cm³),
   X(H\ :sup:`0`), X(H\ :sup:`+`), μ (amu),
   total cooling L (erg/cm³/s),
   L normalised by the chosen density product (erg cm³/s),
   then one column per element giving its normalised cooling contribution,
   then total electron-impact loss, photoelectric gain, and net loss
   — all in erg cm³ s\ :sup:`−1`.

The density normalisation (n\ :sub:`e`\·n\ :sub:`H`, n\ :sub:`H`², n\ :sub:`e`\·n\ :sub:`i`,
etc.) is selected in the setup menu.

``bandSH<pfx>``  (CSV, option H)
----------------------------------

X-ray and UV cooling-band fluxes at each step.  Columns:

   T, n\ :sub:`e`, n\ :sub:`H`, n\ :sub:`i`, X(H\ :sup:`0`), X(H\ :sup:`+`), μ,
   total cooling, cooling length Λ, free-free fraction, then five photon-energy bands:
   0–100 eV, 0.1–0.5 keV, 0.5–1.0 keV, 1.0–2.0 keV, 2.0–10.0 keV, and the total.

``linSH<pfx>`` / ``linPC<pfx>``  (CSV, option L)
--------------------------------------------------

Comma-separated, one row per step.  Columns: step, distance (cm), distance+dr/2 (cm), dr (cm),
T (K), n\ :sub:`e`, n\ :sub:`H`, n\ :sub:`e`\+n\ :sub:`i`, logQH, logUH, logQN (these three are
always 0 — not yet computed for S5), H\ :math:`\beta`, then one column per monitored
line.  Column headers name each line as element symbol + roman numeral ion (e.g. ``O III``);
a second header row gives the line wavelengths.

``ratSH<pfx>`` / ``ratPC<pfx>``  (``.sh5``, option D)
--------------------------------------------------------

Ionisation, recombination, and heating/cooling rates with associated timescales for every ion
at every step.  These are the largest optional output files and can become very large for
multi-element models.

Convergence diagnostics and standard output
============================================

The global shock–precursor iteration loop prints diagnostic information to the terminal that
is **not written to any output file**.  To preserve convergence records, redirect standard
output to a log file::

    ./map52 < model.mv | tee run.log

**Per-iteration convergence test**

After each shock–precursor iteration, ``shock5check`` prints a table comparing six quantities
between the current and previous iterations, plus three further diagnostics:

.. code-block:: text

   ********************************************************
    SHOCK 5  Convergence Test, It.: 02 of 03
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       Psi (Q/v):  1.234e+00   `Psi      :  1.230e+00
       Compress :  4.000e+00   `Compress :  4.001e+00
       T_pre    :  8.500e+03   `T_pre    :  8.450e+03
       T_shock  :  1.200e+06   `T_shock  :  1.198e+06
       ne_pre   :  2.300e+01   `ne_pre   :  2.295e+01
       DelH/He  :  0.031%
       Precursor RMS (Psi,T_pre,ne_pre,DelH/He):  0.010%
       Post-shock RMS (Compress,T_shock)      :  0.003%
       RMS      :  0.008%
       Aitken omega (precursor relaxation)     :  0.732
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     Result: CONVERGED  (precursor RMS  1.0000E-02%, post-shock RMS  3.0000E-03%)
   ********************************************************

The six quantities are the ionisation parameter Ψ (= Q/v), the compression factor, the
precursor temperature, the post-shock temperature, the precursor electron density, and
the H/He ion fraction change.  ``RMS`` is the quadrature sum of all six fractional
differences; convergence is declared when it is below 0.01% (10\ :sup:`−4`).  The
**Precursor RMS** and **Post-shock RMS** lines split that same sum into the four
precursor (pre-shock) terms and the two post-shock terms, so a lingering residual can be
attributed to one side of the shock front or the other — in practice this is almost
always the precursor, and these two values now also directly feed the outcome
classification described below.  **Aitken omega** is the self-tuning relaxation weight
blended into the precursor state that iteration (see :doc:`code_s5`, "How the
precursor↔shock loop actually converges"); it is not fixed, and for a model that is
oscillating rather than converging it is often informative to watch whether omega itself
settles down or keeps cycling.

If convergence is not reached within the requested number of iterations, MAPPINGS does not
stop there: each failed check extends the iteration cap by one and tries again, up to a hard
ceiling of 20 global iterations (``mxshockits`` in ``const.inc``). Whether or not convergence
was ever reached — even after using all 20 — the loop then runs a few more ordinary
iterations followed by one further, final output pass. **There is no failure flag written to
any output file, no non-zero exit code, and no crash** for a run that never converges (with
one exception — see "Known gap: unhandled crash paths" below); the only record is the
"Result:" line on whichever "SHOCK 5 Convergence Test" block was last printed to the
terminal. For a model whose main loop genuinely converges, that last line is now reliable —
a spurious disagreement on the mandatory final independent re-solve no longer overrides an
already-established result. For a model that never converges, the last line classifies *why*
— see :doc:`code_s5`, :ref:`s5_convergence_outcomes`, for the full set of eight possible
outcomes and what each means; the short version is that a model whose only problem is the
characterised, bounded precursor oscillation reports ``NOT CONVERGED (precursor
oscillating)`` and is safe to use as-is, while the other, uncharacterised failure shapes are
reported distinctly and still warrant a look. A run finishing cleanly is still not
automatically the same thing as a run having converged, though — see :doc:`walkthrough_s5`,
"Checking for convergence and other failures", for how to verify this in practice.

**Known gap: unhandled crash paths.**  The outcome classification above only covers models
that run to completion.  Two guards elsewhere in the S5 chain (``cool.f``'s ``Cool out of
range`` check, and ``rankhug.f``'s ``vel0 <= 0`` check) call Fortran ``stop`` directly on an
unphysical state, terminating the whole program with no ``Result:`` line at all — a grid run
just sees a missing or truncated output file, with no graceful record of why.  This is a
pre-existing gap, not something the outcome classification above addresses; fixing it would
mean threading a return/error status up through several call layers instead of a hard
``stop``.

**No-shock outcome**

Some parameter combinations (preshock flow speed below the local Alfvén speed) have no
compressive MHD shock solution at all.  MAPPINGS detects this before attempting the shock
jump and exits that model cleanly rather than crashing; see :doc:`code_s5`, "No shock:
sub-Alfvénic preshock flow", for the exact output.  No ``.sh5``/``.csv`` structure files are
written for a no-shock model.

**Progress heartbeat for slow models**

Some parameter combinations (very low preshock density, in particular) need many thousands
of cooling-zone steps to reach the stopping criterion, which can take a long time even
though the run is working correctly.  Every 100 zone steps, ``compsh5`` prints a progress
line to stdout regardless of the runtime display mode chosen at setup::

    ... SHOCK 5 progress: Global It   4 of   4, Zone Step   601 of 4096, T= 3126. K, Dist= 1.3212E+21 cm

This is flushed immediately, so it is visible in real time even when stdout is redirected to
a file — useful for telling a genuinely slow run apart from one that has actually stopped
making progress.

**Summary of what is in files versus stdout only**

.. list-table::
   :header-rows: 1
   :widths: 52 24 24

   * - Information
     - In ``shck_*.sh5``
     - Stdout only
   * - Per-step structure (T, n\ :sub:`e`, n\ :sub:`H`, B, timlps, …)
     - Yes
     - MINI-mode monitor only
   * - "Model ended" summary (final distance, shock age, T)
     - Yes
     - Also printed
   * - Temperature structural markers
     - Yes
     - Also printed
   * - Integrated emission line list
     - Yes
     - Also printed
   * - Convergence quantities (Ψ, compression, T, n\ :sub:`e`, RMS%,
       precursor/post-shock RMS split, Aitken omega)
     - **No**
     - **Yes — stdout only**
   * - Outcome classification (8 possible ``Result:`` lines, with RMS
       values)
     - **No**
     - **Yes — stdout only**
   * - "No Shock:" summary (Alfvén Mach number)
     - No file written at all
     - **Yes — stdout only**
   * - Zone-step progress heartbeat
     - **No**
     - **Yes — stdout only**
   * - Setup parameters and menu choices
     - **No**
     - **Yes — stdout only**


---------------------------------
Shock Models — S4
---------------------------------

The S4 shock solver (``shock4.f``) is an earlier shock implementation and
produces files with the extension ``.sh4``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Filename pattern
     - Contents
   * - ``shckn<N>.sh4``
     - Main shock zone structure output.
   * - ``dyn<N>.sh4``
     - Dynamical flow variables through the shock.
   * - ``rates<N>.sh4``
     - Reaction and cooling rates.
   * - ``allion<N>.sh4``
     - Complete ion fraction table.
   * - ``bands<N>.sh4``
     - Broadband filter fluxes.
   * - ``spec<N>.csv``
     - Shock emission spectrum.
   * - ``cc<N>.csv``
     - Cooling curve as a function of temperature.


---------------------------------
Non-Equilibrium Cooling
---------------------------------

The non-equilibrium cooling solver (``neqc.f``) computes time-dependent
cooling of a plasma parcel and produces files with the extension ``.neq``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Filename pattern
     - Contents
   * - ``neqcl<N>.neq``
     - Main non-equilibrium cooling output: temperature, ionisation, and
       cooling as a function of time.
   * - ``neqc<N>.csv``
     - Tabulated cooling data.
   * - ``rates<N>.csv``
     - Reaction rates as a function of time.
   * - ``dyn<N>.csv``
     - Dynamical quantities as a function of time.
   * - ``spec<N>.csv``
     - Emission spectrum integrated over the cooling.
   * - ``bands<N>.csv``
     - Broadband fluxes.


---------------------------------
Continuum and Spectrum Files
---------------------------------

These file types are produced across multiple model types depending on which
output options are selected.

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - Extension
     - Contents
   * - ``.lam``
     - Wavelength (Å, air) versus F\ :math:`\lambda`
       (erg s\ :sup:`−1` cm\ :sup:`−2` Å\ :sup:`−1`) — two columns.
       Restricted to 1000–50000 Å, zero-flux bins omitted. **Contains
       both the continuum and (almost) every emission line**, not
       continuum alone — see the caveat below.
   * - ``.nfn``
     - Frequency (Hz) versus :math:`\nu F_\nu`
       (erg s\ :sup:`−1` cm\ :sup:`−2` sr\ :sup:`−1`) — two columns, full
       native energy range (no 1000–50000 Å restriction). Same
       continuum+lines content as ``.lam``, just unrestricted range and
       different units/axes (note the ``/sr`` here — unlike ``.lam``
       this is not multiplied by 4π).
   * - ``.emi``
     - Emission spectrum normalised relative to H\ :math:`\beta`.
   * - ``.dat``
     - X-ray band flux data.
   * - ``.sou``
     - Photon source vector: photon energy (eV) versus J\ :math:`_\nu`
       (erg s\ :sup:`−1` cm\ :sup:`−2` Hz\ :sup:`−1` sr\ :sup:`−1`,
       "1/4π" units) — two columns, whitespace-separated (not comma),
       full native energy range. The rawest, most complete form of this
       same continuum+lines field; can be fed back as input to a
       subsequent model to chain calculations.
   * - ``.txt``
     - Plain-text output from slab geometry models (``slab.f``).

.. admonition:: Caveat — ``.lam``/``.sou``/``.nfn`` are useful for quick
   characterization, not for direct comparison to observations

   These three files are all views of the same ``tphot`` array (built in
   ``totphot2``/``newdif2``, see :doc:`physics_output_spectra`), and that
   array is **not** continuum-only: ``src/localem.f`` injects
   collisionally-excited/forbidden lines (``fmbri``/``febri``/``f3bri``/
   ``fsbri``), He I recombination, and dielectronic/RR satellite lines
   directly into the same per-bin array (``emidif``) as the true
   continuum, and H/He/X-ray-complex lines get their own additional,
   more careful zone-to-zone transport on top of that
   (``src/totphot.f:481-868``). So a ``.lam`` spectrum genuinely shows
   both continuum and lines together — confirmed by matching peak
   wavelengths in real output against the standard optical CEL set
   (Hα, [O III], [O II], [N II], [S II], [O I], [Ne III], [S III], He I,
   etc.).

   However, the wavelength/energy grid these files are written on
   (``photev``, loaded once from ``data/PHOTDAT.txt``) is a **fixed,
   log-uniform grid at a constant resolving power of R ≈ 3500**
   (Δλ/λ ≈ 2.86×10⁻⁴ everywhere, verified directly from real output).
   Every line is simply assigned to whichever native bin its rest
   energy happens to fall in (``mapinit.f``, e.g. ``hbin(line,series)``)
   — no bin is specially narrowed or inserted for it — and the line's
   full luminosity is divided by that bin's native width to produce the
   per-bin flux density. That means:

   - A line's *peak* Flambda/Jν value is not directly comparable between
     lines at different wavelengths, because it is implicitly divided by
     a bin width that scales with wavelength (Δλ = λ/3500). To recover
     an actual, resolution-independent line flux from a peak in these
     files, multiply by the local bin width: ``flux_line ≈
     Flambda_peak(λ) × λ/3500``. Verified against the authoritative
     volume-integrated line list (``specSH*.csv``/``spec<N>.csv``) to
     under 1% agreement across a 2.7× wavelength range.
   - Two real lines separated by less than roughly λ/3500 (≈1.3–2.7 Å
     across the optical) land in the same bin and are indistinguishable.
   - There is no instrumental-resolution convolution, flux calibration,
     or continuum-normalization applied — this is the model's own
     internal radiative-transfer grid, not a synthetic observation.

   **In short**: ``.lam``/``.sou``/``.nfn`` are well suited to a quick
   look at the overall shape and relative strength of a model's
   continuum and lines together, but are not a substitute for a properly
   synthesized, instrument-matched spectrum when comparing to real data
   — use the volume-integrated line fluxes (``specSH*.csv``/``lines<N>.csv``)
   for quantitative line-flux work, and see :doc:`physics_output_spectra`
   for the full mechanism.

``.bln`` — Ionisation Balance File
=====================================

Not a continuum or spectrum file despite living alongside them by extension.
Written unconditionally on the final iteration by ``wbal`` (``src/output.f``)
for every model type — P6, P7, S5, and the slab geometry models (``slab.f``)
all call it. Contents: a plain-text table of the final ionisation state, one
row per ion stage of every element present (atomic number, ion stage, and
population fraction), preceded by a short header naming the run.

---------------------------------
Notes on File Naming
---------------------------------

Sequence numbers
================

Within a single MAPPINGS run, each model increments the sequence counter
so files from successive models do not overwrite one another.  The counter
resets to ``0001`` at the start of each run.

Shock prefix string
====================

For S5 shock models the user is asked at run time::

    Give a prefix for all output files (max 8 chars):

This string (e.g. ``v100s`` for a 100 km/s shock) is embedded in every
output filename for that model, making it straightforward to run and compare
grids of models with varying parameters in the same directory.

MAPPINGS does not insert its own separator between the prefix and the
rest of the filename (element symbol, sequence number, etc.) — it is
simply concatenated.  By convention the prefix itself is typed with a
trailing underscore (e.g. ``v100s_``) so that filenames like
``shck_v100s_0001.sh5`` and ``elSHv100s_O0001.csv`` remain readable; a
prefix without one produces run-together names such as
``elSHv100sO0001.csv``.

Chaining models with ``.sou`` files
=====================================

The ``.sou`` photon source files output by a photoionization or shock model
can be read back by a subsequent model as the incident radiation field.  This
allows, for example, the nebular emission of a model HII region to be used
as the ionising source for a downstream cloud, or a shock precursor to be
iterated self-consistently.
