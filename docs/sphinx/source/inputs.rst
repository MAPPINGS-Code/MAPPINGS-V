.. _inputs:

#####################
Startup Inputs
#####################

When MAPPINGS starts it works through a fixed sequence of initialisation steps
before presenting the model selection menu.  Each step either reads a
configuration file silently or presents an interactive prompt.  This page
describes each step in order.

.. contents:: Contents
   :local:
   :depth: 2

---------------------------------
The Preferences File (map.prefs)
---------------------------------

Before any interactive prompts appear, MAPPINGS reads ``map.prefs``.  This
file is required; MAPPINGS will stop immediately if it cannot be found.
It is searched for in the following order:

1. The current working directory (``./map.prefs``)
2. ``$MAPDATA/map.prefs``
3. ``$MAPDATA/prefs/map.prefs``

The standard installed copy is in ``~/mappings520/lab/map.prefs``.  When
running models in a custom directory, copy ``map.prefs`` into that directory
before starting MAPPINGS.

File format
===========

Lines beginning with ``%`` are comments and are skipped.  The first
non-comment line is a free-text title.  After the title, one line gives the
number of atom types, followed by one line per atom::

    <index>  <symbol>  <Z>  <atomic_weight>  <abundance>  <depletion>  <max_ion>

where:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Meaning
   * - ``index``
     - Sequential atom number (1, 2, 3, …)
   * - ``symbol``
     - Element symbol (e.g. ``H``, ``He``, ``C``)
   * - ``Z``
     - Atomic number
   * - ``atomic_weight``
     - Atomic weight in atomic mass units
   * - ``abundance``
     - Default abundance by number relative to hydrogen (linear scale)
   * - ``depletion``
     - Default dust depletion factor (1.0 = no depletion, <1.0 = depleted)
   * - ``max_ion``
     - Maximum ionisation stage tracked for this element

The default abundances follow Asplund et al. (2009, ARA&A 47, 481) solar
photosphere values.

Optional settings
=================

After the atom table, ``map.prefs`` may contain keyword blocks (also
``%``-commented and preceded by a keyword line).  Currently recognised
keywords:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Keyword
     - Effect
   * - ``PHOTDAT``
     - Override the energy-grid file (default ``PHOTDAT-std.txt`` in the
       data directory).  The replacement file name follows on the next line.
   * - ``Kappa Mode``
     - Set to ``1`` on the following line to enable the kappa electron
       distribution menu at startup (default ``0`` = disabled).

Alternative prefs files
========================

The ``addons/`` directory contains several ready-made ``map.prefs``
variants:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Description
   * - ``mapStd.prefs``
     - Standard 16-element set (H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl,
       Ar, Ca, Fe, Ni).  Safe fallback if ``map.prefs`` is lost.
   * - ``mapFull.prefs``
     - Full 30-element set including all tracked species.
   * - ``mapStdLGC.prefs``
     - Standard set with Large Grain Cooling enabled.
   * - ``mapFullLGC.prefs``
     - Full set with Large Grain Cooling enabled.

Copy the desired file into the working directory and rename it
``map.prefs`` to use it.

Expert switches
===============

A second configuration file, ``data/switches.txt``, controls low-level
physics toggles that are not exposed in the interactive menus.  The
defaults are appropriate for almost all uses; consult the file comments
before changing any switch.  The available switches are:

.. list-table::
   :header-rows: 1
   :widths: 45 30 25

   * - Switch
     - Default
     - Options
   * - Collisional excitation data set
     - 0 (Dere)
     - 1 = A&R, 2 = old Shull & McKee
   * - Photoionisation cross-sections
     - 0 (Verner & Ferland)
     - 1 = old
   * - Verner calculations
     - 0 (effective outer shell)
     - 1 = strict Hartree-Slater
   * - Radiative recombination mode
     - 0 (new)
     - 1 = old, 2 = NORAD for CNO
   * - Dielectronic recombination mode
     - 0 (new)
     - 1 = old, 2 = off
   * - Charge transfer mode
     - 0 (new MV)
     - 1 = A&R85, 2 = ancient, 3 = off
   * - Charge transfer heating
     - 1 (on)
     - 0 = off
   * - Free-free continuum
     - 1 (accurate)
     - 0 = fast
   * - Photoionisation
     - 1 (on)
     - 0 = off
   * - H/He collisional excitation
     - 0 (ABBS02)
     - 1 = GNP89
   * - Lines-only cooling
     - 0 (off)
     - 1 = on
   * - Power-law cooling
     - 0 (off)
     - 1 = on
   * - *(unused, formerly turbulent heating)*
     - 0
     - Ignored; the line is kept so older files still read correctly
   * - Radiation pressure
     - 0 (on)
     - 1 = off


---------------------------------
Step 1 — Elemental Abundances
---------------------------------

After initialisation MAPPINGS displays the current abundance table and
metallicity summary, then asks::

    Changing the global abundances
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    File I/O Abundances:
       Y  :  Select an abundance file
     X,N  :  Exit with no further changes.
    ::

Enter ``Y`` to load an abundance file, or ``N`` (or just press Return) to
continue with the defaults.

If ``Y`` is selected, MAPPINGS prompts for a filename::

    Enter abundance file name :

The file is searched for:

1. As an absolute or relative path as given
2. In the local ``abund/`` subdirectory
3. In ``$MAPDATA/abund/``

A collection of pre-made abundance files is installed in
``~/mappings520/abund/``.

Abundance file format
=====================

Lines beginning with ``%`` are comments.  The structure is:

1. Comment lines (``%`` prefix), skipped
2. A free-text title line (up to ~76 characters)
3. One line: ``<nentries>  <numtype>``
4. ``nentries`` lines, each: ``<Z>  <value>``

The ``numtype`` controls how ``value`` is interpreted:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - numtype
     - Interpretation
   * - 0 (mixed)
     - Negative values are log abundances by number (relative to H = 1);
       positive values are linear abundances by number.
   * - 1
     - All values are log\ :sub:`10` abundances (relative to H = 1)
   * - 2
     - All values are linear abundances by number (relative to H = 1)
   * - 3
     - All values are log\ :sub:`12` (i.e. log\ :sub:`10`\ (n/n\ :sub:`H`) + 12)

Elements with atomic number ``Z`` not in the current ``map.prefs`` atom
list are silently ignored.  After loading, the abundance table is
redisplayed and the prompt repeats, allowing multiple files to be applied
in sequence.  Enter ``N`` or ``X`` to proceed.


---------------------------------
Step 2 — Abundance Offsets
---------------------------------

MAPPINGS then displays the current multiplicative abundance scaling factors
(all 1.0 by default, meaning no change) and asks::

    Change Abundance Offsets (y/N)? :

These scalings are applied *on top of* the base abundances set in Step 1,
allowing targeted per-element adjustments without replacing the entire
abundance table.

Enter ``Y`` to load an offsets file; the format is identical to the
abundance file format described above, but the values are multiplicative
scalings rather than absolute abundances.  A value of 1.0 leaves an
element unchanged; 0.5 halves its abundance.


-----------------------------------------
Step 3 — Electron Energy Distribution
-----------------------------------------

This step only appears if ``Kappa Mode`` is set to ``1`` in ``map.prefs``.
By default it is disabled and this section is skipped.

When enabled, MAPPINGS displays::

    Electron Energy Distributions:
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Maxwellian thermal distributions are being used.
    Use Kappa electron distributions ? (y/N) :

A Maxwellian distribution is assumed by default (equivalent to
κ → ∞).  Enter ``Y`` to use a kappa distribution, which adds a
high-energy power-law tail to the electron energy distribution and can
affect line diagnostics and temperature measurements.

If ``Y`` is selected, MAPPINGS prompts::

    Select a Kappa value
    (2.0 - 1000.0, outside this range disabled)
    (2.0, 3.0, 4.0, 6.0, 10.0, 20.0, 50.0, 100.0 are exact) :

Values at the listed exact grid points are computed directly; other values
in the range 2–1000 are obtained by interpolation.  Values outside this
range disable kappa and revert to a Maxwellian.  Lower kappa values
produce a larger departure from Maxwellian.

See :doc:`adv_kappa` for the physical background and motivation.


---------------------------------
Step 4 — Dust Grain Physics
---------------------------------

MAPPINGS displays the current dust status and asks::

    Dust Physics:
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Dust is currently disabled.
    Include dust calculations? (y/N) :

Enter ``Y`` to enable dust.  If dust is enabled, a further sequence of
prompts follows.  See :doc:`physics_dust` for the physical background
on grain/PAH heating, charging, infrared re-emission, and destruction.

Depletion factors
=================

Enabling dust immediately invokes the depletion setup (``depcha``), which
displays the current depletion factors and asks::

    Change dust depletion (y/N) :

The depletion factor for each element controls what fraction of that
element remains in the gas phase (1.0 = undepleted, 0.5 = half depleted
onto grains).  Entering ``Y`` prompts for a depletion file; its format is
identical to the abundance file format with values representing linear
depletion factors (values ≤ 0 are interpreted as log\ :sub:`10`).

Default depletion values are set in ``map.prefs`` (the sixth column of
each atom line).

Grain destruction
=================

After depletion, MAPPINGS asks::

    Allow grain destruction? (y/N) :

If enabled (``Y``), grains are destroyed above a user-specified radiation
intensity and temperature limit.  MAPPINGS then asks for the limit type::

    Choose the grain photoionisation limit
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Upper ionisation parameter in terms of:
       Q  :  QHDN [nion/nH = ...]
       H  :  QHDH
       U  :  U(H)
    ::

Depending on the choice, the corresponding limit value is requested
(values ≤ 100 are treated as log\ :sub:`10`\ ):

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Option
     - Prompt
   * - Q
     - ``Give QHDN dust limit (<=100 as log [10.500]) :``
   * - H
     - ``Give QHDH dust limit (<=100 as log [10.500]) :``
   * - U
     - ``Give U(H) dust limit (<=100 as log [0.00]) :``

A dust temperature limit is then requested::

    Give dust temperature limit (<=10 as log [5.00]) :

Grain size distribution
=======================

Finally, MAPPINGS asks for the grain size distribution model::

    Choose the grain distribution model
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       P  :  Powerlaw     N(a) = k a^alpha
       M  :  MRN distribution
       S  :  Grain Shattering Profile
      (Note: need full DUSTDATA for P or S)
    ::

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Option
     - Description
   * - M
     - Mathis, Rumpl & Nordsieck (MRN) power-law with α = −3.5,
       carbonaceous grains from 50–2500 Å, silicate grains from 100–2500 Å.
       No further input required.
   * - P
     - User-specified power-law N(a) ∝ a\ :sup:`α`.  MAPPINGS prompts
       for α (MRN default is −3.5) and for the minimum and maximum grain
       radii.
   * - S
     - Grain shattering profile.  Requires the full ``DUSTDATA`` dataset.


---------------------------------
Final Display and Model Selection
---------------------------------

After all the above steps, MAPPINGS prints a summary of the final gas-phase
abundances including the metallicity relative to solar, the total ion-to-H
ratio (n\ :sub:`i`/n\ :sub:`H`), and mean molecular weights for neutral
(μ\ :sub:`neu`), fully ionised (μ\ :sub:`ion`), and hydrogen-only (μ\ :sub:`H`)
plasma.

The model-type selection menu then appears.  The available model types are
described in their own sections of this documentation.  In particular,
each model type prompts separately for its own ending/stopping
condition later in its setup sequence — this is not part of the
startup sequence above.  For shock models (S5), see :doc:`code_s5`,
"Step 6 — Check stopping condition", which also covers the hard zone
cap, the ``terminate`` poll file, and other conditions that can end a
run beyond the chosen stopping criterion.  For photoionization models,
see the equivalent table in :doc:`code_photo`.

Running from a script
=====================

All of the above prompts can be answered non-interactively by piping a
plain-text script to the MAPPINGS executable::

    map52 < mymodel.mv

Each line in the ``.mv`` file corresponds to one prompt response, in
order.  Example scripts for common model types are provided in
``~/mappings520/lab/examples/`` and are also downloadable directly from
:doc:`walkthrough_p6` and :doc:`walkthrough_s5`, alongside a full
line-by-line explanation of each one.
