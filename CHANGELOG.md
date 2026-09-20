# Changelog

A user-facing summary of significant changes in MAPPINGS V since this repository was created (September 2024). This covers changes that affect what users experience — build reliability, correctness of output, and available documentation — not every commit. For the full commit history, see `git log`.

## Getting the code to actually build and run (Sept–Dec 2024)

The repository started as a drop of the existing MAPPINGS V codebase, and the first few months were about making that code reliably buildable and runnable from a fresh checkout rather than adding new capability:

- **A required data directory was missing from the repo** — without it, MAPPINGS could not run at all. Fixed shortly after being noticed.
- **The build failed on a completely fresh clone** for other users; fixed, along with some atomic data files that had gone missing.
- **Apple Silicon (M1/M2/M3) support** — a numerical routine (`absdis`) hung or produced wrong behavior on Apple's ARM chips under gfortran. Fixed, and released as **v5.2.1** specifically so Mac users on newer hardware could be confident the code worked correctly on their machines (v5.2.0 remained the reference version otherwise).
- **Output formatting was made consistent** — CSV files are now comma-separated throughout (some shock-related headers previously weren't).
- Minor portability fix so setup scripts work on modern Linux systems that don't ship `tcsh` by default.

## A long quiet period, then a major documentation effort (Jan 2025 – Jun 2026)

After January 2025 there's a roughly 16-month gap with no repository activity. Work resumed in mid-2026 with the addition of a full **Sphinx documentation site** — before this, understanding how the code worked required reading Fortran source directly. This added structured, browsable documentation covering how the code operates internally and the repository's structure.

## Physics documentation completed (August 12, 2026)

A systematic set of documentation pages was written covering the physics MAPPINGS actually computes: ionization balance, line radiative transfer, dust and PAH physics, heating/cooling channels, the continuum calculation, and non-equilibrium time integration — plus a page directly comparing how shock (S5) and photoionization (P6/P7) models differ internally. This is reference material for understanding *what* the code is doing, not just how to run it.

## Real output-correctness bugs fixed (August 16–19, 2026)

This stretch fixed several bugs that affected the correctness of what users actually see in menus and output files, not just documentation:

- Two menu-handling bugs (multi-line option output, multi-option menu mismatch).
- **Monitor-line wavelength matching was fixed for real** — an earlier attempt hadn't actually solved the underlying matching problem; this pass did, and made the matching fail loudly instead of silently when a line genuinely can't be matched.
- CSV column headers for monitor-line output were corrected to match the data actually in the columns.
- Documentation corrections clarifying what several output file types actually contain (`.bln` is an ionization-balance table, not a continuum spectrum; `.lam`/`.sou`/`.nfn` contain continuum *and* lines, not continuum alone).
- **A real numerical bug (#6): flux values in `.lam` and `flam_<N>.csv` output were inflated** by an extra factor (`evplk`) that shouldn't have been applied. Anyone using those flux values from an affected run had numbers that were quantitatively wrong, not just mislabeled. Fixed.

## Shock model convergence made reliable (August 23–25, 2026, issues #7 & #8)

The most substantial correctness work in the repository's history: the iterative loop that solves a shock model (the post-shock cooling zone and the pre-shock "precursor" gas it ionizes, solved back and forth until they agree) had several independent bugs that caused many models to falsely report as not converged, or in rare cases to hang or crash outright.

- Fixed the convergence check itself, which had been silently miscounting one of its variables.
- Replaced a fixed damping scheme with an adaptive one (Aitken relaxation) that resolves oscillations that were causing legitimate non-convergence in some models.
- Found and fixed the core structural bug: the model's internal solve was only ever accurate to 10%, while the code was comparing that against a 0.01% target — so many models were being judged against a standard they could never have met, iteration after iteration. Fixing this made previously-stuck models converge cleanly, at a modest (~15–20%) increase in run time.
- Fixed a rare crash (`Cool out of range`) affecting cold, dilute, weakly-magnetized shocks, caused by a numerical cancellation in the shock-jump solver (#8).
- Added a clean, correctly-labeled result for shock speeds too low to produce a physical shock at all (previously this could hang or crash instead of just reporting "no shock").
- Fixed a case where the final line of output could misreport a model as unconverged even when it had genuinely converged.

Net effect for a user: shock model runs are now far more likely to converge, no longer hang or crash on certain physically-valid-but-difficult inputs, and the convergence result reported at the end of a run can now actually be trusted. Two known, unusual parameter combinations still don't converge and remain open (tracked in issue #7).

## Ending-condition documentation expanded; a P6/P7 prompt bug found and tracked (September 14, 2026, issue #9)

Documentation for how S5 shock and P6/P7 photoionization models choose when to stop running was significantly expanded in `models.rst`, `code_s5.rst`, and `code_photo.rst`. Previously this information was scattered or missing from the user-facing model-input reference, requiring a trip into the internal "Code Operation" pages to find it. Every ending-condition letter for both model types is now documented where a user actually looks for it, including the exact follow-up prompt text, unit/log conventions, and the always-on safety conditions (hard zone caps, the `terminate` poll file, and P6/P7's automatic recombination floor) that can end a run independently of the chosen condition.

- While auditing P6/P7's `C` (temperature-bounded) ending condition against the source, found that its prompt claims values under 10 are read as `log10(K)` — matching the convention used by S5's equivalent condition and by P6/P7's own `F`/`H` conditions — but no such conversion actually exists in `photo6.f`/`photo7.f`; the value is used directly as Kelvin. A user following the prompt's own stated convention (e.g. entering `4` intending 10,000 K) would silently get a 4 K ending temperature instead. Not yet fixed; tracked in issue #9.

## Read the Docs publishing enabled (September 20, 2026)

Fixed two things that would have blocked the Sphinx documentation from building on Read the Docs' hosted service: `.readthedocs.yaml` pinned a build image RTD has since retired, and `requirements.txt` pinned `docutils<0.18` while leaving `sphinx` unpinned — current Sphinx requires a newer `docutils`, so a clean `pip install` on RTD's infrastructure would have failed even though the same files happened to work locally against an already-installed environment. Also removed leftover cruft in `conf.py` (wrong GitHub repo for issue links, dead `sys.path` entries) inherited from the template the docs were originally cloned from.
