.. MAPPINGS documentation master file, created by
   sphinx-quickstart on Sun Jan 14 18:04:35 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

############
*MAPPINGS V*
############
--------------------------------------------------------------------------------------------------------------
Astrophysical equilibrium and time-dependent photoinisation and steady supersonic shock spectral emission code
--------------------------------------------------------------------------------------------------------------

*MAPPINGS V* is a photoionisation and shock modelling code written in FORTRAN. It allows the compute nebular emission spectra from the far UV-Xrays (~0.1 Å or 100 keV) to the far IR (>1000 μm or 1e-6 eV).

*MAPPINGS* was originally written by Dopita in 1976 and descrbed in a PhD thesis in 1982 by Binette (ANU) and then by PhD thesis by Sutherland in 1993 (ANU) and a cooling paper in 1993, and Sutherland has curated
the code since 1989.  The code has evolved since then and useful early references is are `Sutherland and Dopita 1993` and `Dopita and Sutherland (1996) <https://ui.adsabs.harvard.edu/abs/1996ApJS..102..161D/abstract>`_


*MAPPINGS*
This code calculates ionisation and emission for continuum and line fluxes of recombination and collisionally excited lines
for all elements from hydrogen to zinc, for over 88,000 lines in the current version.

 MAPPINGS began in 1976 as a five-level-atom solver and simple shock model. It has gone through major upgrades since then and is currently at version 5.2.x with 174 multi-level ions and some up to 900 level level ion models with cascades, charge exchange and other processes and full continuum and diffuse field estimation.  It produces output for nebula structures as well as emission spectra both as line lists and f-lambda and f-nu energy and wavelength files to compare with observations of line fluxes and  entire spectra.

It can be used to model equilibrium and time dependent astrophysical plasmas and emission and is used for atomic diagnostics, cooling, model HII regions, Planetary Nebulae, Nova Shells, AGN emission regions and a wide range of Herbig Haro and fast Shockwaves in the interstellar medium.

-------------
Installation
-------------

See :doc:`installation` to build MAPPINGS and run your first model, then
:doc:`walkthrough_p6` or :doc:`walkthrough_s5` for a complete worked
example of the two main model types.

-------------
Issues
-------------

Issues regarding the code and suggestions for improvement the code regarding the should be reported there.  We actively encourage other to make use of the code for their own science. If anyone has questions about whether the code might be useful for a project, we encourage you to contact one of the authors of the code.

-------------
References
-------------

If you make use of *MAPPINGS* in your published research we ask that you reference the following papers:

Add a more limited set of Sutherland and Dopita 1993` and `Dopita and Sutherland (1995) <https://ui.adsabs.harvard.edu/abs/1996ApJS..102..161D/abstract>`_
<https://ui.adsabs.harvard.edu/abs/1996ApJS..102..161D/abstract>`_

-------------
Documentation
-------------

Various documentation exists:

For more information on how this page was generated and how to create more extensive code level documentation for *MAPPINGS*, look at the page for :doc:`documentation on the documentation <meta>`.

-------
Authors
-------
The authors of the *MAPPINGS* code and their institutions are:


Ralph S. Sutherland

Michael A .Dopita
  Research School of Astronomy & Astrophysics
  Australian National University

More recent updates have been makde by:

Yifei Jin
    Westlake University

Knox S.Long
    Space Telescope Science Institute



----------------------------------------

.. toctree::
   :titlesonly:
   :glob:
   :hidden:
   :caption: Documentation

   installation
   repo_structure
   Mappings_guide
   inputs
   models
   popcha
   photsou
   outputs
   walkthrough_p6
   walkthrough_s5
   code_overview
   code_photo
   code_s5
   code_singlezone
   code_testatom
   physics_overview
   physics_ionization
   physics_lines
   physics_shocks
   physics_output_spectra
   physics_dust
   physics_heating_cooling
   physics_continuum
   physics_time_integration
   meta
   *
