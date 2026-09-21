# MAPPINGS V v5.2.1

		Creative Commons v4.0 International
		By Attribution, Share Alike
		CC-BY-SA-4.0Intl https://creativecommons.org
		1976 -- 2022+ Ralph Sutherland,
		Michael Dopita, Luc Binette, Ian Evans,
		Brent Groves, David Nicholls,
		Adam D. Thomas, Yi-Fei Jin, Knox Long



#### References:

When using this code please cite the following papers:
    
	* Sutherland, R., & Dopita, M. A. 1993, ApJS, 88, 253
	* Sutherland, R., & Dopita, M. A. 2017, ApJS, 229, 34S

#### Contact:

	* Ralph Sutherland: Ralph.Sutherland@anu.edu.au          
        * Yifei Jin:        yfjsci@gmail.com
	* Knox Long:        long@stsci.edu


## Documentation

The full set of documentation for Mappings can ge found on ReadtheDocs at https://mappings-v.readthedocs.io/en/latest/index.html

## Quick Start Compiling Install and Run:

		> cd download or git area...
		> make [-j] build

[wait a while.... ][
[for a faster build try make -j build if your make supports it]

MAPPINGS is installed by default into the users home area
~/mappings520

To run

		> cd ~/mappings520/lab
		> ./map52


## Uninstalling:

To remove the MAPPINGS installation simply and leave user files

		> make uninstall
.
Or to remove totally simply delete the `~/mappings520` directory

## The Installed MAPPINGS directory structure:

#### mappings520/
		lab/
		lab/scripts/

		data/
		atmos/
		abund/
		prefs/
		docs/
		tools/

### Key Installed Directories and Files:

* `lab/`  This is where the executable is made and run.  The runtime files such
	as map.prefs are here, the user can create and rename new labs and put map.prefs into them and run
	other models.
* `map52`: The executable.
* `map.prefs`   Essential startup data - must be present.
* `mapStd.prefs`  A standard 16 atom startup in case map.prefs is lost for any reason, can be copied and renamed map.prefs if needed
* `mapFull.prefs`  A full 30 atom startup in case map.prefs is lost for any reason, can be copied and renamed map.prefs if needed
* `scripts/`: A collection of (mostly) useful of UNIX shell and MAPPINGS scripts

* `data/`: Contains the atomic data for MAPPINGS V. Read at runtime in mapinit.f.  Essential and must be present and complete.
* `abund/`:  A set of useful abundance settings that can be read interactively during a run.  Optional.
* `atmos/`:  A set of useful radiation source files and stellar atmosphere models.  Optional.
* `docs/`:  Help, source references and cookbook guides
* `tools/`:  A set of C and other tools for manipilting MAPPINGS output such as reddening

### Key Build Area Directories and Files:

The master copies of the data and other directories above plus:

* `src/` : Contains the code

* `Makefile`   This makefile  controls all the building of MAPPINGS.
It takes one argument to control the operation: build is the main one.

         MAPPINGS V make options:
         -----------------------------------------------------------

         'make, make help'             To see this menu

         'make build'      Build and install, and clean.
                           Creates MAPPINGS in /Users/ralph/mappings520
                           and creates optional environment variable
                           templates for startup scripts.
         'make compile'    to make new '*.o' in the build area only
         'make clean'      to clean up '*.o' files from a compile
         'make distclean'  as clean but also remove map52rss

         'make install'    install afer manual compile
         'make uninstall'  remove installed map52rss

         -----------------------------------------------------------

