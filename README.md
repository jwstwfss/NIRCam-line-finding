# Documentation for NIRCAM-line-finding code - this version is for new file structure

Farhan Hasan, created 2/13/2025 ([fhasan@stsci.edu](mailto:fhasan@stsci.edu))

Last updated 1/7/2026
Added ancillary catalogs functionality - using redshifts from GOODS-N UVCANDELS (Mehta+2024) and JADES (DR4).
These catalogs are in poppies_analysis/anc_data

This software is used to identify line-emitting objects and measure emission line properties in JWST NIRCam WFSS Grism spectra, based on the pure-parallel survey POPPIES (PID#5398).

This code was adapted from the JWST WFSS line-finding code for the PASSAGE survey ([https://github.com/jwstwfss/line-finding](https://github.com/jwstwfss/line-finding)), which itself was adapted from A. Henry's WISP fitting code (with Battisti et al., 2024).

## How it works

A user interactively inspects emission line candidates using 1D and 2D spectra and images.

The user makes a choice between 1) analyzing all objects identified from the photometry, or 2) analyzing the subset of objects that are identified as potential emission line objects using a peak finding algorithm that performs a continuous wavelet transform (CWT) on the spectral data (see SciPys's [find_peaks_cwt](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks_cwt.html)). The code iterates through the potential emission lines, allowing the reviewers to identify lines or reject spurious detections quickly and efficiently.

The user will be asked to enter the number of a parallel field and a username of choice. The fitting is then performed on an object-by-object basis.

----

## How to install and run

Steps for installation are as follows

1. Make sure you have python 3.9 or later installed
2. pip install the following packages:
	- numpy>=1.26.4
	- astropy>=5.3, <7
	- pandas>=1.0, <3
	- scipy==1.14.1 (this one is strict)
	- matplotlib>=3.6, <3.11
3. Git clone [https://github.com/jwstwfss/NIRCam-line-finding.git](https://github.com/jwstwfss/NIRCam-line-finding.git)
4. Install **XPA** from: [https://github.com/ericmandel/xpa](https://github.com/ericmandel/xpa)

a) Git clone from the repo

b) Then, cd into the xpa directory

c) `./configure --prefix=<top_level_install_dir>` # site-specific configuration

`make`  # build the software

`make install` # install it

`make clean` # clean up unneeded temp files

d) Add the path to your bash:

`export PATH=<path/to/software>:$PATH`

e.g., export PATH=/Users/fhasan/Desktop/Software/xpa/bin:$PATH

[See more instructions in *xpa/INSTALL*]

e) Alternatively, you can just download XPA and follow installation instructions from the [SAOImageDS9 site] (https://sites.google.com/cfa.harvard.edu/saoimageds9/download) 

5. Add NIRCAM-line-finding/poppies_analysis to your **$PYTHONPATH**. In your *~/.bash_profile* file, include the following line:

`export PYTHONPATH="<path/to/code>:$PYTHONPATH"`

e.g., export PYTHONPATH= "/Users/fhasan/POPPIES/line-finding/poppies_analysis/:$PYTHONPATH"

(Alternatively, if you use csh, add to ~/.cshrc: `setenv PYTHONPATH <path/to/code>:${PYTHONPATH}`)

6. Install the ***stenv*** Conda environment: [https://stenv.readthedocs.io/en/latest/](https://stenv.readthedocs.io/en/latest/)

[Make sure to activate your stenv, e.g., `conda activate stenv`]

7. Install **DS9** if you haven't already: [https://sites.google.com/cfa.harvard.edu/saoimageds9](https://sites.google.com/cfa.harvard.edu/saoimageds9)


You may also consult the installation steps outlined by Anika Kumar at RIT in this [Google Doc] (https://docs.google.com/document/d/1mjbGBk05_j2NcnEHEbcGD3JI44Kfy0NOJMgJp7XFTW0/edit?tab=t.0)


Steps for running the line-finding tool after installing are as follows:

1. Download the data.
2. Point CODE_DIR, OUTPUT_DIR, and DATA_DIR to the appropriate directories in the [mainPOPPIES.py](http://mainPOPPIES.py) file.
3. Run `python mainPOPPIES.py`
4. Enter the number of the parallel field (e.g., "004").
5. Hit ENTER or 'c' to go through the wavelength decomposition and find peaks, and only go through candidates identified by this method, or 'all' to go through all objects identified photometrically.
6. Enter a username.
7. Inspect objects individually.
8. Choose what to do with each object (accept, reject, etc.). The outputs are saved in the OUTPUT_DIR directory, including line catalogs and 1D spectra figures.

----

## Commands

The following is a list of commands for the interactive inspection.

**OBJECT SPECIFIC OPTIONS:**

a = accept object fit

ac = accept object fit, noting contamination

r = reject object

c = add comment

user = toggle between previously saved fits

contam = specify contamination to line flux and/or continuum

s = print the (in progress) object summary

list = list all the objects in line list

left = list all the objects left to inspect

reset = reset interactive options back to default for this object

ex = prints out existing redshift(s) for this object, or states it doesn't have one

**EMISSION LINE SPECIFIC OPTIONS:**

z = enter a different z guess

w = enter a different emission line wavelength guess

dz = change the allowable redshift difference between lines

n = skip to next brightest line found in this object

2gauss = double gaussian profile for the line being fitted

1gauss = option to go back to 1 gaussian fit after selecting 2 gaussian fit

ha, hb, hg, o31, o32, o2, s2, s31, s32, lya, c4, pa, pb, pg, lya, c4, pah, bra = change strongest emission line

The full list of commands and corresponding lines for the strongest emission line is as follows:

| **Command** | **Line**       | **Vacuum Wavelength (Å)** |
| ----------- | -------------- | ------------------------- |
| lya         | Ly-alpha 1215  | 1215.670                  |
| c4          | CIV 1548       | 1548.203                  |
| o2          | [OII] 3730     | 3729.875                  |
| hg          | H-gamma 4342   | 4341.684                  |
| hb          | H-beta 4863    | 4862.683                  |
| o31         | [OIII] 4959    | 4960.295                  |
| o32         | [OIII] 5007    | 5008.240                  |
| ha          | H-alpha 6563   | 6564.610                  |
| s2          | [SII] 6716     | 6718.290                  |
| s31         | [SIII] 9069    | 9071.100                  |
| s32         | [SIII] 9532    | 9533.200                  |
| he          | HeI 10830      | 10832.86                  |
| pg          | Pa-gamma 10941 | 10941.1                   |
| pb          | Pa-beta 12822  | 12821.6                   |
| pa          | Pa-alpha 18756 | 18756.1                   |
| pah         | PAH-3.3 micron | 32890.0                   |
| bra         | Br-alpha 40523 | 40523.0                   |

**SPECTRUM SPECIFIC OPTIONS:**

fw = change the fwhm guess in pixels

m1, m2, m3, to m8 = mask up to eight discontinuous wavelength regions

nodes = change the wavelengths for the continuum spline nodes

addnodes = add wavelengths for the continuum spline nodes

rmnodes = remove wavelengths from the continuum spline nodes

shiftallnodes = SHIFT ALL nodes used for the continuum spline by some wavelength

bluecut = change the blue cutoff of the spectrum

redcut = change the red cutoff of the spectrum

lincont = fit continuum as a line

polycont = fit continuum as a higher-order polynomial

splinecont = fit continuum as a spline (piecewise) polynomial

grismr = use only Grism-R spectrum for line-fitting (default)

grismrcontam = use only Grism-R spectrum (with contamination) for line-fitting

grismc = use only Grism-C spectrum for line-fitting

grismccontam = use only Grism-C spectrum (with contamination) for line-fitting

**DS9 SPECIFIC OPTIONS:**

lin = linear z-scale

log = logarithmic

dc = recenter direct images

reload = reload direct images and 2D spectra

**SOFTWARE SPECIFIC OPTIONS:**

h = print this message

q = quit

----

In installation, you may also consult **[install.sh](http://install.sh)** located in the top-level directory

You will need the following software and packages:

- *Python >=3.9*
- *numpy (1.26.4 - 2)*
- *astropy (6.0 - 7)*
- *pandas (2.0 - 3)*
- *scipy (1.14.1)*
- *matplotlib (3.6 - 3.11)*
- *SAO Image DS9*
- *XPA* Messaging System
- *Anaconda*
