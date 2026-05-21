# FRET

There are two complementary Jython scripts for [ImageJ](http://imagej.net/Welcome)/[Fiji](https://fiji.sc) version 2.14.0/1.54p (Rueden et al. 2017) to perform [FRET](https://en.wikipedia.org/wiki/Förster_resonance_energy_transfer) based analysis of cell and adhesion dynamics. The first script, FRET_LSM_Timelapse.py, constitutes the core FRET quantification pipeline and is designed for the pixel wise analysis of spectral confocal FRET images from time lapse series or stacks. It enables computation of FRET metrics (e.g., FRET index, A/D, or D/A ratios) over the entire field of view after background subtraction, photobleaching correction, and cell segmentation, and is applicable to a wide range of FRET based assays. The second script, FRET_Wound_Healing.py, is an optional, wound healing specific extension that uses a manually defined wound ROI to generate a 2D mesh of local sampling regions and to compute mean FRET values along the wound edge.

## License
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)<br>
FRET runs under the  [BSD-3 License](https://opensource.org/licenses/BSD-3-Clause)

---
This document describes the workings of the ImageJ/Fiji python scripts **FRET_LSM_Timelapse.py** and **FRET_Wound_Healing.py**, written by Philippe Girard ([email](philippe.girard@ijm.fr)).

