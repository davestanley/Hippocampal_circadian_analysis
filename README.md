Hippocampal_circadian_analysis
==========================

Demo of code for analyzing hippocampal circadian rhythms

Overview
--------

See folder `feature_extraction` for a demo of initial signal processing (PSD and EMD) used to extract features
See folder `circadian_analysis` for a demo analyzing circadian rhythms in the data



Getting started (Mac / Linux)
--------

### Clone main repo

	git clone --recursive git@github.com:davestanley/Hippocampal_circadian_analysis.git
  
### Update submodules
	cd Hippocampal_circadian_analysis
	git checkout master
	git submodule update --init --recursive
	
### Set all submodule to correct branches (optional)
	./pull_submodules.sh

### When pulling, be sure to also update submodules:
	git pull --recurse-submodules	


	
