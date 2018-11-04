Hippocampal_circadian_analysis
==========================

Demo of code for analyzing hippocampal circadian rhythms

Overview
--------

See folder `feature_extraction` for a demo of initial signal processing (PSD and EMD) used to extract features
See folder `circadian_analysis` for a demo analyzing circadian rhythms in the data



Repo set up (Mac / Linux)
--------

### Clone main repo
	git clone --recursive git@github.com:davestanley/Hippocampal_circadian_analysis.git
  
### Update submodules
	cd Hippocampal_circadian_analysis
	git checkout master
	git submodule update --init --recursive
	
### Set all submodule to correct branches (optional)
	./checkout_submodules.sh

### When pulling, be sure to also update submodules
	git pull --recurse-submodules	


Getting started
--------
- `preprocessing/demo_script.m` - A demo of code for preprocessing data and extracting EEG rhythms
- `circadian_EEGrhythms/ratscript_FFT_thetadelta2_arr.m` - Extracts power in EEG frequency bands for later processing
- `circadian_EEGrhythms/run_allrats_ergodic.m` - Analyzes power in EEG frequency bands. Uses ergodicity assumption
