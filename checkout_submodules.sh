#!/bin/bash
#Sets all submodules to appropriate branches

cd submodules
cd SigProc-Plott
git checkout dev

cd ../chronux
git checkout master

cd ../lib_MAScPhD_Matlab
git checkout master

cd ../lib_dav
git checkout master

cd ../..

