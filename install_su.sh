#! /bin/sh
#
# install_su.sh
# Copyright (C) 2021 dlilien <dlilien@hozideh>
#
# Distributed under terms of the MIT license.
#
sudo apt-get install gfortran
export thisdir=$PWD
cd ..
export CWPROOT=$PWD/SeisUnix
git clone https://github.com/JohnWStockwellJr/SeisUnix.git
cd $CWPROOT/src 
mv configs/Makefile.config_Linux_x86_64 Makefile.config
touch LICENSE_44R14_ACCEPTED
touch MAILHOME_44R14
make install <<< y
cd $thisdir
