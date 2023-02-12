#!/bin/sh
COMPILER=e20                 # SELECT COMPILER
VERSION=v09_66_02            # LOOK UP LATEST VERSION ON DUNESW GITHUB
DIRECTORY=larsoft_v09_66_02  # CHOOSE DIRECTORY NAME
USERNAME=username            # USE YOUR FNAL USERNAME
HDIR=/dune/app/users
#HDIR=/build

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
cd ${HDIR}/${USERNAME}
touch ${DIRECTORY}
rm -rf ${DIRECTORY}
mkdir ${DIRECTORY}
cd ${DIRECTORY}
mrb newDev -v ${VERSION} -q ${COMPILER}:prof
source ${HDIR}/${USERNAME}/${DIRECTORY}/localProducts*/setup
mkdir work
cd srcs

# checks out the develop versions of the split repositories
# use the -t <tag> option to check out a specific tag
# you can use mrb g dune_suite to get all the code below plus duneutil

mrb g dunecore
mrb g duneopdet
mrb g dunesim
mrb g dunecalib
mrb g duneprototypes
mrb g dunedataprep
mrb g dunereco
mrb g duneana
mrb g duneexamples
mrb g protoduneana
mrb g dunesw
mrb uc

cd $MRB_BUILDDIR
mrbsetenv

# build the software stack.  Use -j<n> where n is the number of cores on the machine.
# using <n> too large (such as 16 on a dunegpvm machine), will run the computer out of memory
# the dune build nodes have 16 cores and enough memory to run the build with -j16

mrb i -j16
