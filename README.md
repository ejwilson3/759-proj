# discretize_geom
An attempted GPU implementation of a method to know where in a mesh volumes lie.

If the program were complete, you would need to do the following to run it:

Download and/or build hdf5: 

sudo apt-get install libhdf5-dev

or build from https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz


Download and build MOAB:

with the following shell commands on linux (I'm using Ubuntu 16.04)

mkdir moab

cd moab

git clone https://bitbucket.org/fathomteam/moab

cd moab

autoreconf -fi

cd ../

mkdir build

cd build

../moab/configure --enable-shared --enable-dagmc --with-hdf5=<PATH TO HDF5> --prefix=<PATH TO MOAB>

make

make install

export LD_LIBRARY_PATH=<PATH TO MOAB>:$LD_LIBRARY_PATH

export LIBRARY_PATH=<PATH TO MOAB>:$LIBRARY_PATH



Finally install discretize_geom

git clone https://github.com/ejwilson3/759-proj.git

cd 759-proj

mkdir build

cmake ../ -DMOAB_DIR=<PATH TO MOAB>/lib/cmake/MOAB/

make

run from the build directory like so:

./discretize_geom ../tests/unitbox.h5m 10000 0

# NOTE:
I'm not sure that will work. As I don't have access to CUDA on my own machine
and wasn't super keen on installing everything on Euler when I was nowhere near
done with my code, it is untested. I'm certain there will be errors when you try to make it, even if CMAKE goes off without a hitch.

# Files:
The file in the main directory is just an executable wrapper. All of the useful
code is in the src directory.
The tests directory contains two tests that could be used were it working.
The notes directory contains a couple of .txt files that I created while
working through this project.
