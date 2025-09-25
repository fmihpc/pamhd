# Introduction

The following is a sequence of `commands` and their explanations of how to install PAMHD. The commands should be copied and pasted into a terminal. It is assumed that each command succeeds, if not, you probably should **not** continue with the rest of the commands and instead figure out why the command failed. If you don't want to wait for every command to finish to see whether it succeeded you can start a new shell with e.g. `bash -e`, which will exit at the first error, and copy & paste all commands into the shell. You can request help by creating a new issue.

These instruction assume that [git](http://git-scm.com) is installed. A C++20 compiler is also required such as GCC-10 or later. Also make sure that `$HOME/bin` is in your `PATH` environment variable and `$HOME/lib` is in your `LD_LIBRARY_PATH` environment variable (for example in bash run `echo "export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH" >> $HOME/.bashrc`, similarly for PATH, and then log out and log back in).

# Installing prerequisites

As many of the prerequisites are not available in package repositories these instructions install them into your home directory. This also allows better control of the versions to download/install. If a package is available via the system's package manager you can edit paths to the install location of the package in the makefile (default is to use the paths given in makefiles/homedir).

### muparserx
```bash
cd $HOME
git clone https://github.com/beltoforion/muparserx.git
cd muparserx
cmake .
make
mkdir -p $HOME/lib
cp libmuparserx.a $HOME/lib/
mkdir -p $HOME/include
cp parser/*.h $HOME/include/
make clean
```

### Open MPI
If you don't have a supported version of MPI, install one into your home directory:

```bash
cd $HOME
wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.8.tar.bz2
tar xf openmpi-5.0.8.tar.bz2
cd openmpi-5.0.8
./configure --prefix=$HOME --enable-mpi-fortran=no
make
make install
```

### Zoltan
```bash
cd $HOME
wget https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz -O zoltan_distrib_v3.901.tar.gz
tar xf zoltan_distrib_v3.901.tar.gz
mkdir zoltan-build
cd zoltan-build
CC=mpicc CXX=mpic++ ../Zoltan-3.901/configure --prefix=$HOME --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong
make
make install
```

### Boost
```bash
cd $HOME
wget https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.bz2
tar xf boost_1_89_0.tar.bz2
cd boost_1_89_0
./bootstrap.sh
echo "using gcc : $(mpic++ --version|head -1|cut -f3 -d\ ) : mpic++ ;" >> user-config.jam
echo "using mpi : mpic++ : <define>B2MPIJUSTUSEMPI ;" >> user-config.jam
./b2 --user-config=user-config.jam
./b2 --user-config=user-config.jam --prefix=$HOME install
```

# Installing PAMHD

Download PAMHD and it's submodules:
```bash
cd $HOME
git clone --recursive https://github.com/fmihpc/pamhd.git
```

and compile it running (GNU) make from the PAMHD root directory:
```bash
cd pamhd
make
```

Optionally select another environment to compile for:
```bash
make ENVIRONMENT_MAKEFILE=makefiles/macosx_macports_llvm
```

The `makefiles` directory houses all environment dependent makefiles, if none correspond to your environment request a new one by creating a new issue.
