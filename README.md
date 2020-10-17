# casm-utilities
A collection of utilities that make use of the [CASM](https://github.com/prisms-center/CASMcode) libraries. The repository also includes its own library, which wraps around CASM to provide a useable interface, as well as python modules that can call said interface interactively.

## Prerequisites 
This project relies heavily on sections of [CASM](https://github.com/prisms-center/CASMcode), which will be internally compiled.
You will need a compiler with `c++17` support (at least `g++-8`, or `clang-6`).

casm-utilities uses autotools to build and install everything, so if you're installing via git cloning, a few additional packages might need to be installed on your computer:

* autoconf
* automake
* libtool
* autoconf-archive

These should all be readily available via `brew` and `apt-get`.

## Installation
The repository comes with two main components:

* The utilities library (c++)
* Python wrappers for the c++ utilities library

If you want to disable the python module installation, you can use the `--disable-casmutils-python` flag during the configure step. By default, a shared library will be installed. You can force static linkage with the `--disable-shared` flag.

There is a collection of utilities available that use these libraries, which are no longer hosted on this repository.
Each utility exists in an independent module, which can be plugged into the build, as described below.

### Cloning the repository
This repository includes a few submodules that are needed in order to fully compile everything.
If you plan on cloning the repository, be sure to do it recursively:
```bash
git clone --recurse-submodules
```

You should be mindful when switching branches as well, since differnt branches may be on different submodule commits.
When checking out a different branch, remember to follow with
```bash
git submodule update
```
or possibly
```bash
git submodule update --init
```

### Preparing plugins for compilation
When [properly structured](https://github.com/goirijo/casm-utilities-plugin/blob/main/README.md), independent repositories can piggy back onto the build chain, and be linked agains the `casm-utilities` library.
To install a plugin, simply clone its repository inside the `plugins` directory, *before* beginning the compilation steps outlined below.

For example, the [`primify`](https://github.com/goirijo/casm-utilities-primify) plugin can be targeted for compilation with:
```bash
cd plugins
git clone https://github.com/goirijo/casm-utilities-primify
```

If your plugin depends on unstable library features, or is reaching directly into `casm` for calls, it's highly recommended you install without shared linking (`--disable-shared` at the configure step). If you have cloned plugins but want to supress their compilation, you can use the `--disable-plugins` flag during the configure step. 

### Generating the configure script
Once you've cloned the repository, and also the repository of any plugins you want to install, you'll have to generate the `configure` script.
You can do this with:
```bash
./boostrap.sh
```

### Configure the compile environment
Skip this section entirely if you're familiar with the standard `./configure && make && make install` procedure.
If you're unsure what do do, just follow the recommended steps here.

Begin by creating a build directory, this is where all the compilation will take place:
```bash
mkdir build
cd build
```

You are now ready to run the configure script, which will generate a `Makefile`.
This is the step where you get to specify any compilation flags you might want, as well as where to install everything provided by casm-utilities.
You can bypass the need for admin privileges by using the `--prefix` option (recommended): `--prefix=$HOME/.local`

Once you're set on what flags you need, put it all together to run the configure script:
```bash
../configure --prefix=$HOME/.local CXXFLAGS='-Any-flags -You -Might-want'    
```

Note that the python modules have been written for `python3`. If your default python version is still lagging, you'll have to specify `PYTHON_VERSION` during the configure step above. An example of what the configure step could look like:
```bash
CXX=g++-8 PYTHON_VERSION=3.6 ../configure --prefix=$HOME/.local CXXFLAGS='-O3 -DNDEBUG'    
```

For more information on the configure step, see the documentation by calling `../configure -h`.

### Make and install casm-utilities/plugins
This one is easy:
```
make && make install
```

Be sure to have your `PATH` environment variable pointing at whatever value you chose for your `--prefix`, with `/bin` appended to it. You can take advantage of your processors with the `-j` flag.
You can build things insanely fast on the computer cluster with that flag, but you might annoy every other student on the same login node.
