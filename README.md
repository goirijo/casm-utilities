# casm-utilities
A collection of utilities that make use of the [CASM](https://github.com/prisms-center/CASMcode) libraries.
The repository also includes its own library, which wraps around CASM to provide a useable interface, as well
as python modules that can call said interface interactively.

## Getting started
This project relies heavily on [CASM](https://github.com/prisms-center/CASMcode), which you'll need to install
before anything else. Additionally, you'll need to have [pybind11](https://github.com/pybind/pybind11) and
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) headers available on your system.

casm-utilities uses the autotools to build and install everything, so if you're installing via git cloning,
a few additional packages might need to be installed on your computer:

* autoconf
* automake
* libtool (?)
* autoconf-archive

These should all be readily available via ```brew``` and ```apt-get```.

###Installation
The repository comes
