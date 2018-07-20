# casm-utilities
A collection of utilities that make use of the [CASM](https://github.com/prisms-center/CASMcode) libraries. The repository also includes its own library, which wraps around CASM to provide a useable interface, as well as python modules that can call said interface interactively.

## Getting started
This project relies heavily on [CASM](https://github.com/prisms-center/CASMcode), which you'll need to install before anything else. Additionally, you'll need to have [pybind11](https://github.com/pybind/pybind11) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) headers available on your system.

casm-utilities uses the autotools to build and install everything, so if you're installing via git cloning, a few additional packages might need to be installed on your computer:

* autoconf
* automake
* libtool (?)
* autoconf-archive

These should all be readily available via ```brew``` and ```apt-get```.

### Installation
The repository comes with three main components:

* The utilities library (c++)
* A collection of command line utilities
* Python wrappers for the c++ utilities library

Though it is recommended to simply install all three simultaneously, you can toggle any of these components out of the installation process by commenting out the approriate line ```Makemfile.am```.

#### Generate the configure script
If you're trying to install after cloning the git repository, you'll first have to generate the ```configure``` script using the provided script:
```
./boostrap.sh
```

#### Configure the compile environment
Skip this section entirely if you're familiar with the standard ```./configure && make && make install``` procedure. If you're unsure what do do, just follow the recommended steps here.

Begin by creating a build directory, this is where all the compilation will take place:
```
mkdir build
cd build
```

You are now ready to run the configure script, which will generate a ```Makefile```. This is the step where you get to specify any compilation flags you might want, point the compiler at any libraries you have installed in non standard locations (e.g. casm), as well as where to install everything provided by casm-utilities.

* The typical casm user probably has its libraries installed somewhere in their ```$HOME``` directory, typically ```$HOME/.local/lib```. Locate the the path to your ```libcasm.so``` shared library (extension may vary based on your system); you'll it to ```LDFLAGS```: ```LDFLAGS=-L$HOME/.local/lib```.
* If you've installed pybind11 and Eigen in standard locations, you shouldn't need any preprocessor flags. However, since they are header only libraries, installation is not strictly required. You can just point the configure script to their location in a similar matter to ```LDFLAGS```: ```CPPFLAGS="-I/path/to/pybind11/include -I/path/to/Eigen/include"```.
* You can bypass the need for admin privileges by using the ```--prefix``` option (recommended): ```--prefix=$HOME/.local```

Once you're set on what flags you need, put it all together to run the configure script:
```
LDFLAGS=-L$HOME/.local/lib CPPFLAGS="-I/path/to/pybind11/include -I/path/to/Eigen/include" ../configure --prefix=$HOME/.local
```

For more information on the configure step, see the documentation by calling ```../configure -h```.

#### Make and install casm-utilities
This one is easy:
```
make && make install
```

Be sure to have both your ```PATH``` and ```LD_LIBRARY_PATH``` environment variables pointing at whatever value you chose for your ```--prefix```, appending ```/bin``` and ```/lib``` for each variable respectively.
