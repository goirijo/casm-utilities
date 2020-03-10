# List of things to remember
* clang-format
* c++ vs py sections
* .ycm_extra_conf.py
* libraries vs utitilites
* How to document
* When it's OK to invoke CASM

# Naming conventions
* `is_*` for functions/functors returning booleans

# Python modules
* Always use numpy style documentation for every function and class you write (use snippets!).
* Directly binded modules start with leadint `_` (e.g. `_xtal`)
* Methods of directly binded classes get named by the signature (e.g. `_bring_within` and `_bring_within_const`) to dinstinguish between functions that mutate the class and funcitons that return a new instance.
* Once the class is totally binded, create modules in python that rename the methods.
* Usually you'll have to python types per c++ type: mutable vs immutabe. The immutable type is just the name of the class, the mutable one gets a `Mutable` prefix (e.g `Coordinate` vs `MutableCoordinate`).
* Both types have access to all the c++ binded methods, but each is only "supposed" to use the nicely named python functions.
