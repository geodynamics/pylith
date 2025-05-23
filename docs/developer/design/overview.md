(sec-developer-design-overview)=
# Design overview

The PyLith software suite is composed of a C++ library, Python modules, a Python application, and a few Python preprocessing and post-processing utilities.

## Boundary between Python and C++

We use the Pyre framework (written in Python) to collect all user parameters and to launch the MPI application.
As a result, the top-level code is written in Python.
In most cases there is a low-level C++ object of the same name with the low-level implementation of the object.

We limit the Python code to collection of the user parameters, some simple checking of the parameters, and passing the parameters to the corresponding C++ objects.
Everything else is done in C++.
This facilitates debugging (it is easier to track symbols in the C/C++ debugger) and unit testing, and reduces the amount of information that needs to be passed from Python to C++.
The PyLith application and a few other utility functions, like writing the parameter file, are limited to Python.
All other objects have a C++ implementation.
Objects that have user input collect the user input in Python using Pyre and pass it to a corresponding C++ object.
Objects that do not have user input, such as the integrators and constraints, are limited to C++.

:::{important}
Consistent inheritance between C++ and Python is important in order for SWIG to generate a Python interface that is consistent with the C++ interface.
:::
