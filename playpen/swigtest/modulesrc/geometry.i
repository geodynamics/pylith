// SWIG interface
%module(directors="1") geometry
%feature("director"); 

// Header files for module C++ code
%{
#include "Bar.hh"
#include "Sphere.hh"
#include "Scene.hh"
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  }
}

%include "Shape.i"
%include "Bar.i"
%include "Sphere.i"
%include "Scene.i"
