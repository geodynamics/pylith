// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// ----------------------------------------------------------------------
// const Field<SubMesh>::Mesh (SubMesh)
%typemap(in) (const pylith::topology::Field<pylith::topology::SubMesh>::Mesh& mesh)
{
  pylith::topology::Field<pylith::topology::SubMesh>::Mesh* submesh = 0;
  int err = SWIG_ConvertPtr($input, (void**) &submesh, 
			    $descriptor(pylith::topology::Field<pylith::topology::SubMesh>::Mesh*),
				0);
  if (SWIG_IsOK(err))
    $1 = (pylith::topology::Field<pylith::topology::SubMesh>::Mesh*) submesh;
  else {
    err = SWIG_ConvertPtr($input, (void**) &submesh, 
			  $descriptor(pylith::topology::SubMesh*),
			  0);
    if (SWIG_IsOK(err))
      $1 = (pylith::topology::Field<pylith::topology::SubMesh>::Mesh*) submesh;
    else {
      PyErr_SetString(PyExc_TypeError, "Argument must be a submesh.");
      return NULL;
    } // if
  } // if/else
} // typemap(in)


// End of file
