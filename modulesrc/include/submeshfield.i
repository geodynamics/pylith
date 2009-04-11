// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
