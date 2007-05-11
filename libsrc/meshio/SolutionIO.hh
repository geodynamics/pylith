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

#if !defined(pylith_solutionio_solutionio_hh)
#define pylith_solutionio_solutionio_hh

#include "pylith/utils/sievefwd.hh" // USES ALE::Obj, ALE::Mesh

namespace pylith {
  namespace meshio {
    class SolutionIO;
  } // meshio
} // pylith

class pylith::meshio::SolutionIO
{ // SolutionIO

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :
  /// Constructor
  SolutionIO(void);

  /// Destructor
  virtual
  ~SolutionIO(void);

  /** Write solution to file.
   *
   * @param mesh PETSc mesh object
   */
  virtual
  void write(const ALE::Obj<ALE::Mesh>& mesh) = 0;

}; // SolutionIO

#endif // pylith_solutionio_solutionio_hh

// End of file 
