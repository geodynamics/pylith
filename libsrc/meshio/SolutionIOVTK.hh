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

#if !defined(pylith_solutionio_solutioniovtk_hh)
#define pylith_solutionio_solutioniovtk_hh

#include "SolutionIO.hh" // ISA SolutionIO

namespace pylith {
  namespace meshio {
    class SolutionIOVTK;
  } // meshio
} // pylith

class pylith::meshio::SolutionIOVTK : public SolutionIO
{ // SolutionIOVTK

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  SolutionIOVTK(void);

  /// Destructor
  virtual
  ~SolutionIOVTK(void);

  /** Set filename for VTK file.
   *
   * @param filename Name of VTK file.
   */
  void filename(const char* filename);

  /** Write solution to file.
   *
   * @param mesh PETSc mesh object
   */
  void write(const ALE::Obj<ALE::Mesh>& mesh);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
public :

  std::string _filename; ///< Name of VTK file.

}; // SolutionIOVTK

#endif // pylith_solutionio_solutioniovtk_hh

// End of file 
