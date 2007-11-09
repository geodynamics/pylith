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

#include <iosfwd> // HOLDSA std::ofstream

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
  ~SolutionIOVTK(void);

  /** Set filename for VTK file.
   *
   * @param filename Name of VTK file.
   */
  void filename(const char* filename);

  /** Get filename of VTK file.
   *
   * @returns filename Name of VTK file.
   */
  const char* filename(void) const;

  /** Set time format for time stamp in name of VTK file.
   *
   * @param format C style time format for filename.
   */
  void timeFormat(const char* format);

  /** Get time format for time stamp in name of VTK file.
   *
   * @returns C Style time format for filename.
   */
  const char* timeFormat(void) const;

  /** Open output files.
   *
   * @param mesh PETSc mesh object
   */
  void open(const ALE::Obj<ALE::Mesh>& mesh);

  /// Close output files.
  void close(void);

  /** Write solution topology to file.
   *
   * @param mesh PETSc mesh object.
   * @param 
   */
  void writeTopology(const ALE::Obj<ALE::Mesh>& mesh,
		     const spatialdata::geocoords::CoordSys* csMesh);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param field PETSc field over vertices.
   * @param name Name of field.
   * @param mesh PETSc mesh object.
   */
  void writeVertexField(const double t,
			const ALE::Obj<real_section_type>& field,
			const char* name,
			const ALE::Obj<ALE::Mesh>& mesh);

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param field PETSc field over cells.
   * @param name Name of field.
   * @param mesh PETSc mesh object.
   */
  void writeCellField(const double t,
		      const ALE::Obj<real_section_type>& field,
		      const char* name,
		      const ALE::Obj<ALE::Mesh>& mesh);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
public :

  std::string _filename; ///< Name of VTK file.
  std::string _timeFormat; ///< C style time format for time stamp.

  PetscViewer _viewer; ///< Output file

}; // SolutionIOVTK

#include "SolutionIOVTK.icc" // inline methods

#endif // pylith_solutionio_solutioniovtk_hh

// End of file 
