// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/FaultCohesive.hh
 *
 * @brief C++ abstract base class for a fault interface implemented with
 * cohesive cells.
 */

#if !defined(pylith_faults_faultcohesive_hh)
#define pylith_faults_faultcohesive_hh

#include "faultsfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include <string> // HASA std::string

class pylith::faults::FaultCohesive : public pylith::problems::Physics {
    friend class TestFaultCohesive; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesive(void);

    /// Destructor.
    virtual ~FaultCohesive(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set identifier for fault cohesive cells.
     *
     * @param[in] value Fault identifier
     */
    void setInterfaceId(const int value);

    /** Get identifier for fault cohesive cells.
     *
     * @returns Fault identifier
     */
    int getInterfaceId(void) const;

    /** Set label marking surface of interface.
     *
     * @param[in] value Label of surface (from mesh generator).
     */
    void setSurfaceMarkerLabel(const char* value);

    /** Get label marking surface of interface.
     *
     * @returns Label of surface (from mesh generator).
     */
    const char* getSurfaceMarkerLabel(void) const;

    /** Set label marking buried edges of interface surface.
     *
     * @param[in] value Label of buried surface edge (from mesh generator).
     */
    void setBuriedEdgesMarkerLabel(const char* value);

    /** Get label marking buried edges of interface surface.
     *
     * @returns Label of buried surface edge (from mesh generator).
     */
    const char* getBuriedEdgesMarkerLabel(void) const;

    /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void setRefDir1(const PylithReal vec[3]);

    /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void setRefDir2(const PylithReal vec[3]);

    /** Adjust mesh topology for fault implementation.
     *
     * @param mesh[in] PETSc mesh.
     */
    void adjustTopology(pylith::topology::Mesh* const mesh);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _refDir1[3]; ///< First choice reference direction used to compute boundary tangential directions.
    PylithReal _refDir2[3]; ///< Second choice reference direction used to compute boundary tangential directions.

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    int _interfaceId; ///< Identifier for cohesive cells.
    std::string _interfaceLabel; ///< Label identifying vertices associated with fault.
    std::string _buriedEdgesLabel; ///< Label identifying vertices along buried edges of fault.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    FaultCohesive(const FaultCohesive&); ///< Not implemented
    const FaultCohesive& operator=(const FaultCohesive&); ///< Not implemented

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh

// End of file
