# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


subpackage = feassemble
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	PhysicsImplementation.hh \
	DSLabelAccess.hh \
	DSLabelAccess.icc \
	FEKernelKey.hh \
	FEKernelKey.icc \
	Integrator.hh \
	IntegratorDomain.hh \
	IntegratorBoundary.hh \
	IntegratorInterface.hh \
	IntegrationData.hh \
	InterfacePatches.hh \
	JacobianValues.hh \
	UpdateStateVars.hh \
	Constraint.hh \
	ConstraintSpatialDB.hh \
	ConstraintUserFn.hh \
	ConstraintSimple.hh \
	AuxiliaryFactory.hh \
	feassemblefwd.hh

dist_noinst_HEADERS =

dist_noinst_DATA = \
	jacobian2d_lgdeform.wxm \
	jacobian2d_nonsymm_lgdeform.wxm \
	jacobian3d_lgdeform.wxm \
	jacobian3d_nonsymm_lgdeform.wxm \
	tri3_elasticity.wxm \
	tet4_elasticity.wxm


# End of file
