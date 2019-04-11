#!/usr/bin/env python

import numpy
# import pdb
# pdb.set_trace()

logFilePrefs2D = ['drucker_prager_2d','drucker_prager_2d_nonassoc',
				  'drucker_prager_2d_tensile','elastic_2d','gen_max_2d',
				  'maxwell_2d','power_law_2d']
jacobianNames2D = ['dp_2d','dp_2d_nonassoc','dp_2d_tensile','el_2d','gm_2d',
				   'mx_2d','pl_2d']
logFilePrefs3D = ['drucker_prager_3d','drucker_prager_3d_nonassoc',
				  'drucker_prager_3d_tensile','elastic_3d','gen_max_3d',
				  'maxwell_3d','power_law_3d']
jacobianNames3D = ['dp_3d','dp_3d_nonassoc','dp_3d_tensile','el_3d','gm_3d',
				   'mx_3d','pl_3d']

solveSuff = '_solve'
jacobianSuff = '_jacobian'
fdPref = 'finite_diff_'
pyPref = 'pylith_'
diffPref = 'diff_'

nonSymmetricString = 'Jacobian is not symmetric!'

numEquations2D = 16
numEquations3D = 64

tolSymmetric = 1.0e-6
tolZero = 1.0e-10

#------------------------------------------------------------------------------
def parseJacobianFile(fileName, numEquations):
	"""
	Function to extract one or more Jacobians from a text file.
	"""
	f = open(fileName, 'r')
	lines = f.readlines()
	jacobianStartInds = [i for i, x in enumerate(lines) if 'row 0:' in x]
	numJacobians = len(jacobianStartInds)
	jacobians = []
	for jacobianNum in range(numJacobians):
		jacobian = parseJacobianLines(lines[jacobianStartInds[jacobianNum]:],
									  numEquations)
		jacobians.append(jacobian)

	f.close()

	return jacobians


def parseJacobianLines(lines, numEquations):
	"""
	Function to extract Jacobian from lines of text.
	"""
	jacobian = numpy.zeros((numEquations, numEquations), dtype=numpy.float64)

	for rowNum in range(numEquations):
		row = parseRow(lines[rowNum], numEquations)
		jacobian[rowNum] = row

	return jacobian


def parseRow(line, numEquations):
	"""
	Function to parse a row of value pairs enclosed in parentheses.
	"""
	row = numpy.zeros(numEquations, dtype=numpy.float64)
	lineSplit = line.split(' (')
	numSplits = len(lineSplit)
	entries = [lineSplit[i].replace(')','') for i in range(1,numSplits)]
	numEntries = len(entries)
	for entryNum in range(numEntries):
		entrySplit = entries[entryNum].split(',')
		index = int(entrySplit[0])
		val = float(entrySplit[1].strip())
		row[index] = val

	return row


def checkModel(modelNum, modelType, numEquations):
	"""
	Function to check symmetry and equivalence between computed and finite
	difference Jacobian.
	"""
	modelNames = logFilePrefs2D
	modelShortNames = jacobianNames2D
	if (modelType == '3D'):
		modelNames = logFilePrefs3D
		modelShortNames = jacobianNames3D

	logFileName = modelNames[modelNum] + jacobianSuff + '.log'
	pyJacobianName = pyPref + modelShortNames[modelNum] + jacobianSuff + '.txt'
	fdJacobianName = fdPref + modelShortNames[modelNum] + jacobianSuff + '.txt'
	diffJacobianName = diffPref + modelShortNames[modelNum] + jacobianSuff + \
					   '.txt'

	print ""
	print "Model %s:" % modelNames[modelNum]

	# Test symmetry from PETSC.
	f = open(logFileName, 'r')
	lines = f.readlines()
	petscSymmetry = not any(nonSymmetricString in x for x in lines)
	f.close()
	print "  Symmetry from PETSc:                           %s" % petscSymmetry

	# Get Jacobians and differences.
	pyJacobians = parseJacobianFile(pyJacobianName, numEquations)
	fdJacobians = parseJacobianFile(fdJacobianName, numEquations)
	diffJacobians = parseJacobianFile(diffJacobianName, numEquations)
	numJacobians = len(pyJacobians)
	numFdJacobians = len(fdJacobians)
	numDiffJacobians = len(diffJacobians)
	if (numJacobians != numFdJacobians or numJacobians != numDiffJacobians):
		msg = 'Number of Jacobians is different in different files.'
		raise ValueError(msg)

	# Loop over number of Jacobian arrays:
	print "  Number of Jacobian arrays:                     %d" % numJacobians
	for jacobianNum in range(numJacobians):
		print "    Jacobian number                   %d:" % jacobianNum
		(pySymmetry, maxDiff, maxDiffPct,
		 meanDiff, meanDiffPct) = checkJacobian(
			pyJacobians[jacobianNum], diffJacobians[jacobianNum])
		print "      Symmetry from numpy:                       %s" % pySymmetry
		print "      Mean difference from PyLith:               %g" % meanDiff
		print "      Mean percentage difference from PyLith:    %g" % meanDiffPct
		print "      Maximum difference from PyLith:            %g" % maxDiff
		print "      Maximum percentage difference from PyLith: %g" % maxDiffPct
		
	return


def checkJacobian(pyJacobian, diffJacobian):
	"""
	Function to check Jacobian symmetry and find max differences between
	PyLith Jacobian and finite difference Jacobian.
	"""

	numpySymmetry = numpy.allclose(pyJacobian, pyJacobian.transpose(),
								   atol=tolSymmetric)
	absDiff = numpy.absolute(diffJacobian)
	absJac = numpy.absolute(pyJacobian)
	useInds = numpy.where(absJac > tolZero)
	maxDiff = numpy.amax(absDiff)
	pctDiff = 100.0*(absDiff[useInds]/absJac[useInds])
	maxDiffPercent = numpy.amax(pctDiff)
	meanDiff = numpy.mean(absDiff[useInds])
	meanDiffPct = numpy.mean(pctDiff)

	return (numpySymmetry, maxDiff, maxDiffPercent, meanDiff, meanDiffPct)
	
	
#------------------------------------------------------------------------------
num2DModels = len(logFilePrefs2D)
num3DModels = len(logFilePrefs3D)

# Loop over 2D models:
for modelNum in range(num2DModels):
	checkModel(modelNum, '2D', numEquations2D)

# Loop over 3D models:
for modelNum in range(num3DModels):
	checkModel(modelNum, '3D', numEquations3D)
	
