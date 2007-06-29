#!/usr/bin/env python
import os, re, sys

class Verifier:
  def __init__(self):
    self.pylith    = '/Users/brad/tools/cig/gcc-4.0/bin/pylith'
#    self.petscOpts = ['--petsc.ksp_type=preonly', '--petsc.pc_type=lu', '--petsc.mat_type=aijmumps']
    self.petscOpts = ['--petsc.pc_type=asm']
    self.exp       = re.compile(r'SCALARS \w+_verify_t. double 4')
    return

  def readVTK(self, filename):
    f = file(filename)
    found = False
    data = []
    for line in f.readlines():
      if self.exp.match(line):
        found = True
      if found and not line.startswith('LOOKUP') and not line.startswith('SCALARS'):
        data.append(line)
    f.close()
    data.sort()
    return data


  def compare(self, file1, file2):
    def convertLine(line):
      parts = line.split()
      print parts
      return (int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]))
    data1 = self.readVTK(file1)
    data2 = self.readVTK(file2)
    print data1
    print data2
    if not data1 == data2:
      # Do full check
      for line1, line2 in zip(data1, data2):
        v1,x1,y1,z1 = convertLine(line1)
        v2,x2,y2,z2 = convertLine(line2)
        if not v1 == v2:
          sys.exit('ERROR: Nonmatching vertex sets')
        if abs(x1 - x2) > 1.0e-10:
          sys.exit('ERROR: Nonmatching x displacement, vertex %d, %g != %g' % (v1, x1, x2))
        if abs(y1 - y2) > 1.0e-10:
          sys.exit('ERROR: Nonmatching y displacement, vertex %d, %g != %g' % (v1, y1, y2))
        if abs(z1 - z2) > 1.0e-10:
          sys.exit('ERROR: Nonmatching z displacement, vertex %d, %g != %g' % (v1, z1, z2))
    print file1,'matches',file2
    return

  def run(self, cfgFilename):
    outputFiles = []
    for p in [1, 2]:
      basename   = os.path.splitext(os.path.basename(cfgFilename))[0]+'_p'+str(p)
      outputname = basename+'.vtk'
      outputFiles.append(basename+'_t0.vtk')
      cmd = '%s --nodes=%d %s --problem.formulation.output.output.filename=%s %s' % (self.pylith, p, ' '.join(self.petscOpts), outputname, cfgFilename)
      print 'Running',cmd
      os.system(cmd)
    for file1, file2 in zip(outputFiles[:-1], outputFiles[1:]):
      self.compare(file1, file2)
    return

if __name__ == '__main__':
  if not len(sys.argv) == 2:
    print 'Usage: verify.py <cfg file>'
    sys.exit()
  Verifier().run(sys.argv[1])
