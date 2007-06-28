#!/usr/bin/env python
import os, re, sys

class Verifier:
  def __init__(self):
    self.pylith    = '/PETSc3/cig/pylith3d/bin/pylith'
    self.petscOpts = ['--petsc.ksp_type=preonly', '--petsc.pc_type=lu', '--petsc.mat_type=aijmumps']
    return

  def compare(self, file1, file2):
    exp = re.compile(r'SCALARS \w+_verify_t. double 4')
    f1 = file(file1)
    found = False
    data1 = []
    for line in f1.readlines():
      if exp.match(line):
        found = True
      if found:
        data1.append(line)
    f1.close()
    f2 = file(file2)
    found = False
    data2 = []
    for line in f2.readlines():
      if exp.match(line):
        found = True
      if found:
        data2.append(line)
    f2.close()
    data1.sort()
    data2.sort()
    if data1 == data2:
      print file1,'matches',file2
    else:
      sys.exit('ERROR: Files do not match')
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
