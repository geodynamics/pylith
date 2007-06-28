#!/usr/bin/env python
import re, sys


def compare(file1, file2):
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
    print 'Files match'
  else:
    sys.exit('ERROR: Files do not match')
  return

if __name__ == '__main__':
  if not len(sys.argv) == 3:
    print 'Usage: verify.py file1 file2'
    sys.exit()
  compare(sys.argv[1], sys.argv[2])
