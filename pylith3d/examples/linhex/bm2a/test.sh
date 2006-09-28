#!/bin/sh

pylith3dapp.py -pl3dscan.fileRoot=linhex/bm1a/bm1a -pl3dscan.bcInputFile=linhex/bm2a/bm2a.bc -pl3dscan.keywordEqualsValueFile=linhex/bm2a/bm2a.keyval -pl3dscan.timeStepInputFile=linhex/bm2a/bm2a.time -pl3dscan.asciiOutputFile=linhex/bm2a/bm2a.ascii -pl3dscan.ucdOutputRoot=linhex/bm2a/bm2a

# end of file
