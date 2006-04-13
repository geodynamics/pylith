#!/bin/sh

pylith3dapp.py -pl3dscan.fileRoot=lintet/bm1a/bm1a -pl3dscan.bcInputFile=lintet/bm2a/bm2a.bc -pl3dscan.keywordEqualsValueFile=lintet/bm2a/bm2a.keyval -pl3dscan.timeStepInputFile=lintet/bm2a/bm2a.time -pl3dscan.asciiOutputFile=lintet/bm2a/bm2a.ascii -pl3dscan.ucdOutputRoot=lintet/bm2a/bm2a

# end of file
