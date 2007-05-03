time pylith3dapp.py --fileRoot=bm1/bm1 --keywordEqualsValueFile=bm2/bm2.keyval --bcInputFile=bm2/bm2.bc --timeStepInputFile=bm2/bm2.time --asciiOutputFile=bm2/bm2.ascii --plotOutputFile=bm2/bm2.plot
