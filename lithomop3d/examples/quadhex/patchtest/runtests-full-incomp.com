lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt1.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt1-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt2.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt2-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt3.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt3-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt4.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt4-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt5.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt5-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt6.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt6-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt7.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt7-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt8.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt8-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt9.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt9-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt10.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt10-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt11.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt11-full-incomp.ascii
/bin/rm pt1.plot
lithomop3dapp.py --lm3dscan.fileRoot=pt1 --lm3dscan.keywordEqualsValueFile=pt-full.keyval --lm3dscan.bcInputFile=pt12.bc --lm3dscan.materialPropertiesInputFile=pt-incomp.prop --lm3dscan.asciiOutputFile=pt12-full-incomp.ascii
/bin/rm pt1.plot
