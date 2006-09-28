pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt1.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt1-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt2.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt2-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt3.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt3-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt4.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt4-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt5.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt5-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt6.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt6-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt7.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt7-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt8.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt8-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt9.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt9-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt10.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt10-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt11.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt11-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --pl3dscan.fileRoot=pt1 --pl3dscan.keywordEqualsValueFile=pt-bbar.keyval --pl3dscan.bcInputFile=pt12.bc --pl3dscan.materialPropertiesInputFile=pt-incomp.prop --pl3dscan.asciiOutputFile=pt12-bbar-incomp.ascii
/bin/rm pt1.plot
