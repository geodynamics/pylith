pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt1.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt1-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt2.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt2-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt3.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt3-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt4.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt4-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt5.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt5-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt6.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt6-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt7.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt7-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt8.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt8-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt9.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt9-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt10.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt10-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt11.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt11-bbar-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-bbar.keyval --bcInputFile=pt12.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt12-bbar-incomp.ascii
/bin/rm pt1.plot
