pylith3dapp.py --bcInputFile=pt1.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt1-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt2.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt2-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt3.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt3-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt4.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt4-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt5.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt5-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt6.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt6-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt7.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt7-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt8.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt8-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt9.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt9-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt10.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt10-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt11.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt11-full-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt12.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt12-full-incomp.ascii
/bin/rm pt1.plot
