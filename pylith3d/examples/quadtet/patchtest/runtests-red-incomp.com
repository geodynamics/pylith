pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt1.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt1-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt2.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt2-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt3.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt3-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt4.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt4-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt5.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt5-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt6.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt6-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt7.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt7-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt8.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt8-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt9.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt9-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt10.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt10-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt11.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt11-red-incomp.ascii
/bin/rm pt1.plot
pylith3dapp.py --quadratureOrder=Reduced --bcInputFile=pt12.bc --materialPropertiesInputFile=pt-incomp.prop --asciiOutputFile=pt12-red-incomp.ascii
/bin/rm pt1.plot
