import hps
import numpy as np
import cosmology_distance

#this will submit all the cosmos laes with photometry

hpslist = hps.makehpslist(True)
hpslist = hps.fieldtrans(hpslist,'COSMOS','Ly{alpha}') + hps.fieldtrans(hpslist,"HDFN",'Ly{alpha}')

numtofit = len(hpslist)

for i,obj in enumerate(hpslist):
    
    maxage = cosmology_distance.tz(obj.z,0.7,0.7,0.3,0.0) #this is in years
    
    logmaxage = np.log(maxage)
    
    obj.sedparaminfo['Mas']=[16.,9.,35.,0.5]
    obj.sedparaminfo['Age']=[16.,12.,logmaxage,0.5]
    obj.sedparaminfo['Tau']=[-0.100,-5.,15.,0.0]
    obj.sedparaminfo['EBV']=[0.1,0.0,1.0,0.04]
    obj.sedparaminfo['Met']=[-0.7,-2.29,0.45,0.0]
    obj.sedparaminfo['Red']=[obj.z,0.0,1.5,0.0]
    if (i+1) == numtofit:
        obj.sedfit('o2default.ini',4,True)
    else:
        obj.sedfit('o2default.ini',4)
