#this module will SED fit the 3DHST objects in goods-s
#written by alex hagen
#hagen@psu.edu

import numpy as np
import sedtools
from AB2uJy import AB2uJy
#import cosmology_distance

def get_data(filename = "photom_goodss.dat"):
    """
    reads in 3dhst data from filename
    """
    #data is 30 columns -- first is string, rest are floats
    datacols = 30
    
    ids = np.loadtxt(filename, dtype=str, usecols=[0], converters = {0 : lambda x : x.decode()})
    #TODO: need to check IDS to make sure there are no duplicates
    data = np.loadtxt(filename, usecols = range(1,datacols))
    #correcting for some odd negative values in uncertainties
    data[:,3:] = np.abs(data[:,3:])
    
    galaxies = []
    f435wfilter = sedtools.Filter("filtercurves/acsF435W.res")
    f606wfilter = sedtools.Filter("filtercurves/acsF606W.res")
    f775wfilter = sedtools.Filter("filtercurves/acsF775W.res")
    f814wfilter = sedtools.Filter("filtercurves/wfcF814W.res")
    f850lpfilter = sedtools.Filter("filtercurves/acsF850LP.res")
    f105wfilter = sedtools.Filter("filtercurves/wfc3F105W.res")
    f125wfilter = sedtools.Filter("filtercurves/wfc3F125W.res")
    f160wfilter = sedtools.Filter("filtercurves/wfc3F160W.res")
    iracch1filter = sedtools.Filter("filtercurves/iracch1.res")
    iracch2filter = sedtools.Filter("filtercurves/iracch2.res")
    iracch3filter = sedtools.Filter("filtercurves/iracch3.res")
    iracch4filter = sedtools.Filter("filtercurves/iracch4.res")
    
    for i in range(len(ids)):
        fluxlist = []
        name = ids[i]
        line = data[i]
        ra = line[0]
        dec = line[1]
        hbeta = line[2]
        hbetaerr = line[3]
        redshift = line[4]
        f435w = AB2uJy(line[5], line[6])
        fluxlist.append(sedtools.Flux(f435w[0],f435w[1], f435wfilter))
        f606w = AB2uJy(line[7], line[8])
        fluxlist.append(sedtools.Flux(f606w[0],f606w[1], f606wfilter))
        f775w = AB2uJy(line[9], line[10])
        fluxlist.append(sedtools.Flux(f775w[0],f775w[1], f775wfilter))
        f814w = AB2uJy(line[11], line[12])
        fluxlist.append(sedtools.Flux(f814w[0],f814w[1], f814wfilter))
        f850lp = AB2uJy(line[13], line[14])
        fluxlist.append(sedtools.Flux(f850lp[0],f850lp[1], f850lpfilter))
        f105w = AB2uJy(line[15], line[16])
        fluxlist.append(sedtools.Flux(f105w[0],f105w[1], f105wfilter))
        f125w = AB2uJy(line[17], line[18])
        fluxlist.append(sedtools.Flux(f125w[0],f125w[1], f125wfilter))
        f160w = AB2uJy(line[19], line[20])
        fluxlist.append(sedtools.Flux(f160w[0],f160w[1], f160wfilter))
        iracch1 = AB2uJy(line[21], line[22])
        fluxlist.append(sedtools.Flux(iracch1[0],iracch1[1], iracch1filter))
        iracch2 = AB2uJy(line[23], line[24])
        fluxlist.append(sedtools.Flux(iracch2[0],iracch2[1], iracch2filter))
        iracch3 = AB2uJy(line[25], line[26])
        fluxlist.append(sedtools.Flux(iracch3[0],iracch3[1], iracch3filter))
        iracch4 = AB2uJy(line[27], line[28])
        fluxlist.append(sedtools.Flux(iracch4[0],iracch4[1], iracch4filter))
        gal = sedtools.Galaxy(name, fluxlist, redshift)
        gal.ra = ra
        gal.dec = dec
        gal.hbeta = hbeta
        gal.hbetaerr = hbetaerr
        galaxies.append(gal)
        
    return galaxies

def run():
    numchains = 4      
    gals = get_data()
    galmc = sedtools.GalMC("o2default.ini", "SED.zip")
    preamble =  ['#PBS -l nodes=1','#PBS -l walltime=96:00:00', 
                 '#PBS -q lionxf-yuexing','#PBS -j oe', 'cd $PBS_O_WORKDIR',
                 'echo " "' , 'echo "Job started on `hostname` at `date`"']
    postamble = ['echo "Job Ended at `date`"', 'echo " "']
    pbs = sedtools.PBSBackend(preamble, postamble)
    
    for i in gals:
        maxage = cosmology_distance.tz(i.z,0.7,0.7,0.3,0.0) #this is in years
        
        logmaxage = np.log(maxage)
        sedparaminfo = {}
        sedparaminfo['Mas']=[16.,9.,35.,0.5]
        sedparaminfo['Age']=[16.,12.,logmaxage,0.5]
        sedparaminfo['Tau']=[-0.100,-5.,15.,0.0]
        sedparaminfo['EBV']=[0.1,0.0,1.0,0.04]
        sedparaminfo['Met']=[-0.7,-2.29,0.45,0.0]
        sedparaminfo['Red']=[i.z,0.0,1.5,0.0]
        #submit(self, backend, galaxy, numchains, sedparaminfo):
        galmc.submit(pbs, i, numchains, sedparaminfo)


def getsedfits():
    gals = get_data()
    galmc = sedtools.GalMC("o2default.ini","SED.zip")
    getdist = sedtools.GetDist()
    for i in gals:
        galmc.cleanphotometry(i)
        getdist.get_sed_fit(i)
    return gals




