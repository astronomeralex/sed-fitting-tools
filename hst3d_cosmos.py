#this module will SED fit the 3DHST objects in cosmos
#written by alex hagen
#hagen@psu.edu

import numpy as np
import sedtools
#from AB2uJy import AB2uJy
import cosmology_distance

#def prep_data(filename="cosmos_photom_all.dat", outfile = "cosmos_photom_ujy.dat"):
#    """
#    preps the data
#    should run one time only
#    """
#    data = np.loadtxt(filename)
#    mags = data[:,5:]
#    outarr = np.zeros_like(mags)
#    #have 27 filters
#    assert len(mags[0]) == 54
#    for i,line in enumerate(mags):
#        for j in range(27):
#            mag = line[2*j]
#            err = line[2*j+1]
#            if mag == 0 and err == 0:
#                #set limit at 26
#                flux = 0
#                fluxerr = 0.145
#            elif mag == 99.:
#                flux = 0
#                fluxerr = AB2uJy(err,0.1)[0]
#            elif err == -1.:
#                flux = 0
#                fluxerr = AB2uJy(mag,0.1)[0]
#            else:
#                ujy = AB2uJy(mag,err)
#                flux = ujy[0]
#                fluxerr = ujy[1]
#            outarr[i,2*j] = flux
#            outarr[i,2*j+1] = fluxerr
#    np.savetxt(outfile,outarr,fmt="%.6f")

def get_data():
    """
    gets the data and puts it into galaxies
    """
    
    filterlist = []
    filternames = ["CFHTu.res",
    "SubB.res",
    "SubV.res",
    "Subg.res",
    "Subrp.res",
    "Subip.res",
    "Subzp.res",
    "SubNB816.res",
    "SubIB427.res",
    "SubIB464.res",
    "SubIB505.res",
    "SubIB574.res",
    "SubIB709.res",
    "SubIB827.res",
    "SubNB711.res",
    "CFHTK.res",
    "SubIB484.res",
    "SubIB527.res",
    "SubIB624.res",
    "SubIB679.res",
    "SubIB738.res",
    "SubIB767.res",
    "iracch1.res",
    "iracch2.res",
    "hst_acsF814W.res",
    "hst_wfc3F125W.res",
    "hst_wfc3F160W.res"]
    
    for i in filternames:
        filterlist.append(sedtools.Filter("filtercurves/" + i))
        
    gallist = []
    photometry = np.loadtxt("cosmos_photom_ujy.dat")
    fieldnum = np.loadtxt("cosmos_photom_all.dat",usecols=[0,1],dtype=str)
    ids = []
    for i in fieldnum:
        ids.append(i[0] + '_' + i[1])
    
    zs = np.loadtxt("cosmos_photom_all.dat",usecols=[2])
    for i,line in enumerate(photometry):
        fluxlist = []
        for j in range(27):
            fluxlist.append(sedtools.Flux(line[2*j],line[2*j + 1],filterlist[j],"uJy"))
        gal = sedtools.Galaxy(ids[i],fluxlist,zs[i])
        gallist.append(gal)
    return gallist
    
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

    
    
