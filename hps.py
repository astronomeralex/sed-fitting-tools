#this will contain the hps class, in which each object from this class will 
#contain information on an emission line object from the hps survey.

import pickle
import os
import random
import subprocess
import cosmology_distance
import numpy as np
from glob import glob
import extinction_laws


class HpsObj(object):
    def __init__(self,hpsid,silent=False):
        """
        this method will initiate the class by taking in an hps id number and then
        fetching all the pertinent information about it. hps id should be an integer
        between 1 and 479
        """
        self.id = hpsid
        #with the id number, we can now go fetch other information about this emission line
        #detection
        #the tables have been read in and pickled for easy access
        #i will open each of them and extract the relevent info
        table3file = open("hps_table3.pickle")
        table3 = pickle.load(table3file)
        
        for i in table3:
            if i[0]==hpsid:
                table3line = i
                break
        self.rah = table3line[1]
        self.ram = table3line[2]
        self.ras = table3line[3]
        self.decd = table3line[4]
        self.decm = table3line[5]
        self.decs = table3line[6]
        self.wavelength = table3line[7]
        self.fwhm = table3line[8]
        self.snr = table3line[9]
        self.flux = table3line[10]
        self.fluxlowererr = table3line[11]
        self.fluxuppererr = table3line[12]
        self.sfwhm = table3line[13]
        self.sfwhmlowererr = table3line[14]
        self.sfwhmuppererr = table3line[15]
        self.same = table3line[16]
        
        self.raold = self.rah + ':' + self.ram + ':' + self.ras
        self.decolrd = self.decd + ':' + self.decm + ':' + self.decs
        
        self.radecimalold = 15*(float(self.rah) +float(self.ram)/60. + float(self.ras)/3600.)
        self.decdecimalold = float(self.decd) + float(self.decm)/60. + float(self.decs)/3600 
        
        #now the coordinates file
        coordsfile = open("hps_newcoords.pickle")
        coords=pickle.load(coordsfile)
        coordsfile.close()
        
        self.ra = coords[hpsid-1][0]
        self.dec = coords[hpsid-1][1]
        
        #now for table 4
        table4file = open("hps_table4.pickle")
        table4 = pickle.load(table4file)
        for i in table4:
            if i[0] == hpsid:
                table4line = i
                break
        
        self.counterpart = table4line[1]
        self.rmag = table4line[2]
        self.counterprob = table4line[3]
        self.ew = table4line[4]
        self.ewlowererr = table4line[5]
        self.ewuppererr = table4line[6]
        self.restew = table4line[7]
        self.restewlowererr = table4line[8]
        self.restewuppererr = table4line[9]
        self.transition = table4line[10]
        self.z = table4line[11]
        self.lyaprob = table4line[12]
        self.xraycounter = table4line[13]
        
        self.distances = cosmology_distance.distcalc(self.z,0.7,0.7,0.3,0.0)
        self.scale = self.distances['scale']
        self.luminosity = self.flux * 1e-17 * 4 * np.pi * (100*self.distances['dl'])**2 #dl is in meters and is converted to cm here
        self.lumuppererr = (self.fluxuppererr / self.flux) * self.luminosity
        self.lumlowererr = (self.fluxlowererr / self.flux) * self.luminosity
        
        #now i will search for the photometry
        self.photometry = []
        names = ['COSMOS','HDFN','MUNICS','XMM']
        for name in names:
            #print name
            photfile = open(name + '.pickle')
            photlist = pickle.load(photfile)
            photfile.close
            
            filterfile = open(name + '_filters.pickle')
            filters = pickle.load(filterfile)
            filterfile.close()
            
            for i,line in enumerate(photlist):
                #print line[0]
                if line[0] == hpsid:
                    photline = line
                    photfilters = filters
                    break
            if 'photline' in locals():
                break
        #there is a possibility photometry will not be found
        if 'photline' in locals():
            self.field=name
            for i,bandpass in enumerate(photfilters):
                self.photometry.append(Flux(photline[i*2 + 1],photline[i*2 + 2],bandpass,silent))
            
            #if the object is an [OII] emitter, check to see if it has galex photometry
            if self.transition == '[OII]':
                galexfile = open('galex.phot.dict')
                galexdict = pickle.load(galexfile)
                galexfile.close()
                try:
                    galexphot = galexdict[self.id]
                    for i in galexphot:
                        self.photometry.append(Flux(i[1],i[2],i[0]))
                    
                except KeyError:
                    pass
                    
                if self.field == "COSMOS":
                    cosmosiracfile = open("o2.irac.cosmos.pickle")
                    cosmosirac = pickle.load(cosmosiracfile)
                    cosmosiracfile.close()
                    try:
                        cirac = cosmosirac[self.id]
                        for i in cirac:
                            self.photometry.append(Flux(i[1],i[2],i[0]))
                    except KeyError:
                        pass
                
            #if the object is a cosmos LAE, then add additional photometry from CANDELS and SCOSMOS
            if self.field == 'COSMOS' and self.transition == 'Ly{alpha}':
                laedictfile = open('COSMOS-lae-candels.dict.pickle')
                laedict = pickle.load(laedictfile)
                laedictfile.close()
                try:
                    candelsphot = laedict[self.id]
                    for i in candelsphot:
                        self.photometry.append(Flux(i[1],np.sqrt(i[2]**2 + (0.05*i[1]**2)),i[0],silent))
                       
                except KeyError:
                    #if not detected, use limits from Steve
                    #5-sigma "limits"
                    #F814W: 200 nJy
                    #F125W: 110 nJy
                    #F160W: 90 nJy
                    #self.photometry.append(Flux(0.0,0.200,'wfcF814W'))
                    #self.photometry.append(Flux(0.0,0.110,'wfc3F125W'))
                    #self.photometry.append(Flux(0.0,0.090,'wfc3F160W'))
                    pass
                    
                #now for scosmos photometry
                scosmosdictfile = open("cosmos-lya-irac-phot.dict.pickle")
                scosmosdict = pickle.load(scosmosdictfile)
                scosmosdictfile.close()
                scosmosphot = scosmosdict[self.id]
                for i in scosmosphot:
                    self.photometry.append(Flux(i[1],i[2],i[0],silent))
            
            if self.field == 'HDFN' and self.transition == 'Ly{alpha}':
                laedictfile = open('lae_candels_hdfn_photdict.pickle')
                laedict = pickle.load(laedictfile)
                laedictfile.close()
                try:
                    candelsphot = laedict[self.id]
                    for i in candelsphot:
                        self.photometry.append(Flux(i[1],np.sqrt(i[2]**2 + (0.05*i[1]**2)),i[0],silent))
                       
                except KeyError:
                    pass
                    
                #now for goods-n irac photometry
                iracdictfile = open("hdfn-lya-irac-phot.dict.pickle")
                iracdict = pickle.load(iracdictfile)
                iracdictfile.close()
                iracphot = iracdict[self.id]
                for i in iracphot:
                    self.photometry.append(Flux(i[1],i[2],i[0],silent))
            
 
            #now sort the photometry
            self.photometry.sort(key = lambda x: x.filter.central)
            
            #now i need to go through the photometry and make sure that all the photmetry is blue-ward
            #of 3 microns. this is so we avoid contamination from the 3.3 micron PAH feature.
            for i,fluxobj in enumerate(self.photometry):
                filterobj = fluxobj.filter
                if hasattr(filterobj,'long10'):
                    if filterobj.long10 > 30000*(1 + self.z):
                        self.photometry = self.photometry[:i]
                        break
                                            
        else:
            if not silent:
                print "photometry not found for " + str(hpsid)
            self.field=''
                
        #now i will get the half light radius
        hlfile = open('halflight.pickle')
        halflight = pickle.load(hlfile)
        hlfile.close()
        
        try:
            hlpix = halflight[self.id]
            hlarcsec = hlpix * 0.03
            self.halflight = self.scale * hlarcsec #in kpc
            
        except KeyError:
            pass
        
        #this attribute is used for SED fitting. it will be a dictionary
        #that contains the parameters to fit using GalMC    
        self.sedparaminfo={}

    def sedfit(self,paramfilename,numchains=4,email=False):
        """
        this method will SED fit an instance of HpsObj using the GalMC SED fitting
        paramfilename is a string that is the location of a file containing a bunch of GalMC
        parameters
        numchains is an integer that is the number of distinct chains to be run
        """
        #first, check to make sure that the object actually has associated photometry
        if len(self.photometry)==0:
            raise ValueError(" no photometry for HPSID " + str(self.hpsid))
        
        #and check to make sure the parafilename exists
        if not os.path.exists(paramfilename):
            raise SyntaxError(paramfilename + " doesn't exist!")
            
        # and a sanity check on numchains
        if numchains < 0 or numchains > 8:
            raise ValueError(" check the number of chains you want, please")
        
        #see if the sedparaminfo attribute of the object has information
        if len(self.sedparaminfo) == 0:
            raise ValueError("No sedparaminfo!")
            
        #ok! now for the fun part
        #def __init__(self,hpsid,hpsphot,paramfilename,numchains,sedparaminfo)
        sed = GalMC(self.id,self.photometry,paramfilename,numchains,self.sedparaminfo)
        sed.fit(email)
        
    def getsedfit(self,root='.',silent=False):
        """
        this method will fetch the results of the sed fit. root is the directory
        location that contains the folders that each contain the results of the
        sed fit run by sedfit. the results have to be analyzed by getdist prior
        to the running of this method
        this is designed so if the SED results don't exist, it will exit gracefully
        with a message and will not throw an error
        this is so it can be run on a large sample easily without try statements
        """
        
        folder = root + '/' + str(self.id)
        
        #first, check to make sure the folder exists and that files from getdist exist
        if os.path.exists(folder) and os.path.exists(folder + '/hps' + str(self.id) + '.margestats'):
            #read in coverge file -- currently will take R-1 values
            convfile = open(folder + '/hps' + str(self.id) + '.converge')
            convlines = convfile.readlines()
            convfile.close()
            
            #now i need to find the line before R-1 value
            for i,line in enumerate(convlines):
                if 'var(mean)/mean(var) for eigenvalues of covariance of means' in line:
                    linenum = i
                    break
            
            if 'linenum' in locals():        
                #time to get the R values
                Rm1s = []
                for i in xrange(1,50):
                    line = convlines[linenum + i].split()
                    if len(line) < 2:
                        break
                    else:
                        Rm1s.append(float(line[1]))
            
                self.Rm1s = Rm1s        
                #find the largest R value
                self.Rm1 = max(self.Rm1s)
                
            else:
                self.Rm1 = None
                self.Rm1s = None
                
            #now time to read in the marginalized stats
            margefile = open(folder + '/hps' + str(self.id) + '.margestats')
            margelines = margefile.readlines()
            margefile.close()
            
            sedstats = {}
            
            for i,line in enumerate(margelines):
                if i>0:
                    stats = line.split()
                    if len(stats) < 2:
                        break
                    else:
                        paramname = stats[-1]
                        sedstats[paramname] = {}
                        sedstats[paramname]['mean'] = float(stats[1])
                        sedstats[paramname]['stdev'] = float(stats[2])
                        sedstats[paramname]['lower1'] = float(stats[3])
                        sedstats[paramname]['upper1'] = float(stats[4])
                        sedstats[paramname]['lower2'] = float(stats[5])
                        sedstats[paramname]['upper2'] = float(stats[6])
            
            self.sedstats = sedstats
            
            #now get the best fit parameters from the likestats file
            likefile = open(folder + '/hps' + str(self.id) + '.likestats')
            likelines = likefile.readlines()
            likefile.close()
            
            likestats = {}
            
            for i,line in enumerate(likelines):
                stats = line.split()
                if i == 0:
                    self.mloglike = float(stats[-1])
                elif i>2 and len(stats) > 1:
                    paramname = stats[-1]
                    likestats[paramname] = {}
                    likestats[paramname]['best'] = float(stats[1])
                    likestats[paramname]['lower1'] = float(stats[2])
                    likestats[paramname]['upper1'] = float(stats[3])
                    likestats[paramname]['lower2'] = float(stats[4])
                    likestats[paramname]['upper2'] = float(stats[5])
                    
            self.likestats = likestats        
                    
            self.chisq = self.mloglike*2. #note that this is only valid for a flat prior
            numphot = len(self.photometry)
            numsedparams = len(self.sedstats) - 1 #since mass is included twice
            if numphot - 1 > numsedparams:
                self.redchisq = self.chisq / float(numphot - numsedparams - 1)
            else:
                if not silent:
                    print "Cannot calculate reduced chi squared"
            
            if not silent:
                print str(self.id) + ": found"
                
        else:
            if not silent:
                print str(self.id) + ": SED info or folder doesn't exist"
                
    def extinction_correct(self):
        """
        this will extinction correct the photometry out to rest frame 2.2 microns
        this assumes a Calzetti law and that the E(B-V) of the galaxy is already found via SED fitting
        also assumes Rv = 4.05
        """
        #make sure we have E(B-V)
        ebv = self.sedstats['E(B-V)']['mean']
        
        #first, correct ly-alpha luminosity
        self.fluxcor = extinction_laws.calzetti_stellar(self.flux,0.1216,ebv / 0.44) #ebv is over 0.44 from calzetti
        self.luminositycor = self.fluxcor * 1e-17 * 4 * np.pi * (100*self.distances['dl'])**2
        self.lumcor_uppererr = self.luminositycor * self.lumuppererr / self.luminosity
        self.lumcor_lowererr = self.luminositycor * self.lumlowererr / self.luminosity
        
        self.photcor = []
        for i in self.photometry:
            restwave = i.filter.central / ((1+ self.z) * 1e4) #converting to rest wavelength and from angstroms to microns
            if restwave > 0.12 and restwave < 2.2:
                corflux = extinction_laws.calzetti_stellar(i.flux,restwave,ebv)
                if i.flux == 0:
                    corfluxerr = i.err
                else:
                    corfluxerr = corflux * i.err / i.flux
                self.photcor.append(Flux(corflux,corfluxerr,i.filter.name))
            elif restwave > 2.2:
                break
        
    
    def lyasfr(self):
        """
        calculated the 
        """
        #gronwall et al 2007 section 6
        if self.transition == 'Ly{alpha}':
            self.lyasfrcor = self.luminositycor * 9.1e-43
            self.lyasfrcor_uppererr = self.lumcor_uppererr * 9.1e-43
            self.lyasfrcor_lowererr = self.lumcor_lowererr * 9.1e-43
            
            self.lyasfr = self.luminosity * 9.1e-43
            self.lyasfr_uppererr = self.lumuppererr * 9.1e-43
            self.lyasfr_lowererr = self.lumlowererr * 9.1e-43
            
            
    def uvsfr(self):
        #now derive the UV SFR using the average of all bands who's central rest 
        #wavelength is from 1500 - 2800 Angstroms
        uvflux = []
        uverr = []
        for i in self.photometry:
            restcentral = i.filter.central / (1 + self.z)
            if restcentral > 1500 and restcentral < 2800:
                uvflux.append(i.flux)
                uverr.append(i.err)
        if len(uvflux) > 0:
            uvflux = np.array(uvflux).mean()
            uverr = np.array(uverr)
            uverr = np.sqrt(np.sum(uverr**2)) / len(uverr)
            self.uvlumdensity = uvflux * 1e-29 * 4 * np.pi * (100*self.distances['dl'])**2 #this should be in ergs/s/Hz
            self.uverr = uverr * 1e-29 * 4 * np.pi * (100*self.distances['dl'])**2
            self.uvsfr = 1.4e-28 * self.uvlumdensity # in solar masses per year
            self.uvsfrerr = 1.4e-28 * self.uverr
            
        uvfluxcor = []
        uverrcor = []
        for i in self.photcor:
            restcentral = i.filter.central / (1 + self.z)
            if restcentral > 1500 and restcentral < 2800:
                uvfluxcor.append(i.flux)
                uverrcor.append(i.err)
        if len(uvfluxcor) > 0:
            uvfluxcor = np.array(uvfluxcor).mean()
            uverrcor = np.array(uverrcor)
            uverrcor = np.sqrt(np.sum(uverrcor**2)) / len(uverrcor)
            self.uvlumdensitycor = uvfluxcor * 1e-29 * 4 * np.pi * (100*self.distances['dl'])**2 #this should be in ergs/s/Hz
            self.uverrcor = uverrcor * 1e-29 * 4 * np.pi * (100*self.distances['dl'])**2
            self.uvsfrcor = 1.4e-28 * self.uvlumdensitycor # in solar masses per year
            self.uvsfrerrcor = 1.4e-28 * self.uverrcor
            
        
class Flux(object):
    """
    this flux class will contain a single photometric measurement
    """
    def __init__(self,flux,err,bandpass,silent=False):
        self.flux = flux
        self.err = err
        self.filter = Filter(bandpass,silent)
    
class Filter(object):
    """
    the class contains information on the filter, including the name,
    location of the file describing the transmission, and the central wavelength
    """
    def __init__(self,name,silent = False):
        self.name = name
        self.transfile = 'filtercurves/'+self.name + '.res' #the location of the transmission file
        #now check and see if the transmission file exists
        if os.path.exists(self.transfile):
            transfile = open(self.transfile)
            transmission = transfile.readlines()
            transfile.close()
            waves = []
            trans = []
            for i,line in enumerate(transmission):
                transmission[i]=line.split()
                transmission[i][0] = float(transmission[i][0])
                transmission[i][1] = float(transmission[i][1])
                waves.append(transmission[i][0])
                trans.append(transmission[i][1])
            self.central = waves[trans.index(max(trans))]
            trans = np.array(trans)
            waves = np.array(waves)
            #normalize transmission
            trans = trans / trans.max()
            self.transmission = trans
            self.waves = waves
            #now find longest wavelength at 10% transmission
            #this is useful for determining if NIR filters could be contaminated by
            #3.3um PAH feature
            self.long10 = waves[np.where(trans > 0.1)[0][-1]]
            
        else:
            if not silent:
                print "Can't find transmission file at " + self.transfile
        
    
    
class GalMC(object):
    """
    this class will do the heavy lifting for SED fitting. Each instance will
    write all the necessary files and transfer them to the correct places and so on
    """
    def __init__(self,hpsid,hpsphot,paramfilename,numchains,sedparaminfo):
        """
        hpsphot is the photometry attribute from the HpsObj class
        paramfilename is the file name that sets a bunch of GalMC parameters
        numchains in the number of chains to run
        """
        self.hpsid = hpsid
        self.phot = hpsphot
        self.paramfile = paramfilename
        self.numchains = numchains
        self.params = sedparaminfo

    def fit(self, email=False):
        """
        this method will run MCMCfit on the GalMC object
        if email is true, then the code will send an email at the end to notify the recipient
        that the job is done
        this should be used on the last hps object out of a group that needs to be fitted
        """
        #it assumes that it is SED fitting in a new directory and will overwrite previous sed versions
        #make the new directory
        rootdir = os.getcwd()
        hpsdir = str(self.hpsid)
        os.mkdir(hpsdir)
        os.chdir(hpsdir)
        
        #now we'll write the data file
        self.writedata()
        
        #now its time to write out the param files for each of the chains we want
        for i in xrange(self.numchains):
            if i>0:
                self.randparams()
                self.writeparams(i)
            else:
                self.writeparams(i)
        
        #now i will write out the pbs file
        pbsout = ['#PBS -l nodes=1','#PBS -l walltime=96:00:00','#PBS -q lionxf-yuexing','#PBS -j oe','cd $PBS_O_WORKDIR',
        'echo " "','echo "Job started on `hostname` at `date`"','cp ../SED.zip .','unzip SED.zip','rm SED.zip']
        for i in xrange(self.numchains):
            line = './MCMCfit ' + str(self.hpsid) + '_' + str(i) + '.ini'
            pbsout.append(line)
        if email:
            #write line for script to send email
            pbsout.append('cd ..')
            pbsout.append('python hpsemail.py')
            pbsout.append('cd $PBS_O_WORKDIR')
            
        pbsout.append('echo "Job Ended at `date`"')
        pbsout.append('echo " "')
            
        
        pbsfile = open('sed.pbs','w')
        for i in pbsout:
            pbsfile.write(i + '\n')
        pbsfile.close()
        
        proc = ['qsub','sed.pbs']
        
        subprocess.call(proc)
        
        os.chdir(rootdir)
       
    def writedata(self):
        """
        this method will write out the data file in standard GalMC format 
        wavelength flux density error (in microJy)
        """ 
        datafile = open('hps' + str(self.hpsid)+'.dat','w')
        for i in self.phot:
            centralwave = str(i.filter.central)
            flux = str(i.flux)
            err = str(i.err)
            outline = centralwave + ' ' + flux + ' ' + err + '\n'
            datafile.write(outline)
        datafile.close()

        
    def randparams(self):
        """
        ramdomly sets the starting parameters that can be varied uniformly between their 
        lower and upper values
        """     
        for i in self.params:
            if self.params[i][3] > 0: #this checks to see if the parameter has been hardcoded
                self.params[i][0] = random.uniform(self.params[i][1],self.params[i][2])
        

    def writeparams(self,chainnum):
        """
        this writes out the params file
        chainnum is an int that is the number of the chain
        """
        defaultparamsfile = open('../' + self.paramfile)
        defaultparams = defaultparamsfile.read()
        defaultparamsfile.close()
        output= []
        #each entry of output is a line in the params file
        output.append('chain root = ' + 'hps' + str(self.hpsid) + '_' + str(chainnum))
        output.append('Data File = hps' + str(self.hpsid) + '.dat')
        #fitting para
        for i in self.params:
            outline = i + ' = ' + str(self.params[i][0]) + ', ' + str(self.params[i][1]) + ', ' + str(self.params[i][2]) + ', ' + str(self.params[i][3])
            output.append(outline)
        #now for the filters
        for i,fluxobj in enumerate(self.phot):
            output.append('filter('+str(i+1)+') = ../' + fluxobj.filter.transfile)
        
        outputfile = open(str(self.hpsid) + '_' + str(chainnum) + '.ini','w')
        outputfile.write(defaultparams)
        for i in output:
            outputfile.write(i + '\n')
        outputfile.close()

def runGetDist(folders,distparamsfile='distparams.ini'):
    """
    this will run GetDist on the chains outputed from the fit method of the GalMC class. it also will rename the chain
    files so they work with getdist
    
    folders is a list where each element is a string that is the name of a folder where getdist should be run
    it is recommended that glob is used to make this list
    
    distparamsfile is the file name of the parameter file for getdist
    this file cannot have the chain root included in it
    the programs will write that so it works for each folder
    """
    
    distparams = open(distparamsfile).read()
    
    root = os.getcwd()
    
    for folder in folders:
        os.chdir(folder)
        #write out the distparams file
        distparout = open('distparams.ini','w')
        distparout.write('file_root = hps' + folder + '\n')
        distparout.write(distparams)
        distparout.close()
        
        #now its time to rename the chains
        for i in xrange(4):
            proc=['mv','hps' + folder + '_' + str(i) + '_chain.txt','hps' + folder + '_' + str(i) + '.txt']
            subprocess.call(proc)
            
        #now its time to run getdist!
        #make sure its in the .bashrc file
        
        subprocess.call(['getdist','distparams.ini'])

        #and now time to run the python files produced by getdist so we can have plotses

        pyfiles = glob('hps' + folder + '*.py')

        for i in pyfiles:
            execfile(i)
        
        os.chdir(root)
        
def getSEDparam(hpslist,paramname,sigma = 1,tobase10 = False):
    """
    this will take in a list of HpsObj objects, a paramname, which is a string, and a sigma
    sigma can either be 1 or 2 and will specify the type of error bars returned, 1 or 2 sigma
    this will return three numpy arrays, the parameter from the SED fit for each object, 
    lower errors and upper errors
    """

    param = []
    lower = []
    upper = []
    
    for i in hpslist:
        param.append(i.sedstats[paramname]['mean'])
        lower.append(i.sedstats[paramname]['lower' + str(sigma)])
        upper.append(i.sedstats[paramname]['upper' + str(sigma)])
        
    param = np.array(param)
    lower = np.array(lower)
    upper = np.array(upper)
    
    #coverting from values to differences from the mean
    lower = param - lower
    upper = upper - param
    
    if tobase10:
        param = param / np.log(10)
        lower = lower / np.log(10)
        upper = upper / np.log(10)
    
    return param,lower,upper


def getrm1(hpslist):
    """
    this will take in a list of hps objects and return the
    highest R-1 value. In this case that the R-1 value is None,
    the value will be set to 999
    """
    
    rm1s = []
    
    for i in hpslist:
        rm1 = i.Rm1
        if rm1 is None:
            rm1s.append(999)
        else:
            rm1s.append(rm1)
            
    return np.array(rm1s)


def makehpslist(silent=False):
    """
    this will make a list where each element is an emission line detection from the HPS.
    this is writtent to save a bit of time :)
    """
    
    hpslist = []
    for i in xrange(1,480):
        hpslist.append(HpsObj(i,silent))
    
    
    return hpslist
    
def getSEDfits(hpslist,root='.',silent=False):
    """
    this will take in a list of HpsObj objects
    and will retrieve their SED fits
    """
    
    for i in hpslist:
        i.getsedfit(root,silent)
        
    return hpslist
    
def selectseds(hpslist):
    """
    this will make a new list with only the objects that have the
    sedstats attribute, meaning they have SED fit information
    """    
    
    sedlist = []
    
    for i in hpslist:
        if hasattr(i,'sedstats'):
            sedlist.append(i)
                
    return sedlist
    
def selectconverge(hpslist,cutoff = 0.05):
    """
    this will take in a list of hps objects with SED fits
    it will return only the objects who's R-1 value are less
    than the given cutoff
    """
    
    sedlist=[]
    
    for i in hpslist:
        r = i.Rm1
        if r is not None:
            if r < cutoff:
                sedlist.append(i)
    
    return sedlist
    
def hps_start(root='.',cutoff=0.05,silent=True):
    """
    this is another helpful function in interactive sessions
    if will return a list of only the objects with SED fits
    out of the whole HPS sample.
    root is a directory string which is where it looks for SED fits
    cutoff if the maximal R-1 value to be included
    """
    
    hpslist = makehpslist(silent)
    hpslist = getSEDfits(hpslist,root,silent)
    hpslist = selectseds(hpslist)
    hpslist = selectconverge(hpslist,cutoff)
    
    return hpslist
   
    
def fieldtrans(hpslist,field,transition):
    """
    this will take in a list of hps objects, a string field
    and a string transition
    it will return only those objects that are in the field and
    have the transition specified
    """
    
    output = []
    
    for i in hpslist:
        if i.field == field and i.transition == transition:
            output.append(i)
            
    return output
    
    
def radec(hpslist,outputfilename):
    """
    this will take in a list of hps objects and a string that
    is the name of the output file to be created
    
    to the output file will be written the ra,dec of the hpsobjects
    in degrees
    """
    
    output = open(outputfilename,'w')
    for i in hpslist:
        ra = str(i.ra)
        dec = str(i.dec)
        output.write(ra + ' ' + dec + '\n')
        
    output.close()
    
    
def thingstoplot(hpslist):
    """
    gets arrays of things to plot and returns them
    """
    luminosity = []
    z = []
    rmag = []
    halflight = []
    restew = []
    
    for i in hpslist:
        luminosity.append(i.luminosity)
        z.append(i.z)
        rmag.append(i.rmag)
        halflight.append(i.halflight)
        restew.append(i.restew)
        
    luminosity = np.array(luminosity)
    z = np.array(z)
    rmag = np.array(rmag)
    halflight = np.array(halflight)
    restew = np.array(restew)
    
    return luminosity, z, rmag, halflight, restew
    

def checkcomplete(folders):
    """
    this will take in a list of folders and go through each one to see if the
    SED fits have complete. If they have, then the string killed will not be in
    the log file
    also, this assumes that one completed run means the whole directory is good
    so, it assumes that completed runs are run again and that an output file
    that finishes means the whole directory finished
    """

    for i in folders:
        outfiles = glob(i + '/sed.pbs.*')
        killedarr = []
        for filename in outfiles:
            outfile = open(filename)
            output = outfile.read()
            outfile.close()
            if 'killed' in output:
                killedarr.append(True)
            else:
                killedarr.append(False)
        if len(killedarr) == 0:
            print i
        elif False not in killedarr:
            print i

def hasfilter(hpslist,filtername):
    """
    this will only return hpsobjects that have the provided filtername
    """            
    
    newhpslist=[]
    
    for i in hpslist:
        hpsfilters=[]
        for j in i.photometry:
            hpsfilters.append(j.filter.name)
        if filtername in hpsfilters:
            newhpslist.append(i)
    return newhpslist
        
def removefilter(hpslist,filtername):
    """
    this will remove the given filtername from all objects in hpslist
    """

    for obj in hpslist:
        photometry = obj.photometry
        newphot = []
        for i,fluxobj in enumerate(photometry):
            if fluxobj.filter.name != filtername:
                newphot.append(fluxobj)
            else:
                print filtername + " removed from " + str(obj.id)
        obj.photometry = newphot

    return hpslist
    
def removexray(hpslist):
    """
    this will take in a list of hps objects and will return the list with
    all the sources with x-ray counterparts removed
    """
    output = []
    
    for i in hpslist:
        if i.xraycounter == '':
            output.append(i)
        else:
            print "removed " + str(i.id) + ' ' + i.xraycounter
    
    return output
    
    
def lae_start(root='.',cutoff=0.05,silent=True):
    """
    this is another helpful function in interactive sessions
    if will return a list of only the objects with SED fits
    out of the whole HPS sample.
    root is a directory string which is where it looks for SED fits
    cutoff if the maximal R-1 value to be included
    """
    
    hpslist = makehpslist(silent)
    hpslist = getSEDfits(hpslist,root,silent)
    hpslist = selectseds(hpslist)
    hpslist = selectconverge(hpslist,cutoff)
    hpslist = removexray(hpslist)
    
    for i in hpslist:
        i.extinction_correct()
        i.lyasfr()
        i.uvsfr()
    
    return hpslist
    
    
def hpssame(hpslist,hpslist2):
    """
    this will take in two hpslists and return 2 hpslists
    the returned lists will only have objects in both of the lists
    """
    
    same1 = []
    same2 = []
    
    for i in hpslist:
        hpsid = i.id
        for j in hpslist2:
            hpsid2 = j.id
            if hpsid == hpsid2:
                same1.append(i)
                same2.append(j)
                break
                
    return same1, same2
    
    
def binbyz(param,z,step = 0.2,start = 2.0, finish = 3.6):
    """
    this function will bin the data into redshift bins
    and return the median and standard error
    """
    
    zout = []
    paramout = []
    paramstdout = []
    
    while start <= finish:
        bs = z>start
        bf = z<(start + step)
        b = bs*bf
        zout.append(start + step / 2)
        paramb = param[b]
        paramout.append(np.median(paramb))
        paramstdout.append(paramb.std() / np.sqrt(len(paramb)))
        start += step
        
    zout = np.array(zout)
    paramout = np.array(paramout)
    paramstdout = np.array(paramstdout)
    
    return zout, paramout, paramstdout
    
