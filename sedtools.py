# This provides tools for SED fitting many galaxies
# it is designed to be modularized and be useful for many different code
# written by Alex Hagen, a graduate student in astrophysics at penn state
# mr.alex.hagen@gmail.com
#
#Licensed under The Academic License (https://github.com/dfm/license)
#
#Copyright (c) 2014 Alex Hagen
#
#This project includes academic-research code and documents under development.
#You would be a fool to run any of the code.
#Any use of the content requires citation.
#
from abc import ABCMeta, abstractmethod
import numpy as np
import logging
import os
import random
import subprocess
import glob
import shutil

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

################################################################################

class Flux(object):
    """
    this flux class will contain a single photometric measurement
    """
    def __init__(self, flux, err, filterobj, units=None):
        if flux < 0:
            raise ValueError("Input flux is less than zero")
        elif flux == 0:
            logging.warning("The flux is set at zero; are you sure about that?")
        self.flux = flux
        
        if err < 0:
            raise ValueError("Input error is less than zero")
        elif err == 0:
            logging.warning("The error is set at zero; that could cause problems...")
        self.err = err
        
        if type(filterobj) == Filter:
            self.filter = filterobj
        else:
            raise TypeError(str(filterobj) + "isn't a Filter object")
        
        if units is None:
            logging.info("Units for flux object not defined")
        self.units = units
    
################################################################################

class Filter(object):
    """
    the class contains information on the filter, including the name,
    location of the file describing the transmission, and the central wavelength
    """
    def __init__(self,bandpassfile):
        self.transfile = bandpassfile #the location of the transmission file
        #now check and see if the transmission file exists
        if os.path.exists(self.transfile):
            transmission = np.loadtxt(self.transfile)
            assert transmission.shape[-1] == 2
            waves = transmission[:,0]
            trans = transmission[:,1]
            #sanity checks for waves and transmission (>0, sorted)
            if not np.all(waves > 0):
                raise ValueError("One of more wavelengths is negative")
            if not np.all(trans >= 0):
                raise ValueError("One or more transmission values is negative")
            if not np.all( waves == np.sort(waves) ):
                raise ValueError("Transmission file wavelengths are not sorted")
            self._central = np.trapz(trans * waves, x=waves) / np.trapz(trans, x=waves)
            #normalize transmission
            trans = trans / trans.max()
            self._transmission = trans
            self._waves = waves            
            #now find longest wavelength at 10% transmission
            #this is useful for determining if NIR filters could be contaminated by
            #3.3um PAH feature
            self._long10 = waves[np.where(trans > 0.1)[0][-1]]
            self._short10 = waves[np.where(trans > 0.1)[0][0]]
            #find minimum spacing of the waves array
            self.minspacing = np.min( np.diff(waves))
            
        else:
            logging.error("Can't find transmission file at " + self.transfile)
            
    @property
    def central(self):
        """Return the central wavelength"""
        return self._central
        
    @property
    def transmission(self):
        """Return the filter transmission array"""
        return self._transmission
        
    @property
    def waves(self):
        """Return the filter wavelength array"""
        return self._waves
        
    @property
    def long10(self):
        """Returns the wavelength of the red side of 10% transmission"""
        return self._long10
        
    @property
    def short10(self):
        """Returns the wavelength of the blue side of 10% transmission"""
        return self._short10

################################################################################

class Galaxy(object):
    """
    
    """
    def __init__(self, name, fluxlist, redshift):
        """
        
        """
        #sanity checking the inputs
        if len(fluxlist) == 0:
            raise ValueError("Input fluxlist is length zero")
        elif False in [isinstance(i,Flux) for i in fluxlist]:
            raise TypeError("Input fluxlist should only contain instances of Flux")
            
        if type(name) != str:
            try:
                nametype = type(str)
                name = str(name)
                logging.info("Galaxy name cast from " + str(nametype) + "to string: name is " + name)
            except ValueError:
                raise ValueError("Galaxy name input of type " + str(nametype) + " cannot be cast to string")
        #TODO: make sure name is also a valid filename
            
        if type(redshift) != float:
            try: 
                inputztype = type(redshift)
                redshift = float(redshift)
                logging.info("Redshift cast from input " + str(inputztype) + "to float: z = " + str(redshift))
            except ValueError:
                raise ValueError("Galaxy redshift input of type " + str(inputztype) + " cannot be cast to float!")
        if redshift < 0:
            raise ValueError("Redshift is less than zero: z = " + str(redshift))
        
        self.name = name
        self.z = redshift
        self.fluxlist = fluxlist
        self.sedfluxlist = []
        
        #make sure all filter objects have central wavelengths, etc
        for i in fluxlist:
            if not hasattr(i.filter,"central"):
                logging.warning(i.filter.transfile + "doesn't exist -- this filter " + \
                    "won't be included in the SED fit")
            else:
                self.sedfluxlist.append(i)
                
        if len(self.sedfluxlist) == 0:
            logging.critical("Fluxlist for SED fitting is length zero")
        else:
            #sort fluxlist by central filter wavelength
            self.sedfluxlist.sort(key = lambda x: x.filter.central)
        
        self.numphot = len(self.fluxlist)
        self.numsedphot = len(self.sedfluxlist)
            
    def __str__(self):
        return "Galaxy: "+ str(self.name) + " z=" + "%.3f" %self.z
        
    def homogenize_photometry(self,fluxunit):
        """
        
        """
        pass
        #use astropy.units?

################################################################################

class SEDfitter(object, metaclass=ABCMeta):
    @abstractmethod
    def submit(self):
        """this method will submit the galaxy for 
        SED fitting using the given backend"""
        
    def cleanphotometry(self, galaxy):
        """
        
        """
        sedphot = []
        bluelim = self.photlimits[0]
        redlim = self.photlimits[1]
        for i in galaxy.sedfluxlist:
            if (i.filter.short10 > bluelim*(1 + galaxy.z)) and (i.filter.long10 < redlim*(1 + galaxy.z)):
                sedphot.append(i)
            else:
                if i.filter.short10 < bluelim*(1 + galaxy.z):
                    logging.warning("Flux measurement with filter " + i.filter.transfile + " not included -- filter is too blue")
                else:
                    logging.warning("Flux measurement with filter " + i.filter.transfile + " not included -- filter is too red")
        galaxy.cleansedphot = sedphot
        
class GalMC(SEDfitter):
    """
    
    """
    
    def __init__(self, paramfilename, depfile, photlimits = (0,30000)):
        """
        paramfilename is the parameter file
        depfile is the file dependencies needed for the sedfit
        currently assume depfile is a zipped file
        """
        #TODO add path variable to inputs to make the directory structure more
        #flexible
        
        #see if input paramfile exists TODO: make sure paramfile is sane?
        try:
            self.paramfilename = paramfilename
            paramfileobj = open(self.paramfilename)
            self.paramfile = paramfileobj.read()
            paramfileobj.close()
        except IOError:
            raise IOError(paramfilename + " not found")
            
        if os.path.exists(depfile):
            self.depfile = os.path.abspath(depfile)
        else:
            raise IOError(depfile + "doesn't exist")
        try:
            assert len(photlimits) == 2
        except TypeError:
            raise TypeError("Photlimits is not a sequence")
        assert photlimits[0] < photlimits[1]
        assert (photlimits[0] >= 0) and (photlimits[1] >=0)
        
        self.photlimits = photlimits
        
    @staticmethod
    def writedata(galaxy):
        """
        this method will write out the data file in standard GalMC format 
        wavelength flux density error (in microJy)
        """ 
        datafile = open(str(galaxy.name)+'.dat','w')
        for i in galaxy.cleansedphot:
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
        

    def writeparams(self, galaxy, chainnum):
        """
        this writes out the params file
        chainnum is an int that is the number of the chain
        """
        output= []
        #each entry of output is a line in the params file
        output.append('chain root = ' + str(galaxy.name) + '_' + str(chainnum))
        output.append('Data File = ' + str(galaxy.name) + '.dat')
        #fitting params
        for i in self.params:
            outline = i + ' = ' + str(self.params[i][0]) + ', ' + str(self.params[i][1]) + ', ' + str(self.params[i][2]) + ', ' + str(self.params[i][3])
            output.append(outline)
        #now for the filters
        for i,fluxobj in enumerate(galaxy.cleansedphot):
            output.append('filter('+str(i+1)+') = ../' + fluxobj.filter.transfile)
        
        outputfile = open(str(galaxy.name) + '_' + str(chainnum) + '.ini','w')
        outputfile.write(self.paramfile)
        for i in output:
            outputfile.write(i + '\n')
        outputfile.close()
        
    def submit(self, backend, galaxy, numchains, sedparaminfo):
        """
        example sedparaminfo
        maxage = cosmology_distance.tz(obj.z,0.7,0.7,0.3,0.0) #this is in years
        
        logmaxage = np.log(maxage)
        sedparaminfo = {}
        sedparaminfo['Mas']=[16.,9.,35.,0.5]
        sedparaminfo['Age']=[16.,12.,logmaxage,0.5]
        sedparaminfo['Tau']=[-0.100,-5.,15.,0.0]
        sedparaminfo['EBV']=[0.1,0.0,1.0,0.04]
        sedparaminfo['Met']=[-0.7,-2.29,0.45,0.0]
        sedparaminfo['Red']=[obj.z,0.0,1.5,0.0]
        """
        
        self.params = sedparaminfo
        
        #first some sanity checks on the inputs
        assert type(numchains) == int
        assert numchains > 0
        
        #and clean the galaxy photometry & check the filter spacing
        self.cleanphotometry(galaxy)
        for i in galaxy.cleansedphot:
            if i.filter.minspacing < 1.0:
                raise ValueError("Filter spacing for " + str(i.filter.transfile) + " needs to be greater than 1 angstrom")
        
        #it assumes that it is SED fitting in a new directory and will overwrite previous sed versions
        #make the new directory
        rootdir = os.getcwd()
        galaxydir = str(galaxy.name)
        os.mkdir(galaxydir)
        os.chdir(galaxydir)
        
        self.writedata(galaxy)
        
        #now its time to write out the param files for each of the chains we want
        for i in range(numchains):
            if i>0:
                self.randparams()
                self.writeparams(galaxy, i)
            else:
                self.writeparams(galaxy, i)
                
        depcommands = []
        #this split is to get the filename, not the full path
        localname = self.depfile.split('/')[-1]
        depcommands.append('cp ' + self.depfile + " .")
        depcommands.append('unzip ' + localname )
        depcommands.append('rm -rf ' + localname)
        
        commandlist = ['./MCMCfit ' + str(galaxy.name) + '_' + str(i) + '.ini' for i in range(numchains)]
        
        backend.run(depcommands, commandlist)
        
        os.chdir(rootdir)
        

################################################################################

class ComputingBackend(object, metaclass=ABCMeta):
    @abstractmethod
    def run(self):
        """
        
        """
        
class PBSBackend(ComputingBackend):
    """
    
    """
    def __init__(self, preamble, postamble):
        """
        preamble is list -- ['#PBS -l nodes=1','#PBS -l walltime=96:00:00', 
        '#PBS -q lionxf-yuexing','#PBS -j oe', 'cd $PBS_O_WORKDIR',
        'echo " "' , 'echo "Job started on `hostname` at `date`"']
        postamble is list -- ['echo "Job Ended at `date`"', 'echo " "']
        """
        if type(preamble) != list: raise TypeError("Preamble should be a list")
        if type(postamble) != list: raise TypeError("Postamble should be a list")
        for i in preamble:
            if type(i) != str:
                raise TypeError("preamble entries should be strings")
        for i in postamble:
            if type(i) != str:
                raise TypeError("postamble entries should be strings")
        self.preamble = preamble
        self.postamble = postamble
        
    def run(self, depdirections, commandlist, scriptname = "sed.pbs"):
        """
        depdirections and commandlist are lists. each entry is a single command
        """
        #TODO: add way to add or edit preamble or postamble
        
        pbsout = self.preamble + depdirections + commandlist + self.postamble
                
        pbsfile = open(scriptname,'w')
        for i in pbsout:
            pbsfile.write(i + '\n')
        pbsfile.close()
        
        proc = ['qsub','sed.pbs']
        
        subprocess.call(proc)
        
    
class LocalBackend(ComputingBackend):
    """
    
    """
        
    def run(self):
        pass

################################################################################

class SEDAnalysisTool(object, metaclass = ABCMeta):
    """
    
    """
    
    @abstractmethod
    def run(self):
        """
        """


class GetDist(SEDAnalysisTool):
    """
    
    """
    
    def __init__(self, distparamsfile="distparams.ini"):
        """
        
        """
        
        try:
            distparams = open(distparamsfile).read()
        except IOError:
            raise IOError("Distparamsfile not found at " + distparamsfile)
        
        self.distparamsfile = distparamsfile
        self.distparams = distparams
        
    def run(self, folders, numchains=4, makeplots=False):
        """
        this will run GetDist on the chains outputed from the fit method of the GalMC class. it also will rename the chain
        files so they work with getdist
        
        folders is a list where each element is a string that is the name of a folder where getdist should be run
        it is recommended that glob is used to make this list
        
        distparamsfile is the file name of the parameter file for getdist
        this file cannot have the chain root included in it
        the programs will write that so it works for each folder
        """
        
        root = os.getcwd()
        for folder in folders:
            os.chdir(folder)
            #write out the distparams file
            distparout = open('distparams.ini','w')
            distparout.write('file_root = ' + folder + '\n')
            distparout.write(self.distparams)
            distparout.close()
            
            #now its time to rename the chains
            for i in range(numchains):
                proc=['mv',folder + '_' + str(i) + '_chain.txt',folder + '_' + str(i) + '.txt']
                subprocess.call(proc)    
            #now its time to run getdist!
            #make sure its in the .bashrc file
            subprocess.call(['getdist','distparams.ini'])
    
            #and now time to run the python files produced by getdist so we can have plotses
            if makeplots:
                pyfiles = glob.glob(folder + '*.py')
                for i in pyfiles:
                    exec(compile(open(i).read(), i, 'exec'))
            else:
                logging.info("Plots not made")
            
            os.chdir(root)
    
    @staticmethod
    def get_sed_fit(galaxy,root='.'):
        """
        this method will fetch the results of the sed fit. root is the directory
        location that contains the folders that each contain the results of the
        sed fit run by sedfit. the results have to be analyzed by getdist prior
        to the running of this method
        this is designed so if the SED results don't exist, it will exit gracefully
        with a message and will not throw an error
        this is so it can be run on a large sample easily without try statements
        """
        
        folder = root + '/' + str(galaxy.name)
        
        #first, check to make sure the folder exists and that files from getdist exist
        if os.path.exists(folder) and os.path.exists(folder + '/' + galaxy.name + '.margestats'):
            #read in coverge file -- currently will take R-1 values
            convfile = open(folder + '/' + galaxy.name + '.converge')
            convlines = convfile.readlines()
            convfile.close()
            
            #TODO -- see if there's a better way to do this
            #now i need to find the line before R-1 value
            for i,line in enumerate(convlines):
                if 'var(mean)/mean(var) for eigenvalues of covariance of means' in line:
                    linenum = i
                    break
            if 'linenum' in locals():        
                #time to get the R values
                Rm1s = []
                for i in range(1,50):
                    line = convlines[linenum + i].split()
                    if len(line) < 2:
                        break
                    else:
                        Rm1s.append(float(line[1]))
                galaxy.Rm1s = Rm1s        
                #find the largest R value
                galaxy.Rm1 = max(galaxy.Rm1s)
            else:
                galaxy.Rm1 = None
                galaxy.Rm1s = None
                
            #now time to read in the marginalized stats
            margefile = open(folder + '/' + galaxy.name + '.margestats')
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
            
            galaxy.sedstats = sedstats
            
            #now get the best fit parameters from the likestats file
            likefile = open(folder + '/' + galaxy.name + '.likestats')
            likelines = likefile.readlines()
            likefile.close()
            
            likestats = {}
            for i,line in enumerate(likelines):
                stats = line.split()
                if i == 0:
                    galaxy.mloglike = float(stats[-1])
                elif i>2 and len(stats) > 1:
                    paramname = stats[-1]
                    likestats[paramname] = {}
                    likestats[paramname]['best'] = float(stats[1])
                    likestats[paramname]['lower1'] = float(stats[2])
                    likestats[paramname]['upper1'] = float(stats[3])
                    likestats[paramname]['lower2'] = float(stats[4])
                    likestats[paramname]['upper2'] = float(stats[5])
                    
            galaxy.likestats = likestats        
                    
            galaxy.chisq = galaxy.mloglike*2. #note that this is only valid for a flat prior
            numphot = len(galaxy.cleansedphot)
            numsedparams = len(galaxy.sedstats) - 1 #since mass is included twice
            if numphot - 1 > numsedparams:
                galaxy.redchisq = galaxy.chisq / float(numphot - numsedparams - 1)
            else:
                logging.warning("Cannot calculate reduced chi squared")
            
            logging.info(galaxy.name + ": SED fit found")
                
        else:
            logging.warning(galaxy.name + ": SED info or folder doesn't exist")

    @staticmethod
    def get_sed_param(galaxylist,paramname,sigma = 1,tobase10 = False):
        """
        this will take in a list of HpsObj objects, a paramname, which is a string, and a sigma
        sigma can either be 1 or 2 and will specify the type of error bars returned, 1 or 2 sigma
        this will return three numpy arrays, the parameter from the SED fit for each object, 
        lower errors and upper errors
        """
    
        param = []
        lower = []
        upper = []
        
        for i in galaxylist:
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

    @staticmethod
    def add_sfr_to_chain(dirs,numchains=4):
        """
        this will add the SED SFR to all the chains in dirs
        this assumes that runGetDist has been run on all the chains
        """
    
    
        root = os.getcwd()
    
        for i in dirs:
            os.chdir(i)
            chainroot = i
            for j in range(numchains):
                chainfile = chainroot + '_' + str(j) + '.txt'
                #make copy of chain file to save
                chain = np.loadtxt(chainfile)
                shutil.copy2(chainfile,chainfile+'.old')
                #age is 2nd column, galaxymass is 3rd (0-indexed)
                sedsfr = np.exp(chain[:,3]) / np.exp(chain[:,2])
                sedsfr = sedsfr.reshape(len(sedsfr),1) #make column vector for hstack
                chain = np.hstack((chain,sedsfr))
                np.savetxt(chainfile,chain,fmt="%0.6e")
            os.chdir(root)
        
        
