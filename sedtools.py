# This provides tools for SED fitting many galaxies
# it is designed to be modularized and be useful for many different code
# written by Alex Hagen, a graduate student in astrophysics at penn state
# mr.alex.hagen@gmail.com

from abc import ABCMeta, abstractmethod
import numpy as np
import copy
import logging
import os
import pp 
import random
import subprocess

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

################################################################################

class Flux(object):
    """
    this flux class will contain a single photometric measurement
    """
    def __init__(self, flux, err, bandpassfile, units=None):
        self.flux = flux
        if err < 0:
            raise ValueError("Input error is less than zero")
        self.err = err
        self.filter = Filter(bandpassfile)
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
            #TODO: add sanity checks here for waves and transmission (>0, etc)
            
            
            #now find longest wavelength at 10% transmission
            #this is useful for determining if NIR filters could be contaminated by
            #3.3um PAH feature
            self.long10 = waves[np.where(trans > 0.1)[0][-1]]
            self.short10 = waves[np.where(trans > 0.1)[0][0]]
            
        else:
            logging.error("Can't find transmission file at " + self.transfile)

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

class SEDfitter(object):
    __metaclass__ = ABCMeta
#    
#    @abstractmethod
#    def prepare(self):
#        """this method will make any necessary folders and files
#        and the job will be ready to submit"""
#    
    @abstractmethod
    def submit(self):
        """this method will submit the galaxy for 
        SED fitting using the given backend"""
        
    def cleanphotometry(self, galaxy):
        """
        
        """
        sedphot = []
        redlim = self.photlims[0]
        bluelim = self.photlims[1]
        for i in galaxy.sedfluxlist:
            if (i.filter.short10 > bluelim) and (i.filter.long10 < redlim):
                sedphot.append(i)
            else:
                if i.filter.short10 < bluelim:
                    logging.warning("Flux measurement with filter " + i.filter.transfile + "not included -- filter is too blue")
                else:
                    logging.warning("Flux measurement with filter " + i.filter.transfile + "not included -- filter is too red")
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
        #see if input paramfile exists TODO: make sure paramfile is sane?
        try:
            self.paramfilename = paramfilename
            paramfileobj = open(self.paramfilename)
            self.paramfile = paramfileobj.readlines()
            paramfileobj.close()
        except IOError:
            raise IOError(paramfilename + " not found")
            
        if os.path.exists(depfile):
            self.depfile = depfile
        else:
            raise IOError(depfile + "doesn't exist")
        try:
            assert len(photlimits) == 2
        except TypeError:
            raise TypeError("Photlimits is not a sequence")
        assert photlimits[0] < photlimits[1]
        assert (photlimits[0] >= 0) and (photlimits[1] >=0)
        
        self.photlimits = photlimits
        
    
    def writedata(self, galaxy):
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
        
        #and clean the galaxy photometry
        self.cleanphotometry(galaxy)
        
        #it assumes that it is SED fitting in a new directory and will overwrite previous sed versions
        #make the new directory
        rootdir = os.getcwd()
        galaxydir = str(galaxy.name)
        os.mkdir(galaxydir)
        os.chdir(galaxydir)
        
        self.writedata(galaxy)
        
        #now its time to write out the param files for each of the chains we want
        for i in xrange(numchains):
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
        
        commandlist = ['./MCMCfit ' + str(galaxy.name) + '_' + str(i) + '.ini' for i in xrange(numchains)]
        
        backend.run(depcommands, commandlist)
        
        os.chdir(rootdir)
        
        

################################################################################

class ComputingBackend(object):
    __metaclass__ = ABCMeta
    
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
        self.preamble = preamble
        self.postamble = postamble
        
    def run(self, depdirections, commandlist, scriptname = "sed.pbs"):
        """
        depdirections and commandlist are lists. each entry is a single command
        """
        
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
    def __init__(self,**kwargs):
        """
        used pp for local stuffs
        """
        self.ppserver = pp.Server(kwargs)
        
        



