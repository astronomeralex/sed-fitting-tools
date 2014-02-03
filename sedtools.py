# This provides tools for SED fitting many galaxies
# it is designed to be modularized and be useful for many different code
# written by Alex Hagen, a graduate student in astrophysics at penn state
# mr.alex.hagen@gmail.com

from abc import ABCMeta, abstractmethod
import numpy as np
import copy
import logging
import pp
import os

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

################################################################################

class Flux(object):
    """
    this flux class will contain a single photometric measurement
    """
    def __init__(self, flux, err, bandpassfile, units=None):
        self.flux = flux
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
    
    @abstractmethod
    def prepare(self):
        """this method will make any necessary folders and files
        and the job will be ready to submit"""
    
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
            if i.filter.short10 > bluelim and i.filter.long10 < redlim:
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
        depfile is the file or folder of dependencies needed for the sedfit
        """
        #see if input paramfile exists TODO: make sure paramfile is sane?
        try:
            self.paramfilename = paramfilename
            paramfileobj = open(self.paramfilename)
            self.paramfile = paramfileobj.readlines()
            paramfileobj.close()
        except IOError:
            raise IOError(paramfilename + " not found")
            
        #see if dependent file path exists and is a directory or a file
        if os.path.exists(depfile):
            if os.path.isfile(depfile):
                self.depfile = depfile
                self.deptype = 'file'
            elif os.path.isdir(depfile):
                self.depfile = depfile
                self.deptype = "dir"
            else:
                raise IOError(depfile + " isn't a file or directory...something \
                has gone terribly wrong")
        else:
            raise IOError(depfile + "doesn't exist")
        
        assert len(photlimits) == 2
        assert photlimits[0] < photlimits[1]
        self.photlimits = photlimits
        
    def submit(self, backend, galaxy):
        """
        
        """
        pass
        

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
    
    def run(self):
        pass
    
class LocalBackend(ComputingBackend):
    """
    
    """
    def __init__(self,**kwargs):
        """
        used pp for local stuffs
        """
        self.ppserver = pp.Server(kwargs)
        
        



