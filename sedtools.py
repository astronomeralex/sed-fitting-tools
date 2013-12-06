# This provides tools for SED fitting many galaxies
# it is designed to be modularized and be useful for many different code
# written by Alex Hagen, a graduate student in astrophysics at penn state
# mr.alex.hagen@gmail.com

from abc import ABCMeta, abstractmethod
import numpy as np
import copy
import logging

logger = logging.getLogger('sedtools')
logger.setLevel(logging.WARNING)

################################################################################

class Flux(object):
    """
    this flux class will contain a single photometric measurement
    """
    def __init__(self,flux,err,bandpassfile):
        self.flux = flux
        self.err = err
        self.filter = Filter(bandpassfile)
    
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
            
        else:
            logger.error("Can't find transmission file at " + self.transfile)
                
    def 

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
        elif False in [isinstance(Flux,i) for i in fluxlist]:
            raise TypeError("Input fluxlist should only contain instance of Flux")
            
        if type(name) != str:
            try:
                nametype = type(str)
                name = str(name)
                logger.info("Galaxy name cast from " + str(nametype) + "to string: name is " + name)
            except ValueError:
                raise ValueError("Galaxy name input of type " + str(nametype) + " cannot be cast to string")
            
        if type(redshift) != float:
            try: 
                inputztype = type(redshift)
                redshift = float(redshift)
                logger.info("Redshift cast from input " + str(inputztype) + "to float: z = " + str(redshift))
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
                logger.warning(i.filter.transfile + "doesn't exist -- this filter " + \
                    "won't be included in the SED fit")
            else:
                self.sedfluxlist.append(i)
                
        if len(sedfluxlist) == 0:
            logger.critical("Fluxlist for SED fitting is length zero")
        else:
            #sort fluxlist by central filter wavelength
            self.sedfluxlist.sort(key = lambda x: x.filter.central)
        
        self.numphot = len(self.fluxlist)
        self.numsedphot = len(self.sedfluxlist)
            
    def __str__(self):
        return "Galaxy: "+ str(self.name) + " z=" + "%.3f" %self.z
        


################################################################################

class SEDfitter(object):
    __metaclass__ = ABCMeta
    
    #@abstractmethod
    #put methods SED fitters must have here
    
    
class GalMC(SEDfitter):
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
    pass



