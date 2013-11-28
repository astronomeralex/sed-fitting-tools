# This provides tools for SED fitting many galaxies
# it is designed to be modularized and be useful for many different code
# written by Alex Hagen, a graduate student in astrophysics at penn state
# mr.alex.hagen@gmail.com

from abc import ABCMeta, abstractmethod

################################################################################

class Flux(object):
    """
    this flux class will contain a single photometric measurement
    """
    def __init__(self,flux,err,bandpassfile,silent=False):
        self.flux = flux
        self.err = err
        self.filter = Filter(bandpassfile,silent)
    
################################################################################

class Filter(object):
    """
    the class contains information on the filter, including the name,
    location of the file describing the transmission, and the central wavelength
    """
    def __init__(self,bandpassfile,silent = False):
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
            if not silent:
                print "Can't find transmission file at " + self.transfile

################################################################################

class Galaxy(object):
    """
    
    """
    def __init__(self, name, fluxlist, redshift):
        """
        
        """
        self.name = name
        self.z = redshift
        self.photometry = fluxlist

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



