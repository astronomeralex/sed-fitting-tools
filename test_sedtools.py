from sedtools import *
from numpy.testing import assert_almost_equal
import pytest

def test_flux_class():
    fluxobj = Flux(5.0,1.0,Filter('test_data/SubB.res'))
    assert type(fluxobj) == Flux
    assert hasattr(fluxobj, 'flux')
    assert hasattr(fluxobj, 'err')
    assert hasattr(fluxobj, 'filter')
    assert hasattr(fluxobj, 'units')
    assert_almost_equal(fluxobj.flux,5.0)
    assert_almost_equal(fluxobj.err,1.0)

def test_flux_badinputs():
    f = Filter('test_data/SubB.res')
    with pytest.raises(ValueError):
        fluxobj = Flux(5.0, -1.0, f)
    with pytest.raises(ValueError):
        fluxobj = Flux(-5.0, 1.0, f)
    with pytest.raises(TypeError):
        fluxobj = Flux(5.0, 1.0, 'test_data/SubB.res')
    
def test_filter_class():
    filterobj = Filter("test_data/SubB.res")
    assert type(filterobj) == Filter
    assert hasattr(filterobj, 'transfile')
    assert hasattr(filterobj, 'central')
    assert hasattr(filterobj, 'transmission')
    assert hasattr(filterobj, 'waves')
    assert hasattr(filterobj, 'long10')
    assert hasattr(filterobj, 'short10')
    
def test_filter_badinputs():
    with pytest.raises(ValueError):
        f = Filter("test_data/SubB.negtrans.res")
    with pytest.raises(ValueError):
        f = Filter("test_data/SubB.negwave.res")
    with pytest.raises(ValueError):
        f = Filter("test_data/SubB.outoforder.res")    
    
def test_galaxy_class():
    f1 = Flux(5.0,1.0,Filter('test_data/SubB.res'))
    f2 = Flux(10.0,1.0,Filter('test_data/UKIRTJ.res'))
    f3 = Flux(20.0, 5.0,Filter(""))
    fluxlist = [f1,f2,f3]
    name = "test_gal"
    testgal = Galaxy(name,fluxlist,1.0)
    assert type(testgal) == Galaxy
    assert hasattr(testgal, "name")
    assert hasattr(testgal, "fluxlist")
    assert hasattr(testgal, "z")
    assert hasattr(testgal, "sedfluxlist")
    assert len(testgal.fluxlist) == 3
    assert len(testgal.sedfluxlist) == 2
    #test sorting of sedfluxlist
    for i in range(len(testgal.sedfluxlist) - 1):
        assert testgal.sedfluxlist[i].filter.central < testgal.sedfluxlist[i + 1].filter.central
        
def test_galaxy_badinputs():
    f1 = Flux(5.0,1.0,Filter('test_data/SubB.res'))
    f2 = Flux(10.0,1.0,Filter('test_data/UKIRTJ.res'))
    f3 = Flux(20.0, 5.0,Filter(""))
    fluxlist = [f1,f2,f3]
    name = "test_gal"
    #testing fluxlist
    with pytest.raises(TypeError):
        testgal = Galaxy(name, [5,6], 1.0)
    with pytest.raises(ValueError):
        testgal = Galaxy(name, [], 1.0)
    #testing redshift
    with pytest.raises(ValueError):
        testgal = Galaxy(name, fluxlist, -1.0)
    with pytest.raises(TypeError):
        testgal = Galaxy(name, fluxlist, None)

def test_galmc_class():
    paramfile = 'test_data/o2default.ini'
    depfile = "test_data/empty.depfile"
    photlimits = (0,30000)
    testgalmc = GalMC(paramfile, depfile, photlimits)
    assert type(testgalmc) == GalMC
    assert hasattr(testgalmc, "paramfilename")
    assert hasattr(testgalmc, "paramfile")
    assert hasattr(testgalmc, "depfile")
    assert hasattr(testgalmc, "photlimits")
    
def test_galmc_badinputs():
    paramfile = 'test_data/o2default.ini'
    depfile = "test_data/empty.depfile"
    photlimits = (0,30000)
    with pytest.raises(IOError):
        testgalmc = GalMC("", depfile, photlimits)
    with pytest.raises(IOError):
        testgalmc = GalMC(paramfile, "", photlimits)
    with pytest.raises(TypeError):
        testgalmc = GalMC(paramfile, depfile, (0))
    with pytest.raises(AssertionError):
        testgalmc = GalMC(paramfile, depfile, (-1, 2))
    with pytest.raises(AssertionError):
        testgalmc = GalMC(paramfile, depfile, (1,2,3,4))
        
def test_pbsbackend_class():
    preamble =  ['#PBS -l nodes=1','#PBS -l walltime=96:00:00', 
        '#PBS -q lionxf-yuexing','#PBS -j oe', 'cd $PBS_O_WORKDIR',
        'echo " "' , 'echo "Job started on `hostname` at `date`"']
    postamble = ['echo "Job Ended at `date`"', 'echo " "']
    testpbs = PBSBackend(preamble,postamble)
    assert type(testpbs) == PBSBackend
    assert hasattr(testpbs,"preamble")
    assert hasattr(testpbs,"postamble")  
    
def test_pbsbackend_badinputs():
    preamble =  ['#PBS -l nodes=1','#PBS -l walltime=96:00:00', 
        '#PBS -q lionxf-yuexing','#PBS -j oe', 'cd $PBS_O_WORKDIR',
        'echo " "' , 'echo "Job started on `hostname` at `date`"']
    postamble = ['echo "Job Ended at `date`"', 'echo " "']
    with pytest.raises(TypeError):
        testpbs = PBSBackend(5,postamble)
    with pytest.raises(TypeError):
        testpbs = PBSBackend(preamble, "moose")
    with pytest.raises(TypeError):
        testpbs = PBSBackend([5],postamble)
    with pytest.raises(TypeError):
        testpbs = PBSBackend(preamble,[4.567])
        
def test_localbackend_class():
    testlocal = LocalBackend()
    assert type(testlocal) == LocalBackend
    
