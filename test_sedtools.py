from sedtools import *
from numpy.testing import assert_almost_equal
import pytest

def test_flux_class():
    fluxobj = Flux(5.0,1.0,'test_data/SubB.res')
    assert type(fluxobj) == Flux
    assert hasattr(fluxobj, 'flux')
    assert hasattr(fluxobj, 'err')
    assert hasattr(fluxobj, 'filter')
    assert hasattr(fluxobj, 'units')
    assert_almost_equal(fluxobj.flux,5.0)
    assert_almost_equal(fluxobj.err,1.0)

def test_flux_badinputs():
    with pytest.raises(ValueError):
        fluxobj = Flux(5.0, -1.0, 'test_data/SubB.res')
    
def test_filter_class():
    filterobj = Filter("test_data/SubB.res")
    assert type(filterobj) == Filter
    assert hasattr(filterobj, 'transfile')
    assert hasattr(filterobj, 'central')
    assert hasattr(filterobj, 'transmission')
    assert hasattr(filterobj, 'waves')
    assert hasattr(filterobj, 'long10')
    assert hasattr(filterobj, 'short10')
    
def test_galaxy_class():
    f1 = Flux(5.0,1.0,'test_data/SubB.res')
    f2 = Flux(10.0,1.0,'test_data/UKIRTJ.res')
    f3 = Flux(20.0, 5.0,"")
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
    for i in xrange(len(testgal.sedfluxlist) - 1):
        assert testgal.sedfluxlist[i].filter.central < testgal.sedfluxlist[i + 1].filter.central
        
def test_galaxy_badinputs():
    f1 = Flux(5.0,1.0,'test_data/SubB.res')
    f2 = Flux(10.0,1.0,'test_data/UKIRTJ.res')
    f3 = Flux(20.0, 5.0,"")
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

