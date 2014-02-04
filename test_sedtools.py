from sedtools import *
from numpy.testing import assert_almost_equal

def test_flux_class():
    fluxobj = Flux(5.0,1.0,'test_data/SubB.res')
    assert type(fluxobj) == Flux
    assert hasattr(fluxobj, 'flux')
    assert hasattr(fluxobj, 'err')
    assert hasattr(fluxobj, 'filter')
    assert hasattr(fluxobj, 'units')
    
