from pytest import approx
from opendrift.physics.waveperiod import Tm0, Tp, Spectrum

def test_tp():
    t = Tp(12.)
    assert t.tp(Spectrum.JONSWAP) == approx(12.)

def test_tm0():
    ta = Tm0(10.)
    assert ta.tp(Spectrum.JONSWAP) == approx(11.)
