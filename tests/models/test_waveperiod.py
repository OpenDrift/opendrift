from pytest import approx
from opendrift.physics.waveperiod import Tm01, Tp, Spectrum

def test_tp():
    t = Tp(12.)
    assert t.tp(Spectrum.JONSWAP) == approx(12.)

def test_tm0():
    ta = Tm01(10.)
    assert ta.tp(Spectrum.JONSWAP) == approx(11.)

    b = Tm01(10.).into(Tm01, Spectrum.Phillips).into(Tm02, Spectrum.JONSWAP)


    print(b.t)




    ###


    c = b.into(Tp, Spectrum.Phillips)
