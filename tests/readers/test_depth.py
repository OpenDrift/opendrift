"""test depth.py"""

from opendrift.readers.roppy import depth
import numpy as np

import unittest

class TestMyFunction(unittest.TestCase):
    def setUp(self):
        # don't need to test time because always first select time index or
        # is interpolated to time
        self.NZ, self.NY, self.NX = 3, 4, 5
        self.Htot = np.ones((self.NY, self.NX))
        self.hc = 0.1
        self.Cs_r = np.linspace(-1, 0, self.NZ)

    def test_sdepth_zeta_zero_Vtransform1(self):
        """test zero zeta in sdepth, Vtransform 1."""

        Vtransform = 1

        # zeta 0 scalar
        zeta = 0
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        
        expected_output = np.array([-0.98333333, -0.5, -0.01666667])        
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)

        # zeta 0s array
        zeta = np.zeros((self.NY, self.NX))
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        
        expected_output = np.array([-0.98333333, -0.5, -0.01666667])
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)


    def test_sdepth_zeta_zero_Vtransform2(self):
        """test zero zeta in sdepth, Vtransform 2."""

        Vtransform = 2

        # zeta 0 scalar
        zeta = 0
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        expected_output = np.array([-0.98484848, -0.5, -0.01515152])
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)

        # zeta 0s array
        zeta = np.zeros((self.NY, self.NX))
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        expected_output = np.array([-0.98484848, -0.5, -0.01515152])
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)

    def test_sdepth_zeta_nonzero_Vtransform1(self):
        """test nonzero zeta in sdepth, Vtransform 1."""

        Vtransform = 1

        # zeta scalar
        zeta = 1
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        expected_output = np.array([-0.96666667,  0.,  0.96666667])        
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)

        # zeta 1s array
        zeta = np.ones((self.NY, self.NX))
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        expected_output = np.array([-0.96666667,  0.,  0.96666667])
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)


    def test_sdepth_zeta_nonzero_Vtransform2(self):
        """test nonzero zeta in sdepth, Vtransform 2."""

        Vtransform = 2

        # zeta 1 scalar
        zeta = 1
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        expected_output = np.array([-0.96969697,  0.,  0.96969697])
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)

        # zeta 1s array
        zeta = np.ones((self.NY, self.NX))
        z_rho = depth.sdepth(self.Htot, zeta, self.hc, self.Cs_r, Vtransform=Vtransform)
        expected_output = np.array([-0.96969697,  0.,  0.96969697])
        assert z_rho.shape == (self.NZ, self.NY, self.NX)
        assert np.allclose(z_rho[:,0,0], expected_output)



if __name__ == '__main__':
    unittest.main()
