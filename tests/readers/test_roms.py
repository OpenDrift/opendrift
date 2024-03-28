import unittest
from opendrift.readers import reader_ROMS_native
import xarray as xr
import numpy as np
from datetime import datetime
import pytest


class TestROMSReader(unittest.TestCase):
    def setUp(self):
        
        h = np.ones((2, 3))*10
        
        u_eastward = np.array([[[[1.0, 1.1, 1.2], [1.3, 1.4, 1.5]], [[1.6, 1.7, 1.8], [1.9, 2.0, 2.1]], [[2.2, 2.3, 2.4], [2.5, 2.6, 2.7]]], [[[2.8, 2.9, 3.0], [3.1, 3.2, 3.3]], [[3.4, 3.5, 3.6], [3.7, 3.8, 3.9]], [[4.0, 4.1, 4.2], [4.3, 4.4, 4.5]]]])
        
        u = np.array([[[[0.1, 0.2], [0.3, 0.4]], [[0.5, 0.6], [0.7, 0.8]], [[0.9, 1.0], [1.1, 1.2]]], [[[1.3, 1.4], [1.5, 1.6]], [[1.7, 1.8], [1.9, 2.0]], [[2.1, 2.2], [2.3, 2.4]]]])
        v = np.linspace(0, 1, 2*3*1*3).reshape(2, 3, 1, 3)
        Cs_r = np.linspace(-1, 0, 3)
        mask_rho = np.ones((2, 3))
        # use odd mask values so I can check when used in tests
        wetdry_mask_rho = np.ones((2, 2, 3))*2
        mask_u = np.ones((2, 2))
        mask_v = np.ones((1, 3))
        angle = np.zeros((2, 3))
        
        self.ds = xr.Dataset(
            {
                "lon_rho": (["eta_rho", "xi_rho"], [[-152.0, -151.9, -151.8], [-152.0, -151.9, -151.8]]),
                "lat_rho": (["eta_rho", "xi_rho"], [[60.5, 60.5, 60.5], [60.4, 60.4, 60.4]]),
                "lon_u": (["eta_u", "xi_u"], [[-151.95, -151.85], [-151.95, -151.85]]),
                "lat_u": (["eta_u", "xi_u"], [[60.5, 60.5], [60.4, 60.4]]),
                "lon_v": (["eta_v", "xi_v"], [[-152.0, -151.9, -151.8]]),
                "lat_v": (["eta_v", "xi_v"], [[60.45, 60.45, 60.45]]),
                "mask_rho": (["eta_rho", "xi_rho"], mask_rho),
                "h": (["eta_rho", "xi_rho"], h),
                "wetdry_mask_rho": (["ocean_time", "eta_rho", "xi_rho"], wetdry_mask_rho),
                "angle": (["eta_rho", "xi_rho"], angle),
                "mask_u": (["eta_u", "xi_u"], mask_u),
                "mask_v": (["eta_v", "xi_v"], mask_v),
                "ocean_time": ("ocean_time", [0, 1], {"units": "seconds since 1970-01-01"}),
                # "zeta": (["ocean_time", "eta_rho", "xi_rho"], zeta),
                "s_rho": (["s_rho"], Cs_r),
                "Cs_r": (["s_rho"], Cs_r),
                "hc": 16,
                "u_eastward": (("ocean_time", "s_rho", "eta_rho", "xi_rho"), u_eastward),
                "v_northward": (("ocean_time", "s_rho", "eta_rho", "xi_rho"), u_eastward),
                "u": (("ocean_time", "s_rho", "eta_u", "xi_u"), u),
                "v": (("ocean_time", "s_rho", "eta_v", "xi_v"), v),
                }
        )
    
    
    def test_catch_multiple_variables(self):
        """Catch if multiple variables are present."""
        standard_mapping = {"u_eastward": 'x_sea_water_velocity'}
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_mapping)
        with pytest.raises(ValueError):
            reader.get_variables("x_sea_water_velocity", datetime(1970, 1, 1), 0, 0, 0)

    def test_get_variables_u_eastward_wetdry_mask_rho(self):
        """get a variable from the dataset and verify correct mask was used
        
        Since wetdry_mask_rho is available and matches u_eastward dims, it should be used.
        """
        standard_mapping = {"u_eastward": 'x_sea_water_velocity',
                            "v_northward": 'y_sea_water_velocity',
                            'u': 'blah',
                            'v': 'blah'}  # redirect mapping of u/v so these don't overlap
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_mapping)
        var, key = "u_eastward", "x_sea_water_velocity"
        output, masks = reader.get_variables(key, datetime(1970, 1, 1), 0, 0, 0, testing=True)

        # check that variable was handled correctly
        assert key in output
        assert reader.do_not_rotate == ["x_sea_water_velocity", "y_sea_water_velocity"]
        assert np.allclose(output[key].data, reader.Dataset[var][0,-1,0,0:2].values)
        # check that wetdry mask was used (can tell by mask value)
        assert np.allclose(masks["mask_rho"], reader.Dataset["wetdry_mask_rho"][0,0,0:2].values)

    def test_get_variables_u_eastward_mask_rho(self):
        """get a variable from the dataset and verify correct mask was used
        
        wetdry_mask_rho is not available so mask_rho should be used.
        """
        standard_mapping = {"u_eastward": 'x_sea_water_velocity',
                            "v_northward": 'y_sea_water_velocity',
                            'u': 'blah',
                            'v': 'blah'}  # redirect mapping of u/v so these don't overlap
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_mapping)
        # drop wetdry_mask_rho to test that mask_rho is used
        reader.Dataset = reader.Dataset.drop_vars("wetdry_mask_rho")
        var, key = "u_eastward", "x_sea_water_velocity"
        output, masks = reader.get_variables(key, datetime(1970, 1, 1), 0, 0, 0, testing=True)

        # check that variable was handled correctly
        assert key in output
        assert reader.do_not_rotate == ["x_sea_water_velocity", "y_sea_water_velocity"]
        assert np.allclose(output[key].data, reader.Dataset[var][0,-1,0,0:2].values)
        # check that mask was used (can tell by mask value)
        # only mask_rho should have been used since u_eastward on rho grid
        assert "mask_u" not in masks
        assert np.allclose(masks["mask_rho"], reader.Dataset["mask_rho"][0,0:2].values)


    def test_get_variables_u_mask_u(self):
        """get a variable from the dataset and verify correct mask was used
        
        wetdry_mask_u is not available so mask_u should be used.
        """
        standard_mapping = {'u': 'x_sea_water_velocity'} 
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_mapping)
        # drop wetdry_mask_rho to test that mask_rho is used
        reader.Dataset = reader.Dataset.drop_vars("wetdry_mask_rho")
        var, key = "u", "x_sea_water_velocity"
        output, masks = reader.get_variables(key, datetime(1970, 1, 1), 0, 0, 0, testing=True)

        # check that variable was handled correctly
        assert key in output
        assert reader.do_not_rotate == []
        assert np.allclose(output[key].data, reader.Dataset[var][0,-1,0,0:2].values)
        # check that mask was used (can tell by mask value)
        # mask_rho should not have been used since u on u-grid
        assert "mask_rho" not in masks


    def test_get_variables_u_depth_zeta_zero(self):
        """get a variable from the dataset and verify in depth."""
        standard_mapping = {'u': 'x_sea_water_velocity'}
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_mapping)
        # drop wetdry_mask_rho to test that mask_rho is used
        reader.Dataset = reader.Dataset.drop_vars("wetdry_mask_rho")
        var, key = "u", "x_sea_water_velocity"
        zeta = np.zeros((2,2,3))
        reader.Dataset["zeta"] = (["ocean_time", "eta_rho", "xi_rho"], zeta)

        output = reader.get_variables(key, datetime(1970, 1, 1), 0, 0, [-5], testing=False)

        assert key in output
        assert reader.do_not_rotate == []
        # those expected values at 0.5, 0.6 from the input array
        assert np.allclose(output[key].data[1,0,:], self.ds[var][0,1,0,:].values)


    def test_get_variables_u_depth_zeta_nonzero(self):
        """get a variable from the dataset and verify in depth."""
        standard_mapping = {'u': 'x_sea_water_velocity'}
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_mapping)
        # drop wetdry_mask_rho to test that mask_rho is used
        reader.Dataset = reader.Dataset.drop_vars("wetdry_mask_rho")
        var, key = "u", "x_sea_water_velocity"
        zeta = np.ones((2,2,3))*2
        reader.Dataset["zeta"] = (["ocean_time", "eta_rho", "xi_rho"], zeta)

        output = reader.get_variables(key, datetime(1970, 1, 1), 0, 0, [-5], testing=False)

        u_expected = [0.64285714, 0.74285714]
        assert key in output
        assert reader.do_not_rotate == []
        assert np.allclose(output[key].data[1,0,:], u_expected)


class TestROMSReaderRotation(unittest.TestCase):
    def setUp(self):
        
        self.ds = xr.Dataset(
            data_vars={
                "u": (("ocean_time", "Z", "Y", "X"), np.zeros((2, 3, 2, 3))),
                "v": (("ocean_time", "Z", "Y", "X"), np.zeros((2, 3, 2, 3))),
                "Uwind": (("ocean_time", "Y", "X"), np.zeros((2, 2, 3))),
                "Vwind": (("ocean_time", "Y", "X"), np.zeros((2, 2, 3))),
                "Cs_r": (("Z"), np.linspace(-1, 0, 3)),
                "hc": 16,
            },
            coords={
                "ocean_time": ("ocean_time", [0, 1], {"units": "seconds since 1970-01-01"}),
                "s_rho": (("Z"), np.linspace(-1, 0, 3)),
                "lon_rho": (("Y", "X"), np.array([[1, 2, 3], [1, 2, 3]])),
                "lat_rho": (("Y", "X"), np.array([[1, 1, 1], [2, 2, 2]])),
            },
        )

    
    def test_rotate(self):
        """test that correct variables are rotated"""
        
        # currents and winds have normal ROMS names and should be rotated
        reader = reader_ROMS_native.Reader(self.ds)
        assert reader.do_not_rotate == []
        
        # having the "east" or "north" in the standard_name attributes also cause
        # a variable to not be rotated
        # this is to test the code though this standard mapping would not lead to 
        # these standard names actually being used outside of the ROMS reader
        # in opendrift
        standard_name_mapping = {"u": 'eastward_sea_water_velocity', "v": 'northward_sea_water_velocity'}
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_name_mapping)
        assert reader.do_not_rotate == ['eastward_sea_water_velocity', 'northward_sea_water_velocity']
        
        # currents have unusual clearly east/north variable names and should not
        # be rotated
        self.ds = self.ds.rename_vars({"u": "u_eastward", "v": "v_northward"})
        standard_name_mapping = {"u_eastward": 'x_sea_water_velocity', "v_northward": 'y_sea_water_velocity'}
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_name_mapping)
        assert reader.do_not_rotate == ['x_sea_water_velocity', 'y_sea_water_velocity']
        
        # analogous for winds
        standard_name_mapping = {"Uwind": 'eastward_wind', "Vwind": 'northward_wind'}
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_name_mapping)
        assert reader.do_not_rotate == ['eastward_wind', 'northward_wind']

        self.ds = self.ds.rename_vars({"Uwind": "Uwind_eastward", "Vwind": "Vwind_northward"})
        standard_name_mapping = {"Uwind_eastward": 'x_wind', "Vwind_northward": 'y_wind'}
        reader = reader_ROMS_native.Reader(self.ds, standard_name_mapping=standard_name_mapping)
        assert reader.do_not_rotate == ['x_wind', 'y_wind']
        
        
        
        

if __name__ == '__main__':
    unittest.main()
