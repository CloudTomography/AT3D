"""
Generation objects generate the atmospheric optical properties and grid. 
"""
from shdom import float_round
import numpy as np
import argparse
import shdom
import copy
from collections import OrderedDict

class Generator(object):
    """
    TODO 
       
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.
    """
    def __init__(self, args):
        self._args = args
        
    @classmethod
    def update_parser(self, parser): 
        """
        Update the argument parser with parameters relavent to this generation script. 
        
        Parameters
        ----------
        parser: argparse.ArgumentParser()
            The main parser to update.
    
        Returns
        -------
        parser: argparse.ArgumentParser()
            The updated parser.
        """
        return parser
    
    def get_grid(self):
        """
        This method gets the atmospheric grid.
        """
        return None
        
    def get_extinction(self, grid=None):
        """
        This method gets the atmospheric extinction.
        """
        return None
        
    def get_albedo(self, grid=None):
        """
        This method gets the atmospheric albedo.
        """
        self._albedo = albedo
    
    def get_phase(self, grid=None):
        """
        This method gets the atmospheric phase.
        """
        return None  
   
    
    def get_scatterer(self):
        """
        This method gets a shdom.Scatterer object.
        """
        return None
    
    @property
    def args(self):
        return self._args

    
class CloudGenerator(Generator):
    """TODO"""
    def __init__(self, args):
        super(CloudGenerator, self).__init__(args)
        self._mie = OrderedDict()
        
    def add_mie(self, table_path):
        mie = shdom.MiePolydisperse()
        mie.read_table(table_path)
        self._mie[mie.wavelength] = copy.deepcopy(mie)
        
    def get_lwc(self, grid=None):
        """
        This method gets the lwc.
        """
        return None     
      
    def get_reff(self, grid=None):
        """
        This method gets the reff.
        """
        return None          

    def get_veff(self, grid=None):
        """
        This method gets the veff.
        """
        return None
    
    def get_albedo(self, wavelength, grid=None):
        if grid is None:
            grid = self.get_grid()
        reff = self.get_reff(grid)
        veff = self.get_veff(grid)        
        return self.mie[float_round(wavelength)].get_albedo(reff, veff)
    
    def get_phase(self, wavelength, grid=None):
        if grid is None:
            grid = self.get_grid()  
        reff = self.get_reff(grid)
        veff = self.get_veff(grid)          
        return self.mie[float_round(wavelength)].get_phase(reff, veff)
    
    def get_scatterer(self):
        """
        Return a shdom.Scatterer object
        """
        lwc = self.get_lwc()
        reff = self.get_reff()
        veff = self.get_veff()
        if (lwc is not None) and (reff is not None) and (veff is not None):
            scatterer = shdom.MicrophysicalScatterer(lwc, reff, veff)
            for mie in self.mie.values():
                scatterer.add_mie(mie)
        
        else:
            if len(self.mie) == 1:
                wavelength = next(iter(self.mie.keys()))
                scatterer = shdom.OpticalScatterer(
                    wavelength, 
                    self.get_extinction(wavelength), 
                    self.get_albedo(wavelength), 
                    self.get_phase(wavelength)
                )
            
            else:
                scatterer = shdom.MultispectralScatterer()
                for wavelength in self.mie.keys():
                    scatterer.add_scatterer(shdom.OpticalScatterer(
                        wavelength, 
                        self.get_extinction(wavelength), 
                        self.get_albedo(wavelength), 
                        self.get_phase(wavelength)
                    ))
                                 
        return scatterer    

    @property
    def mie(self):
        return self._mie    


class AirGenerator(Generator):
    """
    Generate air optical properties according to rayleigh molecular scattering.
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.
    wavelength: float,
        The wavelength in [microns]
    """
    def __init__(self, args):
        super(AirGenerator, self).__init__(args)
        
    @classmethod
    def update_parser(self, parser): 
        parser.add_argument('--air_max_alt',
                            default=20.0,
                            type=np.float32,
                            help='(default value: %(default)s) Maximum altitude for the air volume.')
        
        parser.add_argument('--air_num_points',
                            default=20,
                            type=int,
                            help='(default value: %(default)s) Number of altitude grid points for the air volume.')     
        return parser        
    
    def set_temperature_profile(self, temperature_profile):
        self._temperature_profile = temperature_profile     
        
    def get_scatterer(self, wavelength):
        """
        Parameters
        -----------
        wavelength: float or list,
            The wavelength or wavelength list in [microns] 
            
        Returns
        -------
        scatterer: shdom.Scatterer or shdom.MultispectralScatterer
        """
        altitudes = shdom.Grid(z=np.linspace(0.0, self.args.air_max_alt, self.args.air_num_points))       
        wavelengths = np.array(wavelength)
        num_bands = wavelengths.size
        
        if  num_bands == 1:
            rayleigh = shdom.Rayleigh(wavelength)
            rayleigh.set_profile(self.temperature_profile.resample(altitudes))
            scatterer = rayleigh.get_scatterer()
        else:
            scatterer = shdom.MultispectralScatterer()
            for i in range(num_bands):
                rayleigh = shdom.Rayleigh(wavelength[i])
                rayleigh.set_profile(self.temperature_profile.resample(altitudes))
                scatterer.add_scatterer(rayleigh.get_scatterer())   
        return scatterer
    
    @property
    def temperature_profile(self):
        return self._temperature_profile    
    
    
class AFGLSummerMidLatAir(AirGenerator):
    """
    Temperature profile is taken from AFGL measurements of summer mid-lat.
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.   
    """
    def __init__(self, args):
        super(AFGLSummerMidLatAir, self).__init__(args)
        temperatures = np.array([292.220, 292.040, 291.860, 291.680, 291.500, 291.320, 291.140, 290.960, 290.780, 
                                 290.600, 290.420, 290.240, 290.060, 289.880, 289.700, 289.920, 290.140, 290.360, 
                                 290.580, 290.800, 291.020, 291.240, 291.460, 291.680, 291.900])
        temp_grid = shdom.Grid(z=np.linspace(0.0, 20.0, len(temperatures)))
        temperature_profile = shdom.GridData(temp_grid, temperatures)
        self.set_temperature_profile(temperature_profile)
       
        
class SingleVoxel(CloudGenerator):
    """
    Define a Medium with a single voxel in center of the grid. 
    It is useful for developing and debugging of derivatives and sensitivity analysis.
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.
    """
    def __init__(self, args):
        super(SingleVoxel, self).__init__(args)

    @classmethod
    def update_parser(self, parser):
        parser.add_argument('--nx', 
                            default=10,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in x (North) direction')
        parser.add_argument('--ny', 
                            default=10,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in y (East) direction')
        parser.add_argument('--nz', 
                            default=10,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in z (Up) direction')
        parser.add_argument('--domain_size', 
                            default=1.0,
                            type=float, 
                            help='(default value: %(default)s) Cubic domain size [km]')
        parser.add_argument('--extinction', 
                            default=10.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Extinction of the center voxel [km^-1]')
        parser.add_argument('--lwc',
                            default=None,
                            type=np.float32, 
                            help='(default value: %(default)s) Liquid water content of the center voxel [g/m^3]. If specified, extinction argument is ignored.')
        parser.add_argument('--reff', 
                            default=10.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective radius of the center voxel [micron]')
        parser.add_argument('--veff', 
                            default=0.1,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective variance of the center voxel')
        parser.add_argument('--mie_table_path',
                            help='Path to a precomputed Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')
        return parser      


    def get_grid(self):
        bb = shdom.BoundingBox(0.0, 0.0, 0.0, self.args.domain_size, self.args.domain_size, self.args.domain_size)
        return shdom.Grid(bounding_box=bb, nx=self.args.nx, ny=self.args.ny, nz=self.args.nz)
            
    def get_extinction(self, wavelength=None, grid=None):
        if grid is None:
            grid = self.get_grid()
            
        if self.args.lwc is None:
            ext_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
            ext_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.extinction
            extinction = shdom.GridData(grid, ext_data)
        else:
            assert wavelength is not None, 'No wavelength provided'
            lwc = self.get_lwc(grid)
            reff = self.get_reff(grid)
            veff = self.get_veff(grid)
            extinction = self.mie[float_round(wavelength)].get_extinction(lwc, reff, veff)
        return extinction

    def get_lwc(self, grid=None):
        if grid is None:
            grid = self.get_grid()
            
        lwc = self.args.lwc
        if lwc is not None:
            lwc_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
            lwc_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.lwc 
            lwc = shdom.GridData(grid, lwc_data)        
        return lwc
        
    def get_reff(self, grid=None):
        if grid is None:
            grid = self.get_grid()        
        reff_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
        reff_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.reff 
        return shdom.GridData(grid, reff_data)
        
    def get_veff(self, grid=None):
        if grid is None:
            grid = self.get_grid()        
        veff_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
        veff_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.veff    
        return shdom.GridData(grid, veff_data)    


class Homogeneous(CloudGenerator):
    """
    Define a homogeneous Medium. 

    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.    
    """
    def __init__(self, args):
        super(Homogeneous, self).__init__(args)
   
    @classmethod
    def update_parser(self, parser):
        parser.add_argument('--nx', 
                            default=10,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in x (North) direction')
        parser.add_argument('--ny', 
                            default=10,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in y (East) direction')
        parser.add_argument('--nz', 
                            default=10,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in z (Up) direction')
        parser.add_argument('--domain_size', 
                            default=1.0,
                            type=float, 
                            help='(default value: %(default)s) Cubic domain size [km]')        
        parser.add_argument('--extinction', 
                            default=1.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Extinction [km^-1]')  
        parser.add_argument('--lwc', 
                            default=None,
                            type=np.float32, 
                            help='(default value: %(default)s) Liquid water content of the center voxel [g/m^3]. If specified, extinction argument is ignored.')
        parser.add_argument('--reff', 
                            default=10.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective radius [micron]')
        parser.add_argument('--veff', 
                            default=0.1,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective variance')        
        parser.add_argument('--mie_table_path', 
                            help='Path to a precomputed polydisperse Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')    
        return parser        

    def get_grid(self):
        bb = shdom.BoundingBox(0.0, 0.0, 0.0, self.args.domain_size, self.args.domain_size, self.args.domain_size)
        return shdom.Grid(bounding_box=bb, nx=self.args.nx, ny=self.args.ny, nz=self.args.nz)
        
    def get_extinction(self, wavelength=None, grid=None):
        if grid is None:
            grid = self.get_grid()
            
        if self.args.lwc is None:
            ext_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=self.args.extinction, dtype=np.float32)  
            extinction = shdom.GridData(grid, ext_data)
        else:
            assert wavelength is not None, 'No wavelength provided'
            lwc = self.get_lwc(grid)
            reff = self.get_reff(grid)
            veff = self.get_veff(grid)
            extinction = self.mie[float_round(wavelength)].get_extinction(lwc, reff, veff)
        return extinction
        
    def get_lwc(self, grid=None):
        if grid is None:
            grid = self.get_grid()
            
        lwc = self.args.lwc
        if lwc is not None:
            lwc_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=lwc, dtype=np.float32)        
            lwc = shdom.GridData(grid, lwc_data)        
        return lwc
        
    def get_reff(self, grid=None):
        if grid is None:
            grid = self.get_grid()     
        reff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=self.args.reff, dtype=np.float32)
        return shdom.GridData(grid, reff_data)
        
    def get_veff(self, grid=None):
        if grid is None:
            grid = self.get_grid()        
        veff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=self.args.veff, dtype=np.float32)    
        return shdom.GridData(grid, veff_data)
        
        
class LesFile(CloudGenerator):
    """
    Generate an optical medium from an LES file which contains a table with: lwc, reff, veff (optional).
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.
    """
    def __init__(self, args):
        super(LesFile, self).__init__(args)
        self._droplets = shdom.MicrophysicalScatterer()
        self._droplets.load_from_csv(args.path, args.veff)
        
    @classmethod
    def update_parser(self, parser): 
        parser.add_argument('--path',
                            help='Path to the LES generated file')
        parser.add_argument('--veff', 
                            default=0.1,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective variance (if not provided as a 3D field)')         
        parser.add_argument('--mie_table_path', 
                            help='Path to a precomputed Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')    
        return parser
    
    def get_grid(self):
        return self._droplets.grid
    
    def get_scatterer(self):
        for mie in self.mie.values():
            self._droplets.add_mie(mie)
        return self._droplets