"""
Generation objects generate the atmospheric optical properties and grid. 
"""

import os 
import numpy as np
import argparse
import shdom

class Generator(object):
    """
    A generation object has the following components:

    1. update_parser(parser) class method used for argparsing
    2. init_grid() method
    3. init_extinction() method
    4. init_albedo() method
    5. init_phase() method
    6. get_medium() method
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.
    """
    def __init__(self, args):
        self._args = args
        self.init_grid()
        self.init_extinction()
        self.init_albedo()
        self.init_phase()

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
    
    def set_grid(self, grid):
        """
        This method sets the atmospheric grid.
        """
        self._grid = grid
        
    def init_grid(self):
        """
        This method generates the atmospheric grid.
        """
        self._grid = None
    
    def set_extinction(self, extinction):
        """
        This method sets the atmospheric extinction.
        """
        self._extinction = extinction    
        
    def init_extinction(self):
        """
        Initialize the optical extinction field.
        """
        self._extinction = None
        
    def set_albedo(self, albedo):
        """
        This method sets the atmospheric albedo.
        """
        self._albedo = albedo
        
    def init_albedo(self):
        """
        Initialize the single scattering albedo field.
        """
        self._albedo = None
    
    def set_phase(self, phase):
        """
        This method sets the atmospheric phase.
        """
        self._phase = phase    
        
    def init_phase(self):
        """
        Initialize the phase function field.
        """
        self._phase = None
    
    def get_medium(self):
        """
        Return a shdom.Medium object with the optical properties on a grid.
        """
        medium = shdom.Medium()
        medium.set_optical_properties(self.extinction, self.albedo, self.phase)
        return medium
    
    @property
    def args(self):
        return self._args
    
    @property
    def grid(self):
        if self._grid is None:
            self.init_grid()
        return self._grid
    
    @property
    def extinction(self):
        if self._extinction is None:
            self.init_extinction()
        return self._extinction
    
    @property
    def albedo(self):
        if self._albedo is None:
            self.init_albedo()
        return self._albedo
    
    @property
    def phase(self):
        if self._phase is None:
            self.init_phase()
        return self._phase    

    
class Air(Generator):
    """
    Generate air optical properties according to rayleigh molecular scattering.
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator. 
    """
    def __init__(self, args):
        super(Air, self).__init__(args)
        
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
        self._rayleigh = shdom.Rayleigh(self.args.wavelength)
        self.rayleigh.init_temperature_profile(temperature_profile)        
        
        
    def get_medium(self, phase_type):
        altitudes = shdom.Grid(z=np.linspace(0.0, self.args.air_max_alt, self.args.air_num_points))
        
        if phase_type == 'Tabulated':
            medium = shdom.AmbientMedium() 
        elif phase_type == 'Grid':
            medium = shdom.Medium()
    
        extinction, albedo, phase = self.rayleigh.get_scattering_field(altitudes, phase_type)
        medium.set_optical_properties(extinction, albedo, phase)        
        return medium

    @property
    def temperature_profile(self):
        return self._temperature_profile    
    
    @property
    def rayleigh(self):
        return self._rayleigh  
    
class AFGLSummerMidLatAir(Air):
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
       
        
class SingleVoxel(Generator):
    """
    Define a Medium with a single voxel in center of the grid. 
    It is useful for developing and debugging of derivatives and gradients and sensitivity analysis.
    
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
                            default=3,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in x (North) direction')
        
        parser.add_argument('--ny', 
                            default=3,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in y (East) direction')
    
        parser.add_argument('--nz', 
                            default=3,
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
    
        parser.add_argument('--albedo', 
                            default=1.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Albedo of the center voxel in range [0, 1]')
    
        parser.add_argument('--reff', 
                            default=10.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective radius of the center voxel [micron]')
    
        parser.add_argument('--mie_table_path', 
                            help='Path to a precomputed Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')
        return parser        


    def init_grid(self):
        bb = shdom.BoundingBox(0.0, 0.0, 0.0, self.args.domain_size, self.args.domain_size, self.args.domain_size)
        self._grid = shdom.Grid(bounding_box=bb, nx=self.args.nx, ny=self.args.ny, nz=self.args.nz) 
        
    def init_extinction(self): 
        ext_data = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32)
        ext_data[self.grid.nx/2, self.grid.ny/2, self.grid.nz/2] = self.args.extinction
        self._extinction = shdom.GridData(self.grid, ext_data)  
    
    def init_albedo(self):   
        alb_data = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32)
        alb_data[self.grid.nx/2, self.grid.ny/2, self.grid.nz/2] = self.args.albedo
        self._albedo = shdom.GridData(self.grid, alb_data)
 
    def init_phase(self):
        mie = shdom.Mie()
        if self.args.mie_table_path:
            mie.read_table(self.args.mie_table_path)
        else:
            mie.set_parameters((self.args.wavelength, self.args.wavelength), 'Water', 'gamma', 7.0)
            mie.compute_table(1, self.args.reff, self.args.reff, 65.0)   
        
        reff_data = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32) 
        reff_data[self.grid.nx/2, self.grid.ny/2, self.grid.nz/2] = self.args.reff
        self._reff = shdom.GridData(self.grid, reff_data)
        self._phase = mie.get_grid_phase(self._reff)
    
    
class Homogeneous(Generator):
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
                            default=3,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in x (North) direction')
        parser.add_argument('--ny', 
                            default=3,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in y (East) direction')
        parser.add_argument('--nz', 
                            default=3,
                            type=int, 
                            help='(default value: %(default)s) Number of grid cell in z (Up) direction')
        parser.add_argument('--domain_size', 
                            default=1.0,
                            type=float, 
                            help='(default value: %(default)s) Cubic domain size [km]')        
        parser.add_argument('--extinction', 
                            default=1.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Extinction of a homogeneous cloud [km^-1]')  
        parser.add_argument('--albedo', 
                            default=1.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Albedo of a homogeneous cloud in range [0, 1]')
        parser.add_argument('--reff', 
                            default=10.0,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective radius of a homogeneous cloud [micron]')
        parser.add_argument('--mie_table_path', 
                            help='Path to a precomputed Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')    
        return parser        

    def init_grid(self):
        bb = shdom.BoundingBox(0.0, 0.0, 0.0, self.args.domain_size, self.args.domain_size, self.args.domain_size)
        self._grid = shdom.Grid(bounding_box=bb, nx=self.args.nx, ny=self.args.ny, nz=self.args.nz)
        
    def init_extinction(self): 
        ext_data = np.full(shape=(self.grid.nx, self.grid.ny, self.grid.nz), fill_value=self.args.extinction, dtype=np.float32)
        self._extinction = shdom.GridData(self.grid, ext_data)
    
    def init_albedo(self):   
        alb_data = np.full(shape=(self.grid.nx, self.grid.ny, self.grid.nz), fill_value=self.args.albedo, dtype=np.float32)
        self._albedo = shdom.GridData(self.grid, alb_data)
    
    def init_phase(self):
        mie = shdom.Mie()
        if self.args.mie_table_path:
            mie.read_table(self.args.mie_table_path)
        else:
            mie.set_parameters((self.args.wavelength, self.args.wavelength), 'Water', 'gamma', 7.0)
            mie.compute_table(1, self.args.reff, self.args.reff, 65.0)   
        
        reff_data = np.full(shape=(self.grid.nx, self.grid.ny, self.grid.nz), fill_value=self.args.reff, dtype=np.float32)        
        self._reff = shdom.GridData(self.grid, reff_data)
        self._phase = mie.get_grid_phase(self._reff)
        
            
class HomogeneousPolarized(Homogeneous):
    """
    Define a homogeneous Medium with Polarized Mie phase function. 
   
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.    
    """
    
    def __init__(self, args):
        super(HomogeneousPolarized, self).__init__(args)
        
    def init_phase(self):
        mie = shdom.MiePolarized()
        if self.args.mie_table_path:
            mie.read_table(self.args.mie_table_path)
        else:
            mie.set_parameters((self.args.wavelength, self.args.wavelength), 'Water', 'gamma', 7.0)
            mie.compute_table(1, self.args.reff, self.args.reff, 65.0)   
        
        reff_data = np.full(shape=(self.grid.nx, self.grid.ny, self.grid.nz), fill_value=self.args.reff, dtype=np.float32)        
        self._reff = shdom.GridData(self.grid, reff_data)
        self._phase = mie.get_tabulated_phase(self._reff)
        
class LesFile(Generator):
    """
    Generate an optical medium from an LES file which contains
    a table with: lwc, reff, veff (optional).
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this generator.
        
    Notes
    -----
    veff in table not supported yet.
    """
    def __init__(self, args):
        mie = shdom.Mie()
        mie.read_table(args.mie_table_path)  
        microphysics = shdom.MicrophysicalMedium()
        microphysics.load_from_csv(args.path)
    
        reff = microphysics.reff.data[microphysics.lwc.data > 0.0]
        if reff.min() < mie.reff.min():
            print('Minimum effective radius in file ({}) is smaller than ' \
                  'minimum effective radius in table ({})'.format(reff.min(), mie.reff.min()))
    
        if reff.max() > mie.reff.max():
            print('Maximum effective radius in file ({}) is larger than ' \
                  'maximum effective radius in table ({}) '.format(reff.max(), mie.reff.max()))    
        
        self._mie = mie
        self._microphysics = microphysics 
        super(LesFile, self).__init__(args)
        
    @classmethod
    def update_parser(self, parser): 
        parser.add_argument('--path', 
                            help='Path to the LES generated file')
        
        parser.add_argument('--mie_table_path', 
                            help='Path to a precomputed Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details') 
        
        parser.add_argument('--phase_type', 
                            choices=['Tabulated', 'Grid'],
                            default='Tabulated',
                            help='Phase function type, see. \
                            notebooks/Forward Rendering.ipynb for more details')     
        return parser
    
    def init_grid(self):
        microphysics = shdom.MicrophysicalMedium()
        grid = microphysics.get_grid(self.args.path)
        return grid
        
    def init_extinction(self):
        self._extinction = self._mie.interpolate_extinction(self._microphysics.lwc, self._microphysics.reff)        
    
    def init_albedo(self):
        self._albedo = self._mie.interpolate_albedo(self._microphysics.reff)   

    def init_phase(self):
        self._phase = self._mie.interpolate_phase(self._microphysics.reff, self.args.phase_type)
            