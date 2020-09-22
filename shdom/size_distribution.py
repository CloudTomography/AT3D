import xarray as xr
import numpy as np
from shdom import core
from collections import OrderedDict

#TODO Add normalization options for size distributions.

def gamma(radii, reff=None,veff=None,alpha=None, particle_density=1.0):
    """
    TODO

    Takes
    """
    assert reff is not None, 'Reff must be specified'
    reff = np.atleast_1d(reff)
    assert (veff is not None) or (alpha is not None), 'One of veff and alpha must be specified.'
    if veff is not None:
        assert len(reff) == len(veff)
        alpha = 1.0/ np.atleast_1d(veff) - 3.0

    elif alpha is not None:
        assert len(reff) == len(alpha)
        alpha = np.atleast_1d(alpha)

    #fortran subroutine requires a 'gamma' variable to be passed.
    gamma = np.zeros(alpha.shape)

    nd = core.make_multi_size_dist(
                distflag='G',
                pardens=particle_density,
                nsize=len(radii),
                radii=radii,
                reff=reff,
                alpha=alpha,
                gamma=gamma,
                ndist=reff.size)
    return nd

def lognormal(radii, reff=None ,veff=None,alpha=None, particle_density=1.0):
    """
    TODO
    """
    assert reff is not None, 'Reff must be specified'
    reff = np.atleast_1d(reff)
    assert (veff is not None) or (alpha is not None), 'One of veff and alpha must be specified.'
    if veff is not None:
        assert len(reff) == len(veff)
        alpha = np.sqrt(np.log(np.atleast_1d(veff) + 1.0))

    elif alpha is not None:
        assert len(reff) == len(alpha)
        alpha = np.atleast_1d(alpha)

    #fortran subroutine requires a 'gamma' attribute to be passed.
    gamma = np.zeros(alpha.shape)

    nd = core.make_multi_size_dist(
                distflag='L',
                pardens=particle_density,
                nsize=len(radii),
                radii=radii,
                reff=reff,
                alpha=alpha,
                gamma=gamma,
                ndist=reff.size)
    return nd



def get_size_distribution_grid(radii, size_distribution_function=gamma,
                               particle_density=1.0, radius_units='micron',
                               **size_distribution_parameters):
    """
    TODO
    A generic interface to get number density of a size distribution
    on a grid of size distribution parameters.

    Each size_distribution_parameter should contain:
       coord_min: float
           The minimum value for that coordinate.
       coord_max: float
           The maximum value for that coordinate.
       n: integer
           The number of points sampling this dimension.
       spacing: string
           The type of spacing. Either 'logarithmic' or 'linear'.
       units: string
           The units of the microphysical dimension.

    """
    coord_list = []
    name_list = []
    size_distribution_list = []

    for name,parameter in size_distribution_parameters.items():
        try:
            coord_min, coord_max,n,spacing,units = parameter
        except:
            coord_min, coord_max,n,spacing = parameter
            units = 'Not specified'
        if spacing == 'logarithmic':
            coord = np.logspace(np.log10(coord_min), np.log10(coord_max), n)
        elif spacing =='linear':
            coord = np.linspace(coord_min, coord_max, n)
        else:
            raise NotImplementedError
        coord_list.append(coord)
        size_distribution_list.append(name + ' [{}, {}, {}, {}, {}]'.format(coord_min, coord_max,n,spacing,units))
        name_list.append(name)

    meshgrid = np.meshgrid(*coord_list,indexing='ij')
    parameter_dict = OrderedDict()
    coord_dict = OrderedDict()
    coord_dict['radius'] = radii
    size_distribution_list.append('radius units [{}]'.format(radius_units))
    for name,grid,coord in zip(name_list,meshgrid,coord_list):

        parameter_dict[name] = grid.ravel()
        coord_dict[name] =coord

    grid_shape = [len(coord) for name,coord in coord_dict.items()]

    number_density_raveled = size_distribution_function(radii, **parameter_dict,
                                              particle_density=particle_density)
    #TODO this can fail silently in some cases (produce nans), add checks.
    number_density = number_density_raveled.reshape(grid_shape)

    size_dist_grid = xr.Dataset(
            data_vars = {
                'number_density': (list(coord_dict.keys()), number_density),
            },
            coords=dict(coord_dict),
            attrs={'size_distribution_inputs': size_distribution_list,
                  'distribution_type': size_distribution_function.__name__,
                  },
                  #TODO add coordmin/max/n and spacing for full traceability.
    )
    return size_dist_grid
