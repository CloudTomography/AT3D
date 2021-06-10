"""
A module containing tools for generating synthetic data.
"""
import numpy as np
import scipy.stats as st
import numpy.fft as fft
import copy


class SurrogateGenerator:
    """
    Generates statistical surrogates of a supplied data field using the
    Iterative Amplitude Adjusted Fourier Transform algorithm (IAAFT).

    This is useful for producing new realizations of cloud data which may
    have been computationally expensive to produce (e.g. from Large Eddy Simulation).

    The IAAFT matches the distribution of the data and attempts to match the
    amplitudes of the power spectrum. This is done by starting from
    a random shuffle of the data and then iteratively adjusting the spectrum and
    the distribution until the power spectrum stops changing. This algorithm
    may reach local minima. Other refinements to the IAAFT
    such as the Wavelet-initialization or Stochastic variants could be implemented
    to further improve this algorithm.

    Parameters
    ----------
    data : np.ndarray
        The array of data to produce surrogates of.

    Reference
    ---------
    VENEMA, V., MEYER, S., GARCÍA, S.G., KNIFFKA, A., SIMMER, C., CREWELL, S.,
    LÖHNERT, U., TRAUTMANN, T. and MACKE, A. (2006),
    Surrogate cloud fields generated with the iterative amplitude adapted Fourier
    transform algorithm. Tellus A, 58: 104-120.
    https://doi.org/10.1111/j.1600-0870.2006.00160.x
    """
    def __init__(self, data):

        if not isinstance(data, np.ndarray):
            raise TypeError("`data` input should be of type {}".format(type(np.ndarray)))

        self._reference_data = data
        self._normalized_data = (data - data.mean())/data.std()
        self._reference_spectrum = fft.fftn(self._normalized_data)
        self._ranked_data = np.sort(self._normalized_data.ravel())
        self._reference_fourier_amplitudes = np.abs(self._reference_spectrum)

    def __call__(self, err_tol=1e-4, maxiter=1000, normalized=False):
        """
        Produces a surrogate through iterative adjustment.

        Notes
        -----
        This algorithm is initialized with a random shuffle so for reproducibility
        the random number generator should be seeded before calls to this method.
        e.g. by running `np.random.seed(1)`.

        Returns
        -------
        surrogate_field : np.ndarray
            The statistical surrogate of the reference data.

        Parameters
        ----------
        err_tol : float
            The termination criterion for the iterative adjustment of the
            surrogate which is the relative reduction in the mean square
            error of the absolute values of the power spectrum.
        maxiter : int
            The maximum number of iterations for adjusting the surrogate to the
            reference data.
        normalized : bool
            If True then the surrogate is returned with unit standard deviation
            and zero mean. If False then these two moments are exactly as in
            the reference data.
        """

        data_to_shuffle = copy.deepcopy(self._normalized_data)
        shuffled_indices = np.meshgrid(*[np.arange(shape) for shape in self._normalized_data.shape], indexing='ij')
        for inds in shuffled_indices:
            np.random.shuffle(inds.ravel())

        surrogate_field = data_to_shuffle[tuple(shuffled_indices)]

        old_cost = None
        for i in range(maxiter): # This could be a while loop but whatever.
            temp_spectra = fft.fftn(surrogate_field)
            cost = np.sqrt(np.mean((self._reference_fourier_amplitudes - np.abs(temp_spectra))**2))
            if old_cost is not None:
                if np.abs(cost - old_cost)/old_cost < err_tol:
                    break

            adjusted_temp_spectra = self._reference_fourier_amplitudes*temp_spectra/np.abs(temp_spectra)
            spectra_adjusted = fft.ifftn(adjusted_temp_spectra).real
            temp_ranks = st.rankdata(spectra_adjusted).astype(np.int) -1
            amplitude_adjusted_data = self._ranked_data[temp_ranks].reshape(self._normalized_data.shape)
            surrogate_field = amplitude_adjusted_data
            old_cost = cost

        if not normalized:
            surrogate_field = surrogate_field*self._reference_data.std() + self._reference_data.mean()

        return surrogate_field


def generate_stochastic_blob(rte_grid, var_profile, var_beta=-5.0/3.0, snr=1.0,
                 var_vertical_correlation=0.0, xy_blob=0.2, z_blob=0.2,
                 cloud_mask_vertical_correlation=0.8, volume_cloud_fraction=0.1,
                 cloud_mask_beta=-3.0, remove_edges=True):
    """
    Generates a stochastic blob variable that is 'cumulus-like'.

    Masks the `rte_grid` using a cloud mask generated from log-normal red noise with
    an isotropic power spectrum with power-law slope `cloud_mask_beta`. The
    power law slope is typically steeper than the variable to make a smooth cloud shell.
    The red noise is linearly mixed in the vertical to induce vertical correlations
    with correlation coefficient `cloud_mask_vertical_correlation` to reduce holes
    internal to the cloud. The log-normal red noise is filtered with a Butter worth filter
    to produce a 'blob' shape. A small fraction of the domain `volume_cloud_fraction`
    is set to cloudy to avoid the rectangular edges of the domain affecting the cloud.

    The cloudy voxels are populated by a similarly generated log-normal red noise field
    however without any ButterWorth filtering and typically with a -5/3 spectral slope.
    The variable has a mean that follows a given vertical profile and with a specified
    relative standard deviation.

    Parameters
    ----------
    rte_grid : xr.Dataset
        A valid SHDOM grid object (see grid.py)
    var_profile : np.ndarray
        A 1D array of floats that describe the mean of the variable as a function of altitude.
        This should be the same length as rte_grid.z
    var_beta: float
        The slope of the isotropic power spectrum of the variable S(k) ~ k^(var_beta).
        Typically -5.0/3.0 for liquid water content.
    snr : float
        The relative standard deviation that determines the magnitude of the variability
        in the output variable.
    var_vertical_correlation, cloud_mask_vertical_correlation : float
        Vertical correlation coefficient for the stochastic white noise for the generated variable
        or cloud mask. In range [0, 1], larger values impose additional vertical smoothness.
    xy_blob, z_blob : float
        Scaling parameters for the Butterworth filter that controls how smooth the cloud mask will be.
        Smaller values increase smoothness. Strictly positive.
    volume_cloud_fraction : float
        The fraction of the domain's voxels that will have cloud.
    cloud_mask_beta : float
        The slope of the isotropic power spectrum of the stochastic variable used
        to determine the volumetric cloud mask S(k) ~ k^(var_beta). A value of -3.0 was decided on
        to make the cloud smooth at the scale of the voxels.
    remove_edges : Boolean
        Sets the edge points to zero in `var` if True.

    Returns
    -------
    var : np.ndarray
        The masked stochastically generated field.
    """
    nx, ny, nz = rte_grid.x.size, rte_grid.y.size, rte_grid.z.size
    domain_x, domain_y, domain_z = rte_grid.x.data.max(), rte_grid.y.data.max(), rte_grid.z.data.max()

    mask_generator = GaussianFieldGenerator(
        nx, ny, nz, cloud_mask_beta, [domain_x, domain_y, domain_z]
    )
    cloud_mask_field_temp = mask_generator.generate_field()

    #add vertical correlations.
    cloud_mask_field = cloud_mask_field_temp.copy()
    for i in range(cloud_mask_field.shape[-1]-1):
        cloud_mask_field[:,:,i+1] = cloud_mask_vertical_correlation*cloud_mask_field[:,:,i] \
        + np.sqrt(1.0 - cloud_mask_vertical_correlation)*cloud_mask_field_temp[:,:,i+1]

    cloud_mask_field = np.exp(cloud_mask_field)

    #apply butterworth filter.
    X,Y,Z = np.meshgrid((rte_grid.x.data - rte_grid.x.data.mean())/domain_x,
                        (rte_grid.y.data - rte_grid.y.data.mean())/domain_y,
                        (rte_grid.z.data - rte_grid.z.data.mean())/domain_z,
                        indexing='ij')
    R = np.sqrt(X**2+Y**2)/xy_blob
    butter_worth_filter = 1.0/np.sqrt(1 + R**8)* 1.0/np.sqrt(1.0 + (Z/z_blob)**8)

    cloud_mask_field *= butter_worth_filter
    threshold = np.percentile(cloud_mask_field,(1.0-volume_cloud_fraction)*100)
    mask = cloud_mask_field > threshold

    var_profile[var_profile <= 0.0] = 0.0
    var_mean = np.repeat(np.repeat(var_profile[np.newaxis, np.newaxis, :], rte_grid.x.size, axis=0),
                         rte_grid.y.size, axis=1)

    var_generator = GaussianFieldGenerator(
        nx, ny, nz, var_beta, [domain_x, domain_y, domain_z]
    )
    field_temp = var_generator.generate_field()

    field = field_temp.copy()
    for i in range(field.shape[-1]-1):
        field[:,:,i+1] = var_vertical_correlation*field[:,:,i] \
        + np.sqrt(1.0 - var_vertical_correlation)*field_temp[:,:,i+1]

    exp_field = np.exp(field[mask])
    exp_field -= exp_field.mean()
    exp_field /= exp_field.std()

    var_values = exp_field*snr*var_mean[mask] + var_mean[mask]
    var = np.zeros(field.shape)
    var[mask] = var_values

    if remove_edges:
    #Set edge points to zero to avoid interaction with open boundary conditions.
        var[0,:,:] = var[-1,:,:] = var[:,0,:] = var[:,-1,:] = var[...,0] = var[...,-1] = 0.0

    return var


class GaussianFieldGenerator:
    """
    Stochastically generates 3D fields in an equispaced cubic domain
    drawn from truncated normal fields with an isotropic power-law
    power spectrum.

    Parameters
    ----------
    nx, ny, nz: int,
        The number of points in the first, second, and thid dimensions of the field.
    beta: float,
        The exponent of the power-law controlling the isotropic spectral slope.
    domain_size: list of floats,
        The size of each of dimension of the domain in physical units.
    field_min, field_max: float,
        The minimum and maximum values of the 'standard' normal distribution at
        which to truncate.
    inner_scale: float,
        The length scale at which fourier components at isotropic frequencies
        above 1.0/inner_scale will be set to zero. Units are the same as domain_size.
    outer_scale: float,
        The length scale at which fourier components at isotropic frequencies
        below 1.0/outer_scale are set to a constant value. Units are the same as
    domain_size.
    seed, int (< 2**32)
        A seed for the random number generator to ensure reproducibility of the
        results.
    """
    def __init__(self, nx, ny, nz, beta, domain_size, field_min=None, field_max=None,
                 inner_scale=None, outer_scale=None, seed=None):

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.beta = beta
        self.domain_size = domain_size


        if field_min is not None:
            self.field_min = field_min
        elif field_max is not None:
            self.field_min = -10**-9
        else:
            self.field_min = None

        if field_max is not None:
            self.field_max = field_max
        elif field_min is not None:
            self.field_max = 10**9
        else:
            self.field_max = None


        if (outer_scale is not None) and (inner_scale is not None):
            assert outer_scale > inner_scale,'power_spectrum outer_scale {} must be \
            larger (greater) than power_spectrum inner_scale {}'.format(outer_scale, inner_scale)

        self.outer_scale = outer_scale
        self.inner_scale = inner_scale

        if seed is not None:
            np.random.seed(seed)
            self.seed = seed
        else:
            self.seed = None

    def _update_args(self, kwargs):
        """
        Updates the attributes of the generator.
        Parameters
        ---------
        kwargs: dictionary
        """
        for name in kwargs:
            setattr(self, name, kwargs[name])

    def generate_field(self, return_spectra=False, **kwargs):
        """
        Generates the field with the specified parameters.
        Any parameters of the generator can be updated by supplying them as kwargs.

        Parameters
        ---------
        return_spectra: Boolean
            Decides whether the spectra and frequencies of the field are
            also returned by the function.
        kwargs: dictionary
            A dictionary of updated parameters for the generator.
            See __init__ for a list of possiblities.
        Returns
        -------
        output_field: float, shape=(nx, ny, nz)
            The generated field.
        if return_spectra is True this function also returns:
        power_spectra: float, shape=(nx, ny, nz)
            The fourier transform of output_field
        X,Y,Z: float, shape=(nx, ny, nz)
            The frequencies in the first, second, and third dimensions
            respectively.
            """

        self._update_args(kwargs)

        # Generate fields, ensuring that they are sufficiently gaussian.
        skew = 3
        kurtosis = 5
        optimal_field = None
        min_skew = 10**9
        min_kurtosis = 10**9
        maxiter = 500
        counter = 0

        while ((np.abs(skew)>0.1) or (np.abs(kurtosis-0.0)>0.5)) and counter<maxiter:
            data_inv = self._initialize_field()
            skew = st.skew(data_inv.ravel())
            kurtosis = st.kurtosis(data_inv.ravel())
            if np.abs(skew) < np.abs(min_skew) and np.abs(kurtosis)<np.abs(min_kurtosis):
                min_skew=skew
                min_kurtosis=kurtosis
                optimal_field = data_inv
            counter += 1

        # Remap normal distribution to truncated dist by rank ordering.
        if self.field_min is not None and self.field_max is not None:
            assert self.field_min < self.field_max, 'field_min {} must be less than field_max {}'.format(
                    self.field_min,self.field_max)
            trunc_field = st.truncnorm.rvs(self.field_min,self.field_max,loc=0.0,scale=1.0,size=(self.nx,self.ny,self.nz))

            sorted_trunc_values = np.sort(trunc_field.ravel())
            sorted_white_args = np.argsort(optimal_field.ravel())

            output_field = np.zeros(sorted_trunc_values.shape)
            output_field[sorted_white_args] = sorted_trunc_values
            output_field = np.reshape(output_field,(self.nx,self.ny,self.nz))

        else:
            output_field = optimal_field

        if return_spectra:
            fourier = fft.fftshift(fft.fftn(output_field))

            x,y,z = fft.fftfreq(self.nx,d=self.domain_size[0]/self.nx),fft.fftfreq(self.ny,d=self.domain_size[1]/self.ny), \
            fft.fftfreq(self.nz,d=self.domain_size[2]/self.nz)
            X,Y,Z = np.meshgrid(fft.fftshift(y),fft.fftshift(x),fft.fftshift(z))
            return output_field,fourier,X,Y,Z
        else:
            return output_field

    def _initialize_field(self):
        """
        Generates a white gaussian field and applies the spectral scaling specified by the generator.

        Returns
        -------
        scaled_field: np.array shape=(nx, ny, nz)
        """

        if self.nx % 2 == 1:
            nx = self.nx+1
        else:
            nx = self.nx
        if self.ny %2 == 1:
            ny = self.ny+1
        else:
            ny = self.ny
        if self.nz % 2 ==1:
            nz = self.nz+1
        else:
            nz=self.nz

        white_field = np.random.normal(loc=0.0,scale=1.0,size=(nx,ny,nz))
        white_spectra = fft.fftshift(fft.fftn(white_field))
        x,y,z = fft.fftfreq(nx,d=self.domain_size[0]/nx),fft.fftfreq(ny,d=self.domain_size[1]/ny), \
        fft.fftfreq(nz,d=self.domain_size[2]/nz)


        X,Y,Z = np.meshgrid(fft.fftshift(x),fft.fftshift(y),fft.fftshift(z),indexing='ij')
        isotropic_wavenumber = np.sqrt(X**2+Y**2+Z**2)
        #np.where(isotropic_wavenumber==0.0)

        #if self.ny % 2 == 1 and self.nx % 2 == 0:
        #    isotropic_wavenumber[self.nx//2,]

        isotropic_wavenumber[nx//2,ny//2,nz//2] = 10**-6
        spectral_scaling = np.abs(isotropic_wavenumber)**(self.beta/2.0)
        spectral_scaling[nx//2,ny//2,nz//2] = 0.0
        scaled_spectra = white_spectra*spectral_scaling

        if self.inner_scale is not None:
            scaled_spectra[np.where(isotropic_wavenumber>1.0/self.inner_scale)] = 0.0

        if self.outer_scale is not None:
        #Power above scale-break is an azimuthally averaged power around the scale break.
            value = np.mean(np.abs(scaled_spectra[np.where(np.abs(isotropic_wavenumber-1.0/self.outer_scale)<0.3)]))
            scaled_spectra[np.where(isotropic_wavenumber<1.0/self.outer_scale)] = value

        scaled_field = np.real(fft.ifftn(fft.fftshift(scaled_spectra)))[:self.nx,:self.ny,:self.nz]
        scaled_field = (scaled_field - scaled_field.mean())/scaled_field.std()

        return scaled_field
