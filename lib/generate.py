"""
A module containing tools for generating synthetic data.
"""

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
    def __init__(self, nx, ny, nz, beta, domain_size, field_min = None,field_max = None,
                 inner_scale = None, outer_scale = None, seed = None):

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
            larger (greater) than power_spectrum inner_scale {}'.format(outer_scale,inner_scale)

        self.outer_scale = outer_scale
        self.inner_scale = inner_scale

        if seed is not None:
            np.random.seed(seed)
            self.seed = seed
        else:
            self.seed = None

    def _update_args(self,kwargs):
        """
        Updates the attributes of the generator.
        Parameters
        ---------
        kwargs: dictionary
        """
        for name in kwargs:
            setattr(self,name,kwargs[name])

    def generate_field(self,return_spectra = False, **kwargs):
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
