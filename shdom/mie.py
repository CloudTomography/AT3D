"""
A python wrapper for make_mie_table.f90 by Aviad Levis Technion Inst. of Technology February 2019.
Source Fortran files were created by Frank Evans University of Colorado May 2003.
For source code documentation see: http://nit.colorado.edu/shdom/shdomdoc/makemietable.html

Description taken from make_mie_table.f90:
! Does Mie computations to create a scattering table as a function of
! effective radius for gamma or lognormal size distributions of spherical
! particles.  The particles may be water or ice (in which case the 
! program provides the index of refraction depending on wavelength) or
! "aerosols" (in which case the index of refraction is user specified).
! For water or ice particles the scattering properties may be averaged 
! over the desired spectral range with Planck function weighting.  
! The phase functions in the output scattering table are represented 
! with Legendre series.   The effective radii in the table may be evenly
! or logarithmically spaced.
"""

import core
import numpy as np


class Mie(object):
    """
    Mie scattering for a particle size distribution. 
    Scattering coefficients are averaged over a range of particle radii and wavelengths. 
    The output is a table where for each effective radius the follwoing is computed:
      1. Extinction-cross section per 1 unit of mass content [g/m^3](liquid water content for water clouds)  
      2. Legendre expansion coefficients of the normalized scattering phase function (first coefficient is always 1.0)
      3. Single scattering albedo, unitless in the range [0, 1].

    Parameters
    ----------
    wavelength_band: (float, float)
        (minimum, maximum) wavelength in microns. 
        This defines the spectral band over which to integrate, if both are equal monochrome quantities are computed. 
    particle_type: string
        Options are 'Water' or 'Aerosol'.
    distribution: string
        Particle size-distribution. Options are 'Gamma' or 'Log-normal'. 
        Gamma:
          n(r) = a r^alpha exp(-b*r).
          r - droplet radius.
          a, b, alpha - gamma distribution parameters. 
        Log-normal:
          n(r) = a/r exp( -[ln(r/r0)]^2 / (2*alpha^2) ).
          r0 - logarithmic mode.
          alpha - standard deviation of the log. 
    alpha: float
        Shape parameter for the size distribution. 
        Gamma:
          N = a Gamma(alpha+1)/ b^(alpha+1) - number concentration.
          r_eff = (alpha+3)/b - effective radius.
          v_eff = 1/(alpha+3) - effective variance.
        Log-normal:
          N = sqrt(2*pi)*alpha*a - number concentration. 
          r_eff = r0*exp(2.5*sigma^2) - effective radius.
          v_eff = exp(sigma^2)-1  - effective variance.
    wavelength_averaging: bool
        True - average scattering properties over the wavelength_band.
        False - scattering properties of the central wavelength. 
    wavelength_resolution: float
        The distance between two wavelength samples in the band. Used only if wavelength_averaging is True.
    
    Notes
    -----
    Aerosol particle type not supported yet.   
    """
    
    def __init__(self,
                 wavelength_band,
                 particle_type, 
                 distribution,
                 alpha, 
                 wavelength_averaging=False,
                 wavelength_resolution=0.001):
        
        # Particle type 
        if particle_type == 'Water':
            self._partype = 'W'
            self._rindex = 1.0
            self._pardens = 1.0
        else:
            raise NotImplementedError('Particle type {} not supported'.format(particle_type))
        
        # Size distribution pdf 
        if distribution == 'Gamma':
            self._distflag = 'G'
        elif distribution == 'Log-normal':
            self._distflag = 'L'
        else:
            raise NotImplementedError('Distribution type {} not supported'.format(distibution))
            
        # Averaging over spectral band or monochrome
        self._wavelen1, self._wavelen2 = wavelength_band
        assert self._wavelen1 <= self._wavelen2, 'Minimum wavelength is smaller than maximum'
        avgflag = 'C'
        if self._wavelen1 == self._wavelen2:
            deltawave = -1
        elif wavelength_averaging:
            avgflag = 'A'
            deltawave = wavelength_resolution
            
        self._avgflag = avgflag  
        self._deltawave = deltawave
        self._alpha = alpha
        
        self._wavelencen = core.get_center_wavelen(
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2)
        
    
    def compute_table(self,
                      num_effective_radii,
                      start_effective_radius,
                      end_effective_radius,
                      max_integration_radius): 
        
        self._nretab = num_effective_radii
        self._sretab = start_effective_radius
        self._eretab = end_effective_radius
        self._maxradius = max_integration_radius
        
        # Calculate the maximum size parameter and the max number of Legendre terms
        if self._avgflag == 'A':
            xmax = 2 * np.pi * max_integration_radius / self._wavelen1
        else:
            xmax = 2 * np.pi * max_integration_radius / self._wavelencen
        self._maxleg = int(np.round(2.0 * (xmax + 4.0*xmax**0.3334 + 2.0)))
        
        self._reff, self._extinct, self._ssalb, self._nleg, self._legcoef = core.get_mie_table(
            nretab=self._nretab, 
            maxleg=self._maxleg,
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2, 
            wavelencen=self._wavelencen,
            deltawave=self._deltawave, 
            pardens=self._pardens, 
            sretab=self._sretab, 
            eretab=self._eretab, 
            alpha=self._alpha, 
            maxradius=self._maxradius, 
            rindex=self._rindex, 
            partype=self._partype, 
            avgflag=self._avgflag, 
            distflag=self._distflag)
        
    
    def write_table(self, file_path): 
        core.write_mie_table(
            mietabfile=file_path,
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2,
            deltawave=self._deltawave,
            partype=self._partype,
            pardens=self._pardens, 
            rindex=self._rindex,
            distflag=self._distflag,
            alpha=self._alpha, 
            nretab=self._nretab,
            sretab=self._sretab, 
            eretab=self._eretab,             
            reff=self._reff,
            extinct=self._extinct,
            ssalb=self._ssalb,
            nleg=self._nleg,
            legcoef=self._legcoef,
            maxleg=self._maxleg)
     