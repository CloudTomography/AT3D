import shdom
from scripts.render_radiance_toa import RenderScript as RadianceRenderScript


class RenderScript(RadianceRenderScript):
    """
    Render: Polarization at Top of the Atmosphere (TOA)
    ---------------------------------------------------
    Forward rendering of an atmospheric medium with at multiple spectral bands with an orthographic sensor measuring
    exiting stokes vector at the top of the the domain. This sensor is an approximation (somewhat crude) for far observing satellites where the rays are parallel.

    As with all `render` scripts a `generator` needs to be specified using the generator flag (--generator).
    The Generator defines the medium parameters: Grid, Liquid Water Content, Effective Radius and Variance with it's own set of command-line flags.

    For example usage see the README.md

    For information about the command line flags see:
      python scripts/render/render_polarization_toa.py --help

    For a tutorial overview of how to operate the forward rendering model see the following notebooks:
     - notebooks/Make Mie Table Polarized.ipynb
     - notebooks/Polarization Rendering [Pure Rayleigh].ipynb
     - notebooks/Polarization Rendering [Pure Mie].ipynb
     - notebooks/Polarization Rendering [Stcu].ipynb
     - notebooks/Polarization Rendering [Cloudbow].ipynb

    Parameters
    ----------
    num_stokes: int
        The number of stokes components. If num_stokes > 1 pyshdom should be compiled with polarization flag.
    scatterer_name: str
        The name of the scatterer that will be optimized.
    """
    def __init__(self, num_stokes=3, scatterer_name='cloud'):
        self.scatterer_name = scatterer_name
        self.num_stokes = num_stokes
        self.sensor = shdom.StokesSensor()

    def rendering_args(self, parser):
        """
        Add common rendering arguments that may be shared across scripts.

        Parameters
        ----------
        parser: argparse.ArgumentParser()
            parser initialized with basic arguments that are common to most rendering scripts.

        Returns
        -------
        parser: argparse.ArgumentParser()
            parser initialized with basic arguments that are common to most rendering scripts.
        """
        parser = super().rendering_args(parser)
        parser.add_argument('--num_stokes',
                            default=3,
                            type=int,
                            help='(default value: %(default)s) The number of stokes components in [1,4]')

        for action in parser._actions:
            if action.dest == 'mie_base_path':
                action.default = 'mie_tables/polydisperse/Water_<wavelength>nm.scatpol'
                action.help = '(default value: %(default)s) Polarized Mie table base file name.'\
                              '<wavelength> will be replaced by the corresponding wavelengths.'
        return parser


if __name__ == "__main__":
    script = RenderScript(num_stokes=3, scatterer_name='cloud')
    script.main()

