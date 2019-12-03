import os, time
import numpy as np
import shdom
from scripts.optimize_extinction_lbfgs import OptimizationScript as ExtinctionOptimizationScript


class OptimizationScript(ExtinctionOptimizationScript):
    """
    Optimize: Micro-physics
    ----------------------
    Estimate micro-physical properties based on multi-spectral radiance/polarization measurements.
    Note that for convergence a fine enough sampling of effective radii and variances should be pre-computed in the
    Mie tables used by the forward model. This is due to the linearization of the phase-function and it's derivatives.

    Measurements are simulated measurements using a forward rendering script
    (e.g. scripts/render_radiance_toa.py).

    For example usage see the README.md

    For information about the command line flags see:
      python scripts/optimize_microphysics_lbfgs.py --help

    Parameters
    ----------
    scatterer_name: str
        The name of the scatterer that will be optimized.
    """
    def __init__(self, scatterer_name='cloud'):
        super().__init__(scatterer_name)

    def medium_args(self, parser):
        """
        Add common medium arguments that may be shared across scripts.

        Parameters
        ----------
        parser: argparse.ArgumentParser()
            parser initialized with basic arguments that are common to most rendering scripts.

        Returns
        -------
        parser: argparse.ArgumentParser()
            parser initialized with basic arguments that are common to most rendering scripts.
        """
        parser.add_argument('--use_forward_lwc',
                            action='store_true',
                            help='Use the ground-truth LWC.')
        parser.add_argument('--use_forward_reff',
                                action='store_true',
                                help='Use the ground-truth effective radius.')
        parser.add_argument('--use_forward_veff',
                            action='store_true',
                            help='Use the ground-truth effective variance.')
        parser.add_argument('--const_lwc',
                            action='store_true',
                            help='Keep liquid water content constant at a specified value (not optimized).')
        parser.add_argument('--const_reff',
                            action='store_true',
                            help='Keep effective radius constant at a specified value (not optimized).')
        parser.add_argument('--const_veff',
                            action='store_true',
                            help='Keep effective variance constant at a specified value (not optimized).')
        parser.add_argument('--radiance_threshold',
                            default=[0.0175],
                            nargs='+',
                            type=np.float32,
                            help='(default value: %(default)s) Threshold for the radiance to create a cloud mask.'
                            'Threshold is either a scalar or a list of length of measurements.')
        parser.add_argument('--lwc_scaling',
                            default=10.0,
                            type=np.float32,
                            help='(default value: %(default)s) Pre-conditioning scale factor for liquid water content estimation')
        parser.add_argument('--reff_scaling',
                            default=1e-1,
                            type=np.float32,
                            help='(default value: %(default)s) Pre-conditioning scale factor for effective radius estimation')
        parser.add_argument('--veff_scaling',
                            default=1.0,
                            type=np.float32,
                            help='(default value: %(default)s) Pre-conditioning scale factor for effective variance estimation')
        return parser

    def get_medium_estimator(self, measurements, ground_truth):
        """
        Generate the medium estimator for optimization.

        Parameters
        ----------
        measurements: shdom.Measurements
            The acquired measurements.
        ground_truth: shdom.Scatterer


        Returns
        -------
        medium_estimator: shdom.MediumEstimator
            A medium estimator object which defines the optimized parameters.
        """
        # Define the grid for reconstruction
        if self.args.use_forward_grid:
            lwc_grid = ground_truth.lwc.grid
            reff_grid = ground_truth.reff.grid
            veff_grid = ground_truth.reff.grid
        else:
            lwc_grid = reff_grid = veff_grid = self.cloud_generator.get_grid()
        grid = lwc_grid + reff_grid + veff_grid

        # Find a cloud mask for non-cloudy grid points
        if self.args.use_forward_mask:
            mask = ground_truth.get_mask(threshold=0.01)
        else:
            carver = shdom.SpaceCarver(measurements)
            mask = carver.carve(grid, agreement=0.9, thresholds=self.args.radiance_threshold)

        # Define micro-physical parameters: either optimize, keep constant at a specified value or use ground-truth
        if self.args.use_forward_lwc:
            lwc = ground_truth.lwc
        elif self.args.const_lwc:
            lwc = self.cloud_generator.get_lwc(lwc_grid)
        else:
            lwc = shdom.GridDataEstimator(self.cloud_generator.get_lwc(lwc_grid),
                                          min_bound=1e-5,
                                          max_bound=2.0,
                                          precondition_scale_factor=self.args.lwc_scaling)
        lwc.apply_mask(mask)

        if self.args.use_forward_reff:
            reff = ground_truth.reff
        elif self.args.const_reff:
            reff = self.cloud_generator.get_reff(reff_grid)
        else:
            reff = shdom.GridDataEstimator(self.cloud_generator.get_reff(reff_grid),
                                           min_bound=ground_truth.min_reff,
                                           max_bound=ground_truth.max_reff,
                                           precondition_scale_factor=self.args.reff_scaling)
        reff.apply_mask(mask)

        if self.args.use_forward_veff:
            veff = ground_truth.veff
        elif self.args.const_veff:
            veff = self.cloud_generator.get_veff(veff_grid)
        else:
            veff = shdom.GridDataEstimator(self.cloud_generator.get_veff(veff_grid),
                                           max_bound=ground_truth.max_veff,
                                           min_bound=ground_truth.min_veff,
                                           precondition_scale_factor=self.args.veff_scaling)
        veff.apply_mask(mask)

        # Define a MicrophysicalScattererEstimator object
        cloud_estimator = shdom.MicrophysicalScattererEstimator(ground_truth.mie, lwc, reff, veff)
        cloud_estimator.set_mask(mask)

        # Create a medium estimator object (optional Rayleigh scattering)
        medium_estimator = shdom.MediumEstimator(
            loss_type=self.args.loss_type,
            stokes_weights=self.args.stokes_weights
        )
        if self.args.add_rayleigh:
            air = self.air_generator.get_scatterer(cloud_estimator.wavelength)
            medium_estimator.set_grid(cloud_estimator.grid + air.grid)
            medium_estimator.add_scatterer(air, 'air')
        else:
            medium_estimator.set_grid(cloud_estimator.grid)

        medium_estimator.add_scatterer(cloud_estimator, self.scatterer_name)
        return medium_estimator

    def load_forward_model(self, input_directory):
        """
        Load the ground-truth medium, rte_solver and measurements which define the forward model

        Parameters
        ----------
        input_directory: str
            The input directory where the forward model is saved

        Returns
        -------
        ground_truth: shdom.OpticalScatterer
            The ground truth scatterer
        rte_solver: shdom.RteSolverArray
            The rte solver with the numerical and scene parameters
        measurements: shdom.Measurements
            The acquired measurements
        """
        # Load forward model and measurements
        medium, rte_solver, measurements = shdom.load_forward_model(input_directory)

        # Get micro-physical medium ground-truth
        ground_truth = medium.get_scatterer(self.scatterer_name)
        return ground_truth, rte_solver, measurements

    def get_summary_writer(self, measurements, ground_truth):
        """
        Define a SummaryWriter object

        Parameters
        ----------
        measurements: shdom.Measurements object
            The acquired measurements.
        ground_truth: shdom.Scatterer
            The ground-truth scatterer for monitoring

        Returns
        -------
        writer: shdom.SummaryWriter object
            A logger for the TensorboardX.
        """
        writer = None
        if self.args.log is not None:
            log_dir = os.path.join(self.args.input_dir, 'logs', self.args.log + '-' + time.strftime("%d-%b-%Y-%H:%M:%S"))
            writer = shdom.SummaryWriter(log_dir)
            writer.save_checkpoints(ckpt_period=20 * 60)
            writer.monitor_loss()

            # Compare estimator to ground-truth
            writer.monitor_scatterer_error(estimator_name=self.scatterer_name, ground_truth=ground_truth)
            writer.monitor_scatter_plot(estimator_name=self.scatterer_name, ground_truth=ground_truth, dilute_percent=0.4, parameters=['lwc'])
            writer.monitor_horizontal_mean(estimator_name=self.scatterer_name, ground_truth=ground_truth, ground_truth_mask=ground_truth.get_mask(threshold=0.01))

        return writer


if __name__ == "__main__":
    script = OptimizationScript(scatterer_name='cloud')
    script.main()



