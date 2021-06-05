"""
This module contains summary writer that can be used for monitoring
the progress of long optimizations.
"""

import tensorboardX

class SummaryWriter(tensorboardX.SummaryWriter):
    """
    A wrapper around the tensorboardX.SummaryWriter where any default
    figures/scalars to monitor can be added to avoid defining them
    at the script level.
    """
    def __init__(self, logdir=None, comment='', purge_step=None, max_queue=10, flush_secs=120,
                 filename_suffix='', write_to_disk=True, **kwargs):
        super().__init__(logdir, comment, purge_step, max_queue, flush_secs,
                         filename_suffix, write_to_disk, **kwargs)

    #TODO add some defaults/examples.
