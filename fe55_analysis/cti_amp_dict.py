"""
Dictionary of [ps]cti values for each amp.
"""
import numpy as np

__all__ = ['pcti_', 'scti']

pcti_ = None
log_cti_min = np.log10(1e-6)
log_cti_max = np.log10(1e-3)
scti_ = {1: 0}
scti_.update(dict([(i + 2, cti) for i, cti in
                   enumerate(np.logspace(log_cti_min, log_cti_max, 15))]))
