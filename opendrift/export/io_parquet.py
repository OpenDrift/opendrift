import sys
import os
from datetime import datetime, timedelta
import logging

logging.captureWarnings(True)
logger = logging.getLogger(__name__)
import string
from shutil import move

import numpy as np
import pandas as pd
from opendrift.models.basemodel import Mode

def init(self, filename):

    self.outfile = filename
    dummy_data = {
        k: pd.Series([], dtype=t) for k, (t, _) in self.history.dtype.fields.items()
    }
    dummy_data["time"] = pd.Series([], dtype="datetime64[ns]")
    df = pd.DataFrame(dummy_data)
    df.to_parquet(self.outfile, engine="fastparquet")


def write_buffer(self):
    num_steps_to_export = self.steps_output - self.steps_exported

    data = {
        k: self.history[k][:, 0:num_steps_to_export][
            ~self.history[k].mask[:, 0:num_steps_to_export]
        ]  # automatically flattens array
        for k in self.history.dtype.fields
    }

    times = [
        self.start_time + n * self.time_step_output
        for n in range(self.steps_exported, self.steps_output)
    ]

    _arr_template = self.history["ID"][:, 0:num_steps_to_export]
    time_arr = np.repeat([times], _arr_template.shape[0], axis=0)
    data["time"] = time_arr[~_arr_template.mask]  # automatically flattens array

    df = pd.DataFrame(data)
    df.to_parquet(self.outfile, engine="fastparquet", append=True)

    logger.info("Wrote %s steps to file %s" % (num_steps_to_export, self.outfile))
    self.history.mask = True  # Reset history array, for new data
    self.steps_exported = self.steps_exported + num_steps_to_export


def close(self):
    logger.warning("`.close` not strictly needed...?")

def import_file(self, filename, times=None, elements=None, load_history=True):
    """Create OpenDrift object from imported file.
    This reimport is potentially very costly anyway 
    """
    logger.info("Skipping reimport")
    return self

def import_file_xarray(self, filename, times=None, elements=None, load_history=True):
    """Create OpenDrift object from file
    Odd if this I/O backend specific feature were required for all of opendrift to run
    """
    raise NotImplementedError("wontfix")
