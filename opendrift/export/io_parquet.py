import os
import logging

import numpy as np
import pandas as pd
import xarray as xr

logging.captureWarnings(True)
logger = logging.getLogger(__name__)


def init(self, filename):
    """Initialise parquet output."""
    self.outfile_name = filename

    if os.path.exists(self.outfile_name):
        logger.warning(f'Deleting existing parquet file {self.outfile_name}')
        os.remove(self.outfile_name)

    self._parquet_initialized = False


def _result_chunk_to_dataframe(self):
    """Flatten current self.result (xarray.Dataset) to tidy DataFrame."""
    ds = self.result

    if "time" not in ds.dims or "trajectory" not in ds.dims:
        raise ValueError(
            "Expected self.result to have 'time' and 'trajectory' dimensions "
            f"but got dims: {list(ds.dims)}"
        )

    df = ds.to_dataframe().reset_index()

    # Drop unseeded/invalid rows (simple proxy: lon is NaN)
    if "lon" in df.columns:
        df = df[~df["lon"].isna()]

    df["time"] = pd.to_datetime(df["time"])

    return df


def write_buffer(self):
    """Append current buffer (self.result) to parquet file."""
    if self.result is None or self.result.sizes.get("time", 0) == 0:
        logger.debug("No timesteps in result; nothing to write to parquet.")
        return

    df = _result_chunk_to_dataframe(self)
    if df.empty:
        logger.debug("Flattened DataFrame is empty; skipping parquet write.")
        return

    if not getattr(self, "_parquet_initialized", False):
        df.to_parquet(self.outfile_name, engine="fastparquet")
        self._parquet_initialized = True
        logger.info(
            "Initialized parquet file %s with %d rows, %d timesteps",
            self.outfile_name,
            len(df),
            self.result.sizes.get("time", 0),
        )
    else:
        df.to_parquet(self.outfile_name, engine="fastparquet", append=True)
        logger.info(
            "Appended %d rows (%d timesteps) to %s",
            len(df),
            self.result.sizes.get("time", 0),
            self.outfile_name,
        )


def close(self):
    """Finalize parquet output and reconstruct self.result as xarray.Dataset."""
    if not getattr(self, "_parquet_initialized", False):
        logger.debug("Parquet file was never written; nothing to close.")
        return

    logger.debug(f"Reopening parquet file {self.outfile_name} to build result")

    df = pd.read_parquet(self.outfile_name, engine="fastparquet")

    # Ensure we have trajectory and time as index
    if not {"trajectory", "time"}.issubset(df.columns):
        raise ValueError(
            "Parquet file must contain 'trajectory' and 'time' columns "
            f"but columns are: {list(df.columns)}"
        )

    df = df.set_index(["trajectory", "time"]).sort_index()

    # Convert back to xarray.Dataset
    ds = df.to_xarray()

    # Make sure time coord is datetime64[ns]
    ds = ds.assign_coords(time=("time", pd.to_datetime(ds.coords["time"].values)))

    self.result = ds
    logger.debug("Parquet close: self.result has been rebuilt from parquet.")


def import_file(self, filename, times=None, elements=None, load_history=True):
    """Parquet import placeholder."""
    logger.info("Parquet import is not implemented; returning self unchanged.")
    return self
