"""
The main entry point of Advent is the :class:`simulation.Simulation` class.
"""

import logging; logger = logging.getLogger(__name__)

from .simulation import Simulation, ExplSimulation
from .readers import Reader, OpendriftReader, ConstantReader

def install_logs(level = logging.DEBUG):
  """
  Set up fanzy, colorful, logging.

  Args:

    level: Log level
  """
  import coloredlogs
  fmt = '%(asctime)s %(levelname)-7s %(name)s: %(message)s'
  fields = coloredlogs.DEFAULT_FIELD_STYLES
  fields['levelname']['color'] = 'magenta'
  coloredlogs.install(level = level, logger = logging.getLogger(), fmt = fmt, datefmt = None, field_styles = fields)


