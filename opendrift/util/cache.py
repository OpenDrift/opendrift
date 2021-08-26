# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2021, Gaute Hope, MET Norway

from functools import partial
import time
import copy
from typing import Optional
from pathlib import Path
import pickle
import os
import xdg
import logging

logger = logging.getLogger(__name__)

def opendrift_cache():
    """
    Return the platform cache directory.
    """
    c = xdg.xdg_cache_home() / "opendrift"
    if not os.path.exists(c):
        os.makedirs(c)

    return c

def file_cache(path, timeout):
    return partial(FileCache, path=path, timeout=timeout)


class FileCache:
    """
    A decorator for caching the output of functions without arguments to file.
    """

    f = None
    result = None
    path: Path

    timeout: Optional[float]
    cache_time: float = 0.0

    def __init__(self, function, path, timeout=None):
        self.f = function
        self.timeout = timeout
        self.path = path

        self.__load_result__()

    def __update_cache__(self):
        self.cache_time = time.time()
        self.result = self.f()

        self.__store_result__()

    def __store_result__(self):
        logging.debug(f'storing result to: {self.path}')
        with open(self.path, 'wb') as fd:
            pickle.dump(self.result, fd)

    def __load_result__(self):
        if os.path.exists(self.path):
            self.cache_time = os.path.getmtime(self.path)
            logging.debug(f'loading result from: {self.path}, mtime: {self.cache_time}')
            with open(self.path, 'rb') as fd:
                self.result = pickle.load(fd)
        else:
            logging.debug(f'{self.path} does not exist, not loading.')

    def __call__(self):
        if self.result is None:
            self.__update_cache__()

        if self.timeout is not None and (time.time() - self.cache_time >
                                         self.timeout):
            self.__update_cache__()

        return copy.deepcopy(self.result)
