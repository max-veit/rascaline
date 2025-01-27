# -*- coding: utf-8 -*-
import json
import ctypes

from .clib import _get_library
from .utils import _call_with_growing_buffer


class Profiler:
    """
    Rascaline uses the `time_graph <https://docs.rs/time-graph/>`_ to
    collect timing information on the calculations. The ``Profiler`` context
    manager provides access to this functionality.

    The profiling code collects the total time spent inside the most important
    functions, as well as the function call graph (which function called which
    other function).

    .. code-block:: python

        import rascaline

        with rascaline.Profiler() as profiler:
            # run some calculations

        print(profiler.as_short_table())
    """

    def __init__(self):
        self._lib = _get_library()

    def __enter__(self):
        self._lib.rascal_profiling_enable(True)
        self._lib.rascal_profiling_clear()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._lib.rascal_profiling_enable(False)

    def as_json(self):
        """Get the current profiling data as a JSON string"""
        return _call_with_growing_buffer(
            lambda b, s: self._lib.rascal_profiling_get("json".encode("utf8"), b, s)
        )

    def as_table(self):
        """Get the current profiling data as a table to be displayed to the user"""
        return _call_with_growing_buffer(
            lambda b, s: self._lib.rascal_profiling_get("table".encode("utf8"), b, s)
        )

    def as_short_table(self):
        """
        Get the current profiling data as a table with short names to be
        displayed to the user
        """
        return _call_with_growing_buffer(
            lambda b, s: self._lib.rascal_profiling_get(
                "short_table".encode("utf8"), b, s
            )
        )
