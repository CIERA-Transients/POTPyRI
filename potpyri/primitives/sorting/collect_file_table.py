"""Public entry point to build or reload the pipeline file list table."""
from __future__ import annotations

from .sorting_primitive import FileSortingPrimitive


def collect_file_table(file_list, paths, tel, incl_bad=False, proc=None,
    no_redo=False, log=None):
    """Discover raw files, classify them, and write (or read) the fixed-width file list.

    If ``no_redo`` and ``file_list`` exists, reads the existing table. Otherwise
    globs raw/data/bad, runs :func:`build_raw_file_table`, and writes the list.

    Parameters
    ----------
    file_list : str
        Path to output (or existing) file list table.
    paths : dict
        Paths dict with 'raw', 'data', 'bad' keys.
    tel : Instrument
        Instrument instance.
    incl_bad : bool, optional
        If True, include bad files in list. Default is False.
    proc : str, optional
        Processor/run identifier for raw_format glob.
    no_redo : bool, optional
        If True and file_list exists, read it instead of regenerating.
    log : ColoredLogger, optional
        Logger for progress.

    Returns
    -------
    astropy.table.Table
        File table with Target, Filter, Type, CalType, File, etc.

    Raises
    ------
    SystemExit
        If no files found or no good files after sorting.
    """
    return FileSortingPrimitive(
        incl_bad=incl_bad,
        proc=proc,
        no_redo=no_redo,
    ).apply(input=file_list, paths=paths, tel=tel, log=log)
