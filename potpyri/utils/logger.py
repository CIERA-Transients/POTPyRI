"""Logging utilities for the POTPyRI pipeline.

Provides a colored console and file logger with UTC timestamps for pipeline
steps. Authors: Kerry Paterson, Charlie Kilpatrick.
"""
from potpyri._version import __version__

import logging
import time
import datetime
import os

# Define some colors for different log levels
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

def get_log(log_dir):
    """Create and return a ColoredLogger writing to the given directory.

    Log filename is generated as log_YYYYMMDD_HHMMSS.log (UTC).

    Parameters
    ----------
    log_dir : str
        Directory path for the log file.

    Returns
    -------
    ColoredLogger
        Logger instance with stream and file handlers.
    """
    datestr = datetime.datetime.now(datetime.UTC).strftime('%Y%m%d_%H%M%S')
    base_logname = f'log_{datestr}.log'
    log_filename = os.path.join(log_dir, base_logname)
    log = ColoredLogger(log_filename)

    return(log)

def formatter_message(message, use_color=True):
    """Replace $RESET and $BOLD placeholders with ANSI codes or empty strings.

    Parameters
    ----------
    message : str
        Format string possibly containing $RESET and $BOLD.
    use_color : bool, optional
        If True, insert color sequences; otherwise strip them.
        Default is True.

    Returns
    -------
    str
        Message with placeholders substituted.
    """
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

# Dictionary of colors for different log levels
COLORS = {
    'WARNING': YELLOW,
    'INFO': GREEN,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}

class ColoredFormatter(logging.Formatter):
    """Formatter that optionally colors the levelname in log records."""

    def __init__(self, msg, use_color=True):
        """
        Parameters
        ----------
        msg : str
            Format string for the formatter.
        use_color : bool, optional
            Whether to colorize level names. Default is True.
        """
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        """Format the log record; optionally colorize levelname."""
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            # Change color of levelname only if it is requested
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        elif not self.use_color and RESET_SEQ in levelname:
            # Remove color of levelname if coloring is not requested
            record.levelname = levelname[7:].replace(RESET_SEQ, "")
        else:
            # No change in levelname
            pass
        return logging.Formatter.format(self, record)
    
class ColoredLogger(logging.Logger):
    """Logger that writes to both console (colored) and a file (UTC, no color)."""

    ST_FMT = "[$BOLD%(filename)s::%(lineno)d$RESET] [%(levelname)s]  %(message)s"
    F_FMT = "[$BOLD%(asctime)s::%(filename)s::%(lineno)d$RESET] [%(levelname)s] %(message)s"
    STREAM_FORMAT = formatter_message(ST_FMT, True)
    FILE_FORMAT = formatter_message(F_FMT, False)
    def __init__(self, filename):
        """Create logger with stream (colored) and file (UTC) handlers.

        Parameters
        ----------
        filename : str
            Path to the log file. File is opened in 'w+' mode.
        """
        # Set logging level
        logging.Logger.__init__(self, None, logging.INFO)                

        # Create a stream handler
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(ColoredFormatter(self.STREAM_FORMAT, use_color=True))

        # Now a file handler
        filehandler = logging.FileHandler(filename, mode='w+')
        filehandler.setLevel(logging.INFO)
        file_formatter = ColoredFormatter(self.FILE_FORMAT, use_color=False)
        file_formatter.converter = time.gmtime #convert time in logger to UTC
        filehandler.setFormatter(ColoredFormatter(self.FILE_FORMAT, use_color=False))

        # Add both handlers
        self.addHandler(streamhandler)
        self.addHandler(filehandler)
        return

    def close(self):
        """Close and remove all handlers.

        Call when the logger is no longer needed (e.g., end of test or script)
        to avoid ResourceWarnings from unclosed file handles.
        """
        for handler in self.handlers[:]:
            handler.close()
            self.removeHandler(handler)

    def shutdown(self):
        """Shut down the logging system and flush all handlers."""
        logging.shutdown()
        return
