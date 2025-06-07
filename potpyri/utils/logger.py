"Function for logging throughout the pipeline."
"Authors: Kerry Paterson, Charlie Kilpatrick"

# Initial version tracking on 09/21/2024
__version__ = "1.0"

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

    datestr = datetime.datetime.now(datetime.UTC).strftime('%Y%m%d_%H%M%S')
    base_logname = f'log_{datestr}.log'
    log_filename = os.path.join(log_dir, base_logname)
    log = ColoredLogger(log_filename)

    return(log)

def logprint(log, message, method='info'):
    if log: 
        if method=='info':
            log.info(message)
        elif method=='error':
            log.error(message)
        elif method=='critical':
            log.critical(message)
        else:
            raise Exception(f'Unknown method {method}')

def formatter_message(message, use_color = True):
    # Bold face some text only when rich text/color is requested.
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
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
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
    
# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    # Define formats for stream and file logging
    ST_FMT = "[$BOLD%(filename)s::%(lineno)d$RESET] [%(levelname)s]  %(message)s"
    F_FMT = "[$BOLD%(asctime)s::%(filename)s::%(lineno)d$RESET] [%(levelname)s] %(message)s"
    STREAM_FORMAT = formatter_message(ST_FMT, True)
    FILE_FORMAT = formatter_message(F_FMT, False)
    # initialize logger
    def __init__(self, filename):
        
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
