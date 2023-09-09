import logging

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

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
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        elif not self.use_color:
            record.levelname = levelname[7:].replace(RESET_SEQ, "")
        return logging.Formatter.format(self, record)
    
# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    ST_FMT = "[$BOLD%(filename)s::%(lineno)d$RESET] [%(levelname)s]  %(message)s"
    F_FMT = "[$BOLD%(asctime)s::%(filename)s::%(lineno)d$RESET] [%(levelname)s] %(message)s"
    STREAM_FORMAT = formatter_message(ST_FMT, True)
    FILE_FORMAT = formatter_message(F_FMT, False)
    def __init__(self, filename):
        
        #logging.Formatter.converter = time.gmtime #convert time in logger to UCT
        logging.Logger.__init__(self, None, logging.INFO)                

        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(ColoredFormatter(self.STREAM_FORMAT, use_color=True))

        filehandler = logging.FileHandler(filename, mode='w+')
        filehandler.setLevel(logging.INFO)
        filehandler.setFormatter(ColoredFormatter(self.FILE_FORMAT, use_color=False))

        self.addHandler(streamhandler)
        self.addHandler(filehandler)
        return