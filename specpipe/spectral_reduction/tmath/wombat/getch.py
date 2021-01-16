def getch():
    """
    get single character
    adapted from https://stackoverflow.com/questions/510357/python-read-a-single-character-from-the-user
    user Louis
    """
    import termios
    import sys, tty
    def _getch():
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch
    return _getch()

