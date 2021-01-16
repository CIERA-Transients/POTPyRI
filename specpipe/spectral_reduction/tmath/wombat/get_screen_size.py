def get_screen_size():
    """
    from tmath.scr import get_monitors
    mon=str(get_monitors())
    index1=mon.find('(')
    index2=mon.find('x')
    index3=mon.find('+')
    xsize=int(mon[index1+1:index2])
    ysize=int(mon[index2+1:index3])
    """
    import tkinter as tk
    root = tk.Tk()
    xsize = root.winfo_screenwidth()
    ysize = root.winfo_screenheight()
    root.withdraw()
    root.destroy()
    return xsize, ysize
