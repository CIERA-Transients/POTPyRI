def wshow():
    import matplotlib.pyplot as plt
#    plt.ioff()
    wm=plt.get_current_fig_manager()
    blah=wm.window.attributes('-topmost',1)
    blah=wm.window.attributes('-topmost',0)
#    plt.ion()
    return

