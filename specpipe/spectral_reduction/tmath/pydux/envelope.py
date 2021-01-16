def envelope(array, window_size):
    import numpy as np
    win_num=len(array)//window_size
    maxarr=np.zeros(win_num)
    minarr=np.zeros(win_num)
    for i in range(0,win_num):
        maxarr[i]=max(array[i*window_size:(i+1)*window_size])
        minarr[i]=min(array[i*window_size:(i+1)*window_size])
    return maxarr, minarr

