def womrebin_tform(rx,xin,arctwo,darcone,in21,acc,rslope,in11,in12,arcone):
    """needed for ashrebin"""
    from tmath.wombat.womply import womply
    n=1
    x=xin
    rl=womply(arctwo,in21,rx)
    while True:
        l=womply(arcone,in11,x)
        dx=(rl-l)*rslope
        x=x+dx
        if (abs(dx) < acc):
            break
        if (n > 100):
            break
        n=n+1
        rslope=1./(womply(darcone,in12,x))
    return x

