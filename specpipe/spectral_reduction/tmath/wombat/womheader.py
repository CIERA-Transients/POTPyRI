def womheader(hop):
    """print header to screen"""
    from tmath.wombat.getch import getch
    if (len(hop[0].header)) < 2:
        print('No header!\n')
        return hop
    for i,_ in enumerate(hop[0].header):
        print(str(hop[0].header.cards[i]))
        if (i > 0) and (i % 19 == 0):
            print('Hit any key to continue')
            any=getch()
    return hop

