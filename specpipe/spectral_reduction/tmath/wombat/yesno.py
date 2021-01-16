def yesno(default):
    from tmath.wombat.getch import getch
    answer=''
    while (answer != 'y') and (answer != 'n'):
        print('(y)es or (n)o? (default/return = {}) '.format(default))
        reply=getch()
        reply=reply.strip()
        if len(reply) == 0:
            reply=default
        answer=reply.lower()[0]
    return answer

