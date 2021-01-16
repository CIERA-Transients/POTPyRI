def inputter_single(statement,ans_string):
    from tmath.wombat.getch import getch
    reply=' '
    while (reply not in ans_string):
        print(statement)
        reply=getch()
        reply=reply.strip()
        reply=reply.lower()
        if len(reply) > 0:
            reply=reply[0]
        else:
            reply=' '
    return reply

