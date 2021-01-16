def inputter_single_mix(statement,ans_string):
    reply=' '
    while (reply not in ans_string):
        reply=input(statement)
        reply=reply.strip()
        reply=reply[0]
    return reply

                
