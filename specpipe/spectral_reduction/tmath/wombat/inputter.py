def inputter(statement, datatype, defaultflag,*default):
    done=False
    while (not done):
        answer=input(statement)
        if (defaultflag) and (answer == ''):
            return default[0]
        elif (datatype == 'string'):
            return answer
        elif (datatype == 'float'):
            try:
                answer=float(answer)
            except ValueError:
                print('Please enter a number. \n')
            else:
                return answer
                done=True
        elif (datatype == 'int'):
            try:
                answer=int(answer)
            except ValueError:
                print('Please enter a number. \n')
            else:
                return answer
                done=True


