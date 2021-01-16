def scaleparse():
    import re
    done=False
    while (not done):
        scaleparse=input('Enter new ymin and ymax values: ')
        if (scaleparse != ''):
            scaleparse=scaleparse.strip()
            scalesplit=re.split('\W+',scaleparse)
            if (len(scalesplit) > 1):
                try:
                    ymin=float(scalesplit[0])
                    ymax=float(scalesplit[1])
                except ValueError:
                    print('Please enter numbers, ### ###,\n')
                    print('###,###,or ###-###\n')
                    print('You entered {}\n'.format(scaleparse))
                else:
                    done=True
            else:
                print('Enter more than one number\n')
    return  ymin, ymax
