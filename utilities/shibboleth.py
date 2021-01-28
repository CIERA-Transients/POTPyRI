"""
    By CDK

    A class for storing username and password data grabbed from a user/group
    read-only file (440).  Usage example:

    shibboleth = shibboleth('/path/to/read/only/file')
    data = {'username': shibboleth.get_username('lcogt'),
            'password': shibboleth.get_password('lcogt')}
    r = requests.post('https://observe.lco.global/api-token-auth/',data=data)
"""

import os, sys
from astropy.io import ascii

class shibboleth(object):
    def __init__(self, shibboleth):

        self.shibboleth = None
        self.filename = shibboleth

        # Verify that only the owner can read this file
        st = os.stat(shibboleth)
        permissions = oct(st.st_mode)[-3:]
        if permissions!='440':
            message = 'ERROR: we only accept owner read-only permissions!'
            print(message)
            sys.exit(1)
        else:
            self.shibboleth = ascii.read(shibboleth,
                names=['type','username','password'])

    def get_shibboleth_row(self, stype):
        mask = self.shibboleth['type']==stype

        if len(self.shibboleth[mask])!=1:
            message = 'ERROR: missing or ambiguous shibboleth type={0} \n'
            message += 'Correct your shibboleth file={1}'
            print(message.format(stype, self.filename))
            return(None)

        return(self.shibboleth[mask][0])

    def get_username(self, stype):
        row = self.get_shibboleth_row(stype)
        if row:
            return(row['username'])
        else:
            return(None)

    def get_password(self, stype):
        row = self.get_shibboleth_row(stype)
        if row:
            return(row['password'])
        else:
            return(None)

