import requests,math,sys,time,datetime
from astroquery.mast import Observations
from datetime import datetime
from astroquery.ned import Ned
from astroquery.simbad import Simbad
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from bs4 import BeautifulSoup
from astropy.time import Time
from datetime import date,timedelta
from astropy.table import Table
import astropy.units as u

class metacheck():
    def __init__(self,ra,dec,mjd):
        # Transient information/search criteria
        self.coord = SkyCoord(ra,dec,unit=(u.hour,u.deg), frame='icrs')
        self.time = Time(mjd,format='mjd').datetime

        # Parameters for searches
        self.tns_radius = 5. # arcseconds
        self.mpc_radius = 1. # arcminutes
        self.ast_radius = 5. # arcseconds
        self.sim_radius = 10 # arcseconds
        self.ned_radius = 3. # arcminutes
        self.mdw_radius = 10 # arcseconds

    # Check the minor planet center's database of asteroids
    def checkMPC(self):
        table_output = []

        # Hack AF
        url = "https://minorplanetcenter.net/cgi-bin/mpcheck.cgi"
        form_data = {"year":"%s" % self.time.strftime("%Y"),
            "month":"%s" % self.time.strftime("%m"),
            "day":"%s" % self.time.strftime("%d"),
            "which":"pos",
            "ra":"%s" %
                self.coord.ra.to_string(unit="hour",pad=True,decimal=False,sep=" "),
            "decl":"%s" %
                self.coord.dec.to_string(unit="deg",pad=True,decimal=False,sep=" "),
            "TextArea":"",
            "radius":"%s" % 10,
            "limit":"%s" % 24,
            "oc":"%s" % 500,
            "sort":"d",
            "mot":"h",
            "tmot":"s",
            "pdes":"u",
            "needed":"f",
            "ps":"n",
            "type":"p"}

        r = requests.post(url, data=form_data)
        soup = BeautifulSoup(r.text,'html5lib')
        pre = soup.find("pre")
        if pre is None:
            return(None)
        else:
            data = []
            for row in [row for row in pre.contents[-1].split("\n")[3:-1]]:
                data.append([row[9:25],row[25:36],row[36:47]])
            table = Table(list(map(list, zip(*data))),names=('name','ra','dec'))
            return(table)

    # Check TNS
    def checkTNS(self):
        url = 'https://wis-tns.weizmann.ac.il/search?'
        url += '&ra={ra}'.format(ra=self.coord.ra.degree)
        url += '&decl={dec}'.format(dec=self.coord.dec.degree)
        url += '&radius={rad}'.format(rad=self.tns_radius)
        url += '&coords_unit=arcsec'

        # Get page info from url
        r = requests.get(url)
        soup = BeautifulSoup(r.text,'html5lib')
        html_table = soup.find("tbody")
        if (html_table) is None:
            return(None)
        tables = html_table.find_all("tr")
        data = []
        for table in tables:
            object_list = table.find_all("td")
            new_row = []
            for element in object_list:
                if ('cell-name' in element['class'] or
                    'cell-ra' in element['class'] or
                    'cell-decl' in element['class'] or
                    'cell-discoverydate' in element['class']):
                    new_row.append(element.get_text())
            data.append(new_row)

        table = Table(list(map(list, zip(*data))),names=('name','ra','dec','discoverydate'))
        return(table)

    # Check OSC
    def checkASTROCATS(self):
        # Get photometry in a 10 arcsecond region around coordinates
        url = 'https://api.astrocats.space/catalog/photometry/time+band+magnitude?'
        url += 'ra={ra}'.format(ra=self.coord.ra.degree)
        url += '&dec={dec}'.format(dec=self.coord.dec.degree)
        url += '&radius={radius}'.format(radius=self.ast_radius)
        url += '&format=csv'

        r = requests.get(url)
        if ('No objects found' in r.text):
            return None
        else:
            data = [l.split(',')
                for l in r.text.strip("event,time,band,magnitude\n").split("\n")]
            table = Table(list(map(list, zip(*data))),
                names=('event','time','band','magnitude'))
        return(table)

    # Check Simbad
    def checkSimbad(self):
        table = Simbad.query_region(self.coord, radius = self.sim_radius * u.arcsec)
        return(table)

    # Check NED
    def checkNED(self):
        table = Ned.query_region(self.coord,
            radius=self.ned_radius * u.arcmin, equinox='J2000.0')

        # Pick only galaxies
        mask = table['Type'] == b'G'
        return(table[mask])

    # Check M dwarf flare catalog
    def checkMDWARF(self):
        try:
            table = heasarc.query_region(self.coord,
                mission='mdwarfasc', radius=self.mdw_radius*u.arcsec)
            return(table)
        except:
            return(None)

if __name__ == '__main__':

    start = time.time()
    ra_list,dec_list = np.loadtxt('test.txt',unpack=True,dtype=str)
    for ra,dec in zip(ra_list,dec_list):
        check = metacheck(ra,dec,57410.1231)
        mpc_table = check.checkMPC()
        if (mpc_table is None):
            length = 0
        else:
            length = len(mpc_table)
        print('Done with MPC, took ',time.time()-start,
            'seconds so far and found',length,'sources')
        tns_table = check.checkTNS()
        if (tns_table is None):
            length = 0
        else:
            length = len(tns_table)
        print('Done with TNS, took ',time.time()-start,
            'seconds so far and found',length,'sources')
        ast_table = check.checkASTROCATS()
        if (ast_table is None):
            length = 0
        else:
            length = len(ast_table)
        print('Done with ASTROCATS, took ',time.time()-start,
            'seconds so far and found',length,'records')
        sim_table = check.checkSimbad()
        if (sim_table is None):
            length = 0
        else:
            length = len(sim_table)
        print('Done with Simbad, took ',time.time()-start,
            'seconds so far and found',length,'sources')
        ned_table = check.checkNED()
        if (ned_table is None):
            length = 0
        else:
            length = len(ned_table)
        print('Done with NED, took ',time.time()-start,
            'seconds so far and found',length,'sources')
        mdw_table = check.checkMDWARF()
        if (mdw_table is None):
            length = 0
        else:
            length = len(mdw_table)
        print('Done with MDWARF, took ',time.time()-start,
            'seconds so far and found',length,'sources')

