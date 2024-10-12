import glob
import os
from setuptools import setup
from setuptools import find_packages

# Begin setup
setup_keywords = dict()
setup_keywords['name'] = 'POTPyRI'
setup_keywords['description'] = 'Pipeline for Optical/Infrared Telescopes in Python for Reducing Images'
setup_keywords['author'] = 'Kerry Paterson, Charlie Kilpatrick, et al.'
setup_keywords['author_email'] = 'ckilpatrick@northwestern.edu'
setup_keywords['license'] = 'GNU General Public License v3.0'
setup_keywords['url'] = 'https://github.com/CIERA-Transients/POTPyRI'
setup_keywords['version'] = '2.0'
# Use README.rst as long_description.
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>3.11.0)']

if os.path.exists('requirements.txt'):
    with open('requirements.txt') as f:
        required = f.read().splitlines()
        setup_keywords['install_requires'] = required

setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = False
setup_keywords['packages'] = find_packages()
setup_keywords['package_data'] = {'potpyri': ['data/staticmasks/*.fz', 
    'config/*', 'data/cal/*/*.fz']}
setup_keywords['include_package_data'] = True
setup_keywords['setup_requires'] = ['pytest-runner']
setup_keywords['tests_require'] = ['pytest']

if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
                                 if not os.path.basename(fname).endswith('.rst')]

setup(**setup_keywords)
