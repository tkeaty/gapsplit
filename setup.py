from setuptools import setup,find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(

    name='gapsplit',
    version='1.0.0',
    description='A sampling algorithm for convex and non-convex metabolic models',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/tkeaty/gapsplit',
    author='Tom Keaty',
    author_email='tkeaty2@gmail.com',
    keywords='sampling biology',
    packages=find_packages(),
    install_requires=['numpy>=1.14.5']

)
