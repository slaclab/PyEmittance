from setuptools import setup, find_packages
from os import path, environ

cur_dir = path.abspath(path.dirname(__file__))

with open(path.join(cur_dir, 'requirements.txt'), 'r') as f:
    requirements = f.read().split()


setup(
    name='pyemittance',
    version='v0.1.1',
    author='Sara Miskovich',
    author_email='smiskov@slac.stanford.edu',
    packages=find_packages(where='pyemittance'),

    package_dir={'':'pyemittance'},
    url='https://github.com/slaclab/PyEmittance',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=requirements,
    include_package_data=True,
    python_requires='>=3.5'
)