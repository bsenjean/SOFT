import io
from setuptools import setup, find_packages

setup(
    name='SOFT',
    author='Bruno Senjean',
    author_email='bruno.senjean@umontpellier.fr',
    url='',
    description=('Site Occupation Functional Theory'),
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
)
