
from setuptools import (find_packages, setup)


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='insilico',
    version='0.0.1',
    description='Insilico filter of molecules based on some properties',
    license='Apache-2.0',
    url='https://github.com/SCM-NV/filterInsilico',
    author=['Felipe Zapata', 'Ivan Infante'],
    author_email='fzapata_esciencecenter.nl',
    keywords='chemistry materials',
    long_description=readme(),
    packages=find_packages(),
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'programming language :: python :: 3',
        'development status :: 4 - Beta',
        'intended audience :: science/research',
        'topic :: scientific/engineering :: chemistry'
    ],
    install_requires=[
        'jsonref', 'jsonschema', 'pandas', 'pubchempy', 'pymongo', 'pyyaml', 'scipy'],
    extras_require={
        'test': ['coverage', 'pytest', 'pytest-cov'],
        'doc': ['sphinx']},
    scripts=[
        'scripts/filter/run_filters.py'],
    include_package_data=True,
    package_data={
        'filterInsilico': ['data/schemas/*json']
    }
)
