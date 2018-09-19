
from setuptools import (find_packages, setup)


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='filterInsilico',
    version='0.0.1',
    description='Insilico filter of molecules based on some properties',
    license='Apache-2.0',
    url='',
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
        'qmflows', 'pymonad', 'jsonref', 'jsonschema'],
    extras_require={
        'test': ['coverage', 'pytest', 'pytest-cov']},
    scripts=[
        'scripts/filter/filter.py'],
    include_package_data=True,
    package_data={
        'filterInsilico': ['data/schemas/*json']
    }
)
