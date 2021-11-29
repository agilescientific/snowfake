from setuptools import setup

CLASSIFIERS = ['Development Status :: 3 - Alpha',
               'Intended Audience :: Science/Research',
               'Natural Language :: English',
               'License :: OSI Approved :: Apache Software License',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Programming Language :: Python :: 3.9',
               'Programming Language :: Python :: 3.10',
              ]

setup(
    name='snowfake',
    version='0.1.0',
    url='https://github.com/agile-geoscience/snowfake',
    author='Matt Hall',
    author_email='matt@agilescientific.com',
    description='Simulate snowflakes with Gravner-Griffeath\'s algorithm.',
    packages=['snowfake'],    
    install_requires=['numpy', 'matplotlib', 'scipy', 'tqdm'],
)
