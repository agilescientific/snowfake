from setuptools import setup
from pathlib import Path

CLASSIFIERS = ['Development Status :: 3 - Alpha',
               'Intended Audience :: Science/Research',
               'Natural Language :: English',
               'License :: OSI Approved :: Apache Software License',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Programming Language :: Python :: 3.8',
               'Programming Language :: Python :: 3.9',
               'Programming Language :: Python :: 3.10',
              ]

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='snowfake',
    version='0.2.1',
    url='https://github.com/agile-geoscience/snowfake',
    author='Matt Hall',
    author_email='matt@agilescientific.com',
    description='Simulate snowflakes with Gravner-Griffeath\'s algorithm.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['snowfake'],
    license_files = ('LICENSE',),
    classifiers=CLASSIFIERS,
    install_requires=['numpy', 'matplotlib', 'scipy', 'tqdm'],
    extras_require = {
        'skimage':  ["scikit-image"],
    },
    tests_require=['scikit-image', 'pytest', 'pytest-cov'],
    test_suite='python run_tests.py'
)
