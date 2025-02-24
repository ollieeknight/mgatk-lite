"""
mgatk: a mitochondrial genome analysis toolkit
"""
from setuptools import find_packages, setup

dependencies = [
    'click', 
    'pysam', 
    'pytest', 
    'snakemake', 
    'biopython', 
    'numpy', 
    'pandas', 
    'optparse-pretty', 
    'regex', 
    'ruamel.yaml',
    'pulp<=2.7.0',
    'matplotlib'
]

setup(
    name='mgatk',
    version='1.0.0',
    url='https://github.com/ollieeknight/mgatk',
    license='MIT',
    author='Oliver Knight',
    author_email='oliver.knight@charite.de',
    description='Mitochondrial genome analysis toolkit lite',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'mgatk = mgatk.cli:main'
        ],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)