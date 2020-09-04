from setuptools import setup

setup(
    name='caiman',
    version='0.0.b1',
    description="""Count Adjustment to Improve the Modeling of Gene Association
    Networks (CAIMAN)""",
    url='https://github.com/dn070017/Caiman',
    author='Ping-Han Hsieh',
    author_email='dn070017@gmail.com',
    license='MIT',
    packages=['caiman'],
    install_requires=[
        'bokeh>=2.2.0',
        'numpy>=1.19.1',
        'pandas>=1.1.1',
        'scipy>=1.5.2',
    ],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    entry_points={
        'console_scripts': [
            'caiman = caiman.__main__:main'
        ],
    },
    zip_safe=False,
    python_requires='>=3.7.7',
)
