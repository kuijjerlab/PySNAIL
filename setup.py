from setuptools import setup

setup(
    name='PySNAIL',
    version='0.4.0',
    description="""Python Implementation for Smooth-quantile Normalization Adaptation for Inference of co-expression Links (PySNAIL)""",
    url='https://github.com/dn070017/Caiman',
    author='Ping-Han Hsieh',
    author_email='dn070017@gmail.com',
    license='MIT',
    packages=['pysnail'],
    install_requires=[
        'bokeh>=2.2.0',
        'numpy>=1.19.1',
        'pandas>=1.1.1',
        'pandarallel>=1.5.4',
        'scipy>=1.5.2',
        'tqdm>=4.62.3'
    ],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={
        'console_scripts': [
            'pysnail = pysnail.__main__:main'
        ],
    },
    zip_safe=False,
    python_requires='>=3.7.7',
)
