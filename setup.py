from setuptools import setup, find_packages

setup(
    cffi_modules=["gmat/remma/_build.py:ffi", "gmat/process_plink/_build.py:ffi"],
    name='gmat',
    version='2021.1.25',
    description='Genomic Multivariate Analysis Tools',
    # long_description=open('README.rst').read(),
    # long_description_content_type="text/markdown",
    author='Chao Ning',
    author_email='ningchao91@gmail.com',
    packages=find_packages(),
    platforms=["all"],
    url="https://github.com/chaoning/GMAT",
    include_package_data = True,
    install_requires=[
        'numpy>=1.16.0',
        'pandas>=0.19.0',
        'scipy>=1.1.1',
        'cffi>=1.12.0',
        'pandas_plink>=2.0.0',
        'pysnptools>=0.4.26',
        'tqdm>=4.43.0',
    ],
)
