import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hifise",
    version="1.0",
    author='Chentao Yang, Guanliang Meng',
    author_email='yangchentao@genomics.cn',
    description="test my hifi-se program",
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires='>=3.5',
    url='https://github.com/comery/HIFI-barcode-SE400;https://github.com/linzhi2013/bold_identification',
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=['bold-identification==0.0.20', 'biopython>=1.5'],

    scripts = [ 'HIFI-SE.py'],

    classifiers=(
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
    ),
)
