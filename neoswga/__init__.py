"""
NeoSWGA - Selective Whole Genome Amplification primer design.

For CLI usage:
    $ neoswga count-kmers -j params.json
    $ neoswga filter -j params.json
    $ neoswga score -j params.json
    $ neoswga optimize -j params.json

Documentation: See README.md
"""

__version__ = "3.0.0"
__author__ = "Andreas Sjodin"
__email__ = "andreas.sjodin@foi.se"
__license__ = "AGPL-3.0-or-later"


def get_version():
    """Get version string."""
    return __version__


def print_info():
    """Print package information."""
    print(f"NeoSWGA version {__version__}")
    print(f"Author: {__author__}")
    print(f"License: {__license__}")
    print()
    print("Enhanced primer design for selective whole genome amplification")
    print("Documentation: https://github.com/FOI-Bioinformatics/neoswga")
