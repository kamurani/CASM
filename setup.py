
from importlib.metadata import entry_points
import io
import os
from setuptools import setup


VERSION = None
HERE = os.path.abspath(os.path.dirname(__file__))
NAME = "CASM"

DESCRIPTION = "Clustering Abstract Representations of Structural Motifs"

# Import the PYPI README and use it as the long-description.
# Note: this will only work if "README.md" is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(HERE, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__))

setup(
    name="CASM",
    version='0.1.0',
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    py_modules=['cli'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'casm = CASM.CLI:main',
        ]
    }
)