"""Setup python package"""

from setuptools import setup

setup(
    name='kimono', 
    version='0.1', 
    packages=[
        "kimono",
    ],
    entry_points={
        "console_scripts": [
            "kimono = kimono.cli:main",
        ],
    },
)