from setuptools import setup, find_packages

setup(
    name="pandamap-color",
    version="2.1",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "matplotlib",
        "biopython",
    ],
    entry_points={
        "console_scripts": [
            "pandamap-color=pandamap_color.cli:main",
        ],
    },
)