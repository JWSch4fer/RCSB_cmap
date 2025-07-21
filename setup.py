# setup.py
import setuptools

setuptools.setup(
    name="RCSB-cmap",                       # your package name
    version="0.1.0",                        # start with 0.1.0, bump on each release
    author="Joseph W Schafer",
    author_email="",
    description="Residue–residue contact map tools",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/JWSch4fer/RCSB_cmap",  # repo or documentation URL
    package_dir={"": "src"},                # tells setuptools to look in src/
    packages=setuptools.find_packages(where="src"),
    install_requires=[
        "biopython>=1.84",
        "numpy",
        "pandas",
        "requests",
        "scipy",
    ],
    entry_points={
        "console_scripts": [
            "rcsb-cmap=project.cli:main",  # makes `rcsb-cmap` command run project.cli.main()
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.11",
)

