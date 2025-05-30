from setuptools import setup

setup(name="syotools",
      version="1.1.13",
      description="Science Yield Optimization Tools (SYOTools)", 
      author="Jason Tumlinson, STScI", 
      author_email="tumlinson@stsci.edu",
      license="BSD",
      keywords=["simulation", "astronomy", "astrophysics"],
      url="https://github.com/spacetelescope/hwo-tools",
      packages=["syotools", "syotools/coronagraph", "syotools","syotools/reference_data", 
                "syotools/reference_data/pysynphot_data","syotools/defaults",  
                "syotools/interface", "syotools/models", "syotools/persistence",  
                "syotools/spectra", "syotools/utils", "syotools/wrappers"],
      package_data={'':['*.yaml', '*.fits', 'fesc/*txt', '*.txt', '*.dat', 
                    'reference_data/pysynphot_data/calspec/*','reference_data/pysynphot_data/comp/*', 
                    'reference_data/pysynphot_data/grid/*', 'reference_data/pysynphot_data/mtab/*', 
                    'reference_data/pysynphot_data/comp/nonhst/*', 'reference_data/pysynphot_data/grid/agn/*', 
                    'reference_data/pysynphot_data/grid/etc_models/*', 'reference_data/pysynphot_data/grid/galactic/*', 
                    'reference_data/pysynphot_data/grid/kc96/*',
                    'reference_data/pysynphot_data/grid/pickles/*', 'reference_data/pysynphot_data/grid/pickles/dat_uvk/*', 
                    'reference_data/pysynphot_data/extinction/*', 'sci_eng_interface/*json']}, 
      include_package_data=True,
      classifiers=[
          "Programming Language :: Python :: 3.9",
      ],
      install_requires=[
            "hwo_sci_eng==0.1.8",
            "setuptools>=61.0", 
            "astropy>6", 
            "pysynphot==2.0.0", 
            "numpy", 
            "scipy", 
            "specutils", 
      ],
)
