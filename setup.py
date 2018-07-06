"""AeroPy: An easy to use aerodynamic tool."""

from setuptools import setup, find_packages

setup(name='aeropy',
      version='0.0',
      description='An easy to use aerodynamic tool.',
      url='NA',
      author='leal26',
      author_email='leal26@tamu.edu',
      license='MIT',
      packages=['aeropy','aeropy.CST_2D','aeropy.CST_3D',
                'aeropy.morphing','aeropy.geometry',
                'aeropy.filehanding'],
      zip_safe=False,
          package_data={
        # If any package contains *.exe and avian files, include them:
        '': ['*.exe', 'avian'],
        # And include any *.exe files found in the 'CST' package, too:
        'CST': ['*.exe', 'avian'],
        # And include any *.exe files found in the 'geometry' package, too:
        'geometry': ['*.exe'],
        }
      )
