# v0.4.0 (2025/02/03)

- Implementation of Dipole as target
- Bugfix for magnets with meshing parameter = 1
- Additional tests for Dipole extension

# v0.3.1 (2024/10/30)

- Improved Mesing to accept simimal scalar value for all classes except Polyline
  and Dipole
- Include warning when meshing is low
- Improve physics tests

# v0.3.0 (2024/10/29)

- Added the current.Circle class to the force computation
- Improved warning msg when meshing parameter not given
- Added more physics tests

# v0.2.0 (2024/10/11)

- Added a warning when no anchor is given to getFT to avoid misuse.
- Added Cylinder computation
- Added Cylinder Segment computation
- improved internals (there is still a lot of room)

# v0.1.10

- improve docstrings
- fix bad sign of torque

# v0.1.9

- fix docstrings
- bugfix remove pint dependency from development

# v0.1.8

- include Sphere magnets
- rework of getFT workflow
- include object rotations
- update Readme

# v0.1.7

- generalize to getFT instead of separate functions
- add excessive FEM testing
- physics tests
- self-consitency tests

# v0.1.6

- solve pypi problems

# v0.1.5

- solve pypi problems

# v0.1.4

- prepare for pypi
- rename into "magpylib-force"

# v0.1.3

- initial release with first working and tested codes
