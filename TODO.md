# DONE SINCE LAST LAST TIME

- magnet-wire force/torque codes completed
    - confirmed by FEM
    - confirmed by self-consistency
    - confirmed by analytical fundamental solutions
- COMSOL confirm ANSYS and bad FE convergence
- We understand the possibilities and limitation (no wire ends) of our force codes


# DONE SINCE LAST TIME

- restructured repository so it can be made installable and used by end-users
- removed unit-dependencies - everything is SI and magpylib 5 compatible
- setting up a testing framework
    - physics tests
    - interface tests
    - FE tests
    - self-consistency tests
- vectorization for multiple inputs
- input checks and auto-formatting
- dynamic wire meshing
- some doc-strings
- add License
- making installable, current version is 0.1.2dev


# TODO CRITICAL

- testing
    - FE tests wire-wire
    - FE tests magnet-wire
    - FE tests finite-sized wires magnet
    - more interface and meshing tests
- reduce computation effort when no anchor is given or proof that this will not help much
- review and clean up the meshers
- implement dynamic magnet meshing ?
- speed-boost: replace integration using newton_cotes
- profiling
    - memory
    - performance
- upscaling computation to SAL servers


# TODO UNCRITICAL

- shape handling
    - include paths
    - avoid sum over all sources if sumup=False
- more doc strings
- make units optional
- include orientations
- allow Collection inputs
- Arc solution
- linting



