# external/

This directory is reserved for optional Git submodules.

Phase 1 intentionally builds without external dependencies.  Suggested optional submodules later:

- Eigen: fixed-size matrix replacement for the in-tree small matrix code.
- CLI11: command line parser replacement for the small parser in apps.
- Catch2: unit-test framework replacement for the minimal test harness.
- CGAL: periodic Delaunay for A-tri, not needed for B-QEM Phase 1.
- OpenMesh: reference decimator only, not used by the Phase 1 core.
