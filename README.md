# pndf-qem-phase1

Phase 1 prototype repository for **normal-only QEM** and **4D QEM** on a periodic position-normal manifold.

This package is intended as a GitHub repository seed for Windows development.  It deliberately avoids mandatory external dependencies so that the first build is simple.  Optional external libraries can be added later under `external/` as Git submodules.

## Goals

- Start from a full periodic 1024x1024 normal texture manifold.
- Preserve periodic boundary conditions by using true torus connectivity.
- Compare:
  - `normal` mode: normal-only attribute QEM.
  - `4d` mode: QEM on the 4D graph `(u,v,n_x,n_y)`.
- Do **not** use p-NDF, Jacobian, footprint, half-vector, or rendering error in Phase 1.

## Non-goals

- No A-tri implementation in this Phase 1 package.
- No p-NDF query implementation yet.
- No GPU hierarchy yet.
- No grid/warped-grid proxy baseline.

## Build on Windows

Install CMake and Ninja.  Visual Studio 2022 is recommended.

```powershell
cmake --preset msvc-release
cmake --build --preset msvc-release
ctest --test-dir build/msvc-release --output-on-failure
```

If you prefer Visual Studio's generator, create another CMake preset or open the folder directly in Visual Studio.

## Input preparation

The C++ core reads a simple binary file:

```text
int32 N
float nx[N*N]
float ny[N*N]
```

Use Python to convert EXR files:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r python\requirements.txt
python python\prepare_nxy.py --map scratched --data-dir C:\path\to\data --out out\scratched_nxy.bin
```

Map conventions:

- `isotropic.exr`: height map -> periodic central differences.
- `scratched.exr`: height map -> periodic central differences.
- `brush.exr`: normal map -> OpenCV BGR channels `[2,1]` are `(n_x,n_y)`.

## Run QEM

```powershell
build\msvc-release\pndf_run_qem.exe --input out\scratched_nxy.bin --output out\scratched_normal_262k.meshbin --target 262144 --mode normal
build\msvc-release\pndf_eval_mesh.exe --input out\scratched_nxy.bin --mesh out\scratched_normal_262k.meshbin --csv out\scratched_normal_262k.vertex_error.csv
```

4D QEM:

```powershell
build\msvc-release\pndf_run_qem.exe --input out\scratched_nxy.bin --output out\scratched_4d_262k.meshbin --target 262144 --mode 4d
```

Budget sweep:

```powershell
.\scripts\run_budget_sweep.ps1 -Preset msvc-release -DataDir C:\path\to\data -OutDir out -Mode normal
.\scripts\run_budget_sweep.ps1 -Preset msvc-release -DataDir C:\path\to\data -OutDir out -Mode 4d
```

## Periodic boundary invariant

The mesh is built with torus connectivity:

```text
(i+1,j) -> ((i+1) mod N, j)
(i,j+1) -> (i, (j+1) mod N)
```

Thus the seam is not an open boundary.  Edge lengths, candidate placement, and face unwrapping use periodic coordinates.

Important tests:

- `test_wrap_delta_midpoint`
- `test_torus_mesh_counts`
- `test_small_torus_collapse_normal_qem`

## Current prototype limitations

1. Evaluation currently reports vertex-sample error only.  A full rasterized piecewise-linear reconstruction evaluator should be added next.
2. 4D QEM placement currently tests midpoint and endpoints rather than solving the unconstrained 4D minimizer.  This is robust but conservative.
3. The collapse queue is a simple priority queue with periodic rebuilds; an indexed heap would be better later.
4. The link condition and non-flip checks are conservative.  More robust topological checks should be added before final experiments.
5. No collapse-history hierarchy is emitted yet.  This should be added for future BVH/k-DOP work.

## Suggested next patches

1. Add rasterized normal-map evaluator for output mesh.
2. Emit collapse history tree.
3. Implement unconstrained / constrained 4D optimal placement.
4. Add true p95/p99 texel-level normal error, not only vertex-sample error.
5. Add AABB vs k-DOP diagnostic for simplified mesh nodes.
