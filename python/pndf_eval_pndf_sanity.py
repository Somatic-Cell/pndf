#!/usr/bin/env python3
"""
Sanity tests for python/pndf_eval_pndf.py.

This script does not change the evaluator and does not implement an accelerated
or analytic PNM triangle-intersection p-NDF. It checks basic invariants of the
current sampled footprint-NDF diagnostic:

  1. GT histogram compared with itself must be exactly zero-error.
  2. Reconstructed-mesh histogram compared with itself must be exactly zero-error.
  3. A constant normal map compared with itself must be exactly zero-error.
  4. Optionally, a texel-centered full mesh reconstructed by the same rasterizer
     should reproduce its source normal map at texel centers.

The last check is a rasterizer sanity test. It uses a texel-centered synthetic
mesh and should not be confused with the current C++ full-torus vertex convention.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

# Import the current evaluator implementation. This script lives in the same
# python/ directory, so normal execution via
#   python python/pndf_eval_pndf_sanity.py ...
# should find pndf_eval_pndf.py without extra PYTHONPATH settings.
try:
    import pndf_eval_pndf as pe
except Exception as exc:  # pragma: no cover - diagnostic path
    raise SystemExit(f"Could not import pndf_eval_pndf.py: {exc}") from exc


ZERO_METRIC_KEYS = [
    "tv_distance",
    "hellinger",
    "js_divergence",
    "l1_error",
    "linf_error",
    "peak_location_error_nxy",
]


def metric_status(metrics: dict[str, float], tol: float) -> str:
    for key in ZERO_METRIC_KEYS:
        if abs(float(metrics.get(key, 0.0))) > tol:
            return "FAIL"
    return "OK"


def row_from_metrics(
    test_name: str,
    query: pe.Query,
    metrics: dict[str, float],
    mass_a: float,
    mass_b: float,
    tol: float,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "test": test_name,
        "query_id": query.query_id,
        "label": query.label,
        "u": query.u,
        "v": query.v,
        "sigma_pixels": query.sigma_pixels,
        "mass_a": mass_a,
        "mass_b": mass_b,
        "mass_diff_abs": abs(mass_a - mass_b),
        "status": metric_status(metrics, tol),
    }
    row.update(metrics)
    if extra:
        row.update(extra)
    return row


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_constant_nxy(height: int, width: int, nx: float, ny: float) -> np.ndarray:
    arr = np.zeros((height, width, 2), dtype=np.float32)
    arr[..., 0] = nx
    arr[..., 1] = ny
    return pe.clamp_nxy(arr)


def build_texel_centered_full_mesh(nxy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Build a synthetic full mesh whose vertices are texel centers.

    This is intentionally a rasterizer identity test for the sampled evaluator.
    It is not claiming to match the C++ full-torus vertex convention.
    """
    height, width = nxy.shape[:2]
    vertices = np.zeros((height * width, 4), dtype=np.float32)
    for y in range(height):
        for x in range(width):
            idx = y * width + x
            vertices[idx, 0] = (x + 0.5) / float(width)
            vertices[idx, 1] = (y + 0.5) / float(height)
            vertices[idx, 2:4] = nxy[y, x]

    faces = np.zeros((2 * height * width, 3), dtype=np.int32)
    f = 0
    for y in range(height):
        y1 = (y + 1) % height
        for x in range(width):
            x1 = (x + 1) % width
            v00 = y * width + x
            v10 = y * width + x1
            v01 = y1 * width + x
            v11 = y1 * width + x1
            faces[f] = (v00, v10, v01)
            faces[f + 1] = (v01, v10, v11)
            f += 2
    return vertices, faces


def max_angular_error_deg(a: np.ndarray, b: np.ndarray) -> float:
    an = pe.nxy_to_unit(a)
    bn = pe.nxy_to_unit(b)
    dot = np.clip(np.sum(an * bn, axis=-1), -1.0, 1.0)
    return float(np.degrees(np.arccos(dot)).max())


def run_histogram_identity_tests(
    name: str,
    nxy_a: np.ndarray,
    nxy_b: np.ndarray,
    queries: list[pe.Query],
    normal_bins: int,
    kernel: str,
    truncate_sigma: float,
    tol: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for q in queries:
        ha, ma = pe.ndf_histogram(nxy_a, q, normal_bins, kernel, truncate_sigma)
        hb, mb = pe.ndf_histogram(nxy_b, q, normal_bins, kernel, truncate_sigma)
        metrics = pe.distribution_metrics(ha, hb)
        rows.append(row_from_metrics(name, q, metrics, ma, mb, tol))
    return rows


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Sanity checks for the sampled footprint-NDF evaluator."
    )
    ap.add_argument("--nxy-bin", required=True, type=Path)
    ap.add_argument("--meshbin", type=Path, default=None)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--queries", type=Path, default=None)
    ap.add_argument("--sigma-pixels", default="8,16,32")
    ap.add_argument("--normal-bins", type=int, default=64)
    ap.add_argument("--kernel", choices=["gaussian", "box"], default="gaussian")
    ap.add_argument("--truncate-sigma", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--random-queries", type=int, default=8)
    ap.add_argument("--feature-queries", type=int, default=4)
    ap.add_argument("--seam-queries", type=int, default=4)
    ap.add_argument("--constant-nx", type=float, default=0.12)
    ap.add_argument("--constant-ny", type=float, default=-0.07)
    ap.add_argument("--tol", type=float, default=1.0e-12)
    ap.add_argument(
        "--check-texel-centered-full-mesh",
        action="store_true",
        help=(
            "also build a synthetic texel-centered full mesh and verify that "
            "the rasterizer reproduces the source normal map at texel centers"
        ),
    )
    ap.add_argument(
        "--full-mesh-max-n",
        type=int,
        default=256,
        help="skip texel-centered full-mesh check if max(width,height) is larger",
    )
    ap.add_argument(
        "--no-fail",
        action="store_true",
        help="write diagnostics but do not return a nonzero exit code on failure",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    gt_nxy = pe.load_nxy_bin(args.nxy_bin)
    height, width = gt_nxy.shape[:2]

    if args.queries is not None:
        queries = pe.load_queries_csv(args.queries)
    else:
        sigmas = pe.parse_float_list(args.sigma_pixels)
        queries = pe.make_auto_queries(
            gt_nxy,
            sigmas,
            args.seed,
            args.random_queries,
            args.feature_queries,
            args.seam_queries,
        )
    pe.write_queries_csv(args.out_dir / "queries_used.csv", queries)

    rows: list[dict[str, Any]] = []

    # 1. GT vs GT: histogram, normalization, and metric identity.
    rows.extend(
        run_histogram_identity_tests(
            "gt_vs_gt",
            gt_nxy,
            gt_nxy,
            queries,
            args.normal_bins,
            args.kernel,
            args.truncate_sigma,
            args.tol,
        )
    )

    # 2. Constant map vs itself: checks a trivial degenerate NDF case.
    const_nxy = make_constant_nxy(height, width, args.constant_nx, args.constant_ny)
    rows.extend(
        run_histogram_identity_tests(
            "constant_vs_constant",
            const_nxy,
            const_nxy,
            queries,
            args.normal_bins,
            args.kernel,
            args.truncate_sigma,
            args.tol,
        )
    )

    # 3. Optional mesh reconstruction self-identity.
    mesh_coverage = None
    if args.meshbin is not None:
        mesh_w, mesh_h, vertices, faces, mesh_format = pe.load_meshbin(args.meshbin)
        rec_nxy, covered = pe.rasterize_mesh_to_grid(vertices, faces, height, width)
        mesh_coverage = float(np.mean(covered))
        mesh_rows = run_histogram_identity_tests(
            "mesh_reconstruction_vs_itself",
            rec_nxy,
            rec_nxy,
            queries,
            args.normal_bins,
            args.kernel,
            args.truncate_sigma,
            args.tol,
        )
        for row in mesh_rows:
            row["meshbin"] = str(args.meshbin)
            row["mesh_format"] = mesh_format
            row["mesh_declared_width"] = mesh_w
            row["mesh_declared_height"] = mesh_h
            row["vertices"] = len(vertices)
            row["faces"] = len(faces)
            row["coverage_ratio"] = mesh_coverage
        rows.extend(mesh_rows)

    # 4. Optional rasterizer identity using a synthetic texel-centered full mesh.
    full_mesh_note = "not_requested"
    full_mesh_max_ang = math.nan
    full_mesh_coverage = math.nan
    if args.check_texel_centered_full_mesh:
        if max(width, height) <= args.full_mesh_max_n:
            full_vertices, full_faces = build_texel_centered_full_mesh(gt_nxy)
            full_rec, full_covered = pe.rasterize_mesh_to_grid(
                full_vertices, full_faces, height, width
            )
            full_mesh_coverage = float(np.mean(full_covered))
            full_mesh_max_ang = max_angular_error_deg(gt_nxy, full_rec)
            full_mesh_note = "ran"
            rows.extend(
                run_histogram_identity_tests(
                    "texel_centered_full_mesh_vs_gt",
                    gt_nxy,
                    full_rec,
                    queries,
                    args.normal_bins,
                    args.kernel,
                    args.truncate_sigma,
                    max(args.tol, 1.0e-7),
                )
            )
        else:
            full_mesh_note = (
                f"skipped_max_dimension_{max(width, height)}_gt_{args.full_mesh_max_n}"
            )

    csv_path = args.out_dir / "sanity_metrics.csv"
    write_csv(csv_path, rows)

    failed = [r for r in rows if r.get("status") != "OK"]
    report_path = args.out_dir / "sanity_report.txt"
    with report_path.open("w", encoding="utf-8") as f:
        f.write("Sampled footprint-NDF evaluator sanity report\n")
        f.write("==============================================\n")
        f.write("This validates the current sampled histogram diagnostic only.\n")
        f.write("It does not validate the analytic continuous triangle-intersection P-NDF.\n\n")
        f.write(f"nxy_bin: {args.nxy_bin}\n")
        f.write(f"meshbin: {args.meshbin}\n")
        f.write(f"resolution: {width} x {height}\n")
        f.write(f"queries: {len(queries)}\n")
        f.write(f"normal_bins: {args.normal_bins}\n")
        f.write(f"kernel: {args.kernel}\n")
        f.write(f"truncate_sigma: {args.truncate_sigma}\n")
        f.write(f"tolerance: {args.tol}\n")
        if mesh_coverage is not None:
            f.write(f"mesh_reconstruction_coverage_ratio: {mesh_coverage:.9g}\n")
        f.write(f"texel_centered_full_mesh_check: {full_mesh_note}\n")
        if full_mesh_note == "ran":
            f.write(f"texel_centered_full_mesh_coverage_ratio: {full_mesh_coverage:.9g}\n")
            f.write(f"texel_centered_full_mesh_max_ang_deg: {full_mesh_max_ang:.9g}\n")
        f.write(f"rows: {len(rows)}\n")
        f.write(f"failed_rows: {len(failed)}\n")
        if failed:
            f.write("\nFirst failed rows:\n")
            for r in failed[:16]:
                f.write(
                    f"  {r['test']} {r['query_id']} sigma={r['sigma_pixels']} "
                    f"TV={float(r.get('tv_distance', 0.0)):.6g} "
                    f"H={float(r.get('hellinger', 0.0)):.6g} "
                    f"JS={float(r.get('js_divergence', 0.0)):.6g}\n"
                )

    print(f"Wrote {csv_path}")
    print(f"Wrote {report_path}")
    if failed:
        print(f"FAILED sanity checks: {len(failed)} rows exceed tolerance {args.tol}")
        return 0 if args.no_fail else 1
    print("All sanity checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
