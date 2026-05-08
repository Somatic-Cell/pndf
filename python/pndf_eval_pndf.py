#!/usr/bin/env python3
"""
Sampled reference p-NDF evaluator and visualizer for PNDF/QEM meshes.

This is intentionally a *reference/evaluation* tool, not an acceleration structure.
It does not use BVH, k-DOP, LOD, CUDA, or stochastic shortcuts.

Definition used in this script
------------------------------
For each query q=(u,v,sigma), the tool estimates a footprint-conditioned NDF by
sampling texel centers under a periodic footprint kernel and histogramming their
projected normals (nx,ny) on a fixed normal-domain grid.  The same estimator is
used for:

  1. GT normal map loaded from nxy.bin
  2. simplified mesh reconstructed to texel centers by periodic PL rasterization

This is a discrete sampled p-NDF diagnostic.  It is designed to avoid misleading
comparisons while staying simple and robust.  It is not the final accelerated
continuous triangle-intersection evaluator used by the 2025 PNM renderer.

Inputs
------
  --nxy-bin : int32 N, float nx[N*N], float ny[N*N]
              or int32 W,H, float nx[W*H], float ny[W*H]
  --meshbin : int32 W,H,V,F, then V*(float u,v,nx,ny), F*(int32 i,j,k)

Outputs
-------
  queries_used.csv
  metrics.csv
  densities/*.npz
  images/*.png

Example
-------
  python python/pndf_eval_pndf.py `
    --nxy-bin out256_normal_rebased/brush_ds4_nxy.bin `
    --meshbin out256_normal_rebased/brush_normal_ds4_v4096.meshbin `
    --out-dir pndf_eval/brush_normal_v4096 `
    --normal-bins 64 `
    --sigma-pixels 4,8,16
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS = 1.0e-12


@dataclass
class Query:
    query_id: str
    u: float
    v: float
    sigma_pixels: float
    label: str


# -----------------------------------------------------------------------------
# Basic normal/mesh I/O
# -----------------------------------------------------------------------------


def clamp_nxy(nxy: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(nxy, axis=-1, keepdims=True)
    return (nxy * np.minimum(1.0, 0.999 / np.maximum(n, EPS))).astype(np.float32)


def nxy_to_unit(nxy: np.ndarray) -> np.ndarray:
    z = np.sqrt(np.maximum(0.0, 1.0 - np.sum(nxy * nxy, axis=-1, keepdims=True)))
    n = np.concatenate([nxy, z], axis=-1)
    n /= np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), EPS)
    return n.astype(np.float32)


def load_nxy_bin(path: Path) -> np.ndarray:
    data = path.read_bytes()
    if len(data) < 4:
        raise ValueError(f"nxy file too small: {path}")
    ints = np.frombuffer(data[:8], dtype=np.int32, count=min(2, len(data) // 4))

    # Current prepare_nxy.py format: int32 N, float nx[N*N], float ny[N*N]
    n = int(ints[0])
    size_single_header = 4 + 2 * n * n * 4
    if n > 0 and len(data) == size_single_header:
        arr = np.frombuffer(data[4:], dtype=np.float32)
        nx = arr[: n * n].reshape(n, n)
        ny = arr[n * n :].reshape(n, n)
        return clamp_nxy(np.stack([nx, ny], axis=-1).copy())

    # Optional format: int32 W, int32 H, float nx[W*H], float ny[W*H]
    if len(ints) >= 2:
        w, h = int(ints[0]), int(ints[1])
        size_two_header = 8 + 2 * w * h * 4
        if w > 0 and h > 0 and len(data) == size_two_header:
            arr = np.frombuffer(data[8:], dtype=np.float32)
            nx = arr[: w * h].reshape(h, w)
            ny = arr[w * h :].reshape(h, w)
            return clamp_nxy(np.stack([nx, ny], axis=-1).copy())

    raise ValueError(f"Could not parse nxy binary: {path}")


def normalize_mesh_uv(vertices: np.ndarray, width: int, height: int) -> np.ndarray:
    uv = vertices[:, :2]
    # Accept old prototypes with pixel-space UVs.
    if np.nanmax(np.abs(uv)) > 2.0:
        vertices[:, 0] /= float(width)
        vertices[:, 1] /= float(height)
    vertices[:, :2] %= 1.0
    vertices[:, 2:4] = clamp_nxy(vertices[:, 2:4])
    return vertices


def load_meshbin(path: Path) -> tuple[int, int, np.ndarray, np.ndarray, str]:
    data = path.read_bytes()
    if len(data) < 12:
        raise ValueError(f"mesh file too small: {path}")

    # Current format: int32 width, int32 height, int32 V, int32 F
    if len(data) >= 16:
        width, height, vcount, fcount = np.frombuffer(data[:16], dtype=np.int32, count=4)
        width, height, vcount, fcount = int(width), int(height), int(vcount), int(fcount)
        expected = 16 + 16 * vcount + 12 * fcount
        if width > 0 and height > 0 and vcount >= 0 and fcount >= 0 and expected == len(data):
            off = 16
            vertices = np.frombuffer(data[off : off + 16 * vcount], dtype=np.float32).reshape(vcount, 4).copy()
            off += 16 * vcount
            faces = np.frombuffer(data[off : off + 12 * fcount], dtype=np.int32).reshape(fcount, 3).copy()
            return width, height, normalize_mesh_uv(vertices, width, height), faces, "width,height,V,F"

    # Legacy: int32 N, int32 V, int32 F
    n, vcount, fcount = np.frombuffer(data[:12], dtype=np.int32, count=3)
    n, vcount, fcount = int(n), int(vcount), int(fcount)
    expected = 12 + 16 * vcount + 12 * fcount
    if n > 0 and vcount >= 0 and fcount >= 0 and expected == len(data):
        off = 12
        vertices = np.frombuffer(data[off : off + 16 * vcount], dtype=np.float32).reshape(vcount, 4).copy()
        off += 16 * vcount
        faces = np.frombuffer(data[off : off + 12 * fcount], dtype=np.int32).reshape(fcount, 3).copy()
        return n, n, normalize_mesh_uv(vertices, n, n), faces, "N,V,F"

    raise ValueError(f"Could not parse mesh binary: {path}")


# -----------------------------------------------------------------------------
# Periodic PL rasterization of simplified mesh to texel centers.
# -----------------------------------------------------------------------------


def wrap_delta(d: np.ndarray) -> np.ndarray:
    return (d + 0.5) % 1.0 - 0.5


def rasterize_mesh_to_grid(vertices: np.ndarray, faces: np.ndarray, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
    """Rasterize a periodic PL normal mesh to texel centers.

    Returns reconstructed nxy and a boolean coverage mask.  Missing pixels, if any,
    are filled from nearest vertex after coverage is recorded.
    """
    uv = vertices[:, :2].astype(np.float64) % 1.0
    vals = vertices[:, 2:4].astype(np.float64)
    out = np.full((height, width, 2), np.nan, dtype=np.float32)
    covered = np.zeros((height, width), dtype=bool)

    for tri in faces:
        ids = tri.astype(np.int64)
        p0 = uv[ids[0]]
        p1 = p0 + wrap_delta(uv[ids[1]] - p0)
        p2 = p0 + wrap_delta(uv[ids[2]] - p0)
        vv = vals[ids]
        base_tri = np.stack([p0, p1, p2], axis=0)

        for sy in (-1, 0, 1):
            for sx in (-1, 0, 1):
                p = base_tri + np.array([sx, sy], dtype=np.float64)
                xmin, xmax = p[:, 0].min(), p[:, 0].max()
                ymin, ymax = p[:, 1].min(), p[:, 1].max()
                if xmax < 0.0 or xmin > 1.0 or ymax < 0.0 or ymin > 1.0:
                    continue

                ix0 = max(0, int(math.floor(xmin * width - 1)))
                ix1 = min(width - 1, int(math.ceil(xmax * width + 1)))
                iy0 = max(0, int(math.floor(ymin * height - 1)))
                iy1 = min(height - 1, int(math.ceil(ymax * height + 1)))
                if ix1 < ix0 or iy1 < iy0:
                    continue

                xs = (np.arange(ix0, ix1 + 1, dtype=np.float64) + 0.5) / width
                ys = (np.arange(iy0, iy1 + 1, dtype=np.float64) + 0.5) / height
                xx, yy = np.meshgrid(xs, ys, indexing="xy")
                q = np.stack([xx.ravel(), yy.ravel()], axis=-1)

                e0 = p[1] - p[0]
                e1 = p[2] - p[0]
                den = e0[0] * e1[1] - e0[1] * e1[0]
                if abs(den) < 1.0e-14:
                    continue

                e2 = q - p[0]
                a = (e2[:, 0] * e1[1] - e2[:, 1] * e1[0]) / den
                b = (e0[0] * e2[:, 1] - e0[1] * e2[:, 0]) / den
                c = 1.0 - a - b
                inside = (a >= -1.0e-7) & (b >= -1.0e-7) & (c >= -1.0e-7)
                if not np.any(inside):
                    continue

                interp = c[:, None] * vv[0] + a[:, None] * vv[1] + b[:, None] * vv[2]
                coords = q[inside]
                px = np.clip(np.floor(coords[:, 0] * width).astype(np.int64), 0, width - 1)
                py = np.clip(np.floor(coords[:, 1] * height).astype(np.int64), 0, height - 1)
                out[py, px] = interp[inside].astype(np.float32)
                covered[py, px] = True

    missing = ~np.isfinite(out).all(axis=-1)
    if np.any(missing):
        # Deterministic slow fallback without scipy.
        pts = uv.copy()
        vals_v = vertices[:, 2:4].astype(np.float32)
        yy, xx = np.nonzero(missing)
        q = np.stack([(xx + 0.5) / width, (yy + 0.5) / height], axis=-1)
        for i, qq in enumerate(q):
            d = np.abs(pts - qq)
            d = np.minimum(d, 1.0 - d)
            nn = int(np.argmin(np.sum(d * d, axis=1)))
            out[yy[i], xx[i]] = vals_v[nn]

    return clamp_nxy(out), covered


# -----------------------------------------------------------------------------
# Query generation / loading.
# -----------------------------------------------------------------------------


def parse_float_list(s: str) -> list[float]:
    vals = [float(x.strip()) for x in s.split(",") if x.strip()]
    if not vals:
        raise ValueError("empty float list")
    return vals


def feature_score(nxy: np.ndarray) -> np.ndarray:
    n = nxy_to_unit(nxy)
    dx = np.linalg.norm(n - np.roll(n, 1, axis=1), axis=-1) + np.linalg.norm(n - np.roll(n, -1, axis=1), axis=-1)
    dy = np.linalg.norm(n - np.roll(n, 1, axis=0), axis=-1) + np.linalg.norm(n - np.roll(n, -1, axis=0), axis=-1)
    return (dx + dy).astype(np.float32)


def make_auto_queries(nxy: np.ndarray, sigmas: list[float], seed: int, random_count: int, feature_count: int, seam_count: int) -> list[Query]:
    height, width = nxy.shape[:2]
    rng = np.random.default_rng(seed)
    positions: list[tuple[str, float, float]] = []

    for i in range(random_count):
        u, v = rng.random(2)
        positions.append((f"random{i:02d}", float(u), float(v)))

    if feature_count > 0:
        score = feature_score(nxy).ravel()
        k = min(max(feature_count * 64, feature_count), score.size)
        idxs = np.argpartition(score, -k)[-k:]
        idxs = idxs[np.argsort(score[idxs])[::-1]]
        chosen: list[int] = []
        min_dist = 0.08
        for idx in idxs:
            y = int(idx // width)
            x = int(idx % width)
            u = (x + 0.5) / width
            v = (y + 0.5) / height
            ok = True
            for old in chosen:
                oy = int(old // width)
                ox = int(old % width)
                ou = (ox + 0.5) / width
                ov = (oy + 0.5) / height
                d = np.abs(np.array([u - ou, v - ov]))
                d = np.minimum(d, 1.0 - d)
                if float(np.linalg.norm(d)) < min_dist:
                    ok = False
                    break
            if ok:
                chosen.append(int(idx))
                positions.append((f"feature{len(chosen)-1:02d}", float(u), float(v)))
                if len(chosen) >= feature_count:
                    break

    seam_positions = [
        ("seam_u0", 0.01, 0.50),
        ("seam_u1", 0.99, 0.50),
        ("seam_v0", 0.50, 0.01),
        ("seam_v1", 0.50, 0.99),
        ("corner00", 0.01, 0.01),
        ("corner11", 0.99, 0.99),
    ]
    for label, u, v in seam_positions[:seam_count]:
        positions.append((label, u, v))

    queries: list[Query] = []
    qid = 0
    for label, u, v in positions:
        for sigma in sigmas:
            queries.append(Query(f"q{qid:04d}", u % 1.0, v % 1.0, float(sigma), label))
            qid += 1
    return queries


def load_queries_csv(path: Path) -> list[Query]:
    queries: list[Query] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"u", "v", "sigma_pixels"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"query CSV missing columns: {sorted(missing)}")
        for i, row in enumerate(reader):
            qid = row.get("query_id") or f"q{i:04d}"
            label = row.get("label") or qid
            queries.append(Query(qid, float(row["u"]) % 1.0, float(row["v"]) % 1.0, float(row["sigma_pixels"]), label))
    return queries


def write_queries_csv(path: Path, queries: list[Query]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["query_id", "label", "u", "v", "sigma_pixels"])
        writer.writeheader()
        for q in queries:
            writer.writerow({"query_id": q.query_id, "label": q.label, "u": q.u, "v": q.v, "sigma_pixels": q.sigma_pixels})


# -----------------------------------------------------------------------------
# Sampled p-NDF histogram / metrics.
# -----------------------------------------------------------------------------


def footprint_samples(nxy: np.ndarray, q: Query, kernel: str, truncate_sigma: float) -> tuple[np.ndarray, np.ndarray]:
    height, width = nxy.shape[:2]
    sigma = max(float(q.sigma_pixels), 1.0e-6)
    if kernel == "gaussian":
        radius = max(1, int(math.ceil(truncate_sigma * sigma)))
    elif kernel == "box":
        radius = max(1, int(math.ceil(sigma)))
    else:
        raise ValueError(f"unknown kernel: {kernel}")

    cx = q.u * width - 0.5
    cy = q.v * height - 0.5
    ix_center = int(round(cx))
    iy_center = int(round(cy))
    xs = np.arange(ix_center - radius, ix_center + radius + 1, dtype=np.int64)
    ys = np.arange(iy_center - radius, iy_center + radius + 1, dtype=np.int64)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")

    u = (xx.astype(np.float64) + 0.5) / width
    v = (yy.astype(np.float64) + 0.5) / height
    du = ((u - q.u + 0.5) % 1.0 - 0.5) * width
    dv = ((v - q.v + 0.5) % 1.0 - 0.5) * height
    r2 = du * du + dv * dv

    if kernel == "gaussian":
        w = np.exp(-0.5 * r2 / (sigma * sigma)).astype(np.float64)
        if truncate_sigma > 0:
            w[r2 > (truncate_sigma * sigma) ** 2] = 0.0
    else:
        w = ((np.abs(du) <= sigma) & (np.abs(dv) <= sigma)).astype(np.float64)

    vals = nxy[yy % height, xx % width].reshape(-1, 2)
    weights = w.reshape(-1)
    valid = weights > 0
    return vals[valid], weights[valid]


def ndf_histogram(nxy: np.ndarray, q: Query, normal_bins: int, kernel: str, truncate_sigma: float) -> tuple[np.ndarray, float]:
    vals, weights = footprint_samples(nxy, q, kernel, truncate_sigma)
    hist = np.zeros((normal_bins, normal_bins), dtype=np.float64)
    if len(vals) == 0:
        return hist, 0.0
    ix = np.floor((vals[:, 0] + 1.0) * 0.5 * normal_bins).astype(np.int64)
    iy = np.floor((vals[:, 1] + 1.0) * 0.5 * normal_bins).astype(np.int64)
    valid = (ix >= 0) & (ix < normal_bins) & (iy >= 0) & (iy < normal_bins)
    np.add.at(hist, (iy[valid], ix[valid]), weights[valid])
    return hist, float(np.sum(weights[valid]))


def normalize_density(h: np.ndarray) -> np.ndarray:
    s = float(np.sum(h))
    if s <= 0:
        return np.zeros_like(h, dtype=np.float64)
    return h.astype(np.float64) / s


def distribution_metrics(gt_hist: np.ndarray, mesh_hist: np.ndarray) -> dict[str, float]:
    p = normalize_density(gt_hist)
    q = normalize_density(mesh_hist)
    tv = 0.5 * float(np.sum(np.abs(p - q)))
    overlap = float(np.sum(np.minimum(p, q)))
    hellinger = float(math.sqrt(max(0.0, 0.5 * np.sum((np.sqrt(p) - np.sqrt(q)) ** 2))))
    m = 0.5 * (p + q)
    mask_p = p > 0
    mask_q = q > 0
    js = 0.0
    if np.any(mask_p):
        js += 0.5 * float(np.sum(p[mask_p] * np.log(p[mask_p] / np.maximum(m[mask_p], EPS))))
    if np.any(mask_q):
        js += 0.5 * float(np.sum(q[mask_q] * np.log(q[mask_q] / np.maximum(m[mask_q], EPS))))

    gt_peak = np.unravel_index(int(np.argmax(p)), p.shape)
    mesh_peak = np.unravel_index(int(np.argmax(q)), q.shape)
    bins = p.shape[0]
    # Convert bin-index distance to projected-normal coordinate distance.
    peak_dist = float(np.linalg.norm((np.array(gt_peak[::-1]) - np.array(mesh_peak[::-1])) * (2.0 / bins)))
    peak_ratio = float(np.max(q) / max(float(np.max(p)), EPS))

    return {
        "tv_distance": tv,
        "overlap": overlap,
        "hellinger": hellinger,
        "js_divergence": js,
        "l1_error": float(np.sum(np.abs(p - q))),
        "linf_error": float(np.max(np.abs(p - q))),
        "peak_location_error_nxy": peak_dist,
        "peak_value_ratio": peak_ratio,
        "gt_peak_bin_x": int(gt_peak[1]),
        "gt_peak_bin_y": int(gt_peak[0]),
        "mesh_peak_bin_x": int(mesh_peak[1]),
        "mesh_peak_bin_y": int(mesh_peak[0]),
    }


# -----------------------------------------------------------------------------
# Visualization.
# -----------------------------------------------------------------------------


def safe_name(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s.strip("_") or "query"


def save_density_comparison(path: Path, gt_hist: np.ndarray, mesh_hist: np.ndarray, q: Query, title_prefix: str) -> None:
    p = normalize_density(gt_hist)
    r = normalize_density(mesh_hist)
    diff = np.abs(p - r)
    vmax = max(float(np.percentile(p, 99.8)), float(np.percentile(r, 99.8)), EPS)
    diff_vmax = max(float(np.percentile(diff, 99.8)), EPS)

    fig, ax = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    extent = (-1.0, 1.0, -1.0, 1.0)
    im0 = ax[0].imshow(p, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=vmax)
    ax[0].set_title("GT sampled p-NDF")
    im1 = ax[1].imshow(r, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=vmax)
    ax[1].set_title("mesh sampled p-NDF")
    im2 = ax[2].imshow(diff, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=diff_vmax)
    ax[2].set_title("abs difference")
    for a in ax:
        a.set_xlabel("n_x")
        a.set_ylabel("n_y")
    fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.02)
    fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.02)
    fig.colorbar(im2, ax=ax[2], fraction=0.046, pad=0.02)
    fig.suptitle(f"{title_prefix}: {q.query_id} {q.label}, u={q.u:.4f}, v={q.v:.4f}, sigma={q.sigma_pixels:g}px")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description="Sampled p-NDF reference evaluator and visualizer.")
    ap.add_argument("--nxy-bin", required=True, type=Path)
    ap.add_argument("--meshbin", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--queries", type=Path, default=None, help="optional query CSV with columns query_id,u,v,sigma_pixels,label")
    ap.add_argument("--sigma-pixels", default="4,8,16", help="comma-separated sigmas for auto queries")
    ap.add_argument("--normal-bins", type=int, default=64)
    ap.add_argument("--kernel", choices=["gaussian", "box"], default="gaussian")
    ap.add_argument("--truncate-sigma", type=float, default=3.0, help="Gaussian support radius in sigmas")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--random-queries", type=int, default=4)
    ap.add_argument("--feature-queries", type=int, default=4)
    ap.add_argument("--seam-queries", type=int, default=4)
    ap.add_argument("--dump-densities", action="store_true", help="write per-query .npz density dumps")
    ap.add_argument("--no-images", action="store_true", help="skip PNG comparison images")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    image_dir = args.out_dir / "images"
    density_dir = args.out_dir / "densities"
    if not args.no_images:
        image_dir.mkdir(parents=True, exist_ok=True)
    if args.dump_densities:
        density_dir.mkdir(parents=True, exist_ok=True)

    gt_nxy = load_nxy_bin(args.nxy_bin)
    height, width = gt_nxy.shape[:2]
    mesh_width, mesh_height, vertices, faces, mesh_format = load_meshbin(args.meshbin)
    rec_nxy, covered = rasterize_mesh_to_grid(vertices, faces, height, width)
    coverage_ratio = float(np.mean(covered))

    if args.queries is not None:
        queries = load_queries_csv(args.queries)
    else:
        sigmas = parse_float_list(args.sigma_pixels)
        queries = make_auto_queries(gt_nxy, sigmas, args.seed, args.random_queries, args.feature_queries, args.seam_queries)
    write_queries_csv(args.out_dir / "queries_used.csv", queries)

    metrics_path = args.out_dir / "metrics.csv"
    rows: list[dict[str, object]] = []
    title_prefix = f"{Path(args.meshbin).stem}"

    for qi, query in enumerate(queries):
        gt_hist, mass_gt = ndf_histogram(gt_nxy, query, args.normal_bins, args.kernel, args.truncate_sigma)
        mesh_hist, mass_mesh = ndf_histogram(rec_nxy, query, args.normal_bins, args.kernel, args.truncate_sigma)
        m = distribution_metrics(gt_hist, mesh_hist)
        row: dict[str, object] = {
            "query_id": query.query_id,
            "label": query.label,
            "u": query.u,
            "v": query.v,
            "sigma_pixels": query.sigma_pixels,
            "kernel": args.kernel,
            "normal_bins": args.normal_bins,
            "mass_gt": mass_gt,
            "mass_mesh": mass_mesh,
            "mass_ratio": mass_mesh / max(mass_gt, EPS),
            "coverage_ratio": coverage_ratio,
            "nxy_bin": str(args.nxy_bin),
            "meshbin": str(args.meshbin),
            "mesh_format": mesh_format,
            "nxy_width": width,
            "nxy_height": height,
            "mesh_declared_width": mesh_width,
            "mesh_declared_height": mesh_height,
            "vertices": len(vertices),
            "faces": len(faces),
        }
        row.update(m)
        rows.append(row)

        name = f"{query.query_id}_{safe_name(query.label)}_s{query.sigma_pixels:g}"
        if args.dump_densities:
            np.savez_compressed(
                density_dir / f"{name}.npz",
                gt_hist=gt_hist.astype(np.float32),
                mesh_hist=mesh_hist.astype(np.float32),
                gt_prob=normalize_density(gt_hist).astype(np.float32),
                mesh_prob=normalize_density(mesh_hist).astype(np.float32),
                u=np.float32(query.u),
                v=np.float32(query.v),
                sigma_pixels=np.float32(query.sigma_pixels),
            )
        if not args.no_images:
            save_density_comparison(image_dir / f"{name}.png", gt_hist, mesh_hist, query, title_prefix)

        print(
            f"[{qi+1:03d}/{len(queries):03d}] {query.query_id} {query.label} "
            f"sigma={query.sigma_pixels:g} TV={m['tv_distance']:.5f} "
            f"H={m['hellinger']:.5f} peak={m['peak_location_error_nxy']:.5f}",
            flush=True,
        )

    with metrics_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames: list[str] = []
        for r in rows:
            for k in r.keys():
                if k not in fieldnames:
                    fieldnames.append(k)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Compact summary grouped by sigma for quick inspection.
    summary_path = args.out_dir / "summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write(f"nxy_bin: {args.nxy_bin}\n")
        f.write(f"meshbin: {args.meshbin}\n")
        f.write(f"resolution: {width} x {height}\n")
        f.write(f"vertices: {len(vertices)}\n")
        f.write(f"faces: {len(faces)}\n")
        f.write(f"coverage_ratio: {coverage_ratio:.8g}\n")
        f.write(f"queries: {len(queries)}\n")
        for sigma in sorted({float(r["sigma_pixels"]) for r in rows}):
            sub = [r for r in rows if float(r["sigma_pixels"]) == sigma]
            f.write(
                f"sigma={sigma:g}: "
                f"mean_TV={np.mean([float(r['tv_distance']) for r in sub]):.6g}, "
                f"mean_H={np.mean([float(r['hellinger']) for r in sub]):.6g}, "
                f"mean_peak_err={np.mean([float(r['peak_location_error_nxy']) for r in sub]):.6g}\n"
            )

    print(f"Wrote {metrics_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
