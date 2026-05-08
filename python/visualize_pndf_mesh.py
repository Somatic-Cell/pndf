#!/usr/bin/env python3
"""
Visualize PNDF/QEM normal texture meshes.

This is a standalone helper script.  It reads:

  - nxy.bin produced by python/prepare_nxy.py
      format A: int32 N, float nx[N*N], float ny[N*N]
      format B: int32 W, int32 H, float nx[W*H], float ny[W*H]

  - meshbin produced by write_mesh_binary
      current format: int32 width, int32 height, int32 V, int32 F,
                      V * (float u, v, nx, ny), F * (int32 i0, i1, i2)
      legacy format:  int32 N, int32 V, int32 F,
                      V * (float u, v, nx, ny), F * (int32 i0, i1, i2)

The mesh UVs are expected in [0,1).  If they appear to be pixel coordinates,
the script normalizes them using mesh width/height.

Outputs:
  01_gt_normal_full.png
  02_gt_with_mesh_overlay_full.png
  03_gt_normal_crop.png
  04_gt_with_mesh_overlay_crop.png
  optional reconstruction and error images.

Examples:

  python python/visualize_pndf_mesh.py `
    --nxy-bin out256_normal/scratched_ds4_nxy.bin `
    --meshbin out256_normal/scratched_normal_v16384.meshbin `
    --out-dir vis/scratched_v16384 `
    --crop auto `
    --crop-size 0.15 `
    --crop-pixels 2048

  python python/visualize_pndf_mesh.py `
    --nxy-bin out256_normal/scratched_ds4_nxy.bin `
    --meshbin out256_normal/scratched_normal_v16384.meshbin `
    --out-dir vis/scratched_v16384 `
    --crop auto `
    --crop-size 0.15 `
    --crop-pixels 2048 `
    --render-reconstruction `
    --error-heatmap
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

EPS = 1e-12


def clamp_nxy(nxy: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(nxy, axis=-1, keepdims=True)
    return (nxy * np.minimum(1.0, 0.999 / np.maximum(n, EPS))).astype(np.float32)


def nxy_to_unit(nxy: np.ndarray) -> np.ndarray:
    z = np.sqrt(np.maximum(0.0, 1.0 - np.sum(nxy * nxy, axis=-1, keepdims=True)))
    n = np.concatenate([nxy, z], axis=-1)
    n /= np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), EPS)
    return n.astype(np.float32)


def normal_rgb(nxy: np.ndarray) -> np.ndarray:
    return np.clip(0.5 * (nxy_to_unit(nxy) + 1.0), 0.0, 1.0)


def angular_error_deg(gt: np.ndarray, rec: np.ndarray) -> np.ndarray:
    d = np.sum(nxy_to_unit(gt) * nxy_to_unit(rec), axis=-1).clip(-1.0, 1.0)
    return np.degrees(np.arccos(d)).astype(np.float32)


def wrap_delta(d: np.ndarray) -> np.ndarray:
    """Map coordinate differences to [-0.5, 0.5)."""
    return (d + 0.5) % 1.0 - 0.5


def load_nxy_bin(path: Path) -> np.ndarray:
    data = path.read_bytes()
    if len(data) < 4:
        raise ValueError(f"nxy file too small: {path}")
    ints = np.frombuffer(data[:8], dtype=np.int32, count=min(2, len(data)//4))

    # Current prepare_nxy.py format: int32 N, float nx[N*N], float ny[N*N]
    N = int(ints[0])
    size_single_header = 4 + 2 * N * N * 4
    if len(data) == size_single_header:
        arr = np.frombuffer(data[4:], dtype=np.float32)
        nx = arr[:N*N].reshape(N, N)
        ny = arr[N*N:].reshape(N, N)
        return clamp_nxy(np.stack([nx, ny], axis=-1).copy())

    # Optional format: int32 W, int32 H, float nx[W*H], float ny[W*H]
    if len(ints) >= 2:
        W, H = int(ints[0]), int(ints[1])
        size_two_header = 8 + 2 * W * H * 4
        if W > 0 and H > 0 and len(data) == size_two_header:
            arr = np.frombuffer(data[8:], dtype=np.float32)
            nx = arr[:W*H].reshape(H, W)
            ny = arr[W*H:].reshape(H, W)
            return clamp_nxy(np.stack([nx, ny], axis=-1).copy())

    raise ValueError(
        f"Could not parse nxy binary {path}. "
        f"First int={N}, file bytes={len(data)}."
    )


def load_meshbin(path: Path) -> tuple[int, int, np.ndarray, np.ndarray, str]:
    data = path.read_bytes()
    if len(data) < 12:
        raise ValueError(f"mesh file too small: {path}")

    # Prefer current format: int32 width, int32 height, int32 V, int32 F
    if len(data) >= 16:
        width, height, V, F = np.frombuffer(data[:16], dtype=np.int32, count=4)
        width, height, V, F = int(width), int(height), int(V), int(F)
        expected = 16 + 16 * V + 12 * F
        if width > 0 and height > 0 and V >= 0 and F >= 0 and expected == len(data):
            off = 16
            vertices = np.frombuffer(data[off:off+16*V], dtype=np.float32).reshape(V, 4).copy()
            off += 16 * V
            faces = np.frombuffer(data[off:off+12*F], dtype=np.int32).reshape(F, 3).copy()
            vertices = normalize_mesh_uv(vertices, width, height)
            return width, height, vertices, faces, "width,height,V,F"

    # Legacy format: int32 N, int32 V, int32 F
    N, V, F = np.frombuffer(data[:12], dtype=np.int32, count=3)
    N, V, F = int(N), int(V), int(F)
    expected = 12 + 16 * V + 12 * F
    if N > 0 and V >= 0 and F >= 0 and expected == len(data):
        off = 12
        vertices = np.frombuffer(data[off:off+16*V], dtype=np.float32).reshape(V, 4).copy()
        off += 16 * V
        faces = np.frombuffer(data[off:off+12*F], dtype=np.int32).reshape(F, 3).copy()
        vertices = normalize_mesh_uv(vertices, N, N)
        return N, N, vertices, faces, "N,V,F"

    raise ValueError(
        f"Could not parse mesh binary {path}. file bytes={len(data)}. "
        f"Tried current and legacy layouts."
    )


def normalize_mesh_uv(vertices: np.ndarray, width: int, height: int) -> np.ndarray:
    """Accept uv in [0,1) or pixel coordinates; normalize to [0,1)."""
    uv = vertices[:, :2]
    # Heuristic: current TorusMesh stores [0,1). Some earlier prototypes stored [0,N).
    if np.nanmax(np.abs(uv)) > 2.0:
        vertices[:, 0] = vertices[:, 0] / float(width)
        vertices[:, 1] = vertices[:, 1] / float(height)
    vertices[:, 0:2] = vertices[:, 0:2] % 1.0
    vertices[:, 2:4] = clamp_nxy(vertices[:, 2:4])
    return vertices


def sample_periodic_bilinear(nxy: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    H, W = nxy.shape[:2]
    x = (u % 1.0) * W - 0.5
    y = (v % 1.0) * H - 0.5
    x0 = np.floor(x).astype(np.int64)
    y0 = np.floor(y).astype(np.int64)
    tx = (x - x0).astype(np.float32)
    ty = (y - y0).astype(np.float32)
    x0 %= W
    y0 %= H
    x1 = (x0 + 1) % W
    y1 = (y0 + 1) % H
    v00 = nxy[y0, x0]
    v10 = nxy[y0, x1]
    v01 = nxy[y1, x0]
    v11 = nxy[y1, x1]
    out = ((1-tx)[..., None]*(1-ty)[..., None]*v00 +
           tx[..., None]*(1-ty)[..., None]*v10 +
           (1-tx)[..., None]*ty[..., None]*v01 +
           tx[..., None]*ty[..., None]*v11)
    return clamp_nxy(out)


def crop_normal_image(nxy: np.ndarray, u0: float, v0: float, size: float, pixels: int) -> np.ndarray:
    xs = u0 + (np.arange(pixels, dtype=np.float32) + 0.5) / pixels * size
    ys = v0 + (np.arange(pixels, dtype=np.float32) + 0.5) / pixels * size
    uu, vv = np.meshgrid(xs, ys, indexing="xy")
    return sample_periodic_bilinear(nxy, uu, vv)


def choose_auto_crop(nxy: np.ndarray, size: float, block: int = 16) -> tuple[float, float]:
    """Choose a crop around strong normal variation. Returns top-left u0,v0."""
    n = nxy_to_unit(nxy)
    # Local angular-ish variation proxy.
    dx = np.linalg.norm(n - np.roll(n, 1, axis=1), axis=-1) + np.linalg.norm(n - np.roll(n, -1, axis=1), axis=-1)
    dy = np.linalg.norm(n - np.roll(n, 1, axis=0), axis=-1) + np.linalg.norm(n - np.roll(n, -1, axis=0), axis=-1)
    score = dx + dy
    H, W = score.shape
    # Periodic box smoothing using roll accumulation over a modest window.
    win = max(3, int(round(size * min(W, H) / 4)))
    sm = np.zeros_like(score, dtype=np.float32)
    count = 0
    for oy in range(-win, win + 1, max(1, win // 8)):
        for ox in range(-win, win + 1, max(1, win // 8)):
            sm += np.roll(np.roll(score, oy, axis=0), ox, axis=1)
            count += 1
    sm /= max(count, 1)
    y, x = np.unravel_index(np.argmax(sm), sm.shape)
    u_center = (x + 0.5) / W
    v_center = (y + 0.5) / H
    return (u_center - 0.5 * size) % 1.0, (v_center - 0.5 * size) % 1.0


def parse_crop_arg(crop: str, nxy: np.ndarray, crop_size: float) -> tuple[float, float, float]:
    if crop.lower() == "auto":
        u0, v0 = choose_auto_crop(nxy, crop_size)
        return u0, v0, crop_size
    parts = [float(x) for x in crop.split(",")]
    if len(parts) != 3:
        raise ValueError("--crop must be 'auto' or 'u0,v0,size'")
    u0, v0, size = parts
    if size <= 0 or size > 1:
        raise ValueError("crop size must be in (0,1]")
    return u0 % 1.0, v0 % 1.0, size


def unique_edges(faces: np.ndarray) -> np.ndarray:
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0)
    e = np.sort(e.astype(np.int64), axis=1)
    e = np.unique(e, axis=0)
    return e


def edge_segments_full(uv: np.ndarray, faces: np.ndarray, max_edges: int, seed: int = 1) -> list[np.ndarray]:
    edges = unique_edges(faces)
    if max_edges > 0 and len(edges) > max_edges:
        rng = np.random.default_rng(seed)
        edges = edges[rng.choice(len(edges), size=max_edges, replace=False)]
    segs = []
    for a, b in edges:
        p0 = uv[a]
        d = wrap_delta(uv[b] - p0)
        q0 = p0
        q1 = p0 + d
        # Draw shifted copies that intersect [0,1]^2 so seam-crossing edges are visible without long lines.
        for sy in (-1, 0, 1):
            for sx in (-1, 0, 1):
                s0 = q0 + np.array([sx, sy], dtype=np.float32)
                s1 = q1 + np.array([sx, sy], dtype=np.float32)
                xmin, xmax = min(s0[0], s1[0]), max(s0[0], s1[0])
                ymin, ymax = min(s0[1], s1[1]), max(s0[1], s1[1])
                if xmax >= 0 and xmin <= 1 and ymax >= 0 and ymin <= 1:
                    segs.append(np.stack([s0, s1], axis=0))
    return segs


def edge_segments_crop(uv: np.ndarray, faces: np.ndarray, u0: float, v0: float, size: float, max_edges: int, seed: int = 1) -> list[np.ndarray]:
    edges = unique_edges(faces)
    if max_edges > 0 and len(edges) > max_edges:
        rng = np.random.default_rng(seed)
        edges = edges[rng.choice(len(edges), size=max_edges, replace=False)]
    center = np.array([(u0 + 0.5 * size) % 1.0, (v0 + 0.5 * size) % 1.0], dtype=np.float32)
    half = 0.5 * size
    segs = []
    for a, b in edges:
        pa = uv[a]
        pb = uv[b]
        p0 = wrap_delta(pa - center)          # centered local coords in [-0.5,0.5)
        d = wrap_delta(pb - pa)
        p1 = p0 + d
        xmin, xmax = min(p0[0], p1[0]), max(p0[0], p1[0])
        ymin, ymax = min(p0[1], p1[1]), max(p0[1], p1[1])
        if xmax >= -half and xmin <= half and ymax >= -half and ymin <= half:
            # local crop coordinates [0,1]
            s0 = np.array([(p0[0] + half) / size, (p0[1] + half) / size], dtype=np.float32)
            s1 = np.array([(p1[0] + half) / size, (p1[1] + half) / size], dtype=np.float32)
            segs.append(np.stack([s0, s1], axis=0))
    return segs

def parse_uv_shift_arg(s: str) -> np.ndarray:
    parts = [float(x) for x in s.split(",")]
    if len(parts) != 2:
        raise ValueError("--uv-shift must be 'du,dv'")
    return np.array([parts[0], parts[1]], dtype=np.float32)


def view_shift_uv(uv: np.ndarray, shift: np.ndarray) -> np.ndarray:
    """Change only the visualization cut.

    If shift=(0.5,0.5), the displayed [0,1)^2 window corresponds to
    original coordinates [0.5,1.5)^2 modulo 1.
    """
    return (uv - shift[None, :]) % 1.0


def shifted_nxy_for_view(nxy: np.ndarray, shift: np.ndarray) -> np.ndarray:
    """Return normal image for the shifted visualization cut.

    Pixel at displayed coordinate q samples original texture at q + shift.
    """
    H, W = nxy.shape[:2]
    xs = (np.arange(W, dtype=np.float32) + 0.5) / W + float(shift[0])
    ys = (np.arange(H, dtype=np.float32) + 0.5) / H + float(shift[1])
    uu, vv = np.meshgrid(xs, ys, indexing="xy")
    return sample_periodic_bilinear(nxy, uu, vv)


def save_vertex_density(
    path: Path,
    uv: np.ndarray,
    title: str,
    bins: int = 128,
    dpi: int = 180,
):
    """Visualize vertex density in the displayed fundamental domain."""
    uv = uv % 1.0
    H, xedges, yedges = np.histogram2d(
        uv[:, 1],
        uv[:, 0],
        bins=bins,
        range=[[0.0, 1.0], [0.0, 1.0]],
    )

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    im = ax.imshow(
        np.log1p(H),
        origin="upper",
        extent=(0, 1, 1, 0),
        cmap="magma",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="log(1 + vertex count)")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def seam_density_stats(
    uv: np.ndarray,
    faces: np.ndarray,
    width: int,
    height: int,
    band_pixels: int,
    shift: np.ndarray,
) -> dict:
    """Measure whether vertices/edges are concentrated near the displayed cut.

    This does not change the mesh. It changes only the fundamental-domain cut.
    """
    view_uv = view_shift_uv(uv, shift)
    bw = float(band_pixels) / float(width)
    bh = float(band_pixels) / float(height)

    # Clamp for pathological tiny resolutions.
    bw = min(max(bw, 0.0), 0.5)
    bh = min(max(bh, 0.0), 0.5)

    in_x_band = (view_uv[:, 0] < bw) | (view_uv[:, 0] >= 1.0 - bw)
    in_y_band = (view_uv[:, 1] < bh) | (view_uv[:, 1] >= 1.0 - bh)
    in_band = in_x_band | in_y_band

    vertex_band_count = int(np.count_nonzero(in_band))
    vertex_count = int(len(view_uv))
    vertex_band_ratio = vertex_band_count / max(vertex_count, 1)

    # Area fraction of the union of x/y seam bands.
    expected_area_ratio = 1.0 - max(0.0, 1.0 - 2.0 * bw) * max(0.0, 1.0 - 2.0 * bh)

    edges = unique_edges(faces)
    if len(edges) > 0:
        p0 = view_uv[edges[:, 0]]
        p1 = view_uv[edges[:, 1]]
        raw_d = np.abs(p1 - p0)
        seam_crossing = (raw_d[:, 0] > 0.5) | (raw_d[:, 1] > 0.5)

        e_in_band = in_band[edges[:, 0]] | in_band[edges[:, 1]]
        edge_band_count = int(np.count_nonzero(e_in_band))
        seam_crossing_count = int(np.count_nonzero(seam_crossing))
    else:
        edge_band_count = 0
        seam_crossing_count = 0

    edge_count = int(len(edges))
    return {
        "shift_u": float(shift[0]),
        "shift_v": float(shift[1]),
        "band_pixels": int(band_pixels),
        "band_width_u": float(bw),
        "band_width_v": float(bh),
        "vertex_count": vertex_count,
        "vertex_band_count": vertex_band_count,
        "vertex_band_ratio": float(vertex_band_ratio),
        "expected_band_area_ratio": float(expected_area_ratio),
        "vertex_band_over_expected": float(vertex_band_ratio / max(expected_area_ratio, 1e-12)),
        "edge_count": edge_count,
        "edge_band_count": edge_band_count,
        "edge_band_ratio": float(edge_band_count / max(edge_count, 1)),
        "seam_crossing_edge_count": seam_crossing_count,
        "seam_crossing_edge_ratio": float(seam_crossing_count / max(edge_count, 1)),
    }


def format_seam_stats(prefix: str, stats: dict) -> list[str]:
    lines = [f"{prefix}:"]
    for k, v in stats.items():
        if isinstance(v, float):
            lines.append(f"{prefix}.{k}: {v:.8g}")
        else:
            lines.append(f"{prefix}.{k}: {v}")
    return lines


def save_image(path: Path, image: np.ndarray, title: Optional[str] = None, dpi: int = 180):
    h, w = image.shape[:2]
    fig_w = min(10.0, max(4.0, w / 160.0))
    fig_h = fig_w * h / w
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    ax.imshow(image, origin="upper")
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def save_overlay(path: Path, image: np.ndarray, segs: list[np.ndarray], title: str, linewidth: float, alpha: float, dpi: int = 200):
    h, w = image.shape[:2]
    fig_w = min(12.0, max(5.0, w / 180.0))
    fig_h = fig_w * h / w
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    ax.imshow(image, origin="upper", extent=(0, 1, 1, 0))
    if segs:
        lc = LineCollection(segs, colors=(0, 0, 0, alpha), linewidths=linewidth)
        ax.add_collection(lc)
    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def save_error_heatmap(path: Path, err: np.ndarray, title: str, dpi: int = 180, vmax: Optional[float] = None):
    if vmax is None:
        vmax = max(1.0, float(np.percentile(err, 99.5)))
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    im = ax.imshow(err, origin="upper", cmap="magma", vmin=0.0, vmax=vmax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="angular error [deg]")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def rasterize_mesh_to_grid(vertices: np.ndarray, faces: np.ndarray, H: int, W: int) -> np.ndarray:
    """Slow but robust CPU rasterization of periodic PL mesh to texel centers.

    Intended for visualization/evaluation at N=128/256/512.  For large 1024 meshes with
    many triangles, use with care.
    """
    uv = vertices[:, :2].astype(np.float64) % 1.0
    vals = vertices[:, 2:4].astype(np.float64)
    out = np.full((H, W, 2), np.nan, dtype=np.float32)

    for tri in faces:
        ids = tri.astype(np.int64)
        p0 = uv[ids[0]]
        p1 = p0 + wrap_delta(uv[ids[1]] - p0)
        p2 = p0 + wrap_delta(uv[ids[2]] - p0)
        vv = vals[ids]
        base_tri = np.stack([p0, p1, p2], axis=0)
        # Draw periodic copies that touch [0,1]^2.
        for sy in (-1, 0, 1):
            for sx in (-1, 0, 1):
                p = base_tri + np.array([sx, sy], dtype=np.float64)
                xmin, xmax = p[:, 0].min(), p[:, 0].max()
                ymin, ymax = p[:, 1].min(), p[:, 1].max()
                if xmax < 0 or xmin > 1 or ymax < 0 or ymin > 1:
                    continue
                ix0 = max(0, int(math.floor(xmin * W - 1)))
                ix1 = min(W - 1, int(math.ceil(xmax * W + 1)))
                iy0 = max(0, int(math.floor(ymin * H - 1)))
                iy1 = min(H - 1, int(math.ceil(ymax * H + 1)))
                if ix1 < ix0 or iy1 < iy0:
                    continue
                xs = (np.arange(ix0, ix1 + 1, dtype=np.float64) + 0.5) / W
                ys = (np.arange(iy0, iy1 + 1, dtype=np.float64) + 0.5) / H
                xx, yy = np.meshgrid(xs, ys, indexing="xy")
                q = np.stack([xx.ravel(), yy.ravel()], axis=-1)
                e0 = p[1] - p[0]
                e1 = p[2] - p[0]
                den = e0[0] * e1[1] - e0[1] * e1[0]
                if abs(den) < 1e-14:
                    continue
                e2 = q - p[0]
                a = (e2[:, 0] * e1[1] - e2[:, 1] * e1[0]) / den
                b = (e0[0] * e2[:, 1] - e0[1] * e2[:, 0]) / den
                c = 1.0 - a - b
                inside = (a >= -1e-7) & (b >= -1e-7) & (c >= -1e-7)
                if not np.any(inside):
                    continue
                interp = c[:, None] * vv[0] + a[:, None] * vv[1] + b[:, None] * vv[2]
                coords = q[inside]
                px = np.clip(np.floor(coords[:, 0] * W).astype(np.int64), 0, W - 1)
                py = np.clip(np.floor(coords[:, 1] * H).astype(np.int64), 0, H - 1)
                out[py, px] = interp[inside].astype(np.float32)

    missing = ~np.isfinite(out).all(axis=-1)
    if np.any(missing):
        # Fallback nearest vertex for any tiny holes.
        from scipy.spatial import cKDTree  # optional; only needed for fallback
        pts = uv.copy()
        # periodic replicated KD fallback
        shifts = np.array([(i, j) for j in (-1, 0, 1) for i in (-1, 0, 1)], dtype=np.float64)
        pts_aug = (pts[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
        vals_aug = np.tile(vertices[:, 2:4], (9, 1))
        tree = cKDTree(pts_aug)
        yy, xx = np.nonzero(missing)
        q = np.stack([(xx + 0.5) / W, (yy + 0.5) / H], axis=-1)
        _, nn = tree.query(q)
        out[yy, xx] = vals_aug[nn]
    return clamp_nxy(out)


def main():
    ap = argparse.ArgumentParser(description="Visualize PNDF/QEM mesh normal maps and connectivity.")
    ap.add_argument("--nxy-bin", required=True, type=Path)
    ap.add_argument("--meshbin", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--crop", default="auto", help="auto or u0,v0,size in normalized coordinates")
    ap.add_argument("--crop-size", type=float, default=0.15, help="crop size for --crop auto, in [0,1]")
    ap.add_argument("--crop-pixels", type=int, default=2048, help="output crop resolution")
    ap.add_argument("--full-max-edges", type=int, default=120000, help="max edges drawn on full overlay; <=0 draws all")
    ap.add_argument("--crop-max-edges", type=int, default=400000, help="max edges considered for crop overlay; <=0 draws all")
    ap.add_argument("--line-width", type=float, default=0.25)
    ap.add_argument("--crop-line-width", type=float, default=0.35)
    ap.add_argument("--line-alpha", type=float, default=0.55)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--render-reconstruction", action="store_true")
    ap.add_argument("--error-heatmap", action="store_true")
    ap.add_argument("--skip-full-overlay", action="store_true", help="skip full mesh overlay if mesh is huge")
    ap.add_argument(
        "--diagnose-seam",
        action="store_true",
        help="write seam/cut diagnostics: shifted-cut overlay, vertex-density maps, and seam-band stats",
    )
    ap.add_argument(
        "--uv-shift",
        default="0.5,0.5",
        help="view-space shift for seam diagnostics, e.g. 0.5,0.5. "
             "This changes only the visualization cut, not the mesh.",
    )
    ap.add_argument(
        "--seam-band-pixels",
        type=int,
        default=4,
        help="band width in texels used for seam density diagnostics",
    )
    ap.add_argument(
        "--density-bins",
        type=int,
        default=128,
        help="resolution of vertex-density diagnostic heatmap",
    )
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    nxy = load_nxy_bin(args.nxy_bin)
    H, W = nxy.shape[:2]
    width, height, vertices, faces, mesh_format = load_meshbin(args.meshbin)
    uv = vertices[:, :2]

    diag_shift = parse_uv_shift_arg(args.uv_shift)
    diag_summary_lines: list[str] = []

    u0, v0, crop_size = parse_crop_arg(args.crop, nxy, args.crop_size)
    gt_full = normal_rgb(nxy)
    gt_crop_nxy = crop_normal_image(nxy, u0, v0, crop_size, args.crop_pixels)
    gt_crop = normal_rgb(gt_crop_nxy)

    save_image(args.out_dir / "01_gt_normal_full.png", gt_full, title="GT normal map")

    if not args.skip_full_overlay:
        segs_full = edge_segments_full(uv, faces, max_edges=args.full_max_edges, seed=args.seed)
        save_overlay(
            args.out_dir / "02_gt_with_mesh_overlay_full.png",
            gt_full,
            segs_full,
            title=f"GT + mesh overlay (V={len(vertices)}, F={len(faces)}, edges drawn={len(segs_full)})",
            linewidth=args.line_width,
            alpha=args.line_alpha,
            dpi=220,
        )
    
    if args.diagnose_seam:
        # 1) Original cut diagnostics.
        stats0 = seam_density_stats(
            uv=uv,
            faces=faces,
            width=W,
            height=H,
            band_pixels=args.seam_band_pixels,
            shift=np.array([0.0, 0.0], dtype=np.float32),
        )
        diag_summary_lines.extend(format_seam_stats("seam_diag_original_cut", stats0))

        save_vertex_density(
            args.out_dir / "11_vertex_density_original_cut.png",
            uv,
            title="Vertex density, original cut",
            bins=args.density_bins,
            dpi=200,
        )

        # 2) Shifted cut diagnostics.
        view_uv = view_shift_uv(uv, diag_shift)
        shifted_nxy = shifted_nxy_for_view(nxy, diag_shift)
        shifted_rgb = normal_rgb(shifted_nxy)

        stats_shift = seam_density_stats(
            uv=uv,
            faces=faces,
            width=W,
            height=H,
            band_pixels=args.seam_band_pixels,
            shift=diag_shift,
        )
        diag_summary_lines.extend(format_seam_stats("seam_diag_shifted_cut", stats_shift))

        save_vertex_density(
            args.out_dir / "12_vertex_density_shifted_cut.png",
            view_uv,
            title=f"Vertex density, shifted cut du={diag_shift[0]:.3f}, dv={diag_shift[1]:.3f}",
            bins=args.density_bins,
            dpi=200,
        )

        if not args.skip_full_overlay:
            segs_shifted = edge_segments_full(
                view_uv,
                faces,
                max_edges=args.full_max_edges,
                seed=args.seed,
            )
            save_overlay(
                args.out_dir / "13_gt_with_mesh_overlay_shifted_cut.png",
                shifted_rgb,
                segs_shifted,
                title=(
                    f"GT + mesh overlay, shifted cut "
                    f"du={diag_shift[0]:.3f}, dv={diag_shift[1]:.3f} "
                    f"(V={len(vertices)}, F={len(faces)}, edges drawn={len(segs_shifted)})"
                ),
                linewidth=args.line_width,
                alpha=args.line_alpha,
                dpi=220,
            )

    save_image(args.out_dir / "03_gt_normal_crop.png", gt_crop, title=f"GT crop u0={u0:.4f}, v0={v0:.4f}, size={crop_size:.4f}", dpi=220)

    segs_crop = edge_segments_crop(uv, faces, u0, v0, crop_size, max_edges=args.crop_max_edges, seed=args.seed)
    save_overlay(
        args.out_dir / "04_gt_with_mesh_overlay_crop.png",
        gt_crop,
        segs_crop,
        title=f"GT crop + mesh overlay (edges drawn={len(segs_crop)})",
        linewidth=args.crop_line_width,
        alpha=args.line_alpha,
        dpi=260,
    )

    rec = None
    if args.render_reconstruction or args.error_heatmap:
        rec = rasterize_mesh_to_grid(vertices, faces, H, W)
        rec_full = normal_rgb(rec)
        rec_crop_nxy = crop_normal_image(rec, u0, v0, crop_size, args.crop_pixels)
        rec_crop = normal_rgb(rec_crop_nxy)
        save_image(args.out_dir / "05_reconstructed_normal_full.png", rec_full, title="Reconstructed normal map")
        if not args.skip_full_overlay:
            segs_full = edge_segments_full(uv, faces, max_edges=args.full_max_edges, seed=args.seed)
            save_overlay(
                args.out_dir / "06_reconstructed_with_mesh_overlay_full.png",
                rec_full,
                segs_full,
                title="Reconstructed + mesh overlay",
                linewidth=args.line_width,
                alpha=args.line_alpha,
                dpi=220,
            )
        save_image(args.out_dir / "07_reconstructed_normal_crop.png", rec_crop, title="Reconstructed crop", dpi=220)
        save_overlay(
            args.out_dir / "08_reconstructed_with_mesh_overlay_crop.png",
            rec_crop,
            segs_crop,
            title="Reconstructed crop + mesh overlay",
            linewidth=args.crop_line_width,
            alpha=args.line_alpha,
            dpi=260,
        )

    if args.error_heatmap:
        if rec is None:
            rec = rasterize_mesh_to_grid(vertices, faces, H, W)
        err = angular_error_deg(nxy, rec)
        save_error_heatmap(args.out_dir / "09_angular_error_heatmap_full.png", err, title="Angular error full")
        err_crop = angular_error_deg(gt_crop_nxy, crop_normal_image(rec, u0, v0, crop_size, args.crop_pixels))
        save_error_heatmap(args.out_dir / "10_angular_error_heatmap_crop.png", err_crop, title="Angular error crop")
        err_stats = {
            "mean_deg": float(err.mean()),
            "p95_deg": float(np.percentile(err, 95)),
            "p99_deg": float(np.percentile(err, 99)),
            "max_deg": float(err.max()),
        }
    else:
        err_stats = {}

    summary = []
    summary.append(f"nxy_bin: {args.nxy_bin}")
    summary.append(f"meshbin: {args.meshbin}")
    summary.append(f"mesh_format: {mesh_format}")
    summary.append(f"nxy_resolution: {W} x {H}")
    summary.append(f"mesh_declared_resolution: {width} x {height}")
    summary.append(f"vertices: {len(vertices)}")
    summary.append(f"faces: {len(faces)}")
    summary.append(f"crop: u0={u0:.8f}, v0={v0:.8f}, size={crop_size:.8f}, pixels={args.crop_pixels}")
    summary.append(f"full_overlay_max_edges: {args.full_max_edges}")
    summary.append(f"crop_overlay_edges_drawn: {len(segs_crop)}")
    if args.diagnose_seam:
        summary.extend(diag_summary_lines)
    for k, v in err_stats.items():
        summary.append(f"{k}: {v}")
    (args.out_dir / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    print(f"wrote outputs to {args.out_dir}")


if __name__ == "__main__":
    main()
