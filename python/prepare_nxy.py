#!/usr/bin/env python3
"""Prepare C++ nxy binary input from the demo EXR maps.

Binary layout:
  int32 N
  float nx[N*N]
  float ny[N*N]
"""
import argparse
import os
from pathlib import Path

os.environ.setdefault("OPENCV_IO_ENABLE_OPENEXR", "1")
import cv2
import numpy as np

EPS = 1e-12
KIND = {"isotropic": "height", "scratched": "height", "brush": "normal"}


def clamp_nxy(nxy: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(nxy, axis=-1, keepdims=True)
    return (nxy * np.minimum(1.0, 0.999 / np.maximum(norm, EPS))).astype(np.float32)


def periodic_block_average(nxy: np.ndarray, factor: int) -> np.ndarray:
    if factor <= 1:
        return clamp_nxy(nxy)
    h, w = nxy.shape[:2]
    if h % factor != 0 or w % factor != 0:
        raise ValueError(f"downsample factor {factor} must divide map size {w}x{h}")
    # Full-texture block average preserves the torus domain.  This is not a crop:
    # no new boundary is introduced, and the resulting nxy.bin remains periodic.
    out = nxy.reshape(h // factor, factor, w // factor, factor, 2).mean(axis=(1, 3))
    return clamp_nxy(out)


def load_nxy(map_name: str, data_dir: Path) -> np.ndarray:
    path = data_dir / f"{map_name}.exr"
    img = cv2.imread(str(path), cv2.IMREAD_UNCHANGED)
    if img is None:
        raise FileNotFoundError(path)
    img = img.astype(np.float32)
    if KIND[map_name] == "height":
        h = img[..., 0] if img.ndim == 3 else img
        nx = 0.5 * (np.roll(h, -1, axis=1) - np.roll(h, 1, axis=1))
        ny = 0.5 * (np.roll(h, -1, axis=0) - np.roll(h, 1, axis=0))
        nxy = np.stack([nx, ny], axis=-1)
    else:
        # OpenCV loads BGR(A); projected normal is RGB.RG = BGR indices [2, 1].
        nxy = img[..., [2, 1]].copy()
    return clamp_nxy(nxy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", choices=sorted(KIND), required=True)
    ap.add_argument("--data-dir", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--downsample", type=int, default=1,
                    help="Periodically downsample the full projected-normal map by this integer factor. "
                         "This averages the whole torus, not a crop.")
    args = ap.parse_args()
    if args.downsample < 1:
        raise ValueError("--downsample must be >= 1")
    nxy_full = load_nxy(args.map, Path(args.data_dir))
    nxy = periodic_block_average(nxy_full, args.downsample)
    h, w = nxy.shape[:2]
    if h != w:
        raise ValueError("Phase 1 binary writer expects square maps")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as f:
        f.write(np.int32(w).tobytes())
        f.write(np.ascontiguousarray(nxy[..., 0]).astype(np.float32).tobytes())
        f.write(np.ascontiguousarray(nxy[..., 1]).astype(np.float32).tobytes())
    print(f"wrote {out} N={w} vertices={w*w} triangles={2*w*w} downsample={args.downsample}")


if __name__ == "__main__":
    main()
