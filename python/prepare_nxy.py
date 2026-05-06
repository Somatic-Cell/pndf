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
    norm = np.linalg.norm(nxy, axis=-1, keepdims=True)
    nxy = nxy * np.minimum(1.0, 0.999 / np.maximum(norm, EPS))
    return nxy.astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", choices=sorted(KIND), required=True)
    ap.add_argument("--data-dir", default=".")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    nxy = load_nxy(args.map, Path(args.data_dir))
    h, w = nxy.shape[:2]
    if h != w:
        raise ValueError("Phase 1 binary writer expects square maps")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as f:
        f.write(np.int32(w).tobytes())
        f.write(np.ascontiguousarray(nxy[..., 0]).astype(np.float32).tobytes())
        f.write(np.ascontiguousarray(nxy[..., 1]).astype(np.float32).tobytes())
    print(f"wrote {out} N={w} vertices={w*w} triangles={2*w*w}")


if __name__ == "__main__":
    main()
