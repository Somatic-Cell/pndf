#!/usr/bin/env python3
"""
Analyze and visualize localized UV embedding failures in PNDF meshbin files.

Expected meshbin layout:
  int32 width, int32 height, int32 V, int32 F
  V records of float32: u, v, nx, ny
  F records of int32: i0, i1, i2

This script is intentionally conservative and diagnostic-oriented. It localizes:
  - negative / tiny UV-area faces
  - non-incident edge intersections on the periodic [0,1)^2 torus
  - non-incident vertex-inside-triangle events on the periodic domain

It also outputs overview and zoom PNGs.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np


Vec = np.ndarray


def read_meshbin(path: Path):
    with open(path, "rb") as f:
        header = np.fromfile(f, dtype=np.int32, count=4)
        if header.size != 4:
            raise ValueError(f"Failed to read meshbin header: {path}")
        width, height, V, F = map(int, header)
        verts = np.fromfile(f, dtype=np.float32, count=V * 4)
        if verts.size != V * 4:
            raise ValueError(f"Failed to read {V} vertices from {path}")
        verts = verts.reshape((V, 4))
        faces = np.fromfile(f, dtype=np.int32, count=F * 3)
        if faces.size != F * 3:
            raise ValueError(f"Failed to read {F} faces from {path}")
        faces = faces.reshape((F, 3))
    return width, height, verts, faces


def unwrap_near(p: Vec, ref: Vec) -> Vec:
    q = np.array(p, dtype=np.float64)
    r = np.array(ref, dtype=np.float64)
    d = q - r
    d -= np.round(d)
    return r + d


def cross2(a: Vec, b: Vec) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def orient2d(a: Vec, b: Vec, c: Vec) -> float:
    return cross2(b - a, c - a)


def face_unwrapped_uv(uv: np.ndarray, tri: Sequence[int]) -> np.ndarray:
    p0 = np.array(uv[tri[0]], dtype=np.float64)
    p1 = unwrap_near(uv[tri[1]], p0)
    p2 = unwrap_near(uv[tri[2]], p0)
    return np.stack([p0, p1, p2], axis=0)


def edge_unwrapped_uv(uv: np.ndarray, i: int, j: int) -> Tuple[np.ndarray, np.ndarray]:
    p0 = np.array(uv[i], dtype=np.float64)
    p1 = unwrap_near(uv[j], p0)
    return p0, p1


def mod01(p: Vec) -> np.ndarray:
    return np.mod(p, 1.0)


@dataclass(frozen=True)
class Edge:
    eid: int
    i: int
    j: int


@dataclass
class EdgeCopy:
    eid: int
    i: int
    j: int
    p0: np.ndarray
    p1: np.ndarray
    shift: Tuple[int, int]


@dataclass
class FaceCopy:
    fid: int
    tri: Tuple[int, int, int]
    pts: np.ndarray
    shift: Tuple[int, int]


@dataclass
class EdgeIntersection:
    edge0: Tuple[int, int]
    edge1: Tuple[int, int]
    point: Tuple[float, float]
    raw_point: Tuple[float, float]
    eid0: int
    eid1: int


@dataclass
class VertexInside:
    vertex: int
    face: int
    point: Tuple[float, float]
    bary: Tuple[float, float, float]


def build_edges(faces: np.ndarray) -> List[Edge]:
    edge_map: Dict[Tuple[int, int], int] = {}
    edges: List[Edge] = []
    for tri in faces:
        a, b, c = map(int, tri)
        for x, y in ((a, b), (b, c), (c, a)):
            if x == y:
                continue
            key = (min(x, y), max(x, y))
            if key not in edge_map:
                eid = len(edges)
                edge_map[key] = eid
                edges.append(Edge(eid, key[0], key[1]))
    return edges


def bbox_overlaps_unit(pts: np.ndarray, margin: float = 1e-12) -> bool:
    mn = np.min(pts, axis=0)
    mx = np.max(pts, axis=0)
    return bool(mx[0] >= -margin and mx[1] >= -margin and mn[0] <= 1.0 + margin and mn[1] <= 1.0 + margin)


def make_edge_copies(uv: np.ndarray, edges: List[Edge]) -> List[EdgeCopy]:
    copies: List[EdgeCopy] = []
    shifts = [(sx, sy) for sx in (-1, 0, 1) for sy in (-1, 0, 1)]
    for e in edges:
        p0, p1 = edge_unwrapped_uv(uv, e.i, e.j)
        base = np.stack([p0, p1], axis=0)
        for sx, sy in shifts:
            s = np.array([sx, sy], dtype=np.float64)
            pts = base + s
            if bbox_overlaps_unit(pts):
                copies.append(EdgeCopy(e.eid, e.i, e.j, pts[0], pts[1], (sx, sy)))
    return copies


def make_face_copies(uv: np.ndarray, faces: np.ndarray) -> List[FaceCopy]:
    copies: List[FaceCopy] = []
    shifts = [(sx, sy) for sx in (-1, 0, 1) for sy in (-1, 0, 1)]
    for fid, tri_arr in enumerate(faces):
        tri = tuple(map(int, tri_arr))
        pts0 = face_unwrapped_uv(uv, tri)
        for sx, sy in shifts:
            s = np.array([sx, sy], dtype=np.float64)
            pts = pts0 + s
            if bbox_overlaps_unit(pts):
                copies.append(FaceCopy(fid, tri, pts, (sx, sy)))
    return copies


def cells_for_bbox(mn: np.ndarray, mx: np.ndarray, grid_n: int) -> Iterable[Tuple[int, int]]:
    mn = np.maximum(mn, 0.0)
    mx = np.minimum(mx, 1.0)
    if mx[0] < 0.0 or mx[1] < 0.0 or mn[0] > 1.0 or mn[1] > 1.0:
        return []
    ix0 = int(max(0, min(grid_n - 1, math.floor(mn[0] * grid_n))))
    iy0 = int(max(0, min(grid_n - 1, math.floor(mn[1] * grid_n))))
    ix1 = int(max(0, min(grid_n - 1, math.floor(mx[0] * grid_n))))
    iy1 = int(max(0, min(grid_n - 1, math.floor(mx[1] * grid_n))))
    return ((ix, iy) for ix in range(ix0, ix1 + 1) for iy in range(iy0, iy1 + 1))


def segment_grid(edge_copies: List[EdgeCopy], grid_n: int):
    grid: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for idx, ec in enumerate(edge_copies):
        pts = np.stack([ec.p0, ec.p1], axis=0)
        mn = np.min(pts, axis=0) - 1e-12
        mx = np.max(pts, axis=0) + 1e-12
        for cell in cells_for_bbox(mn, mx, grid_n):
            grid[cell].append(idx)
    return grid


def face_grid(face_copies: List[FaceCopy], grid_n: int):
    grid: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for idx, fc in enumerate(face_copies):
        mn = np.min(fc.pts, axis=0) - 1e-12
        mx = np.max(fc.pts, axis=0) + 1e-12
        for cell in cells_for_bbox(mn, mx, grid_n):
            grid[cell].append(idx)
    return grid


def strict_seg_intersect(a: Vec, b: Vec, c: Vec, d: Vec, eps: float = 1e-14) -> bool:
    o1 = orient2d(a, b, c)
    o2 = orient2d(a, b, d)
    o3 = orient2d(c, d, a)
    o4 = orient2d(c, d, b)
    return ((o1 > eps and o2 < -eps) or (o1 < -eps and o2 > eps)) and \
           ((o3 > eps and o4 < -eps) or (o3 < -eps and o4 > eps))


def line_intersection(a: Vec, b: Vec, c: Vec, d: Vec) -> np.ndarray:
    r = b - a
    s = d - c
    den = cross2(r, s)
    if abs(den) < 1e-30:
        return (a + b + c + d) * 0.25
    t = cross2(c - a, s) / den
    return a + t * r


def find_edge_intersections(edge_copies: List[EdgeCopy], grid_n: int, max_samples: int) -> List[EdgeIntersection]:
    grid = segment_grid(edge_copies, grid_n)
    seen_pairs: Set[Tuple[int, int]] = set()
    out: List[EdgeIntersection] = []
    for cell, indices in grid.items():
        n = len(indices)
        for aa in range(n):
            eca = edge_copies[indices[aa]]
            for bb in range(aa + 1, n):
                ecb = edge_copies[indices[bb]]
                if eca.eid == ecb.eid:
                    continue
                # Incident edges are allowed to meet at shared endpoints.
                if eca.i in (ecb.i, ecb.j) or eca.j in (ecb.i, ecb.j):
                    continue
                key = (min(eca.eid, ecb.eid), max(eca.eid, ecb.eid))
                if key in seen_pairs:
                    continue
                if strict_seg_intersect(eca.p0, eca.p1, ecb.p0, ecb.p1):
                    seen_pairs.add(key)
                    p = line_intersection(eca.p0, eca.p1, ecb.p0, ecb.p1)
                    pm = mod01(p)
                    out.append(EdgeIntersection(
                        edge0=(eca.i, eca.j), edge1=(ecb.i, ecb.j),
                        point=(float(pm[0]), float(pm[1])),
                        raw_point=(float(p[0]), float(p[1])),
                        eid0=eca.eid, eid1=ecb.eid,
                    ))
                    if max_samples > 0 and len(out) >= max_samples:
                        return out
    return out


def barycentric_strict(p: Vec, tri_pts: np.ndarray, eps: float = 1e-14) -> Optional[Tuple[float, float, float]]:
    a, b, c = tri_pts
    area = orient2d(a, b, c)
    if abs(area) < eps:
        return None
    w0 = orient2d(b, c, p) / area
    w1 = orient2d(c, a, p) / area
    w2 = orient2d(a, b, p) / area
    if w0 > eps and w1 > eps and w2 > eps:
        return (float(w0), float(w1), float(w2))
    return None


def find_vertex_inside(uv: np.ndarray, faces: np.ndarray, face_copies: List[FaceCopy], grid_n: int, max_samples: int) -> List[VertexInside]:
    grid = face_grid(face_copies, grid_n)
    out: List[VertexInside] = []
    seen: Set[Tuple[int, int]] = set()
    for vid, p in enumerate(uv):
        cell = (int(max(0, min(grid_n - 1, math.floor(float(p[0]) * grid_n)))),
                int(max(0, min(grid_n - 1, math.floor(float(p[1]) * grid_n)))))
        for idx in grid.get(cell, []):
            fc = face_copies[idx]
            if vid in fc.tri:
                continue
            key = (vid, fc.fid)
            if key in seen:
                continue
            bary = barycentric_strict(np.array(p, dtype=np.float64), fc.pts)
            if bary is not None:
                seen.add(key)
                out.append(VertexInside(vertex=vid, face=fc.fid,
                                        point=(float(p[0]), float(p[1])), bary=bary))
                if max_samples > 0 and len(out) >= max_samples:
                    return out
    return out


def analyze_bad_faces(uv: np.ndarray, faces: np.ndarray, area_eps: float):
    rows = []
    for fid, tri in enumerate(faces):
        pts = face_unwrapped_uv(uv, tuple(map(int, tri)))
        area2 = orient2d(pts[0], pts[1], pts[2])
        edges = [np.linalg.norm(pts[(k + 1) % 3] - pts[k]) for k in range(3)]
        rows.append({
            "face": int(fid),
            "i0": int(tri[0]), "i1": int(tri[1]), "i2": int(tri[2]),
            "area2": float(area2),
            "abs_area2": float(abs(area2)),
            "negative": bool(area2 < -area_eps),
            "tiny": bool(abs(area2) <= area_eps),
            "max_edge": float(max(edges)),
            "cx": float(np.mean(mod01(pts), axis=0)[0] % 1.0),
            "cy": float(np.mean(mod01(pts), axis=0)[1] % 1.0),
        })
    return rows


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def draw_overview(out_path: Path, uv: np.ndarray, faces: np.ndarray, edges: List[Edge],
                  intersections: List[EdgeIntersection], inside: List[VertexInside],
                  bad_faces: List[dict], max_draw_edges: int = 200000):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PolyCollection

    fig, ax = plt.subplots(figsize=(10, 10), dpi=180)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("UV embedding failure overview")

    edge_copies = make_edge_copies(uv, edges[:max_draw_edges])
    segs = []
    for ec in edge_copies:
        # Clip only by drawing all visible copies; Matplotlib will clip to axes.
        segs.append([ec.p0, ec.p1])
    if segs:
        lc = LineCollection(segs, linewidths=0.18, colors=[(0, 0, 0, 0.22)])
        ax.add_collection(lc)

    bad_polys = []
    for r in bad_faces:
        if r["negative"] or r["tiny"]:
            pts = face_unwrapped_uv(uv, faces[int(r["face"])])
            for sx in (-1, 0, 1):
                for sy in (-1, 0, 1):
                    pp = pts + np.array([sx, sy])
                    if bbox_overlaps_unit(pp):
                        bad_polys.append(pp)
    if bad_polys:
        pc = PolyCollection(bad_polys, facecolors=[(1.0, 0.65, 0.0, 0.35)], edgecolors=[(1.0, 0.35, 0.0, 0.9)], linewidths=0.6)
        ax.add_collection(pc)

    # Highlight intersecting edges.
    hi_edges = set()
    for it in intersections:
        hi_edges.add(tuple(sorted(it.edge0)))
        hi_edges.add(tuple(sorted(it.edge1)))
    hi_segs = []
    for e in edges:
        if tuple(sorted((e.i, e.j))) not in hi_edges:
            continue
        p0, p1 = edge_unwrapped_uv(uv, e.i, e.j)
        base = np.stack([p0, p1])
        for sx in (-1, 0, 1):
            for sy in (-1, 0, 1):
                pp = base + np.array([sx, sy])
                if bbox_overlaps_unit(pp):
                    hi_segs.append([pp[0], pp[1]])
    if hi_segs:
        lc2 = LineCollection(hi_segs, linewidths=1.1, colors=[(1, 0, 0, 0.9)])
        ax.add_collection(lc2)

    if intersections:
        pts = np.array([it.point for it in intersections], dtype=np.float64)
        ax.scatter(pts[:, 0], pts[:, 1], s=18, marker="x", color="red", label="edge intersections")
    if inside:
        pts = np.array([vi.point for vi in inside], dtype=np.float64)
        ax.scatter(pts[:, 0], pts[:, 1], s=16, marker="o", facecolors="none", edgecolors="magenta", label="vertex inside")
    if intersections or inside:
        ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def draw_zoom(out_path: Path, uv: np.ndarray, faces: np.ndarray, edges: List[Edge], center: Tuple[float, float], radius: float,
              title: str, involved_edges: Optional[Set[Tuple[int, int]]] = None, involved_faces: Optional[Set[int]] = None):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PolyCollection

    cx, cy = center
    fig, ax = plt.subplots(figsize=(6, 6), dpi=220)
    ax.set_aspect("equal")
    ax.set_xlim(cx - radius, cx + radius)
    ax.set_ylim(cy - radius, cy + radius)
    ax.set_title(title, fontsize=9)

    segs = []
    hi = []
    for e in edges:
        p0, p1 = edge_unwrapped_uv(uv, e.i, e.j)
        base = np.stack([p0, p1])
        for sx in (-1, 0, 1):
            for sy in (-1, 0, 1):
                pp = base + np.array([sx, sy])
                mn = np.min(pp, axis=0)
                mx = np.max(pp, axis=0)
                if mx[0] >= cx - radius and mn[0] <= cx + radius and mx[1] >= cy - radius and mn[1] <= cy + radius:
                    key = tuple(sorted((e.i, e.j)))
                    if involved_edges and key in involved_edges:
                        hi.append([pp[0], pp[1]])
                    else:
                        segs.append([pp[0], pp[1]])
    if segs:
        ax.add_collection(LineCollection(segs, linewidths=0.35, colors=[(0, 0, 0, 0.25)]))
    if involved_faces:
        polys = []
        for fid in involved_faces:
            pts = face_unwrapped_uv(uv, faces[fid])
            for sx in (-1, 0, 1):
                for sy in (-1, 0, 1):
                    pp = pts + np.array([sx, sy])
                    mn = np.min(pp, axis=0)
                    mx = np.max(pp, axis=0)
                    if mx[0] >= cx - radius and mn[0] <= cx + radius and mx[1] >= cy - radius and mn[1] <= cy + radius:
                        polys.append(pp)
        if polys:
            ax.add_collection(PolyCollection(polys, facecolors=[(1.0, 0.65, 0.0, 0.25)], edgecolors=[(1.0, 0.35, 0.0, 0.9)], linewidths=0.8))
    if hi:
        ax.add_collection(LineCollection(hi, linewidths=1.5, colors=[(1, 0, 0, 0.95)]))
    ax.scatter([cx], [cy], s=25, marker="x", color="blue")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--meshbin", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--grid", type=int, default=256)
    ap.add_argument("--area-eps", type=float, default=1e-18)
    ap.add_argument("--max-samples", type=int, default=200, help="max samples for edge intersections / vertex-inside; 0 = unlimited")
    ap.add_argument("--zoom-radius", type=float, default=0.035)
    ap.add_argument("--max-zooms", type=int, default=32)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    width, height, verts, faces = read_meshbin(args.meshbin)
    uv = verts[:, :2].astype(np.float64)
    edges = build_edges(faces)
    edge_copies = make_edge_copies(uv, edges)
    face_copies = make_face_copies(uv, faces)

    bad_rows = analyze_bad_faces(uv, faces, args.area_eps)
    intersections = find_edge_intersections(edge_copies, args.grid, args.max_samples)
    inside = find_vertex_inside(uv, faces, face_copies, args.grid, args.max_samples)

    bad_csv_rows = []
    for r in bad_rows:
        if r["negative"] or r["tiny"]:
            rr = dict(r)
            rr["negative"] = int(rr["negative"])
            rr["tiny"] = int(rr["tiny"])
            bad_csv_rows.append(rr)
    write_csv(args.out_dir / "bad_faces.csv", bad_csv_rows,
              ["face", "i0", "i1", "i2", "area2", "abs_area2", "negative", "tiny", "max_edge", "cx", "cy"])

    write_csv(args.out_dir / "edge_intersections.csv",
              [{"eid0": x.eid0, "edge0_i": x.edge0[0], "edge0_j": x.edge0[1],
                "eid1": x.eid1, "edge1_i": x.edge1[0], "edge1_j": x.edge1[1],
                "u": x.point[0], "v": x.point[1], "raw_u": x.raw_point[0], "raw_v": x.raw_point[1]}
               for x in intersections],
              ["eid0", "edge0_i", "edge0_j", "eid1", "edge1_i", "edge1_j", "u", "v", "raw_u", "raw_v"])

    write_csv(args.out_dir / "vertex_inside.csv",
              [{"vertex": x.vertex, "face": x.face, "u": x.point[0], "v": x.point[1],
                "b0": x.bary[0], "b1": x.bary[1], "b2": x.bary[2]}
               for x in inside],
              ["vertex", "face", "u", "v", "b0", "b1", "b2"])

    draw_overview(args.out_dir / "embedding_failures_overview.png", uv, faces, edges, intersections, inside, bad_rows)

    z = 0
    for it in intersections[:args.max_zooms]:
        involved = {tuple(sorted(it.edge0)), tuple(sorted(it.edge1))}
        draw_zoom(args.out_dir / f"zoom_{z:03d}_edge_intersection.png", uv, faces, edges, it.point,
                  args.zoom_radius, f"edge intersection: {it.edge0} x {it.edge1}", involved_edges=involved)
        z += 1
    for vi in inside[:max(0, args.max_zooms - z)]:
        draw_zoom(args.out_dir / f"zoom_{z:03d}_vertex_inside.png", uv, faces, edges, vi.point,
                  args.zoom_radius, f"vertex {vi.vertex} inside face {vi.face}", involved_faces={vi.face})
        z += 1
    for r in bad_csv_rows[:max(0, args.max_zooms - z)]:
        center = (float(r["cx"]), float(r["cy"]))
        draw_zoom(args.out_dir / f"zoom_{z:03d}_bad_face.png", uv, faces, edges, center,
                  args.zoom_radius, f"bad face {r['face']} area2={float(r['area2']):.3e}", involved_faces={int(r["face"])})
        z += 1

    summary = {
        "meshbin": str(args.meshbin),
        "width": width,
        "height": height,
        "vertices": int(verts.shape[0]),
        "faces": int(faces.shape[0]),
        "unique_edges": len(edges),
        "negative_faces": sum(1 for r in bad_rows if r["negative"]),
        "tiny_faces": sum(1 for r in bad_rows if r["tiny"]),
        "edge_intersections_sampled_or_total": len(intersections),
        "vertex_inside_sampled_or_total": len(inside),
        "grid": args.grid,
        "max_samples": args.max_samples,
    }
    write_csv(args.out_dir / "summary.csv", [summary], list(summary.keys()))

    md = args.out_dir / "report.md"
    with open(md, "w", encoding="utf-8") as f:
        f.write("# Embedding failure localization report\n\n")
        for k, v in summary.items():
            f.write(f"- {k}: {v}\n")
        f.write("\n## Files\n\n")
        f.write("- `embedding_failures_overview.png`\n")
        f.write("- `bad_faces.csv`\n")
        f.write("- `edge_intersections.csv`\n")
        f.write("- `vertex_inside.csv`\n")
        f.write("- `zoom_*.png`\n")
        f.write("\n## Interpretation note\n\n")
        f.write("This script reports localized periodic-UV defects. It is diagnostic, not a replacement for the pipeline audit.\n")
    print(f"Wrote {args.out_dir}")
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
