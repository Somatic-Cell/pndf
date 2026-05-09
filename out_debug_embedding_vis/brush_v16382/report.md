# Embedding failure localization report

- meshbin: out_quality_phase1_guard_v4_quick\meshes\brush_pnm_he_normal_ds4_v16384.meshbin
- width: 256
- height: 256
- vertices: 16384
- faces: 32768
- unique_edges: 49152
- negative_faces: 0
- tiny_faces: 1
- edge_intersections_sampled_or_total: 0
- vertex_inside_sampled_or_total: 0
- grid: 256
- max_samples: 500

## Files

- `embedding_failures_overview.png`
- `bad_faces.csv`
- `edge_intersections.csv`
- `vertex_inside.csv`
- `zoom_*.png`

## Interpretation note

This script reports localized periodic-UV defects. It is diagnostic, not a replacement for the pipeline audit.
