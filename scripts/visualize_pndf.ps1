# python python\pndf_eval_pndf.py `
#   --nxy-bin out_quality_local_guard_eps1e-9\nxy\brush_ds4_nxy.bin `
#   --meshbin out_quality_local_guard_eps1e-9\meshes\brush_pnm_he_4d_ds4_v8192.meshbin `
#   --out-dir out_quality_local_guard_eps1e-9\metrics_sampled_pndf\brush_pnm_he_4d_ds4_v8192 `
#   --normal-bins 64 `
#   --sigma-pixels 64 `
#   --random-queries 16 `
#   --feature-queries 16 `
#   --seam-queries 8

python python\pndf_eval_pndf.py `
  --nxy-bin out_quality_local_guard_eps1e-9\nxy\isotropic_ds4_nxy.bin `
  --meshbin out_quality_local_guard_eps1e-9\meshes\isotropic_pnm_he_4d_ds4_v8192.meshbin `
  --out-dir out_quality_local_guard_eps1e-9\metrics_sampled_pndf\isotropic_pnm_he_4d_ds4_v8192 `
  --normal-bins 64 `
  --sigma-pixels 64 `
  --random-queries 16 `
  --feature-queries 16 `
  --seam-queries 8

python python\pndf_eval_pndf.py `
  --nxy-bin out_quality_local_guard_eps1e-9\nxy\scratched_ds4_nxy.bin `
  --meshbin out_quality_local_guard_eps1e-9\meshes\scratched_pnm_he_4d_ds4_v8192.meshbin `
  --out-dir out_quality_local_guard_eps1e-9\metrics_sampled_pndf\scratched_pnm_he_4d_ds4_v8192 `
  --normal-bins 64 `
  --sigma-pixels 64 `
  --random-queries 16 `
  --feature-queries 16 `
  --seam-queries 8
