python python\pndf_eval_pndf.py `
  --nxy-bin out256_normal_rebased\brush_ds4_nxy.bin `
  --meshbin out256_normal_rebased\brush_normal_ds4_v16384.meshbin `
  --out-dir pndf_eval\brush_normal_v16384 `
  --normal-bins 64 `
  --sigma-pixels 8,16,32 `
  --random-queries 32 `
  --feature-queries 16 `
  --seam-queries 8

python python\pndf_eval_pndf.py `
  --nxy-bin out256_normal_rebased\isotropic_ds4_nxy.bin `
  --meshbin out256_normal_rebased\isotropic_normal_ds4_v16384.meshbin `
  --out-dir pndf_eval\isotropic_normal_v16384 `
  --normal-bins 64 `
  --sigma-pixels 8,16,32 `
  --random-queries 32 `
  --feature-queries 16 `
  --seam-queries 8

python python\pndf_eval_pndf.py `
  --nxy-bin out256_normal_rebased\scratched_ds4_nxy.bin `
  --meshbin out256_normal_rebased\scratched_normal_ds4_v16384.meshbin `
  --out-dir pndf_eval\scratched_normal_v16384 `
  --normal-bins 64 `
  --sigma-pixels 8,16,32 `
  --random-queries 32 `
  --feature-queries 16 `
  --seam-queries 8

python python\pndf_eval_pndf.py `
  --nxy-bin out256_4d_rebased\brush_ds4_nxy.bin `
  --meshbin out256_4d_rebased\brush_4d_ds4_v16384.meshbin `
  --out-dir pndf_eval\brush_4d_v16384 `
  --normal-bins 64 `
  --sigma-pixels 8,16,32 `
  --random-queries 32 `
  --feature-queries 16 `
  --seam-queries 8

python python\pndf_eval_pndf.py `
  --nxy-bin out256_4d_rebased\isotropic_ds4_nxy.bin `
  --meshbin out256_4d_rebased\isotropic_4d_ds4_v16384.meshbin `
  --out-dir pndf_eval\isotropic_4d_v16384 `
  --normal-bins 64 `
  --sigma-pixels 8,16,32 `
  --random-queries 32 `
  --feature-queries 16 `
  --seam-queries 8

python python\pndf_eval_pndf.py `
  --nxy-bin out256_4d_rebased\scratched_ds4_nxy.bin `
  --meshbin out256_4d_rebased\scratched_4d_ds4_v16384.meshbin `
  --out-dir pndf_eval\scratched_4d_v16384 `
  --normal-bins 64 `
  --sigma-pixels 8,16,32 `
  --random-queries 32 `
  --feature-queries 16 `
  --seam-queries 8
