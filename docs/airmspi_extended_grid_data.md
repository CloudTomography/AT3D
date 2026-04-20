# AirMSPI/AT3D 扩展 grid data 方案（多 mode + 可选 BRDF/BPDF + 每格折射率）

本文给出当前仓库可执行的落地方案。

## 可行性结论

1. **一个 grid data 文件包含多个颗粒 mode：可行。**  
   `at3d.util.load_from_csv` 会读取除 `x,y,z` 外的任意列，因此可以直接加入 `mode1_*`, `mode2_*` ... 列，再由上层脚本决定如何使用。

2. **在 grid data 中加入 surface BRDF/BPDF 参数：可行（建议“可选列”设计）。**  
   若列不存在，可退化为 `lambertian`。若存在 `surface_model` 与配套参数，可在场景构建阶段映射到 `at3d.surface.*`。

3. **每个 grid 不同折射率：数据层面可行，正演层面需要额外策略。**  
   数据文件可以直接存 `mode*_mr/mode*_mi`。但 AT3D 光学性质通常由 Mie 表驱动，严格逐格折射率需要“分箱/聚类+多表混合”或“逐格查表”，计算量显著增加。建议先实现数据写出与读取，再在求解脚本里分阶段接入。

## 推荐扩展列

- 基础列（保持兼容）
  - `x,y,z,cv,reff,veff,lat,lon`
- 多 mode 列（每个 mode 一组）
  - `mode{n}_fraction,mode{n}_reff,mode{n}_veff,mode{n}_mr,mode{n}_mi`
- 可选 surface 列
  - `surface_model`
  - `surface_albedo`
  - `surface_wind_speed`
  - `surface_pigmentation`
  - `surface_refractive_index_real`
  - `surface_refractive_index_imag`
  - `surface_parameter_1/2/3`

## surface_model 到 at3d 的建议映射

- `lambertian` -> `at3d.surface.lambertian`
- `wave_fresnel` -> `at3d.surface.wave_fresnel`
- `diner` -> `at3d.surface.diner`
- `ocean_unpolarized` -> `at3d.surface.ocean_unpolarized`
- `rpv_unpolarized` -> `at3d.surface.rpv_unpolarized`

## 新增工具函数

本仓库新增 `scripts/grid_data_builder.py`，用于：

- 统一构建扩展 DataFrame；
- 校验 mode 分数和 surface model；
- 按 AT3D 兼容头格式写出 CSV。

该函数与 MATLAB 反演结果输出解耦，适合作为“新的专用函数”规范。
