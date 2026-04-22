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

## 两种输入模式的公共函数建议

建议把输入统一为两条路径，然后都产出同一个“扩展 CSV”规范：

- 路径 A（retrieval-1D netCDF）: 用 `build_from_retrieval_1d_netcdf(...)`
- 路径 B（LES/WRF arrays）: 用 `build_from_les_arrays(...)`

这两个函数都在 `scripts/grid_data_builder.py` 中，并且都返回：

- `df`（可直接写 CSV）
- `geometry`（`nx,ny,nz,dx,dy,z_levels`）
- `options`（mode/surface 等校验配置）

## 你的示例：从 retrieval_1d 生成 grid data（Spyder 方式）

如果你的文件在：

- `data/retrieval_1d/2019_0806_1839_N_Pxl25_3_3.nc`

请打开 `scripts/grid_data_builder.py`，编辑 `if __name__ == "__main__":` 下方变量：

- `input_nc`
- `output_csv`
- `dx_km / dy_km`
- `z_levels_km`
- `mode_count`
- `wavelength_index`

然后在 Spyder 里直接 **Run File**。

脚本里也提供了可直接调用的函数：

- `run_retrieval_case(...)`

- `wavelength_index`（重要）：当折射率等变量是 3D（例如 `(band, y, x)`）时，选择使用哪个波段切片。

> `dx/dy` 需要你按反演网格实际大小填写（文档中的 0.16 km 只是例子）。

## airmspi_image_simulation.py 需要做的调整

本仓库已做兼容更新：

1. 不再强制 `density='cv'`，改为自动识别 `cv/lwc/density`。
2. 若输入没有 `reff/veff`（LES 常见），自动补默认值 `reff=12`, `veff=0.10`。

因此在配置层通常只需要改：

- `scene.input_path` 指向你新生成的 CSV。

其余流程可保持不变，先跑通数据链路，再逐步接入 per-grid 折射率到光学性质模块。


## 为什么很多项目仍保留 CLI（即使你现在不用）

- 批处理方便（一次跑很多 case）。
- 可重复性强（命令即参数快照，便于论文/复现实验）。
- 更容易接入调度系统（crontab/slurm/CI）。
- 不依赖 IDE，服务器无图形界面时更稳。

你当前以 Spyder 为主是完全合理的：开发调试阶段效率更高。


## reff/veff 限定范围在哪里改

当前 `airmspi_image_simulation.py` 会在构建光学插值网格前对 `reff/veff` 做 sanitize + clip。
这些阈值不再硬编码，改为从配置读取：

- `aerosol.reff_default`
- `aerosol.veff_default`
- `aerosol.reff_clip_min`
- `aerosol.reff_clip_max`
- `aerosol.veff_clip_min`
- `aerosol.veff_clip_max`

建议直接在 `scripts/config_v5b.yaml` 的 `aerosol:` 段修改。


## 关于 reff/veff 与 mode1_reff/mode1_veff 重复

当前 builder 默认采用扩展列（`mode1_reff`, `mode1_veff`）作为主定义，不再强制写重复的标量列。
在仿真读取阶段，如果没有标量 `reff/veff`，会自动根据 `mode*_fraction` 对 `mode*_reff/veff` 做加权合成，生成 AT3D 所需的标量 `reff/veff`。

## resample 后异常值处理

`build_scene_and_sensors_single_band` 在 `resample_onto_grid` 之后会：

- 对 `reff/veff` 的 NaN/Inf/<=0 做替换；
- 对清空气柱（`density<=0`）回填默认 `reff/veff`，避免插值伪负值影响；
- clip 到配置区间后再用于光学插值。

## 多 mode 下的 reff/veff 范围

`reff_min/max` 与 `veff_min/max` 现在会同时考虑：

- 标量 `reff/veff`；
- 所有存在且 `mode*_fraction>0` 的 `mode*_reff/veff`。


## lat/lon 写入保证

`build_from_retrieval_1d_netcdf(...)` 现在会优先读取文件中的 `lat/lon`，若缺失或全零，自动用
`fallback_lat0/fallback_lon0 + dx/dy` 生成规则网格经纬度；默认 fallback 为 `(35.0, -112.0)`。

## MAXNBC / 内存压力调参位置

若出现 `BOUNDARY_PNTS: MAXNBC exceeded` 或内存相关 warning，可在 `config_v5b.yaml` 中调：

- `solver.split_accuracy`（增大可减少自适应细分）
- `solver.adapt_grid_factor`
- `solver.cell_to_point_ratio`
- `solver.max_total_mb`
- `aerosol.density_floor`（把极小密度直接置零，减少无意义细分）


## 关于 z 方向插值策略

为贴合你给的 MATLAB 版本语义，`build_scene_and_sensors_single_band` 现在只对 `density`
做 `resample_onto_grid`，其余字段（如 `lat/lon`, `mode*_fraction`, `mode*_reff/veff`, `mode*_mr/mi`）
按“z 不变”处理：取首层并沿 z 广播，避免在 z 方向被插值改写。


## Cv_total 到竖直质量浓度的写入方式

`build_from_retrieval_1d_netcdf(...)` 不再把 `Cv_total` 直接复制到每个 z 层。

现在默认按 MATLAB 逻辑把列积分量分配到竖直：

- 默认 `vertical_distribution="from_nc"`：使用 `Hmean_aerosol` 与 `Saerosol` 作为高斯中心/宽度；
- 可选 `vertical_distribution="fixed"`：使用 `fixed_h_km` / `fixed_sigma_km`；
- 可选 `vertical_distribution="uniform"`：仅用于回退。

质量浓度写入遵循：

- `mass_density_z = rho_p[g/cm^3] * Cv_total[um] * w(z) * 1e-3`  (g/m^3)

并保证 `sum(w * dz_km)=1`。


## ENU 到 NEU 的坐标映射

NC 网格通常可理解为 ENU（行/列近似 North/East）。AT3D 使用 NEU，且 `x=North`, `y=East`。
因此 builder 在写 CSV 时会把：

- `x` 轴映射到原网格的 North 方向（原二维数组行方向）；
- `y` 轴映射到原网格的 East 方向（原二维数组列方向）。

## dx/dy 自动估算（来自 lat/lon）

当 `dx_km` 或 `dy_km` 传入为 `None` 时，builder 会根据 nc 中经纬度估算（支持变量名 `latitude/longitude` 或 `lat/lon`，并自动忽略 0 与 NaN）：

- `dx_km`：相邻“行”点的地理距离中位数（North 方向）；
- `dy_km`：相邻“列”点的地理距离中位数（East 方向）。

这样 `dx` 与 `dy` 可不相等，并且通常接近你说的 ~0.25 km 量级（取决于具体场景分辨率）。

若 nc 缺少 lat/lon 且你又未显式给 `dx_km/dy_km`，builder 会先用 0.25 km 作为 seed 生成规则经纬网，再据此更新估算值。
