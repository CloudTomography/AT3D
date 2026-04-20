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

然后在 Spyder 里直接 **Run File**。

脚本里也提供了可直接调用的函数：

- `run_retrieval_case(...)`

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
