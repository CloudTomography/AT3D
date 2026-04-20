# Cross-track 模式下 `v_out_map` 的来源与 VAA 计算说明

> 结论先行：
>
> - 普通模式走 `perspective_projection`；cross-track 现改为**直接构造扫描射线**（不再用 projective pinhole 逐 view 建相机）。
> - 两者最终都通过 `cam_mu/cam_phi` 还原 `v_out_map`，所以 `v_out_map` 差异来自几何输入（扫描轨迹与扫描角）而不是角度计算函数本身。

---

## 1. 总体链路

1. 先根据轨迹模式计算每个 view 的外参：
   - 相机位置 `position_vector`
   - 注视点 `lookat_vector`
   - 相机上方向 `up_vector`
2. 普通模式：调用 `at3d.sensor.perspective_projection(...)`。  
   cross-track 模式：调用 `cross_track_scan_projection(...)`，直接生成每个像元的 `(cam_x, cam_y, cam_z, cam_mu, cam_phi)`。
3. `compute_vout_map_from_sensor(sensor_ds)` 统一读取 `cam_mu/cam_phi`，重建每个像元在世界坐标中的 `v_out_map`。

也就是说：`compute_vout_map_from_sensor` 是“解码器”，模式差异主要由前面的“编码器”（轨迹 + 投影参数）决定。

---

## 2. cross-track 与普通模式到底共用了什么？

## 共用部分

- 都由 `compute_vout_map_from_sensor(...)` 从 `cam_mu/cam_phi` 统一反推出 `v_out_map`。  
- 都使用 `_compute_angle_maps_from_sensor(...)` 统一算 `VZA/VAA/RAA/Scattering Angle`。

## 不同部分（决定 `v_out_map` 差异的核心）

- 普通模式：`lookat/up` + `fov/resolution` 决定像元视线。  
- `cross_track` 模式：
  - 用 `cross_track_x1..z2` + `spacing` 生成 along-track 扫描行；
  - 用 `scan1/scan2/delscan` 生成跨轨扫描列；
  - 可选为每条 along-track 扫描线指定独立 pitch（线性插值或手动列表）；
  - 逐像元直接写入 `cam_mu/cam_phi`（一行一个扫描位置，一列一个扫描角）。

因此，即使两种模式最终都走同一角度反演流程，几何输入不同也会让 `v_out_map` 不同。

---

## 3. `v_out_map` 的数学与代码对应

`compute_vout_map_from_sensor` 的关键步骤：

1. 读取像元方向变量 `cam_mu/cam_phi` 与图像尺寸 `x_resolution/y_resolution`。
2. 根据球坐标关系恢复方向向量：
   - `sin(theta)=sqrt(1-mu^2)`
   - `vx=sin(theta)*cos(phi)`
   - `vy=sin(theta)*sin(phi)`
   - `vz=mu`
3. 组合得到 `v_out_map(y, x, 3)`；其定义是“**场景 -> 传感器**”传播方向。

---

## 4. VAA 是如何从 `v_out_map` 得到的

在 `_compute_angle_maps_from_sensor(...)` 中：

1. `v_out_map = compute_vout_map_from_sensor(sensor_ds)`
2. `VZA = arccos(v_out[...,2])`
3. `VAA = atan2(v_out[...,1], v_out[...,0])` 映射到 `[0, 360)`
4. 然后额外做 `VAA = (VAA + 180) % 360`（图像方向约定）
5. `RAA = (VAA - SAA) % 360`

所以你看到的 VAA 图，最终由三类因素共同决定：

- 像元方向（`cam_mu/cam_phi`）
- 像元坐标对应的相机位置（`cam_x/cam_y/cam_z`，影响地面投影）
- 角度约定后处理（尤其 `+180` 这一项）

---

## 5. 关键变量清单：哪些会让 cross-track 和普通模式产生不同 `v_out_map`

### 直接影响（强相关）

- `cam_mu`
- `cam_phi`
- `cam_x`, `cam_y`, `cam_z`
- `x_resolution`, `y_resolution`

### cross-track 特有输入（通过上面三个外参间接影响）

- `cross_track_x1/y1/z1`
- `cross_track_x2/y2/z2`
- `cross_track_spacing`
- `cross_track_scan1_deg`
- `cross_track_scan2_deg`
- `cross_track_delscan_deg`
- `cross_track_pitch_start_deg` / `cross_track_pitch_end_deg`（按扫描线线性插值）
- `cross_track_pitch_list_deg`（手动逐扫描线 pitch 列表）

### 当前实现里仅做元数据记录（不直接进 `v_out_map` 计算）

- `cross_track_nbytes`
- `cross_track_scale`

---

## 6. 一个易混点：为什么“共用 projection”，结果还是不一样？

因为最终角度产品只依赖每个像元的方向向量。

- 普通模式通过 pinhole 投影生成 `cam_mu/cam_phi`；
- cross-track 模式直接生成 `cam_mu/cam_phi`；
- 两者进入同一 VZA/VAA/RAA 计算函数，自然可比较且结果一致。

---

## 7. 现阶段实现边界（为避免误解）

当前 `cross_track` 采用“单个 sensor 承载整幅扫描图”的方式：

- `y` 方向表示 along-track 的扫描行；
- `x` 方向表示 cross-track 的扫描列；
- 因而不会像“每个 sample 一个 view”那样造成大量视角对象，速度/内存更稳定。
- 传感器数据集通过 `shdom_cross_track_sensor_wrapper(...)` 直接调用
  `at3d.sensor.make_sensor_dataset(...)` 生成（AT3D 对 SHDOM 的封装路径）。

如果后续要与 SHDOM 的 `V` 模式逐像元严格对齐，通常还需要进一步把扫描几何映射到
每个像元列（或时间序列）层面，而不只是 view 级采样。
