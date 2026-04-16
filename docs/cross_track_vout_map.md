# Cross-track 模式下 `v_out_map` 的来源与 VAA 计算说明

> 结论先行：
>
> - `cross_track` 和普通 `projection`（当前脚本里是 `perspective_projection`）**共用同一套成像与反演方向向量流程**。
> - 两者产生不同 `v_out_map` 的根因，不在 `compute_vout_map_from_sensor` 本身，而在**构建 sensor 时传入的几何参数不同**（主要是 `position_vector / lookat_vector / up_vector`）。

---

## 1. 总体链路

1. 先根据轨迹模式计算每个 view 的外参：
   - 相机位置 `position_vector`
   - 注视点 `lookat_vector`
   - 相机上方向 `up_vector`
2. 用这些外参调用 `at3d.sensor.perspective_projection(...)` 生成 sensor dataset。
3. 在该 dataset 的 `attrs` 中写入：
   - `rotation_matrix`
   - `sensor_to_camera_transform_matrix`（内参矩阵 `K`）
   - `x_resolution`, `y_resolution`, `fov` 等
4. `compute_vout_map_from_sensor(sensor_ds)` 读取这些 attrs，重建每个像元在世界坐标中的光线方向，得到 `v_out_map`。

也就是说：`compute_vout_map_from_sensor` 是“解码器”，模式差异主要由前面的“编码器”（轨迹 + 投影参数）决定。

---

## 2. cross-track 与普通模式到底共用了什么？

## 共用部分

- 都走 `at3d.sensor.perspective_projection(...)`。  
- 都通过同一套公式生成 `rotation_matrix` 与 `K`。  
- 都由 `compute_vout_map_from_sensor(...)` 统一反推出 `v_out_map`。  
- 都使用 `_compute_angle_maps_from_sensor(...)` 统一算 `VZA/VAA/RAA/Scattering Angle`。

## 不同部分（决定 `v_out_map` 差异的核心）

- 普通模式（当前脚本）：`lookat` 通常用固定中心点，或由 aircraft 姿态间接确定。  
- `cross_track` 模式：
  - 用 `cross_track_x1..z2` + `spacing` 生成 along-track 相机位置序列；
  - 用 `scan1/scan2/delscan` 对基准下视方向做绕航迹轴旋转，生成每个 sample 的 `lookat`；
  - 再基于 `look_dir` 与航迹轴求 `up_vector`。

因此，虽然投影函数完全一样，但传入的外参不同，最终 `rotation_matrix` 不同，`v_out_map` 也就不同。

---

## 3. `v_out_map` 的数学与代码对应

`compute_vout_map_from_sensor` 的关键步骤：

1. 从 `sensor_ds.attrs` 读取：
   - `rotation_matrix`（相机坐标到世界坐标旋转）
   - `sensor_to_camera_transform_matrix`（内参矩阵 `K`）
   - `x_resolution`, `y_resolution`
2. 在归一化像平面构造像元中心网格 `(x_s, y_s, 1)`。
3. 通过 `inv(K)` 得到相机坐标系下射线 `rays_cam`。
4. 用 `rotation_matrix @ rays_cam` 转到世界坐标 `rays_world`。
5. 输出 `v_out_map = -rays_world`：定义为“**场景 -> 传感器**”传播方向。

注意最后的负号：
- `rays_world` 是“传感器看向场景”的方向；
- `v_out` 被定义为反方向（场景到传感器），用于后续几何角度约定。

---

## 4. VAA 是如何从 `v_out_map` 得到的

在 `_compute_angle_maps_from_sensor(...)` 中：

1. `v_out_map = compute_vout_map_from_sensor(sensor_ds)`
2. `VZA = arccos(v_out[...,2])`
3. `VAA = atan2(v_out[...,1], v_out[...,0])` 映射到 `[0, 360)`
4. 然后额外做 `VAA = (VAA + 180) % 360`（图像方向约定）
5. `RAA = (VAA - SAA) % 360`

所以你看到的 VAA 图，最终由三类因素共同决定：

- 外参（`position/lookat/up`）
- 内参（`fov`, `x_resolution`, `y_resolution`）
- 角度约定后处理（尤其 `+180` 这一项）

---

## 5. 关键变量清单：哪些会让 cross-track 和普通模式产生不同 `v_out_map`

### 直接影响（强相关）

- `position_vector`
- `lookat_vector`
- `up_vector`
- `fov`
- `x_resolution`, `y_resolution`

### cross-track 特有输入（通过上面三个外参间接影响）

- `cross_track_x1/y1/z1`
- `cross_track_x2/y2/z2`
- `cross_track_spacing`
- `cross_track_scan1_deg`
- `cross_track_scan2_deg`
- `cross_track_delscan_deg`

### 当前实现里仅做元数据记录（不直接进 `v_out_map` 计算）

- `cross_track_nbytes`
- `cross_track_scale`

---

## 6. 一个易混点：为什么“共用 projection”，结果还是不一样？

因为 `perspective_projection` 只是一个通用成像算子；
它像 `f(camera_extrinsic, camera_intrinsic)`。

- 算子本身不变；
- 输入（尤其外参）一变，`rotation_matrix` 就变；
- 同样代码也会生成不同 `v_out_map`。

这正是 cross-track 与普通模式“共用代码但几何结果不同”的本质原因。

---

## 7. 现阶段实现边界（为避免误解）

当前 `cross_track` 实现是“按 sample 生成离散 view 外参”并交给 pinhole 相机渲染，
并非严格复刻 SHDOM `V` 模式里“扫描线/像元级时序采样器”的完整测量过程。

补充（当前代码行为）：

- cross-track 轨迹函数会先生成**全部** `(scan_position, scan_angle)` 组合样本。  
- 若配置里 `sensor.views.names` 只有 1 个基名，会自动扩展为 `base_ct_0000`, `base_ct_0001`, ...，保证每个样本都会建一个 view。  
- 若手工提供了多个 `sensor.views.names`，其数量必须与生成样本数一致，否则会报错。  

如果后续要与 SHDOM 的 `V` 模式逐像元严格对齐，通常还需要进一步把扫描几何映射到
每个像元列（或时间序列）层面，而不只是 view 级采样。
