# AirMSPI 几何参数对齐调参指南（VZA / VAA / RAA / Scattering Angle）

本项目当前几何链路里，和你问题最相关的可调项是：

- `scene.lookat_center_km`：控制“相机看向哪里”，会直接移动 `VZA` 最小值（近似主点）在图上的位置。
- `sensor.views.zenith_deg`：控制视线离天顶角，决定 VZA 场整体幅度和梯度。
- `sensor.views.azimuth_deg`：控制观测方位主方向，影响 VAA / RAA 方向分布。
- `sensor.trajectory.manual_flight_azimuth_deg`（当 `mode=manual_azimuth`）：控制图像坐标系绕视线方向的“旋转感”，通常是让 VAA 等值线与目标区域边界平行的首要参数。
- `solar.saa_deg` / `solar.sza_deg`：主要影响 RAA 与散射角，不直接决定 VZA 主点位置。

---

## 推荐调参顺序（先几何、后太阳）

1. **先固定太阳参数**（`sza/saa` 不动），只调相机几何。
2. **用 `lookat_center_km` 把 VZA 低值中心移到与 AirMSPI 一致的位置**。  
   - 如果模拟的 VZA 中心偏右，就减小 `lookat_center_km[0]`（x）；  
   - 偏上就减小 `lookat_center_km[1]`（y）。
3. **用 `manual_flight_azimuth_deg` 对齐 VAA 主方向**，让 VAA 梯度方向或等值线走向与目标区域边界一致。
4. **微调 `views.azimuth_deg`**，把 VAA 的绝对角度偏差（整体加减）对齐。
5. **微调 `views.zenith_deg`**，让 VZA 数值范围（例如 min/max）接近 AirMSPI。
6. **最后调 `solar.saa_deg/sza_deg`**，使 RAA 和 scattering angle 的空间梯度与数值范围一致。

---

## 经验规则（针对你描述的现象）

- “**VAA 要和区域某个边界平行**”：优先调 `manual_flight_azimuth_deg`。  
- “**VZA 中心不在区域中心**”：优先调 `lookat_center_km`，这在斜视几何下是正常现象，不要求一定在几何中心。  
- “**散射角梯度方向不对**”：优先检查 `solar.saa_deg` 与 `views.azimuth_deg` 是否同一坐标约定（都按顺时针、从 +x 方向起算）。

---

## 建议的自动化搜索（小范围网格）

可在当前初值附近做三阶段网格搜索（每次只保留前若干最优组合）：

1. 粗搜：`manual_flight_azimuth_deg ± 20°`、`lookat_center_km` 每轴 ±1 km  
2. 中搜：在粗搜最优点附近，角度步长 2°，平移步长 0.2 km  
3. 细搜：角度步长 0.5°，平移步长 0.05 km

目标函数可定义为：

`J = w1*RMSE(VZA) + w2*RMSE(VAA_wrapped) + w3*RMSE(RAA_wrapped) + w4*RMSE(SCA)`

其中 VAA/RAA 误差应使用环向角度差（-180°~180° 包裹差值），避免 0/360 跳变。
