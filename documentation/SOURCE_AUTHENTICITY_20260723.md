# 真实数据来源核验

核验日期：2026-07-23

## 发布者与服务

- 美国 Data.gov/NOAA 目录将该系列标记为 MARACOOS 发布、Maryland Department of Natural Resources 连续监测计划数据，访问级别为 public，许可为 CC0-1.0：
  https://catalog.data.gov/dataset/mddnr-station-otter-point-creek
- 目录说明 YSI 6600 V2/EXO2 数据记录器通常以 15 分钟间隔采样，参数包括溶解氧浓度、氧饱和度、水温等。
- MARACOOS ERDDAP 当日仍可公开访问。例如：
  - https://erddap.maracoos.org/erddap/tabledap/mddnr_Bishopville_Prong.html
  - https://erddap.maracoos.org/erddap/tabledap/mddnr_Greys_Creek.html
  - https://erddap.maracoos.org/erddap/tabledap/mddnr_Harris_Creek_Downstream.html

ERDDAP 页面明确给出机构为 MDDNR，溶解氧变量为 `mass_concentration_of_oxygen_in_sea_water`，单位 mg/L，并提供 CSV/JSON/NetCDF 子集访问。

## 本地冻结

- 选择规则在下载前冻结，不使用 DO 值选站：按 dataset ID 排序，排除明示 Bottom/Surface 成对标题，从 2025 向 2017 逆序查找最近且至少 5,000 行的暖季，取前 8 站。
- 冻结目录：`04_真实数据/MARACOOS_MDDNR_discovery_v3`。
- `download_manifest.json` SHA-256：`c06cb797bed4541135b0726ccb07ea04df1eb51929f60645cfce60d530733d2d`。
- ERDDAP 搜索响应和 8 站的 CSV/info 共 17/17 个被 manifest 列出的哈希已通过。

## 能够与不能声称的事

- 可以声称这些是可公开访问、有发布者元数据和本地字节哈希的真实监测数据。
- 不能把 15 分钟折线当作未观测连续 DO 真值。
- 不能把校准站的经验斜率分位数或最大值当作物理 Lipschitz 上界。
- 被 DP 删除的记录只是一个数学 replacement 见证，不证明该记录是传感器错误。
