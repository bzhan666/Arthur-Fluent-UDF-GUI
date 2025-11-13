**文件名:** `fluent_c3m_udf.inp`  
**位置:** 放在Fluent工作目录下

```ini
! ========================================
! Fluent C3M UDF Input File
! 文件名必须是: fluent_c3m_udf.inp
! 放置在Fluent工作目录下
! ========================================

! ========================================
! 1. 煤质分析 (收到基, %)
! ========================================
fc_ar 54.1
vm_ar 41.8
ash_ar 1.5
moist_ar 2.6
f_ep_a 0.5

! ========================================
! 2. 脱挥发分模型开关
! 注意: 只能开启一个脱挥发分模型
! ========================================
pccl_devol false
mgas_devol true
cpd_devol false
fgdvc_devol false
hptr_devol false

! ========================================
! 3. 水分释放模型
! ========================================
pccl_moisture false
mgas_moisture true

! ========================================
! 4. 二次热解和焦油裂解
! ========================================
pccl_2nd_pyro false
mgas_tarcracking true
pccl_tarcracking false

! ========================================
! 5. 气化模型
! ========================================
mgas_gasif true
pccl_gasif false
pccl_soot_gasif false

! ========================================
! 6. 燃烧模型
! ========================================
mgas_char_combustion true
pccl_char_combustion false
pccl_soot_oxidation false
tar_oxidation true
mgas_gas_phase_oxidation true

! ========================================
! 7. 水煤气变换反应
! ========================================
mgas_wgs true

! ========================================
! 8. 水分释放动力学参数
! ========================================
a_moisture_release 2.0e5
e_moisture_release 7.4e7
moisture_flux 0.001
wg3 0.014

! ========================================
! 9. 脱挥发分动力学参数 (双步反应)
! ========================================
a1_devolatilization 2.0e5
e1_devolatilization 7.4e7
a2_devolatilization 1.3e7
e2_devolatilization 1.4e8

! ========================================
! 10. 焦油裂解动力学参数
! Tar -> Light gases
! ========================================
a_tar_cracking 1.0e8
e_tar_cracking 1.0e8

! ========================================
! 11. 水蒸气气化 (C + H2O -> CO + H2)
! ========================================
a_steam_gasification 4.93e3
e_steam_gasification 3.18e4
k_steam_gasification 0.0
n_steam_gasification 1.0
annealing_steam_gasification 1.0

! ========================================
! 12. CO2气化 (C + CO2 -> 2CO)
! ========================================
a_co2_gasification 4.4e3
e_co2_gasification 3.63e4
k_co2_gasification 0.0
n_co2_gasification 1.0
annealing_co2_gasification 1.0

! ========================================
! 13. H2气化 (C + 2H2 -> CH4)
! ========================================
a_h2_gasification 1.0e3
e_h2_gasification 2.5e4
n_h2_gasification 1.0
annealing_h2_gasification 1.0

! ========================================
! 14. 焦炭燃烧 (C + O2 -> CO/CO2)
! ========================================
a_c_combustion 8710.0
e_c_combustion 27000.0
n_c_combustion 0.0
annealing_c_combustion 1.0

! ========================================
! 15. 碳烟-水蒸气气化
! ========================================
a_soot_steam_gasification 1.0e3
e_soot_steam_gasification 3.0e4
k_soot_steam_gasification 0.0
n_soot_steam_gasification 1.0
annealing_soot_steam_gasification 1.0

! ========================================
! 16. 碳烟-CO2气化
! ========================================
a_soot_co2_gasification 1.0e3
e_soot_co2_gasification 3.5e4
k_soot_co2_gasification 0.0
n_soot_co2_gasification 1.0
annealing_soot_co2_gasification 1.0

! ========================================
! 17. 碳烟-H2气化
! ========================================
a_soot_h2_gasification 5.0e2
e_soot_h2_gasification 2.0e4
n_soot_h2_gasification 1.0
annealing_soot_h2_gasification 1.0

! ========================================
! 18. 碳烟燃烧
! ========================================
a_soot_combustion 1.0e5
e_soot_combustion 2.0e4

! ========================================
! END OF INPUT FILE
! ========================================
```

## 使用说明

### 文件保存

1. 将上述内容复制到文本编辑器
2. 另存为 `fluent_c3m_udf.inp`
3. 放在Fluent的工作目录（与.cas/.dat文件同目录）

### 参数修改指南

#### 煤质分析参数

- **fc_ar**: 固定碳含量 (%)
- **vm_ar**: 挥发分含量 (%)
- **ash_ar**: 灰分含量 (%)
- **moist_ar**: 水分含量 (%)
- 四者之和应为100%

#### 动力学参数

- **a_xxx**: 指前因子 (单位各异)
- **e_xxx**: 活化能 (cal/mol)
- **k_xxx**: 抑制系数
- **n_xxx**: 反应级数
- **annealing_xxx**: 退火因子

#### 模型开关

- 使用 `true` 或 `false`
- 注意只能开启一个脱挥发分模型

### 验证配置

在Fluent中加载UDF后,查看控制台输出:

```
start of read_c3m_data
fopen OK
PCCL_Devol = 0
MGAS_Devol = 1
fc_ar = 54.100000
vm_ar = 41.800000
...
read_c3m_data completed successfully
```

如果看到 "ERROR: Could not open fluent_c3m_udf.inp file!",检查:

1. 文件名拼写
2. 文件位置
3. 文件权限