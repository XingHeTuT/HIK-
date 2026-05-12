! export_psf_defocus.ZPL
! ============================================================================
! 自动化惠更斯PSF离焦扫描导出
!
! 用途: 对D50.4f200+ML1系统, 在指定离焦位置批量导出PSF强度数据
!
! 使用方法:
!   1. 在Zemax OpticStudio中打开目标系统文件 (.zos)
!   2. 菜单: Programming -> Edit/Run ZPL -> 打开本文件 -> Run
!   3. 修改下方"用户配置区"的参数后运行
!
! 原理:
!   修改像面前最后一个面的厚度来添加离焦, 然后运行惠更斯PSF分析并导出
!
! 输出: Zemax_txt/ 目录下生成 PSF_{dz}mm.txt 格式的文件
! ============================================================================

! ==================== 用户配置区 (请修改) ====================

! 输出目录 (PSF txt文件保存位置)
output_dir$ = "D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量\validate_theory\zemax_export\alpha_0.50_dz_0.181mm\Zemax_txt\"

! 系统文件路径 (仅用于记录, 实际运行时请已在Zemax中打开)
system_file$ = "D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量\D50.4f50.4AS+ML1_test_v1\D50.4f200+ML1.zos"

! 离焦位置列表 (单位: mm, 以最佳聚焦面为原点)
! 可以在ZPL中直接编辑, 或从外部文件读取
!
! 格式: dz_values(N) = 离焦量 [mm]
! 注意: ZPL数组需要预先声明大小

N_PLANES = 9   ! 离焦面数量, 必须与下面的赋值一致

DECLARE dz_values, DOUBLE, 1, N_PLANES

! ---- 编辑此处的离焦值 (mm) ----
dz_values(1) = -0.725
dz_values(2) = -0.544
dz_values(3) = -0.363
dz_values(4) = -0.181
dz_values(5) = 0.000
dz_values(6) = 0.181
dz_values(7) = 0.363
dz_values(8) = 0.544
dz_values(9) = 0.725
! ---------------------------------

! ==================== 宏开始 ====================

! 创建输出目录 (如果不存在)
! 注意: Zemax ZPL没有直接创建目录的命令, 需手动确保目录存在
! 这里使用SYSTEM函数尝试创建
SYSTEM "mkdir " + output_dir$

PRINT "======================================"
PRINT "  离焦PSF自动导出"
PRINT "  输出目录: ", output_dir$
PRINT "  离焦面数: ", N_PLANES
PRINT "======================================"

! 获取系统表面信息
img_surf = NSUR()
last_surf = img_surf - 1

PRINT "图像面: ", img_surf
PRINT "最后光学面: ", last_surf

! 保存原始像面距离
orig_thickness = THIC(last_surf)
PRINT "原始后截距: ", orig_thickness, " mm"

! 报告原始斯特列尔比 (可选验证)
PRINT ""
PRINT "开始扫描..."

! 主循环
FOR i, 1, N_PLANES, 1

    dz = dz_values(i)

    PRINT "--- [", i, "/", N_PLANES, "] 离焦 = ", dz, " mm ---"

    ! 设置新的像面距离: BFL + dz
    new_thickness = orig_thickness + dz
    SETTHICKNESS last_surf, new_thickness

    ! 更新系统
    UPDATE ALL

    ! 打开惠更斯PSF分析窗口
    ! 注意: OpticStudio不同版本中, "HuygensPSF" 的分析编号可能不同
    ! 常见编号: 34 (OpticStudio 17+), 或使用名称
    ! 如果下面命令不工作, 请尝试替换为 ANALYSIS 34
    ANALYSIS "HuygensPSF"

    ! 等待计算完成 (0.5秒, 复杂系统可能需要更长)
    DELAY 0.5

    ! 构造输出文件名
    ! PSF_{dz}mm.txt - 注意处理负号和科学计数法
    IF (dz >= 0)
        dz_str$ = $STR(dz)
    ELSE
        dz_str$ = $STR(dz)
    ENDIF

    filename$ = output_dir$ + "PSF_" + dz_str$ + "mm.txt"

    ! 导出文本
    GETTEXTFILE filename$

    PRINT "  已保存: ", filename$

NEXT

! 恢复原始像面距离
SETTHICKNESS last_surf, orig_thickness
UPDATE ALL

PRINT ""
PRINT "======================================"
PRINT "  扫描完成! 共导出 ", N_PLANES, " 个文件"
PRINT "  已恢复原始系统配置"
PRINT "======================================"
