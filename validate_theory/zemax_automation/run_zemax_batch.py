"""
run_zemax_batch.py  (v3: 9 alpha, full K, Z_span=+-5mm)
=========================================================
Python ZOS-API 批量PSF离焦扫描导出 — 系统A (D50.4f200+ML1)

功能:
  连接Zemax OpticStudio, 按9组alpha方案自动批量导出惠更斯PSF数据.
  总计590个PSF文件, 预计耗时约30分钟.

前置条件:
  1. Zemax OpticStudio已安装 (标准路径或自定义路径)
  2. Preferences -> ZOS-API -> Enable 已勾选
  3. pip install pythonnet (在Zemax所在机器上)
  4. 编辑下方 ZOS_FILE_PATH 和 BASE_OUTPUT_DIR

使用方法:
  python run_zemax_batch.py

输出:
  zemax_export/
    alpha_0.10_dz_0.036mm/Zemax_txt/PSF_*.txt  (277 files)
    alpha_0.20_dz_0.073mm/Zemax_txt/PSF_*.txt  (138 files)
    ...
    alpha_3.00_dz_1.088mm/Zemax_txt/PSF_*.txt  (11 files)
"""

import os
import sys
import time

# ==================== 用户配置 ====================

ZOS_FILE_PATH = (
    r"D:\LiLinhan_2026\CC-多距离强度相位复原"
    r"\HIK相位测量\D50.4f200+ML1_test_v1\D50.4f200+ML1.zos"
)

BASE_OUTPUT_DIR = (
    r"D:\LiLinhan_2026\CC-多距离强度相位复原"
    r"\HIK相位测量\validate_theory\zemax_export"
)

# ==================== alpha配置 (由 generate_alpha_configs.m 生成) ====================
# Z_span = 10mm fixed, symmetric around best focus [-5.0, +5.0] mm
# Total: 590 PSF files

alpha_configs = [
    {
        "alpha": 0.100, "dz_mm": 0.036, "output_dir": "alpha_0.10_dz_0.036mm",
        "dz_positions_mm": [
            -5.000,-4.964,-4.927,-4.891,-4.855,-4.819,-4.782,-4.746,-4.710,
            -4.674,-4.637,-4.601,-4.565,-4.529,-4.492,-4.456,-4.420,-4.384,
            -4.347,-4.311,-4.275,-4.239,-4.202,-4.166,-4.130,-4.094,-4.057,
            -4.021,-3.985,-3.948,-3.912,-3.876,-3.840,-3.803,-3.767,-3.731,
            -3.695,-3.658,-3.622,-3.586,-3.550,-3.513,-3.477,-3.441,-3.405,
            -3.368,-3.332,-3.296,-3.260,-3.223,-3.187,-3.151,-3.115,-3.078,
            -3.042,-3.006,-2.969,-2.933,-2.897,-2.861,-2.824,-2.788,-2.752,
            -2.716,-2.679,-2.643,-2.607,-2.571,-2.534,-2.498,-2.462,-2.426,
            -2.389,-2.353,-2.317,-2.281,-2.244,-2.208,-2.172,-2.136,-2.099,
            -2.063,-2.027,-1.990,-1.954,-1.918,-1.882,-1.845,-1.809,-1.773,
            -1.737,-1.700,-1.664,-1.628,-1.592,-1.555,-1.519,-1.483,-1.447,
            -1.410,-1.374,-1.338,-1.302,-1.265,-1.229,-1.193,-1.157,-1.120,
            -1.084,-1.048,-1.011,-0.975,-0.939,-0.903,-0.866,-0.830,-0.794,
            -0.758,-0.721,-0.685,-0.649,-0.613,-0.576,-0.540,-0.504,-0.468,
            -0.431,-0.395,-0.359,-0.323,-0.286,-0.250,-0.214,-0.178,-0.141,
            -0.105,-0.069,-0.032,0.000,0.004,0.040,0.076,0.113,0.149,0.185,
            0.221,0.258,0.294,0.330,0.366,0.403,0.439,0.475,0.511,0.548,
            0.584,0.620,0.656,0.693,0.729,0.765,0.801,0.838,0.874,0.910,
            0.947,0.983,1.019,1.055,1.092,1.128,1.164,1.200,1.237,1.273,
            1.309,1.345,1.382,1.418,1.454,1.490,1.527,1.563,1.599,1.635,
            1.672,1.708,1.744,1.780,1.817,1.853,1.889,1.926,1.962,1.998,
            2.034,2.071,2.107,2.143,2.179,2.216,2.252,2.288,2.324,2.361,
            2.397,2.433,2.469,2.506,2.542,2.578,2.614,2.651,2.687,2.723,
            2.759,2.796,2.832,2.868,2.905,2.941,2.977,3.013,3.050,3.086,
            3.122,3.158,3.195,3.231,3.267,3.303,3.340,3.376,3.412,3.448,
            3.485,3.521,3.557,3.593,3.630,3.666,3.702,3.738,3.775,3.811,
            3.847,3.884,3.920,3.956,3.992,4.029,4.065,4.101,4.137,4.174,
            4.210,4.246,4.282,4.319,4.355,4.391,4.427,4.464,4.500,4.536,
            4.572,4.609,4.645,4.681,4.717,4.754,4.790,4.826,4.863,4.899,
            4.935,4.971
        ],
    },
    {
        "alpha": 0.200, "dz_mm": 0.073, "output_dir": "alpha_0.20_dz_0.073mm",
        "dz_positions_mm": [
            -5.000,-4.927,-4.855,-4.782,-4.710,-4.637,-4.565,-4.492,-4.420,
            -4.347,-4.275,-4.202,-4.130,-4.057,-3.985,-3.912,-3.840,-3.767,
            -3.695,-3.622,-3.550,-3.477,-3.405,-3.332,-3.260,-3.187,-3.115,
            -3.042,-2.969,-2.897,-2.824,-2.752,-2.679,-2.607,-2.534,-2.462,
            -2.389,-2.317,-2.244,-2.172,-2.099,-2.027,-1.954,-1.882,-1.809,
            -1.737,-1.664,-1.592,-1.519,-1.447,-1.374,-1.302,-1.229,-1.157,
            -1.084,-1.011,-0.939,-0.866,-0.794,-0.721,-0.649,-0.576,-0.504,
            -0.431,-0.359,-0.286,-0.214,-0.141,-0.069,0.004,0.076,0.149,
            0.221,0.294,0.366,0.439,0.511,0.584,0.656,0.729,0.801,0.874,
            0.947,1.019,1.092,1.164,1.237,1.309,1.382,1.454,1.527,1.599,
            1.672,1.744,1.817,1.889,1.962,2.034,2.107,2.179,2.252,2.324,
            2.397,2.469,2.542,2.614,2.687,2.759,2.832,2.905,2.977,3.050,
            3.122,3.195,3.267,3.340,3.412,3.485,3.557,3.630,3.702,3.775,
            3.847,3.920,3.992,4.065,4.137,4.210,4.282,4.355,4.427,4.500,
            4.572,4.645,4.717,4.790,4.863,4.935,
        ],
    },
    {
        "alpha": 0.500, "dz_mm": 0.181, "output_dir": "alpha_0.50_dz_0.181mm",
        "dz_positions_mm": [
            -5.000,-4.819,-4.637,-4.456,-4.275,-4.094,-3.912,-3.731,-3.550,
            -3.368,-3.187,-3.006,-2.824,-2.643,-2.462,-2.281,-2.099,-1.918,
            -1.737,-1.555,-1.374,-1.193,-1.011,-0.830,-0.649,-0.468,-0.286,
            -0.105,0.000,0.076,0.258,0.439,0.620,0.801,0.983,1.164,1.345,
            1.527,1.708,1.889,2.071,2.252,2.433,2.614,2.796,2.977,3.158,
            3.340,3.521,3.702,3.884,4.065,4.246,4.427,4.609,4.790,4.971,
        ],
    },
    {
        "alpha": 0.800, "dz_mm": 0.290, "output_dir": "alpha_0.80_dz_0.290mm",
        "dz_positions_mm": [
            -5.000,-4.710,-4.420,-4.130,-3.840,-3.550,-3.260,-2.969,-2.679,
            -2.389,-2.099,-1.809,-1.519,-1.229,-0.939,-0.649,-0.359,-0.069,
            0.000,0.221,0.511,0.801,1.092,1.382,1.672,1.962,2.252,2.542,
            2.832,3.122,3.412,3.702,3.992,4.282,4.572,4.863,
        ],
    },
    {
        "alpha": 1.200, "dz_mm": 0.435, "output_dir": "alpha_1.20_dz_0.435mm",
        "dz_positions_mm": [
            -5.000,-4.565,-4.130,-3.695,-3.260,-2.824,-2.389,-1.954,-1.519,
            -1.084,-0.649,-0.214,0.000,0.221,0.656,1.092,1.527,1.962,2.397,
            2.832,3.267,3.702,4.137,4.572,
        ],
    },
    {
        "alpha": 1.571, "dz_mm": 0.570, "output_dir": "alpha_1.57_dz_0.570mm",
        "dz_positions_mm": [
            -5.000,-4.430,-3.861,-3.291,-2.722,-2.152,-1.583,-1.013,-0.444,
            0.000,0.126,0.696,1.265,1.835,2.404,2.974,3.543,4.113,4.682,
        ],
    },
    {
        "alpha": 2.000, "dz_mm": 0.725, "output_dir": "alpha_2.00_dz_0.725mm",
        "dz_positions_mm": [
            -5.000,-4.275,-3.550,-2.824,-2.099,-1.374,-0.649,0.000,0.076,
            0.801,1.527,2.252,2.977,3.702,4.427,
        ],
    },
    {
        "alpha": 2.500, "dz_mm": 0.906, "output_dir": "alpha_2.50_dz_0.906mm",
        "dz_positions_mm": [
            -5.000,-4.094,-3.187,-2.281,-1.374,-0.468,0.000,0.439,1.345,
            2.252,3.158,4.065,4.971,
        ],
    },
    {
        "alpha": 3.000, "dz_mm": 1.088, "output_dir": "alpha_3.00_dz_1.088mm",
        "dz_positions_mm": [
            -5.000,-3.912,-2.824,-1.737,-0.649,0.000,0.439,1.527,2.614,
            3.702,4.790,
        ],
    },
]

# ==================== ZOS-API 连接 (ConnectAsExtension) ====================

def connect_zemax():
    """通过ConnectAsExtension连接到已打开的Zemax OpticStudio实例"""
    import clr
    import winreg

    # 从注册表找到Zemax安装目录
    aKey = winreg.OpenKey(winreg.ConnectRegistry(None, winreg.HKEY_CURRENT_USER),
                          r"Software\Zemax", 0, winreg.KEY_READ)
    zemaxData = winreg.QueryValueEx(aKey, 'ZemaxRoot')
    zemax_root = zemaxData[0]
    winreg.CloseKey(aKey)
    print(f"  Zemax root: {zemax_root}")

    # 加载 NetHelper
    net_helper = os.path.join(zemax_root, r'ZOS-API\Libraries\ZOSAPI_NetHelper.dll')
    clr.AddReference(net_helper)
    import ZOSAPI_NetHelper

    success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize()
    if not success:
        raise RuntimeError("Cannot find OpticStudio")

    zemax_dir = ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory()
    print(f"  Zemax dir: {zemax_dir}")

    # 加载 ZOSAPI 程序集
    clr.AddReference(os.path.join(zemax_dir, r'ZOSAPI.dll'))
    clr.AddReference(os.path.join(zemax_dir, r'ZOSAPI_Interfaces.dll'))
    import ZOSAPI

    connection = ZOSAPI.ZOSAPI_Connection()
    if connection is None:
        raise RuntimeError("Unable to initialize NET connection to ZOSAPI")

    app = connection.ConnectAsExtension(0)
    if app is None:
        raise RuntimeError(
            "Unable to connect to Zemax.\n"
            "请确认: 1) Zemax已启动  2) Programming -> Interactive Extension 已启用"
        )

    if not app.IsValidLicenseForAPI:
        raise RuntimeError("License is not valid for ZOSAPI use")

    print("  Connected to OpticStudio!")
    print(f"  Serial: {app.SerialCode}")

    return app


def format_filename(dz_mm):
    """格式化离焦值为文件名"""
    return f"PSF_{dz_mm:+.3f}mm.txt"


def run_psf_export(app, config, base_output_dir):
    """对单个alpha方案导出所有离焦面的PSF"""
    the_system = app.PrimarySystem

    alpha = config["alpha"]
    dz_positions = config["dz_positions_mm"]
    output_subdir = config["output_dir"]

    output_path = os.path.join(base_output_dir, output_subdir, "Zemax_txt")
    os.makedirs(output_path, exist_ok=True)

    n = len(dz_positions)
    print(f"\n{'='*60}")
    print(f"  alpha = {alpha:.2f}, {n} defocus planes")
    print(f"  Output: {output_path}")
    print(f"{'='*60}")

    # 获取图像面信息
    img_surf = the_system.LDE.NumberOfSurfaces
    last_surf = img_surf - 1
    orig_thickness = the_system.LDE.GetSurfaceAt(last_surf).Thickness
    print(f"  Image surface: {img_surf}, last surf thickness: {orig_thickness:.4f} mm")

    # 设置PSF采样参数 (匹配现有代码: 512x512, 0.425um)
    psf_check = the_system.Tools.OpenHuygensPsf()
    psf_settings = psf_check.GetSettings()
    print(f"  PSF ImageDelta: {psf_settings.ImageDelta} um")
    print(f"  PSF ImageSampling: {psf_settings.ImageSampling}")
    psf_check.Close()

    errors = 0
    for idx, dz in enumerate(dz_positions):
        filename = format_filename(dz)
        filepath = os.path.join(output_path, filename)

        if os.path.exists(filepath):
            print(f"  [{idx+1}/{n}] {filename} (skip)")
            continue

        print(f"  [{idx+1}/{n}] dz={dz:+.3f}mm ...", end=" ", flush=True)

        try:
            # 方法1: 通过PSF设置直接设离焦 (推荐, 不修改系统)
            psf_tool = the_system.Tools.OpenHuygensPsf()
            psf_cfg = psf_tool.GetSettings()
            # Defocus单位是透镜单位 (mm)
            psf_cfg.Defocus = dz
            psf_tool.ApplyAndWaitForCompletion()

            results = psf_tool.GetResults()
            results.GetTextFile(filepath)
            psf_tool.Close()
            print("OK")

        except Exception as e1:
            # 方法2 (回退): 修改像面厚度
            try:
                if abs(orig_thickness) < 1e-10 or orig_thickness > 1e10:
                    # 厚度不合理, 先获取近轴像距
                    the_system.LDE.GetSurfaceAt(last_surf).Thickness = 231.1605 + dz
                else:
                    the_system.LDE.GetSurfaceAt(last_surf).Thickness = orig_thickness + dz

                psf_tool = the_system.Tools.OpenHuygensPsf()
                psf_tool.ApplyAndWaitForCompletion()
                results = psf_tool.GetResults()
                results.GetTextFile(filepath)
                psf_tool.Close()
                print("OK(method2)")
            except Exception as e2:
                print(f"FAILED: {e2}")
                errors += 1

    if errors > 0:
        print(f"  WARNING: {errors}/{n} files failed")
    print(f"  Done: {n - errors} files exported to {output_path}")


def main():
    print("=" * 60)
    print("  Zemax PSF Defocus Scan — Batch Export")
    print("  System: D50.4f200+ML1, Z_span=+-5mm")
    print("=" * 60)

    total_files = sum(len(c["dz_positions_mm"]) for c in alpha_configs)
    print(f"  Configurations: {len(alpha_configs)}")
    print(f"  Total PSF files: {total_files}")
    print(f"  Estimated time: ~{total_files * 3 / 60:.0f} minutes")
    print()

    # 连接Zemax (ConnectAsExtension - 用户已打开Zemax并加载系统)
    app = connect_zemax()
    the_system = app.PrimarySystem

    if the_system is None:
        raise RuntimeError("Unable to acquire Primary system")

    print(f"\n  Currently loaded: {the_system.SystemFile}")
    print()

    try:
        os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

        # 逐个方案导出
        total_errors = 0
        total_start = time.time()

        for i, config in enumerate(alpha_configs):
            print(f"\n>>> Configuration [{i+1}/{len(alpha_configs)}] <<<")
            t0 = time.time()

            try:
                run_psf_export(app, config, BASE_OUTPUT_DIR)
            except Exception as e:
                print(f"  *** Configuration failed: {e}")
                total_errors += 1

            elapsed = time.time() - t0
            print(f"  Time: {elapsed:.1f}s")

        total_elapsed = time.time() - total_start
        print(f"\n{'='*60}")
        print(f"  All done! Total time: {total_elapsed/60:.1f} minutes")
        if total_errors > 0:
            print(f"  WARNING: {total_errors} configurations had errors")
        print(f"{'='*60}")

    finally:
        print("  Disconnecting from Zemax...")
        # 不关闭Zemax, 让用户检查结果
        # connection.TheApplication.CloseApplication()


if __name__ == "__main__":
    main()
