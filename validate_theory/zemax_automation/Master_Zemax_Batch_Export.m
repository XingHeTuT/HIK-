%% Master_Zemax_Batch_Export.m
% 基于用户已有 ZOS-API 代码, 批量导出系统A的离焦PSF数据
% 参考: D:\LiLinhan_2026\相位反演\离焦还原\Auto\Master_Auto_Workflow.m
%
% 前置条件:
%   1. 关闭所有 Zemax 窗口
%   2. 打开一个 Zemax, 加载 D50.4f200+ML1.zos
%   3. Programming -> Interactive Extension, 确认 Instance: 1
%
% 输出: 9组alpha方案, 共590个PSF txt文件

clc; clear classes;

%% ===== 1. 全局配置 =====
export_root = ['D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量\', ...
               'validate_theory\zemax_export'];

% Zemax 安装路径 (与旧代码一致)
installPath = 'D:\Program Files\Ansys Zemax OpticStudio 2023 R1.00';

% === 离焦面索引 ===
% 修改像面前最后一个光学面的厚度来实现离焦
% 在脚本运行时会自动从 NSurf-1 获取

% === 扫描方案: 9组alpha, Z_span=±5mm (由 generate_alpha_configs.m 生成) ===
load(fullfile(fileparts(mfilename('fullpath')), 'alpha_configs.mat'), ...
     'configs', 'alpha_list', 'z_c_mm', 'Z_span_mm');

fprintf('===== Zemax PSF 批量导出 =====\n');
fprintf('系统: D50.4f200+ML1\n');
fprintf('离焦范围: ±%.0f mm\n', Z_span_mm/2);
fprintf('Alpha方案数: %d\n', length(configs));
fprintf('总PSF文件数: %d\n', sum(arrayfun(@(c) c.K, [configs{:}])));
fprintf('预计耗时: ~%.0f 分钟\n\n', sum(arrayfun(@(c) c.K, [configs{:}])) * 3 / 60);

%% ===== 2. 初始化 ZOS-API =====
% 检查是否已有连接
if ~exist('TheSystem', 'var') || isempty(TheSystem)
    fprintf('正在初始化 ZOS-API 环境...\n');

    try
        domain = System.AppDomain.CurrentDomain;
        check_asm = domain.GetAssemblies;
        loaded = false;
        for i = 0:check_asm.Length-1
            if contains(char(check_asm.Get(i).FullName), 'ZOSAPI_Interfaces')
                loaded = true; break;
            end
        end

        if ~loaded
            NET.addAssembly(fullfile(installPath, 'ZOSAPI_NetHelper.dll'));
            NET.addAssembly(fullfile(installPath, 'ZOSAPI_Interfaces.dll'));
            NET.addAssembly(fullfile(installPath, 'ZOSAPI.dll'));
        end

        import ZOSAPI_NetHelper.*;
        ZOSAPI_Initializer.Initialize(installPath);
    catch ME
        error('环境加载失败。路径: %s\n错误: %s', installPath, ME.message);
    end

    import ZOSAPI.*;
    fprintf('正在连接 Zemax (Instance 1)...\n');

    try
        TheConnection = ZOSAPI_Connection();
        TheApplication = TheConnection.ConnectToApplication();
    catch ME
        error(['连接失败!\n%s\n' ...
               '请: 1) 关闭所有Zemax窗口\n' ...
               '   2) 重新打开Zemax -> Interactive Extension -> 确认Instance: 1'], ME.message);
    end

    if isempty(TheApplication)
        error('连接返回为空。请检查 Interactive Extension 状态。');
    end

    TheSystem = TheApplication.PrimarySystem;
    if isempty(TheSystem)
        error('已连接但无法获取 PrimarySystem。');
    end

    fprintf('连接成功! 当前文件: %s\n', char(TheSystem.SystemFile));
else
    try
        TheSystem.SystemID;
        fprintf('使用现有 Zemax 连接...\n');
    catch
        error('现有连接已失效。请执行 "clear classes" 后重试。');
    end
end

%% ===== 3. 获取表面信息 =====
TheLDE = TheSystem.LDE;
NSurf = TheLDE.NumberOfSurfaces;
fprintf('\n系统共 %d 个面 (像面 = %d)\n', NSurf, NSurf);

% 显示所有表面信息, 帮助确认离焦面
fprintf('表面列表:\n');
for s = 1:NSurf
    surfObj = TheLDE.GetSurfaceAt(s);
    fprintf('  面%d: 厚度=%.4f mm, 材料=%s, 类型=%s\n', ...
        s, surfObj.Thickness, char(surfObj.Material), char(surfObj.TypeName));
end

% 离焦面 = 像面前最后一个光学面
TargetSurface_Index = NSurf - 1;
Surf_Defocus = TheLDE.GetSurfaceAt(TargetSurface_Index);
Original_Thickness = Surf_Defocus.Thickness;
fprintf('\n离焦控制面: Surface %d (原始厚度=%.4f mm)\n', ...
    TargetSurface_Index, Original_Thickness);

if abs(Original_Thickness) < 1e-6 || Original_Thickness > 1e8
    warning('厚度值异常 (inf或0), 将尝试自动修正为BFL...');
    % 如果厚度异常, 设为标称BFL
    Original_Thickness = 231.1605;
    Surf_Defocus.Thickness = Original_Thickness;
    fprintf('  已修正为: %.4f mm\n', Original_Thickness);
end

%% ===== 4. 创建并配置 Huygens PSF 分析 =====
fprintf('\n创建 Huygens PSF 分析...\n');
HuygensPSF = TheSystem.Analyses.New_HuygensPsf();
if isempty(HuygensPSF)
    error('无法创建 Huygens PSF 分析窗口。');
end

H_Settings = HuygensPSF.GetSettings();
fprintf('配置 PSF 参数...\n');

try
    % 光瞳采样: 512x512
    H_Settings.PupilSampleSize = ZOSAPI.Analysis.SampleSizes.S_512x512;

    % 像面采样: 512x512
    H_Settings.ImageSampleSize = ZOSAPI.Analysis.SampleSizes.S_512x512;

    % 像面采样间距: 0.425 um (匹配现有代码 phi_re_Gemini_final_v5forHIK_v5.m)
    H_Settings.ImageDelta = 0.425;

    % 波长: 主波长
    H_Settings.Wavelength.SetWavelengthNumber(1);

    % 类型: Linear
    try
        H_Settings.Type = ZOSAPI.Analysis.Settings.Psf.HuygensPsfTypes.Linear;
    catch
        H_Settings.Type = ZOSAPI.Analysis.Settings.HuygensPsfTypes.Linear;
    end

    % 归一化: False
    H_Settings.Normalize = false;

    % 视场: On-Axis
    H_Settings.Field.SetFieldNumber(1);

    fprintf('  PSF配置: Pupil=512x512, Image=512x512, Delta=0.425um\n');

catch ME
    HuygensPSF.Close();
    error('PSF 参数设置失败: %s', ME.message);
end

%% ===== 5. 循环导出 (9组alpha方案) =====
total_exported = 0;
total_errors = 0;
t_total = tic;

for cfg_idx = 1:length(configs)
    cfg = configs{cfg_idx};
    dz_positions = cfg.dz_positions_mm;
    n_pos = length(dz_positions);

    % 输出目录
    out_dir = fullfile(export_root, cfg.output_dir, 'Zemax_txt');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    fprintf('\n>>> [%d/%d] alpha=%.3f, %d面, 输出: %s\n', ...
        cfg_idx, length(configs), cfg.alpha, n_pos, out_dir);

    t_cfg = tic;
    n_ok = 0;

    for i = 1:n_pos
        dz = dz_positions(i);
        fileName = sprintf('PSF_%+.4fmm.txt', dz);
        filePath = fullfile(out_dir, fileName);

        % 断点续传: 跳过已存在文件
        if exist(filePath, 'file')
            n_ok = n_ok + 1;
            if mod(i, 50) == 0
                fprintf('  [%d/%d] (skip existing)\n', i, n_pos);
            end
            continue;
        end

        % 设置离焦
        Surf_Defocus.Thickness = Original_Thickness + dz;

        % 运行分析
        HuygensPSF.ApplyAndWaitForCompletion();

        % 获取结果并保存
        HuygensPSF_Results = HuygensPSF.GetResults();
        success = HuygensPSF_Results.GetTextFile(filePath);

        % 重试
        if ~success
            pause(0.3);
            HuygensPSF.ApplyAndWaitForCompletion();
            HuygensPSF_Results = HuygensPSF.GetResults();
            success = HuygensPSF_Results.GetTextFile(filePath);
        end

        if success
            n_ok = n_ok + 1;
        else
            fprintf('  [%d/%d] FAILED: %s\n', i, n_pos, fileName);
            total_errors = total_errors + 1;
        end

        % 进度报告
        if mod(i, 20) == 0
            elapsed = toc(t_cfg);
            fprintf('  [%d/%d] %.1fmin (latest: %s)\n', i, n_pos, elapsed/60, fileName);
        end
    end

    cfg_time = toc(t_cfg);
    total_exported = total_exported + n_ok;
    fprintf('  完成: %d/%d 文件, 耗时 %.1f 分钟\n', n_ok, n_pos, cfg_time/60);
end

% 恢复原始厚度
Surf_Defocus.Thickness = Original_Thickness;
HuygensPSF.Close();

total_time = toc(t_total);
fprintf('\n===== 全部完成! =====\n');
fprintf('总导出: %d 文件\n', total_exported);
fprintf('总错误: %d\n', total_errors);
fprintf('总耗时: %.1f 分钟\n', total_time/60);
