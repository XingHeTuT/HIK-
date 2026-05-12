% compute_theoretical_params.m
% 基于论文 Section 2.2-2.3，计算各系统的理论特征参数
% 核心公式:
%   z_c = lambda / (pi * NA^2)     -- 特征离焦尺度 (Eq 7)
%   alpha = Delta_z / z_c           -- 无量纲步长 (Eq 9)
%   alpha_opt = pi/2                -- 高频最优步长 (Eq 18)
%   SBP ~ D * NA / lambda           -- 空间带宽积

clc; fprintf('===== 理论特征参数计算 =====\n\n');

%% 系统 A: D50.4f200+ML1 (主实验, 密集采样)
sysA.name    = '系统A: D50.4f200+ML1 (密集0.05mm)';
sysA.lambda  = 10.6e-6;
sysA.D_ml    = 44.5977e-3;    % 超透镜有效口径 [m]
sysA.BFL     = 231.1605e-3;   % 后截距 [m]
sysA.Delta_z = 0.05e-3;       % 当前物理步长 [m]
sysA.N_planes = 59;           % 离焦面数
sysA.z_range = [-1.70, 0.90] * 1e-3; % [z_min, z_max] in m

%% 系统 B: 更改焦平面步长
sysB.name    = '系统B: 更改焦平面步长 (1.0mm)';
sysB.lambda  = 10.6e-6;
sysB.D_ml    = 44.5832e-3;
sysB.BFL     = 231.1605e-3;
sysB.Delta_z = 1.0e-3;
sysB.N_planes = 38;
sysB.z_range = [-11, 25] * 1e-3;

%% 系统 C: 大口径/硫系玻璃
sysC.name    = '系统C: 大口径硫系玻璃 (0.2mm)';
sysC.lambda  = 10.6e-6;
sysC.D_ml    = 44.5832e-3;
sysC.BFL     = 231.1605e-3;
sysC.Delta_z = 0.2e-3;
sysC.N_planes = 17;
sysC.z_range = [-1.6, 1.6] * 1e-3;

%% 系统 D: 代码正确性测试 (10mm球超)
sysD.name    = '系统D: 10mm球超测试 (0.02mm, lambda=10um)';
sysD.lambda  = 10.0e-6;
sysD.D_ml    = 44.5832e-3;    % MATLAB代码中使用相同口径，但Zemax模型是10mm
sysD.BFL     = 256.7808e-3;
sysD.Delta_z = 0.02e-3;
sysD.N_planes = 25;
sysD.z_range = [-0.24, 0.24] * 1e-3;

systems = [sysA, sysB, sysC, sysD];

fprintf('%-50s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    '系统', 'NA', 'z_c[mm]', 'alpha', 'a_opt', 'dZopt[mm]', 'SBP');
fprintf('%s\n', repmat('-', 1, 110));

for s = systems
    NA      = (s.D_ml/2) / s.BFL;
    f_c     = NA / s.lambda;
    z_c     = s.lambda / (pi * NA^2);
    alpha   = s.Delta_z / z_c;
    a_opt   = pi/2;  % 理论最优无量纲步长 (Eq 18)
    dZ_opt  = a_opt * z_c;
    SBP     = s.D_ml * NA / s.lambda;
    z_span  = s.z_range(2) - s.z_range(1);
    beta    = z_span / z_c;

    fprintf('%-50s | %8.4f | %8.3f | %8.3f | %8.3f | %8.3f | %8.1f\n', ...
        s.name, NA, z_c*1e3, alpha, a_opt, dZ_opt*1e3, SBP);
    fprintf('  -> 当前扫描范围: %.1f mm (beta=%.1f), 面数: %d\n', ...
        z_span*1e3, beta, s.N_planes);

    if alpha < 0.1
        fprintf('  *** 警告: alpha << 1, 当前步长远小于特征离焦尺度!\n');
        fprintf('  *** 预测: 高频信息按 alpha^2=%.4f 塌缩, 相邻面近重复编码\n', alpha^2);
    end
end

fprintf('\n===== 关键理论预测 =====\n');
fprintf('1. 高频最优步长: alpha_opt = pi/2 = %.4f (式18)\n', pi/2);
a_current = systems(1).Delta_z / (systems(1).lambda / (pi * ((systems(1).D_ml/2)/systems(1).BFL)^2));
fprintf('2. 系统A当前 alpha=%.3f << pi/2, 信息塌缩因子 = %.4f\n', a_current, a_current^2);
fprintf('3. 系统A最优物理步长建议: %.2f mm (约 %.0f 倍当前步长)\n', ...
    systems(1).Delta_z * (pi/2) / a_current, (pi/2)/a_current);
fprintf('4. 低频扫描范围 beta_min 应正比于 SBP^2\n');
