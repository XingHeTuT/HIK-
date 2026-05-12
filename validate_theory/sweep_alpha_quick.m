% sweep_alpha_quick.m
% 快速验证: 仅测试4个关键alpha值, 确认理论与实验趋势一致
% 系统A: D50.4f200+ML1, z_c=0.363mm
%
% 测试的alpha值:
%   alpha=0.14  (当前步长0.05mm, 信息塌缩区)
%   alpha=0.55  (过渡区, dz=0.2mm)
%   alpha=1.57  (理论最优, dz=0.57mm, pi/2)
%   alpha=2.76  (过最优区, dz=1.0mm)

clc; clear all; %#ok<CLALL>

%% ===== 系统A参数 =====
sys = struct();
sys.lambda_um       = 10.6;
sys.n_lens          = 2.40266;
sys.n_substrate     = 3.422;
sys.R1_mm           = 2.805400000000000E+002;
sys.Tc_lens_mm      = 3.500000000000000E+000;
sys.d_air_gap_mm    = 2.305480752288749E+001;
sys.d_substrate_mm  = 7.250000000000000E-001;
sys.BFL_mm          = 2.311604993761137E+002;
sys.D_stop1_mm      = 2.540000000000000E+001 * 2;
sys.D_metalens_mm   = 2.229883854002862E+001 * 2;
sys.Grid_Window_mm  = 46;
sys.N_sim           = 2*2048;
sys.N_det           = 512;
sys.det_pitch_um    = 0.425;
sys.A_coeffs_ideal  = [5.5968453070e-01, -1.5760996841e-03, 4.5486150095e-06, ...
                       -3.4017277167e-09, 0, 0, 0, 0, 0];
sys.radius_norm_mm  = 1.0;
sys.dataFolder      = ['D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量\', ...
                       'D50.4f200+ML1_test_v1\Zemax_xlsx'];

%% ===== 理论参数 =====
lambda  = sys.lambda_um * 1e-6;
D_ml    = sys.D_metalens_mm * 1e-3;
BFL     = sys.BFL_mm * 1e-3;
NA_eff  = (D_ml/2) / BFL;
z_c     = lambda / (pi * NA_eff^2);
a_opt   = pi/2;

fprintf('===== 系统A理论参数 =====\n');
fprintf('NA=%.4f, z_c=%.3fmm, alpha_opt=%.3f, Dz_opt=%.3fmm\n', ...
    NA_eff, z_c*1e3, a_opt, a_opt*z_c*1e3);

%% ===== 4个关键alpha测试点 =====
alpha_list = [0.138, 0.55, 1.57, 2.76];
dz_list_mm = alpha_list * z_c * 1e3;
labels = {'塌缩区(当前)', '过渡区', '理论最优(pi/2)', '过最优区'};

%% ===== 算法参数 (快速版: 15次迭代) =====
alg = struct();
alg.max_iter    = 15;   % 快速验证用15次迭代
alg.alpha_start = 1.0;
alg.alpha_end   = 0.8;
alg.beta_amp    = 0.8;

% 固定扫描范围: 覆盖约1.5mm范围, 确保每个方案有足够平面
z_range_mm = [-1.0, 0.5];

fprintf('\n===== 快速验证 (%d个alpha值, %d次迭代) =====\n', ...
    length(alpha_list), alg.max_iter);

results = cell(length(alpha_list), 1);
t_start = tic;

for i = 1:length(alpha_list)
    fprintf('\n>>> [%d/%d] alpha=%.3f (dz=%.3fmm) -- %s\n', ...
        i, length(alpha_list), alpha_list(i), dz_list_mm(i), labels{i});

    samp = struct();
    samp.delta_z_mm  = dz_list_mm(i);
    samp.z_range_mm  = z_range_mm;
    samp.max_iter    = alg.max_iter;
    samp.alpha_start = alg.alpha_start;
    samp.alpha_end   = alg.alpha_end;
    samp.beta_amp    = alg.beta_amp;

    try
        r = run_mdpr_validation(sys, samp);
        results{i} = r;
        fprintf('  -> RMSE_2D=%.4f, Err_high=%.4f, 面数=%d\n', ...
            r.RMSE_2D, r.coeff_err_mean_high, r.nPlanes);
    catch ME
        fprintf('  *** 出错: %s\n', ME.message);
        results{i} = struct('RMSE_2D', NaN, 'coeff_err_mean_high', NaN, 'nPlanes', 0);
    end
end

elapsed = toc(t_start);
fprintf('\n总耗时: %.1f 分钟\n', elapsed/60);

%% ===== 结果汇总 =====
fprintf('\n===== 验证结果 =====\n');
fprintf('%-12s | %-8s | %-12s | %-12s\n', 'alpha', 'dz[mm]', 'RMSE_2D', 'Err_high');
fprintf('%s\n', repmat('-', 1, 48));

rmse_vals = zeros(length(alpha_list), 1);
err_high  = zeros(length(alpha_list), 1);
n_planes  = zeros(length(alpha_list), 1);

for i = 1:length(alpha_list)
    rmse_vals(i) = results{i}.RMSE_2D;
    err_high(i)  = results{i}.coeff_err_mean_high;
    n_planes(i)  = results{i}.nPlanes;
    fprintf('%-12.3f | %-8.3f | %-12.4f | %-12.4f\n', ...
        alpha_list(i), dz_list_mm(i), rmse_vals(i), err_high(i));
end

% 找最优
[~, idx_best] = min(rmse_vals);
fprintf('\n最优alpha: %.3f (dz=%.3fmm), 理论: %.3f (dz=%.3fmm)\n', ...
    alpha_list(idx_best), dz_list_mm(idx_best), a_opt, a_opt*z_c*1e3);
fprintf('偏差: %.1f%%\n', abs(alpha_list(idx_best)-a_opt)/a_opt*100);

%% ===== 快速绘图 =====
figure('Color','w','Position',[100,100,900,400]);

yyaxis left;
plot(alpha_list, rmse_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 10);
ylabel('RMSE_{2D} (rad)');
hold on;

yyaxis right;
plot(alpha_list, err_high, 'rs-', 'LineWidth', 2, 'MarkerSize', 10);
ylabel('高阶系数相对误差');

xline(a_opt, 'r--', 'LineWidth', 2);
xlabel('\alpha = \Deltaz / z_c');
grid on;
title(sprintf('System A: Phase Retrieval Error vs \\alpha (z_c=%.2fmm)', z_c*1e3));
legend('RMSE_{2D}', 'Err_{high}', '\alpha_{opt}=\pi/2', 'Location', 'best');

% 保存
save_dir = fileparts(mfilename('fullpath'));
save(fullfile(save_dir, 'sweep_alpha_quick_results.mat'), ...
    'alpha_list', 'dz_list_mm', 'rmse_vals', 'err_high', 'n_planes', ...
    'z_c', 'a_opt', 'results');
fprintf('\n结果已保存\n');
