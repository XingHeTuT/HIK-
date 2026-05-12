% sweep_alpha_validation.m
% 高频最优步长验证: 扫描 alpha 值, 测试不同采样步长下的相位恢复质量
%
% 理论预测 (论文 Eq 17-18):
%   - 高频 Fisher 信息增量: Delta_F_high ∝ sin²(alpha)
%   - 最优无量纲步长: alpha_opt = pi/2 ≈ 1.57
%   - 小步长塌缩: alpha << 1 时信息按 alpha² 衰减
%
% 验证策略:
%   从系统A的密集数据(0.05mm步长)子采样构造不同的alpha方案,
%   对每种方案运行 MDPR, 比较高阶系数的恢复精度

clc; clear all; %#ok<CLALL>

%% ===== 系统配置 (系统A: D50.4f200+ML1, 密集采样) =====
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
sys.A_coeffs_ideal  = [5.5968453070e-01, -1.5760996841e-03, 4.5486150095e-06, -3.4017277167e-09, ...
                       0, 0, 0, 0, 0];
sys.radius_norm_mm  = 1.0;
sys.dataFolder      = ['D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量\', ...
                       'D50.4f200+ML1_test_v1\Zemax_xlsx'];

%% ===== 计算理论特征参数 =====
lambda  = sys.lambda_um * 1e-6;
D_ml    = sys.D_metalens_mm * 1e-3;
BFL     = sys.BFL_mm * 1e-3;
NA_eff  = (D_ml/2) / BFL;
z_c     = lambda / (pi * NA_eff^2);  % 特征离焦尺度 [m]
alpha_opt_theory = pi/2;              % 理论最优
fprintf('===== 系统特征参数 =====\n');
fprintf('lambda = %.1f um, D_ml = %.2f mm, BFL = %.1f mm\n', sys.lambda_um, sys.D_metalens_mm, sys.BFL_mm);
fprintf('NA_eff = %.4f\n', NA_eff);
fprintf('z_c = %.4f mm (特征离焦尺度)\n', z_c * 1e3);
fprintf('alpha_opt (理论) = pi/2 = %.4f\n', alpha_opt_theory);
fprintf('Delta_z_opt (理论) = %.3f mm\n', alpha_opt_theory * z_c * 1e3);

%% ===== alpha 扫描范围设计 =====
% 基于系统A的原始数据 (0.05mm步长, 范围-1.70到+0.90mm)
% alpha = Delta_z / z_c

% 定义测试的 alpha 值
% 覆盖: 塌缩区 (alpha << 1) -> 过渡区 -> 最优区 -> 过采样区
alpha_list = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0];
% 对应的物理步长
dz_list_mm = alpha_list * z_c * 1e3;

fprintf('\n===== 扫描计划 =====\n');
fprintf('测试 %d 个 alpha 值\n', length(alpha_list));
fprintf('alpha 范围: [%.2f, %.2f]\n', alpha_list(1), alpha_list(end));
fprintf('物理步长范围: [%.3f, %.3f] mm\n', dz_list_mm(1), dz_list_mm(end));
fprintf('理论最优 alpha = %.2f (步长 %.3f mm)\n', alpha_opt_theory, alpha_opt_theory * z_c * 1e3);

%% ===== 算法参数 (固定) =====
alg = struct();
alg.max_iter    = 30;
alg.alpha_start = 1.0;
alg.alpha_end   = 0.8;
alg.beta_amp    = 0.8;

% 固定扫描范围: -1.0 ~ +0.5 mm (确保所有alpha方案覆盖相同物理范围)
z_range_mm = [-1.0, 0.5];

%% ===== 执行扫描 =====
fprintf('\n===== 开始 alpha 扫描 =====\n');

n_alpha = length(alpha_list);
results = cell(n_alpha, 1);

for i = 1:n_alpha
    fprintf('\n--- [%d/%d] alpha = %.3f, Delta_z = %.3f mm ---\n', ...
        i, n_alpha, alpha_list(i), dz_list_mm(i));

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
        fprintf('  -> RMSE_2D = %.4f, 高阶系数误差 = %.4f, 面数 = %d\n', ...
            r.RMSE_2D, r.coeff_err_mean_high, r.nPlanes);
    catch ME
        fprintf('  *** 出错: %s\n', ME.message);
        results{i} = struct('RMSE_2D', NaN, 'coeff_err_mean_high', NaN, 'nPlanes', 0);
    end
end

%% ===== 结果汇总与绘图 =====
fprintf('\n===== 汇总结果 =====\n');

% 提取指标
rmse_vals   = cellfun(@(r) r.RMSE_2D, results);
err_high    = cellfun(@(r) r.coeff_err_mean_high, results);
n_planes    = cellfun(@(r) r.nPlanes, results);

valid = ~isnan(rmse_vals);

fprintf('%-8s | %-10s | %-12s | %-12s | %-6s\n', ...
    'alpha', 'dz[mm]', 'RMSE_2D', 'Err_high', 'K');
fprintf('%s\n', repmat('-', 1, 56));
for i = 1:n_alpha
    fprintf('%-8.3f | %-10.3f | %-12.4f | %-12.4f | %-6d\n', ...
        alpha_list(i), dz_list_mm(i), rmse_vals(i), err_high(i), n_planes(i));
end

% 找到最优 alpha
if any(valid)
    [~, idx_best_rmse] = min(rmse_vals(valid));
    [~, idx_best_high] = min(err_high(valid));
    valid_idx = find(valid);
    fprintf('\n最优 alpha (按 RMSE_2D): %.3f (步长 %.3f mm)\n', ...
        alpha_list(valid_idx(idx_best_rmse)), dz_list_mm(valid_idx(idx_best_rmse)));
    fprintf('最优 alpha (按高阶误差): %.3f (步长 %.3f mm)\n', ...
        alpha_list(valid_idx(idx_best_high)), dz_list_mm(valid_idx(idx_best_high)));
    fprintf('理论最优: alpha = %.3f (步长 %.3f mm)\n', ...
        alpha_opt_theory, alpha_opt_theory * z_c * 1e3);
end

% 绘图
figure('Color', 'w', 'Position', [100, 100, 1400, 500]);

% 子图1: RMSE vs alpha
subplot(1, 2, 1);
yyaxis left;
plot(alpha_list, rmse_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'RMSE_{2D} (整体)');
xlabel('无量纲步长 \alpha');
ylabel('RMSE_{2D} (rad)');
hold on;

% 标注理论最优
xline(alpha_opt_theory, 'r--', 'LineWidth', 2, 'DisplayName', '\alpha_{opt}=\pi/2');
grid on; legend('Location', 'best');
title('恢复误差 vs 无量纲步长');

% 子图2: 高频系数误差 vs alpha (半对数)
subplot(1, 2, 2);
semilogy(alpha_list, err_high, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', '高阶系数相对误差');
hold on;
xline(alpha_opt_theory, 'r--', 'LineWidth', 2, 'DisplayName', '\alpha_{opt}=\pi/2');

% 叠画理论曲线: sin²(alpha) 的倒数比例
alpha_dense = linspace(0.01, pi, 200);
theory_curve = 1 ./ (sin(alpha_dense).^2 + 0.01);  % 避免除零
theory_curve = theory_curve / max(theory_curve) * max(err_high(valid));
plot(alpha_dense, theory_curve, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', '理论: 1/sin^2(\alpha)');
xlabel('无量纲步长 \alpha');
ylabel('高阶系数相对误差');
grid on; legend('Location', 'best');
title('高频恢复误差 vs 无量纲步长 (对数坐标)');
xlim([0, max(alpha_list)]);

sgtitle(['高频最优步长验证: \alpha_{opt}^{theory} = \pi/2 = ', ...
    num2str(alpha_opt_theory, '%.2f'), ...
    ', z_c = ', num2str(z_c*1e3, '%.2f'), ' mm'], ...
    'FontSize', 14, 'FontWeight', 'bold');

% 保存结果
save_dir = fileparts(mfilename('fullpath'));
save(fullfile(save_dir, 'sweep_alpha_results.mat'), ...
    'alpha_list', 'dz_list_mm', 'rmse_vals', 'err_high', 'n_planes', ...
    'z_c', 'alpha_opt_theory', 'sys', 'alg', 'results');
fprintf('\n结果已保存至: %s\n', fullfile(save_dir, 'sweep_alpha_results.mat'));
