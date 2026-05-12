% sweep_alpha_v2.m
% 改进版高频最优步长验证
%
% 改进点:
%   1. 固定面数K=9 (而非固定扫描范围), 消除面数差异的混淆
%   2. 使用30次迭代, 保证充分收敛
%   3. 高阶系数误差作为主要指标 (对应论文高频Fisher信息分析)
%   4. 扫描范围随步长自动缩放: Z_span = K * Delta_z
%
% 理论预测:
%   对固定K, 总Fisher信息 F_total ∝ Σ_m sin²(m*alpha)
%   当K较大时, 最优alpha趋近于pi/2

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
fprintf('NA=%.4f, z_c=%.4fmm, alpha_opt=%.3f, Dz_opt=%.3fmm\n', ...
    NA_eff, z_c*1e3, a_opt, a_opt*z_c*1e3);

%% ===== 实验设计: 固定K=9 =====
K_fixed = 9;

% 测试6个alpha值: 覆盖塌缩区→过渡区→最优区
alpha_list = [0.1, 0.2, 0.5, 0.8, 1.2, 1.57];
dz_list_mm = alpha_list * z_c * 1e3;
z_span_list_mm = K_fixed * dz_list_mm;  % 扫描范围随步长缩放

fprintf('\n===== 实验设计 (K=%d固定) =====\n', K_fixed);
fprintf('%-10s | %-10s | %-12s | %-8s\n', 'alpha', 'dz[mm]', 'Z_span[mm]', 'planes');
fprintf('%s\n', repmat('-', 1, 48));
for i = 1:length(alpha_list)
    fprintf('%-10.3f | %-10.3f | %-12.2f | %-8d\n', ...
        alpha_list(i), dz_list_mm(i), z_span_list_mm(i), K_fixed);
end

% 理论Fisher信息预测 (对K个等步长面, 以最佳聚焦面为中心)
fprintf('\n理论预测 (K=%d, 总Fisher ∝ Σ sin²(m*alpha)):\n', K_fixed);
M = (K_fixed - 1) / 2;  % 半范围
F_theory = zeros(size(alpha_list));
for i = 1:length(alpha_list)
    a = alpha_list(i);
    F_theory(i) = 0;
    for m = -M:M
        F_theory(i) = F_theory(i) + sin(m * a)^2;
    end
    fprintf('  alpha=%.2f: F_total=%.3f\n', a, F_theory(i));
end
[~, idx_opt_theory] = max(F_theory);
fprintf('  理论最优: alpha=%.2f\n', alpha_list(idx_opt_theory));

%% ===== 算法参数 =====
alg = struct();
alg.max_iter    = 30;
alg.alpha_start = 1.0;
alg.alpha_end   = 0.8;
alg.beta_amp    = 0.8;

%% ===== 执行扫描 =====
fprintf('\n===== 开始K固定扫描 =====\n');

results = cell(length(alpha_list), 1);
t_start = tic;

for i = 1:length(alpha_list)
    fprintf('\n>>> [%d/%d] alpha=%.3f (dz=%.3fmm, span=%.1fmm, K=%d)\n', ...
        i, length(alpha_list), alpha_list(i), dz_list_mm(i), z_span_list_mm(i), K_fixed);

    % 扫描范围: 以最佳聚焦面为中心
    half_span = z_span_list_mm(i) / 2;
    z_range_mm = [-half_span, half_span];

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
        fprintf('  -> RMSE_2D=%.4f, Err_high=%.4f, err_low=%.4f, 面数=%d\n', ...
            r.RMSE_2D, r.coeff_err_mean_high, ...
            mean(r.coeff_err_low(~isnan(r.coeff_err_low))), r.nPlanes);
    catch ME
        fprintf('  *** 出错: %s\n', ME.message);
        results{i} = struct('RMSE_2D', NaN, 'coeff_err_mean_high', NaN, ...
            'coeff_err_low', [NaN, NaN], 'nPlanes', 0);
    end
end

elapsed = toc(t_start);
fprintf('\n总耗时: %.1f 分钟\n', elapsed/60);

%% ===== 结果汇总 =====
fprintf('\n===== 验证结果 =====\n');
fprintf('%-10s | %-10s | %-12s | %-12s | %-12s | %-6s\n', ...
    'alpha', 'dz[mm]', 'RMSE_2D', 'Err_high', 'Err_low', 'K');
fprintf('%s\n', repmat('-', 1, 70));

rmse_vals = zeros(length(alpha_list), 1);
err_high  = zeros(length(alpha_list), 1);
err_low   = zeros(length(alpha_list), 1);
n_planes  = zeros(length(alpha_list), 1);

for i = 1:length(alpha_list)
    rmse_vals(i) = results{i}.RMSE_2D;
    err_high(i)  = results{i}.coeff_err_mean_high;
    el = results{i}.coeff_err_low;
    err_low(i)   = mean(el(isfinite(el)));
    n_planes(i)  = results{i}.nPlanes;
    fprintf('%-10.3f | %-10.3f | %-12.4f | %-12.4f | %-12.4f | %-6d\n', ...
        alpha_list(i), dz_list_mm(i), rmse_vals(i), err_high(i), err_low(i), n_planes(i));
end

% 找到实验最优
valid = ~isnan(rmse_vals);
[~, idx_best_rmse] = min(rmse_vals(valid));
[~, idx_best_high] = min(err_high(valid));
valid_idx = find(valid);

fprintf('\n--- 最优值对比 ---\n');
fprintf('按RMSE_2D:    alpha=%.3f (dz=%.3fmm)\n', ...
    alpha_list(valid_idx(idx_best_rmse)), dz_list_mm(valid_idx(idx_best_rmse)));
fprintf('按Err_high:   alpha=%.3f (dz=%.3fmm)\n', ...
    alpha_list(valid_idx(idx_best_high)), dz_list_mm(valid_idx(idx_best_high)));
fprintf('理论预测:     alpha=%.3f (dz=%.3fmm)\n', a_opt, a_opt*z_c*1e3);
fprintf('理论(K=%d):   alpha=%.3f (dz=%.3fmm)\n', K_fixed, ...
    alpha_list(idx_opt_theory), dz_list_mm(idx_opt_theory));

%% ===== 绘图 =====
figure('Color','w','Position',[50,50,1400,550]);

% 子图1: RMSE vs alpha
subplot(1,3,1);
yyaxis left;
h1 = plot(alpha_list, rmse_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('RMSE_{2D} (rad)');
yyaxis right;
h2 = plot(alpha_list, err_high, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('高阶系数相对误差');
xline(a_opt, 'r--', 'LineWidth', 1.5);
xlabel('\alpha'); grid on;
title('恢复误差 vs \alpha (K=9固定)');
legend([h1, h2], {'RMSE_{2D}', 'Err_{high}'}, 'Location', 'best');

% 子图2: 低阶 vs 高阶 对比
subplot(1,3,2);
plot(alpha_list, err_low, 'g^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '低阶(r^2,r^4)');
hold on;
plot(alpha_list, err_high, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '高阶(r^6+)');
xline(a_opt, 'r--', 'LineWidth', 1.5);
xlabel('\alpha'); ylabel('相对误差'); grid on;
title('低阶 vs 高阶恢复误差');
legend('Location', 'best');

% 子图3: 理论 vs 实验
subplot(1,3,3);
% 归一化理论曲线
F_norm = F_theory / max(F_theory);
% 归一化实验(高阶误差的倒数, 越大越好)
E_norm = (1 ./ (err_high + 0.01));
E_norm = E_norm / max(E_norm(valid));
plot(alpha_list, F_norm, 'k-', 'LineWidth', 2, 'DisplayName', '理论 F_{total} (归一化)');
hold on;
plot(alpha_list(valid), E_norm(valid), 'rs-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', '实验 1/Err_{high} (归一化)');
xline(a_opt, 'r--', 'LineWidth', 1.5);
xlabel('\alpha'); ylabel('归一化指标'); grid on;
title('理论预测 vs 实验结果');
legend('Location', 'best');

sgtitle(sprintf(['高频最优步长验证 (系统A, K=%d): ', ...
    'z_c=%.2fmm, \\alpha_{opt}=\\pi/2=%.2f'], ...
    K_fixed, z_c*1e3, a_opt), 'FontSize', 14);

% 保存
save_dir = fileparts(mfilename('fullpath'));
save(fullfile(save_dir, 'sweep_alpha_v2_results.mat'), ...
    'alpha_list', 'dz_list_mm', 'rmse_vals', 'err_high', 'err_low', ...
    'n_planes', 'z_c', 'a_opt', 'K_fixed', 'F_theory', 'results');
fprintf('\n结果已保存\n');
