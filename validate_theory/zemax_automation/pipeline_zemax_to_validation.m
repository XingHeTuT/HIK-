% pipeline_zemax_to_validation.m  (v3: 9 alpha, full K, Z_span=+-5mm)
% 完整管道: Zemax导出txt -> xlsx转换 -> MDPR验证 -> 结果汇总
%
% 前置: run_zemax_batch.py已完成, PSF txt文件在 zemax_export/ 下

clc; clear all; %#ok<CLALL>

%% ===== 配置 =====
export_root = ['D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量\', ...
               'validate_theory\zemax_export'];

% 系统A参数
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

% alpha配置 (与 run_zemax_batch.py 一致)
alpha_list     = [0.10, 0.20, 0.50, 0.80, 1.20, pi/2, 2.00, 2.50, 3.00];
dz_list_mm     = [0.036, 0.073, 0.181, 0.290, 0.435, 0.570, 0.725, 0.906, 1.088];
output_subdirs = {'alpha_0.10_dz_0.036mm', 'alpha_0.20_dz_0.073mm', ...
                  'alpha_0.50_dz_0.181mm', 'alpha_0.80_dz_0.290mm', ...
                  'alpha_1.20_dz_0.435mm', 'alpha_1.57_dz_0.570mm', ...
                  'alpha_2.00_dz_0.725mm', 'alpha_2.50_dz_0.906mm', ...
                  'alpha_3.00_dz_1.088mm'};
K_expected     = [277, 138, 57, 36, 24, 19, 15, 13, 11];

% 理论参数
lambda  = sys.lambda_um * 1e-6;
D_ml    = sys.D_metalens_mm * 1e-3;
BFL     = sys.BFL_mm * 1e-3;
NA_eff  = (D_ml/2) / BFL;
z_c_mm  = lambda * 1e3 / (pi * NA_eff^2);
a_opt   = pi/2;

fprintf('===== 管道: Zemax导出 -> 验证 =====\n');
fprintf('z_c = %.4f mm, alpha_opt = %.3f (%.3f mm)\n\n', z_c_mm, a_opt, a_opt*z_c_mm);

%% ===== 第1步: txt -> xlsx 批量转换 =====
fprintf('>>> 第1步: txt -> xlsx 转换\n');
fprintf('注意: 等待时间长, 请勿中断...\n\n');

for i = 1:length(output_subdirs)
    txt_dir = fullfile(export_root, output_subdirs{i}, 'Zemax_txt');
    xlsx_dir = fullfile(export_root, output_subdirs{i}, 'Zemax_xlsx');

    if ~exist(txt_dir, 'dir')
        fprintf('  [%d/%d] 跳过 (txt目录不存在): %s\n', i, length(output_subdirs), txt_dir);
        continue;
    end

    if ~exist(xlsx_dir, 'dir'), mkdir(xlsx_dir); end

    files = dir(fullfile(txt_dir, '*.txt'));
    n_files = numel(files);
    fprintf('  [%d/%d] %s: %d txt文件 (预期%d)...', ...
        i, length(output_subdirs), output_subdirs{i}, n_files, K_expected(i));

    if n_files < K_expected(i)
        fprintf(' 警告: 文件数少于预期!');
    end

    n_ok = 0;
    t_start = tic;

    for j = 1:n_files
        txt_path = fullfile(txt_dir, files(j).name);
        [~, base, ~] = fileparts(files(j).name);
        xlsx_path = fullfile(xlsx_dir, [base, '.xlsx']);

        try
            M = readmatrix(txt_path, 'FileType', 'text');
            [rows, cols] = size(M);

            if ~isempty(M)
                % 清洗逻辑: 去除序号列/行
                if size(M,2) > 1
                    col1 = M(:,1);
                    d = diff(col1);
                    if (~isempty(d) && all(abs(d(~isnan(d))-1) < 1e-4)) || ...
                       all(abs(col1(~isnan(col1))) < 1e-4) || all(isnan(col1))
                        M(:,1) = [];
                    end
                end
                if size(M,1) > 1
                    row1 = M(1,:);
                    d = diff(row1);
                    if (~isempty(d) && all(abs(d(~isnan(d))-1) < 1e-4)) || ...
                       all(abs(row1(~isnan(row1))) < 1e-4) || all(isnan(row1))
                        M(1,:) = [];
                    end
                end
            end

            % 验证尺寸 (期望512x512)
            [r, c] = size(M);
            if r ~= 512 || c ~= 512
                fprintf('\n    警告: %s 尺寸为 %dx%d (期望512x512)\n', files(j).name, r, c);
            end

            if exist(xlsx_path, 'file'), delete(xlsx_path); end
            writematrix(M, xlsx_path);
            n_ok = n_ok + 1;

        catch ME
            fprintf('\n    错误: %s 转换失败: %s\n', files(j).name, ME.message);
        end

        if mod(j, 50) == 0
            fprintf('\n    进度: %d/%d', j, n_files);
        end
    end

    elapsed = toc(t_start);
    fprintf(' 完成 (%d OK, %.1fs)\n', n_ok, elapsed);
end

%% ===== 第2步: MDPR验证扫描 =====
fprintf('\n>>> 第2步: MDPR验证扫描 (预计6-8小时)\n');

alg = struct();
alg.max_iter    = 30;
alg.alpha_start = 1.0;
alg.alpha_end   = 0.8;
alg.beta_amp    = 0.8;

results = cell(length(alpha_list), 1);
t_total = tic;

for i = 1:length(alpha_list)
    dataFolder = fullfile(export_root, output_subdirs{i}, 'Zemax_xlsx');

    if ~exist(dataFolder, 'dir')
        fprintf('  [%d/%d] 跳过: xlsx目录不存在\n', i, length(alpha_list));
        continue;
    end

    % 验证xlsx文件数
    xlsx_files = dir(fullfile(dataFolder, '*.xlsx'));
    n_xlsx = numel(xlsx_files);
    fprintf('\n  [%d/%d] alpha=%.3f, dz=%.3fmm, K=%d (预期%d)\n', ...
        i, length(alpha_list), alpha_list(i), dz_list_mm(i), n_xlsx, K_expected(i));

    sys.dataFolder = dataFolder;

    % 扫描范围: 使用全范围 [-5, 5] mm
    z_range_mm = [-5.0, 5.0];

    samp = struct();
    samp.delta_z_mm  = dz_list_mm(i);
    samp.z_range_mm  = z_range_mm;
    samp.max_iter    = alg.max_iter;
    samp.alpha_start = alg.alpha_start;
    samp.alpha_end   = alg.alpha_end;
    samp.beta_amp    = alg.beta_amp;

    t_run = tic;
    try
        r = run_mdpr_validation(sys, samp);
        results{i} = r;
        elapsed = toc(t_run);
        fprintf('  -> RMSE=%.4f, Err_high=%.4f, Err_low=%.4f, K=%d, time=%.0fs\n', ...
            r.RMSE_2D, r.coeff_err_mean_high, mean(r.coeff_err_low,'omitnan'), ...
            r.nPlanes, elapsed);
    catch ME
        fprintf('  *** 出错: %s\n', ME.message);
        results{i} = struct('RMSE_2D', NaN, 'coeff_err_mean_high', NaN, ...
            'coeff_err_low', [NaN,NaN], 'nPlanes', 0);
    end

    % 中间保存 (防止崩溃丢失数据)
    save_dir = fileparts(mfilename('fullpath'));
    save(fullfile(save_dir, 'pipeline_intermediate.mat'), ...
        'results', 'alpha_list', 'i');
end

total_elapsed = toc(t_total);
fprintf('\n总耗时: %.1f 小时\n', total_elapsed/3600);

%% ===== 第3步: 结果汇总与绘图 =====
fprintf('\n===== 最终结果 =====\n');

rmse_vals = zeros(length(alpha_list), 1);
err_high  = zeros(length(alpha_list), 1);
err_low   = zeros(length(alpha_list), 1);
k_actual  = zeros(length(alpha_list), 1);

for i = 1:length(alpha_list)
    rmse_vals(i) = results{i}.RMSE_2D;
    err_high(i)  = results{i}.coeff_err_mean_high;
    el = results{i}.coeff_err_low;
    err_low(i)   = mean(el(isfinite(el)));
    k_actual(i)  = results{i}.nPlanes;
end

valid = ~isnan(rmse_vals);

fprintf('%-10s | %-10s | %-6s | %-12s | %-12s | %-12s\n', ...
    'alpha', 'dz[mm]', 'K', 'RMSE_2D', 'Err_high', 'Err_low');
fprintf('%s\n', repmat('-', 1, 72));
for i = 1:length(alpha_list)
    fprintf('%-10.3f | %-10.3f | %-6d | %-12.4f | %-12.4f | %-12.4f\n', ...
        alpha_list(i), dz_list_mm(i), k_actual(i), rmse_vals(i), err_high(i), err_low(i));
end

% 最优对比
valid_idx = find(valid);
[~, idx_best] = min(err_high(valid));
fprintf('\n--- 最优值 ---\n');
fprintf('实验最优 (Err_high): alpha=%.3f (dz=%.3fmm)\n', ...
    alpha_list(valid_idx(idx_best)), dz_list_mm(valid_idx(idx_best)));
fprintf('理论预测:            alpha=%.3f (dz=%.3fmm)\n', a_opt, a_opt*z_c_mm);

% 计算理论Fisher信息
F_theory = zeros(size(alpha_list));
for i = 1:length(alpha_list)
    K = k_actual(i);
    if K == 0, continue; end
    M_2 = (K - 1) / 2;
    for m = -M_2:M_2
        F_theory(i) = F_theory(i) + sin(m * alpha_list(i))^2;
    end
end

%% ===== 绘图 =====
figure('Color','w','Position',[50,50,1600,500]);

% 子图1: RMSE vs alpha
subplot(1,3,1);
yyaxis left;
plot(alpha_list, rmse_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('RMSE_{2D} (rad)');
yyaxis right;
plot(alpha_list, err_high, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Err_{high}');
xline(a_opt, 'r--', 'LineWidth', 1.5, 'Label', '\alpha_{opt}=\pi/2');
xlabel('\alpha'); grid on;
title('Recovery Error vs \alpha (Z_{span}=10mm)');

% 子图2: 低阶 vs 高阶
subplot(1,3,2);
plot(alpha_list, err_low, 'g^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Low-order');
hold on;
plot(alpha_list, err_high, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'High-order');
xline(a_opt, 'r--', 'LineWidth', 1.5);
xlabel('\alpha'); ylabel('Relative Error'); grid on;
title('Low vs High Order Error');
legend('Location', 'best');

% 子图3: 理论 vs 实验
subplot(1,3,3);
F_norm = F_theory / max(F_theory(valid));
E_norm = 1 ./ (err_high + 0.01);
E_norm = E_norm / max(E_norm(valid));
plot(alpha_list, F_norm, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory: F_{total}');
hold on;
plot(alpha_list(valid), E_norm(valid), 'rs-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Experiment: 1/Err_{high}');
xline(a_opt, 'r--', 'LineWidth', 1.5);
xlabel('\alpha'); ylabel('Normalized'); grid on;
title('Theory vs Experiment');
legend('Location', 'best');

sgtitle(sprintf(['Phase Retrieval Error vs \\alpha  ', ...
    '(System A, Z_{span}=10mm, z_c=%.2fmm)'], z_c_mm), 'FontSize', 14);

% 保存
save_dir = fileparts(mfilename('fullpath'));
save(fullfile(save_dir, 'pipeline_final_results.mat'), ...
    'alpha_list', 'dz_list_mm', 'rmse_vals', 'err_high', 'err_low', ...
    'k_actual', 'z_c_mm', 'a_opt', 'F_theory', 'results');
fprintf('\n最终结果已保存\n');
