% generate_alpha_configs.m
% 系统A (D50.4f200+ML1): 生成9组alpha方案的离焦位置
%
% v3: 固定Z_span=±5mm, K不压缩, alpha覆盖最优值两侧

clc;

%% 系统A参数
D_ml_mm   = 44.5977;
BFL_mm    = 231.1605;
lambda_um = 10.6;
NA_eff    = (D_ml_mm/2) / BFL_mm;
z_c_mm    = lambda_um * 1e-3 / (pi * NA_eff^2);

fprintf('系统A: z_c = %.4f mm, NA = %.4f\n', z_c_mm, NA_eff);
fprintf('理论最优: alpha_opt = pi/2 = %.4f, dz_opt = %.4f mm\n\n', pi/2, pi/2*z_c_mm);

%% 实验设计: 固定Z_span=10mm, 9个alpha值
Z_span_mm = 10;
z_min = -Z_span_mm/2;
z_max = +Z_span_mm/2;

alpha_list = [0.10, 0.20, 0.50, 0.80, 1.20, pi/2, 2.00, 2.50, 3.00];

fprintf('固定扫描范围: [%.1f, %.1f] mm (Z_span=%.0fmm, beta=%.1f)\n\n', ...
    z_min, z_max, Z_span_mm, Z_span_mm/z_c_mm);

fprintf('%-10s | %-10s | %-6s\n', 'alpha', 'dz[mm]', 'K');
fprintf('%s\n', repmat('-', 1, 30));

configs = cell(length(alpha_list), 1);

for i = 1:length(alpha_list)
    a = alpha_list(i);
    dz_mm = a * z_c_mm;

    % 生成离焦位置: 从z_min到z_max, 步长dz_mm
    dz_positions_mm = z_min : dz_mm : z_max;
    K = length(dz_positions_mm);

    % 确保包含0 (最佳聚焦面)
    if ~any(abs(dz_positions_mm) < dz_mm/10)
        dz_positions_mm = sort(unique([dz_positions_mm, 0]));
        K = length(dz_positions_mm);
    end

    cfg = struct();
    cfg.alpha = a;
    cfg.dz_mm = dz_mm;
    cfg.z_span_mm = Z_span_mm;
    cfg.K = K;
    cfg.dz_positions_mm = dz_positions_mm;
    cfg.output_dir = sprintf('alpha_%.2f_dz_%.3fmm', a, dz_mm);

    configs{i} = cfg;

    fprintf('%-10.3f | %-10.3f | %-6d\n', a, dz_mm, K);

    % 打印前5个和最后5个位置
    if K <= 10
        fprintf('  位置: ');
        fprintf('%.3f ', dz_positions_mm);
        fprintf('\n');
    else
        fprintf('  位置: ');
        fprintf('%.3f ', dz_positions_mm(1:min(3,K)));
        fprintf(' ... ');
        fprintf('%.3f ', dz_positions_mm(max(1,K-2):end));
        fprintf('\n');
    end
end

%% 理论Fisher信息预测
fprintf('\n===== 理论预测 =====\n');
fprintf('%-10s | %-6s | %-12s\n', 'alpha', 'K', 'F_total');
fprintf('%s\n', repmat('-', 1, 32));

for i = 1:length(alpha_list)
    a = alpha_list(i);
    K = configs{i}.K;
    M = (K - 1) / 2;

    % F_total ∝ Σ sin²(m*alpha) for symmetric planes about focus
    F_total = 0;
    for m = -M:M
        F_total = F_total + sin(m * a)^2;
    end

    fprintf('%-10.3f | %-6d | %-12.3f\n', a, K, F_total);
end

%% 保存
save_dir = fileparts(mfilename('fullpath'));
save(fullfile(save_dir, 'alpha_configs.mat'), 'configs', 'alpha_list', 'z_c_mm', 'Z_span_mm');
fprintf('\n配置已保存\n');

%% 打印Python脚本可用的格式
fprintf('\n===== Python可用格式 =====\n');
fprintf('alpha_configs = [\n');
for i = 1:length(alpha_list)
    c = configs{i};
    fprintf('    {  # alpha=%.3f, K=%d\n', c.alpha, c.K);
    fprintf('        "alpha": %.4f,\n', c.alpha);
    fprintf('        "dz_mm": %.3f,\n', c.dz_mm);
    fprintf('        "output_dir": "%s",\n', c.output_dir);
    fprintf('        "dz_positions_mm": [\n');
    for j = 1:c.K
        fprintf('            %.3f', c.dz_positions_mm(j));
        if j < c.K, fprintf(','); end
        fprintf('\n');
    end
    fprintf('        ],\n');
    fprintf('    },\n');
end
fprintf(']\n');

fprintf('\n总导出文件数: %d\n', sum(cellfun(@(c) c.K, configs)));
