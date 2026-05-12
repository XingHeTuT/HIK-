% =========================================================================
% 双面超透镜理想合成相位计算 (Double-sided Metalens Ideal Phase Calculation)
% 基于角谱衍射理论 (Angular Spectrum Method)
% =========================================================================
clc; clear all; close all;

%% 1. 参数设置 (Parameters)
% --- 物理参数 (请确保与你的主代码完全一致) ---
lambda_um = 10.6; % 假设工作波长为 10.6 um，请根据实际修改
lambda = lambda_um * 1e-6;

% 衬底参数 (以你的主代码中的参数为准)
n_substrate = 3.4177615; 
d_substrate_mm = 0.725;    % 衬底厚度 0.5 mm
d_sub = d_substrate_mm * 1e-3;

% --- 仿真网格 ---
N_sim = 4*2048;
Grid_Window_mm = 44.64; 
dx_sim = (Grid_Window_mm * 1e-3) / N_sim; 

x_sim = ((1:N_sim) - N_sim/2 - 1) * dx_sim; 
[X_sim, Y_sim] = meshgrid(x_sim, x_sim);
Rho_sim = sqrt(X_sim.^2 + Y_sim.^2); 

% --- 孔径与光阑 ---
D_stop_mm = 44.64; 
radius_stop = (D_stop_mm * 1e-3) / 2; 

% 软边光阑 (防止衍射计算中的吉布斯振铃效应)
Soft_Mask = exp(- (Rho_sim / radius_stop).^60); 
Soft_Mask(Rho_sim > radius_stop * 1.01) = 0;

%% 2. 正反面 Zemax 系数定义 (需要你手动替换)
radius_norm_mm = 1.0; 
r_norm_z = double(Rho_sim) ./ (radius_norm_mm * 1e-3); 
r_sq_z = r_norm_z.^2;

% 【请在此替换为你真实的 Zemax Binary 2 面系数】
% 此处我使用两组测试系数，仅为演示代码逻辑
A_coeffs_S1 = [ -4.758394571634000E+000,  2.194016497211000E-002, -5.340807755416000E-005,  3.323255468294000E-008, 0, 0, 0, 0, 0]; % 正面系数 (S1)
A_coeffs_S2 = [ 5.318079102338000E+000,  -2.351626465626000E-002, 5.795669256364000E-005,  -3.663428239962000E-008, 0, 0, 0, 0, 0]; % 反面系数 (S2)

% 生成 S1 和 S2 的独立理论相位
Phi_1 = zeros(N_sim, N_sim);
Phi_2 = zeros(N_sim, N_sim);
for m = 1:length(A_coeffs_S1)
    Phi_1 = Phi_1 + A_coeffs_S1(m) * (r_sq_z.^m);
end
for m = 1:length(A_coeffs_S2)
    Phi_2 = Phi_2 + A_coeffs_S2(m) * (r_sq_z.^m);
end

%% 3. 核心计算：物理衍射过程模拟
fprintf('>>> 开始计算双面超透镜合成相位...\n');

% (1) 理想平面波垂直入射到 S1
U_in = Soft_Mask .* ones(N_sim, N_sim); 

% (2) 经过正面超透镜 (S1) 调制
U_after_S1 = U_in .* exp(1i * Phi_1);

% (3) 在衬底中自由传播
% 【核心物理考量】：光在介质中传播时，波长会缩短，必须除以折射率！
lambda_sub = lambda / n_substrate; 
fprintf('    正在模拟光在厚度 %.2f mm 的衬底中的传播...\n', d_substrate_mm);
U_before_S2 = ASM_Propagate(U_after_S1, d_sub, lambda_sub, dx_sim);

% (4) 经过反面超透镜 (S2) 调制
U_after_S2 = U_before_S2 .* exp(1i * Phi_2);

% (5) 提取最终理想合相位
Phi_total_ideal = angle(U_after_S2);

% 为了与测量对比，通常将参考相位平移，使其中心值为0
center_idx = floor(N_sim/2) + 1;
Phi_total_ideal = angle( exp(1i * (Phi_total_ideal - Phi_total_ideal(center_idx, center_idx))) );

fprintf('>>> 计算完成！\n\n');

%% 4. 可视化对比：直接相加 vs 物理传播 (平移归零修正版)
fprintf('>>> 正在生成相位对比图...\n');

% 提取中心剖面的物理坐标
x_mm = x_sim * 1e3;
center_idx = floor(N_sim/2) + 1;

% 直接代数相加的结果 (错误做法的演示)
Phi_algebraic_sum = Phi_1 + Phi_2;

% 提取 1D 剖面并解包裹
prof_algebraic = unwrap(Phi_algebraic_sum(center_idx, :));
prof_physical  = unwrap(Phi_total_ideal(center_idx, :));

% 【核心修正】：强制解包裹后的 1D 曲线在 r=0 (即 center_idx) 处严格为 0
prof_algebraic = prof_algebraic - prof_algebraic(center_idx);
prof_physical  = prof_physical - prof_physical(center_idx);

% 取有效光阑区域内的点进行绘制
valid_mask = abs(x_mm) <= (radius_stop * 1e3 * 0.95);

figure('Color', 'w', 'Position', [200, 200, 1000, 500]);

% --- 左图：绝对相位对比 ---
subplot(1,2,1);
plot(x_mm(valid_mask), prof_algebraic(valid_mask), 'b--', 'LineWidth', 2, 'DisplayName', '\Phi_1 + \Phi_2 (直接代数相加)');
hold on;
plot(x_mm(valid_mask), prof_physical(valid_mask), 'r-', 'LineWidth', 2, 'DisplayName', '物理传播合成相位');
grid on;
xlabel('Radius (mm)', 'FontWeight', 'bold');
ylabel('Phase (rad)', 'FontWeight', 'bold');
title('1D Phase Profile Comparison');
legend('Location', 'best');

% --- 右图：相位误差图 ---
subplot(1,2,2);
% 此时两者中心都严格为0，直接作差即可反映真实的曲率差异
phase_error = prof_physical - prof_algebraic; 
plot(x_mm(valid_mask), phase_error(valid_mask), 'k-', 'LineWidth', 2);
grid on;
xlabel('Radius (mm)', 'FontWeight', 'bold');
ylabel('\Delta Phase (rad)', 'FontWeight', 'bold');
title('衬底厚度引入的本征相位误差');
subtitle('中心点已强制归零对齐');

%% 5. 终极合相位等效 Zemax 系数拟合 (修正数值爆炸版)
fprintf('\n>>> 开始拟合物理传播合相位的等效 Zemax 系数 (修正数值稳定性)...\n');

% 1. 提取坐标与相位
r_fit_mm = abs(x_mm(valid_mask))';     % 物理半径 (mm, 列向量)
Z_target = prof_physical(valid_mask)'; % 展开相位 (列向量)

% [采纳你的建议] 强制 r=0 处的相位为 0
[~, center_idx_fit] = min(r_fit_mm);
Z_target = Z_target - Z_target(center_idx_fit); 

% 2. [核心修复] 坐标归一化 (防止 r^20 导致矩阵数值爆炸)
R_max = max(r_fit_mm);
rho_stable = r_fit_mm / R_max; % 将所有坐标映射到 [0, 1] 区间

% 构建稳定的拟合矩阵 M
N_terms = length(A_coeffs_S1);
M_stable = ones(length(rho_stable), N_terms + 1); 
for k = 1:N_terms
    M_stable(:, k+1) = rho_stable.^(2*k);
end

% 3. 稳定线性最小二乘拟合
C_stable = M_stable \ Z_target;

% 4. [系数还原] 将 [0,1] 区间的系数换算回 Zemax R_norm = 1.0 的标准
Piston_opt = C_stable(1);
A_coeffs_total_opt = zeros(1, N_terms);

for k = 1:N_terms
    % 换算公式: A_zemax = C_stable / (R_max / R_norm)^(2k)
    % 此处 R_norm = 1.0
    A_coeffs_total_opt(k) = C_stable(k+1) / (R_max^(2*k));
end

% 代数直接相加的系数 (作为参考)
A_coeffs_algebraic = A_coeffs_S1 + A_coeffs_S2;

% --- 打印高精度对比表格 ---
fprintf('========================================================================================\n');
fprintf('  双面超透镜等效合相位系数 (归一化半径 R_norm = %.1f mm)\n', radius_norm_mm);
fprintf('========================================================================================\n');
fprintf('  %-6s | %-22s | %-22s | %-20s\n', 'Term', '代数相加 (S1+S2)', '物理传播拟合 (True)', '偏差 (Difference)');
fprintf('----------------------------------------------------------------------------------------\n');

labels = {'r^2', 'r^4', 'r^6', 'r^8', 'r^10', 'r^12', 'r^14', 'r^16', 'r^18', 'r^20'};
for k = 1:N_terms
    if k <= length(labels)
        lb = labels{k}; 
    else
        lb = sprintf('r^%d', 2*k); 
    end
    val_alg = A_coeffs_algebraic(k);
    val_opt = A_coeffs_total_opt(k);
    diff_val = val_opt - val_alg;
    
    fprintf('  %-6s | %22.10e | %22.10e | %20.10e\n', lb, val_alg, val_opt, diff_val);
end
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('  Piston : %22.4f rad\n', Piston_opt);
fprintf('========================================================================================\n\n');

% 验证稳定拟合的精度
Z_fit_check = M_stable * C_stable; % 使用稳定矩阵回算
RMSE_fit = sqrt(mean((Z_target - Z_fit_check).^2));
fprintf('>>> 多项式拟合自身残差 RMSE: %.4e rad (越小说明偶次多项式模型越适用)\n', RMSE_fit);

%% --- 辅助函数：角谱传播 (Angular Spectrum Method) ---
function U_out = ASM_Propagate(U_in, z, lambda, dx)
    [Ny, Nx] = size(U_in); 
    df = 1 / (Nx * dx); 
    f = ((1:Nx) - Nx/2 - 1) * df;
    [FX, FY] = meshgrid(f, f); 
    
    term = 1 - lambda^2 * (FX.^2 + FY.^2); 
    term(term < 0) = 0; % 滤除倏逝波
    
    H = exp(1i * 2 * pi * z / lambda * sqrt(term));
    
    U_f = fftshift(fft2(ifftshift(U_in))); 
    U_out = fftshift(ifft2(ifftshift(U_f .* H)));
end