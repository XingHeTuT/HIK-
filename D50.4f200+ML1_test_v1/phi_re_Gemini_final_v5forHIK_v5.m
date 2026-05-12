clc; clear all; close all;

%% 1. 参数设置 (Parameters)
% --- 物理参数 ---
lambda_um = 10.6;
lambda = lambda_um * 1e-6;

% --- 系统几何参数 (Zemax 提取值) ---
n_lens = 2.40266; 
n_substrate = 3.422; 

% 1. 折射透镜
R1_mm = 2.805400000000000E+002;  % 前表面 (凸)
R2_mm = Inf;     % 后表面 (平)
Tc_lens_mm = 3.500000000000000E+000; % 厚度

% 2. 间隙与基底
d_air_gap_mm = 2.305480752288749E+001; % 紧贴
d_substrate_mm = 7.250000000000000E-001;

% 3. 后截距 (BFL)
BFL_mm = 2.311604993761137E+002; 
BFL = BFL_mm * 1e-3;

% --- [关键] 孔径定义 ---
% 物理含义拆分：
% 1) 第一面就是系统 STOP，其口径 = 折射透镜口径
% 2) 第四面是超透镜面，其通光范围 = 超透镜有效口径
%
% 后续约定：
%   - 第一面发光线 / stop 裁剪：统一使用 radius_stop1
%   - 第四面超透镜有效区域 / 反演口径 / 结果评价：统一使用 radius_metalens

% 第一面（折射透镜 / STOP）口径
D_stop1_mm   = 2.540000000000000E+001 * 2;   % 50.8 mm
radius_stop1 = (D_stop1_mm * 1e-3) / 2;      % m

% 第四面（超透镜）有效口径
D_metalens_mm   = 2.229883854002862E+001 * 2; % 44.59767708 mm
radius_metalens = (D_metalens_mm * 1e-3) / 2; % m

% --- 仿真网格 ---
N_sim = 2*2048;
Grid_Window_mm = 46; 
Grid_Window = Grid_Window_mm * 1e-3; 
dx_sim = Grid_Window / N_sim; 

% --- 探测器参数 (0.425um 分辨率) ---
N_det = 512; 
det_pitch_um = 0.425; 
L_det = N_det * det_pitch_um * 1e-6; % 探测器物理总宽度 (0.3072 mm)

% --- Zemax 理想系数 ---
radius_norm_mm = 1.0; 
A_coeffs_ideal = [5.4342456759e-01, -1.2748303966e-03, 2.0384613231e-06, 7.9430544489e-09, ...
    -2.8816204931e-11, 3.7742704706e-14, -2.7815489016e-17, 1.6979040962e-20, -7.6949989033e-24];

max_iter = 30;
% 请修改为您存放 xlsx 文件的文件夹路径
dataFolder = 'D:\LiLinhan_2026\相位反演\HIK相位测量\D50.4f200+ML1_test_v1\Zemax_xlsx'; 

%% 2. 正向仿真：第一面为 STOP，第四面为超透镜有效口径
fprintf('>>> 启动光线追踪 (构建物理入射场, STOP at Surface 1)...\n');

% -------------------------------------------------------------------------
% 1. 第一面发射光线：按第一面 STOP 口径发射
% -------------------------------------------------------------------------
r_launch = radius_stop1 * 1.05;
N_rays_linear = 5000;
ln = linspace(-r_launch, r_launch, N_rays_linear);
[Ray_X, Ray_Y] = meshgrid(ln, ln);

% 第一面 STOP 裁剪
valid_rays_stop1 = (Ray_X.^2 + Ray_Y.^2) <= radius_stop1^2;
rx = Ray_X(valid_rays_stop1);
ry = Ray_Y(valid_rays_stop1);
rz = zeros(size(rx));

dir_x = zeros(size(rx));
dir_y = zeros(size(rx));
dir_z = ones(size(rx));
OPL   = zeros(size(rx));

% -------------------------------------------------------------------------
% 2. 光线追迹
% -------------------------------------------------------------------------
% Surface 1: 球面 (凸)
R1 = R1_mm * 1e-3;
Center1 = [0, 0, R1];
[P1, t1] = intersect_sphere(rx, ry, rz, dir_x, dir_y, dir_z, Center1, R1);
OPL = OPL + t1 * 1.0;
[dir_x, dir_y, dir_z] = refract_sphere(P1, dir_x, dir_y, dir_z, Center1, R1, 1.0, n_lens, -1);

% Surface 2: 平面
Tc = Tc_lens_mm * 1e-3;
[P2, t2] = intersect_plane(P1(:,1), P1(:,2), P1(:,3), dir_x, dir_y, dir_z, Tc);
OPL = OPL + t2 * n_lens;
[dir_x, dir_y, dir_z] = refract_plane(dir_x, dir_y, dir_z, n_lens, 1.0);

% Surface 3: 超透镜前表面 / 基底前表面
z3 = Tc + d_air_gap_mm * 1e-3;
[P3, t3] = intersect_plane(P2(:,1), P2(:,2), P2(:,3), dir_x, dir_y, dir_z, z3);
OPL = OPL + t3 * 1.0;
[dir_x, dir_y, dir_z] = refract_plane(dir_x, dir_y, dir_z, 1.0, n_substrate);

% Surface 4: 超透镜后表面（相位调制后出射面所在重建面）
z4 = z3 + d_substrate_mm * 1e-3;
[P4, t4] = intersect_plane(P3(:,1), P3(:,2), P3(:,3), dir_x, dir_y, dir_z, z4);
OPL = OPL + t4 * n_substrate;

% -------------------------------------------------------------------------
% 3. 第四面按“超透镜有效口径”筛选
% -------------------------------------------------------------------------
R4 = sqrt(P4(:,1).^2 + P4(:,2).^2);
valid_rays_ml = (R4 <= radius_metalens);

P4_valid  = P4(valid_rays_ml, :);
OPL_valid = OPL(valid_rays_ml);

fprintf(' - 第一面 STOP 半径          = %.4f mm\n', radius_stop1 * 1e3);
fprintf(' - 第四面超透镜有效半径     = %.4f mm\n', radius_metalens * 1e3);
fprintf(' - 参与第四面波前重建的光线数 = %d\n', nnz(valid_rays_ml));

if nnz(valid_rays_ml) < 100
    error('第四面落入超透镜有效口径的有效光线过少，请检查口径或光路参数。');
end

% -------------------------------------------------------------------------
% 4. 在第四面重建波前
% -------------------------------------------------------------------------
fprintf(' - 在第四面重建波前并施加超透镜有效口径软边约束...\n');

x_sim = ((1:N_sim) - N_sim/2 - 1) * dx_sim;
[X_sim, Y_sim] = meshgrid(x_sim, x_sim);
Rho_sim = sqrt(X_sim.^2 + Y_sim.^2);

OPL_rel = OPL_valid - min(OPL_valid);
F_interp = scatteredInterpolant(P4_valid(:,1), P4_valid(:,2), OPL_rel, 'natural', 'linear');
OPD_grid = F_interp(X_sim, Y_sim);
OPD_grid(~isfinite(OPD_grid)) = 0;

Phi_Inc = (2 * pi / lambda) * OPD_grid;

% 第四面使用“超透镜有效口径”作为软边包络
Soft_Mask = exp(- (Rho_sim / radius_metalens).^60);
Soft_Mask(Rho_sim > radius_metalens * 1.01) = 0;

U_inc = Soft_Mask .* exp(1i * Phi_Inc);

% 归一化入射场能量
U_inc = U_inc / sqrt(sum(abs(U_inc(:)).^2));

%% 3. 数据读取与处理 (读取原始 unnormalized xlsx)
files = dir(fullfile(dataFolder, '*.xlsx'));
if isempty(files), error('未找到 .xlsx 文件，请检查路径。'); end

stack_amp = zeros(N_sim, N_sim, numel(files), 'single'); 
dz_list = zeros(numel(files), 1);

% 物理尺寸映射：探测器覆盖的仿真像素数
N_roi = round(L_det / dx_sim); 
if mod(N_roi, 2) ~= 0, N_roi = N_roi + 1; end 

sim_center = floor(N_sim/2) + 1;
valid_count = 0;

fprintf('>>> 读取原始强度数据 (Raw Intensity)...\n');

for i = 1:numel(files)
    fname = files(i).name;
    tok = regexp(fname, 'PSF_([-+0-9.eE]+)mm', 'tokens', 'once');
    if isempty(tok), continue; end
    
    try
        full_path = fullfile(dataFolder, fname);
        I_raw = readmatrix(full_path); 
    catch
        warning('文件读取失败: %s', fname); continue;
    end
    
    % --- 尺寸检查与调整 ---
    [rows, cols] = size(I_raw);
    if rows ~= N_det || cols ~= N_det
        % 简单裁切或缩放以匹配 N_det
        if rows > N_det && cols > N_det
             cr = floor(rows/2)+1; cc = floor(cols/2)+1;
             I_raw = I_raw(cr-256:cr+255, cc-256:cc+255);
        else
             I_raw = imresize(I_raw, [N_det, N_det]);
        end
    end
    
    % --- [关键] 原始数据清洗 ---
    I_raw = double(I_raw);
    I_raw(isnan(I_raw)) = 0;
    
    % 1. 去除负值 (物理上不可能)
    I_raw(I_raw < 0) = 0;
    
    % 2. 背景底噪抑制 (Background Denoising)
    % 对于未归一化数据，这很重要。假设小于峰值 0.1% 的是底噪。
    bg_thresh = 0.001 * max(I_raw(:));
    I_raw(I_raw < bg_thresh) = 0;
    
    % --- 质心对齐 ---
    % 即使数据未归一化，质心计算依然有效
    sum_I = sum(I_raw(:));
    if sum_I > 0
        [rs, cs] = ndgrid(1:size(I_raw,1), 1:size(I_raw,2));
        cy = sum(rs(:).*I_raw(:)) / sum_I;
        cx = sum(cs(:).*I_raw(:)) / sum_I;
        sy = (size(I_raw,1)/2 + 0.5) - cy; 
        sx = (size(I_raw,2)/2 + 0.5) - cx;
        I_raw = imtranslate(I_raw, [sx, sy], 'bilinear'); % 使用双线性插值平移
    end
    
    % --- 映射到仿真网格 ---
    % 注意：这里不再除以 max()，保留原始强度的相对大小
    % 将 'bicubic' 改为 'bilinear'
    % I_resized = imresize(I_raw, [N_roi, N_roi], 'bicubic');
    I_resized = imresize(I_raw, [N_roi, N_roi], 'bilinear');
    I_resized(I_resized < 0) = 0; % 插值可能产生负值，需去除
    
    % 填入全场矩阵
    Full_Grid = zeros(N_sim, N_sim, 'single');
    half = floor(N_roi/2);
    idx = (sim_center - half) : (sim_center - half + N_roi - 1);
    Full_Grid(idx, idx) = single(I_resized);
    
    valid_count = valid_count + 1;
    stack_amp(:,:,valid_count) = sqrt(Full_Grid); % 存储振幅 (Intensity开根号)
    dz_list(valid_count) = str2double(tok{1}) * 1e-3; 
    
    if mod(i, 5) == 0, fprintf('  已读取: %s (Peak Intensity: %.1f)\n', fname, max(I_raw(:))); end
end

% 排序与整理
stack_amp = stack_amp(:,:,1:valid_count);
dz_list = dz_list(1:valid_count);
[dz_list, idx] = sort(dz_list);
stack_amp = stack_amp(:,:,idx);
prop_dists = BFL + dz_list; 
nPlanes = valid_count;

fprintf('成功加载 %d 个平面。\n', nPlanes);

%% 4. 相位反演迭代 (Energy Constrained MDPR)
fprintf('使用 Zemax 理想相位热启动...\n');

% 理想相位初始猜测
r_norm_z = double(Rho_sim) ./ (radius_norm_mm * 1e-3); 
phi_init = zeros(N_sim, N_sim);
r_sq_z = r_norm_z.^2;
for m = 1:9, phi_init = phi_init + A_coeffs_ideal(m) * (r_sq_z.^m); end

U_total_est = U_inc .* exp(1i * phi_init);

% 权重策略
alpha_start = 1.0;
alpha_end = 0.8; % 稍微放松一点，允许算法微调振幅以匹配物理衍射

fprintf('开始迭代 (Energy Constrained)...\n');

for iter = 1:max_iter
    alpha = alpha_start - (alpha_start - alpha_end) * (iter / max_iter);
    U_accum = zeros(N_sim, N_sim);
    
    % 获取当前估计场的总能量 (理论上应该是守恒的，但数值计算会有微小波动)
    E_source = sum(abs(U_inc(:)).^2); 
    
    for k = 1:nPlanes
        dist = prop_dists(k);
        meas_amp_raw = double(stack_amp(:,:,k));
        
        % 1. 前向传播
        U_prop = ASM_Propagate(U_total_est, dist, lambda, dx_sim);
        
        % 2. [核心能量约束] Energy Normalization
        % 计算当前传播场的总能量
        E_calc = sum(abs(U_prop(:)).^2);
        % 计算原始测量数据的总能量
        E_meas = sum(meas_amp_raw(:).^2);
        
        % 计算缩放因子：将测量数据强制拉伸到与物理仿真能量一致
        if E_meas > 0
            Scale_Factor = sqrt(E_calc / E_meas);
        else
            Scale_Factor = 0;
        end
        meas_amp_scaled = meas_amp_raw * Scale_Factor;
        
        % 3. 加权幅度替换
        % 此时 meas_amp_scaled 和 abs(U_prop) 处于同一能量量级，可以直接加权
        Amp_updated = alpha * meas_amp_scaled + (1 - alpha) * abs(U_prop);
        U_corrected = Amp_updated .* exp(1i * angle(U_prop));
        
        % 4. 反向传播
        U_back = ASM_Propagate(U_corrected, -dist, lambda, dx_sim);
        U_accum = U_accum + U_back;
    end
    
    % 平均更新
    U_total_temp = U_accum / nPlanes;
    
    % % 物体域约束 (振幅重置为 U_inc，仅保留相位更新)
    % U_total_new = abs(U_inc) .* exp(1i * angle(U_total_temp));
    
    % --- [修改版] 放松约束策略 (Relaxed Constraint) ---
    % 1. 计算反演回来的振幅
    Amp_back = abs(U_total_temp);
    
    % 2. 混合策略：在光阑内，允许振幅根据实验数据微调；在光阑外，强制为0
    %    Soft_Mask 已经在 U_inc 里了，我们利用它
    Mask_Binary = double(Rho_sim <= radius_metalens);    
    % 定义由于实验光束可能是高斯，而仿真由于RayTracing可能是平顶，
    % 我们混合一下：80% 信任理论模型，20% 信任反演结果（为了消化高斯分布差异）
    beta = 0.8; 
    Amp_Updated = beta * abs(U_inc) + (1 - beta) * Amp_back;
    
    % 3. 施加光阑约束 (光阑外必须无光)
    Amp_Final = Amp_Updated .* Mask_Binary; 
    
    % 4. 更新
    U_total_new = Amp_Final .* exp(1i * angle(U_total_temp));
    
    % 能量归一化保持不变
    U_total_new = U_total_new / sqrt(sum(abs(U_total_new(:)).^2)) * sqrt(E_source);
    
    err = norm(U_total_new(:) - U_total_est(:), 'fro') / norm(U_total_est(:), 'fro');
    U_total_est = U_total_new;
    
    if mod(iter, 10) == 0
        fprintf('Iter %d: Err = %.4f\n', iter, err); 
    end
end

%% 5. 结果验证 (改进版：解包裹 + 鲁棒的区域对齐)
fprintf('\n>>> Step 5: 结果验证与绘图...\n');

% 1. 提取原始反演相位
U_metalens_extracted = U_total_est .* conj(U_inc);
% --- [新增] 相位平滑 (去除中心非物理的高频毛刺) ---
% 使用一个小的高斯滤波器对复振幅进行平滑，比直接平滑相位更安全
H_smooth = fspecial('gaussian', [5 5], 1.0); % 5x5窗口，sigma=1
U_metalens_extracted = imfilter(U_metalens_extracted, H_smooth, 'replicate');
% 理想超透镜复振幅
U_ideal = Soft_Mask .* exp(1i * phi_init);

% --- 2. 系统误差校正 (拟合去除 Piston, Tilt, Defocus) ---
[Xg, Yg] = meshgrid(x_sim, x_sim);
U_diff = U_metalens_extracted .* conj(U_ideal);
phase_diff_map = angle(U_diff);

% 掩膜
Valid_Mask = (Rho_sim <= 0.98 * radius_metalens);
Norm_Radius = radius_metalens;

x_val = Xg(Valid_Mask);
y_val = Yg(Valid_Mask);
z_val = phase_diff_map(Valid_Mask);

% 归一化坐标
x_n = x_val / Norm_Radius;
y_n = y_val / Norm_Radius;
r2_n = x_n.^2 + y_n.^2;

% [拟合矩阵] 仅去除: Piston, Tilt, Defocus
A_fit = [ones(size(x_n)), x_n, y_n, r2_n]; 
coeffs = A_fit \ z_val;

% 构建校正相位
X_norm = Xg / Norm_Radius;
Y_norm = Yg / Norm_Radius;
R2_norm = X_norm.^2 + Y_norm.^2;

Correction_Phase = coeffs(1) + ...             % Piston
                   coeffs(2)*X_norm + ...      % Tilt X
                   coeffs(3)*Y_norm + ...      % Tilt Y
                   coeffs(4)*R2_norm;          % Defocus

% 得到修正后的复振幅 (包含真实球差)
U_final = U_metalens_extracted .* exp(-1i * Correction_Phase);

% --- 3. 准备绘图数据 ---
mid = floor(N_sim/2) + 1;
x_mm = double(x_sim) * 1e3;

% (A) 提取 1D 剖面
phase_ideal_1d = angle(U_ideal(mid, :));
phase_final_1d = angle(U_final(mid, :));

% (B) 解包裹 (Unwrap)
p_ideal_unwrap = unwrap(phase_ideal_1d);
p_final_unwrap = unwrap(phase_final_1d);

% (C) [核心改进] 鲁棒的中心区域对齐 (Robust Center Alignment)
% 定义中心对齐区域：例如半径 0.2mm 内的区域
align_radius = 0.2; % mm
center_mask_1d = (abs(x_mm) <= align_radius);

if sum(center_mask_1d) > 0
    % 计算该区域内的相位均值
    mean_ideal = mean(p_ideal_unwrap(center_mask_1d));
    mean_final = mean(p_final_unwrap(center_mask_1d));
    
    % 计算偏移量
    offset = mean_final - mean_ideal;
    
    fprintf('采用区域对齐 (半径 %.1fmm, 包含 %d 个点)\n', align_radius, sum(center_mask_1d));
else
    % 如果像素太少，回退到单点对齐
    offset = p_final_unwrap(mid) - p_ideal_unwrap(mid);
    fprintf('采用单点对齐 (中心区域点数不足)\n');
end

% 应用对齐
p_final_aligned = p_final_unwrap - offset;

% (D) 计算残差 RMSE (剔除解包裹常数漂移的影响)
% 使用 exp(1i*diff) 计算误差，不仅准确，而且不受 2pi 跳变影响
complex_error = U_final .* conj(U_ideal); 
% 这里的 angle 会自动处理到 [-pi, pi]，实际上就是最小相位差
phase_error_val = angle(complex_error); 

% 计算 RMSE
rmse = sqrt(mean(phase_error_val(Valid_Mask).^2));
fprintf('最终残差 RMSE = %.4f rad\n', rmse);

% --- 4. 绘图展示 ---
figure('Color','w','Position',[100,100,1400,600]);

% [左图] 解包裹后的相位对比
subplot(1, 2, 1);
plot(x_mm, p_ideal_unwrap, 'k', 'LineWidth', 2, 'DisplayName', 'Ideal Design'); hold on;
plot(x_mm, p_final_aligned, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Extracted (Robust Aligned)');

% 标出对齐区域示意
x_fill = [-align_radius, align_radius, align_radius, -align_radius];
y_lims = ylim;
y_fill = [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];
patch(x_fill, y_fill, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Alignment Region');

xlim([-D_metalens_mm/2*1.1, D_metalens_mm/2*1.1]);
grid on; legend('Location', 'best');
xlabel('Position (mm)'); ylabel('Phase (rad)');
title(sprintf('Unwrapped Profile (RMSE = %.3f rad)', rmse));
subtitle('绿色区域为自动对齐基准区，消除中心震荡影响');

% [右图] 2D 残差图
subplot(1, 2, 2);
Err_Map = zeros(N_sim, N_sim);
Err_Map(Valid_Mask) = phase_error_val(Valid_Mask); 
Err_Map(~Valid_Mask) = NaN;

h = imagesc(x_mm, x_mm, Err_Map);
set(h, 'AlphaData', ~isnan(Err_Map));
axis image; colormap jet; colorbar; 
caxis([-1, 1]); 
xlabel('x (mm)'); ylabel('y (mm)');
title('Residual Error Map');

% 打印信息
fprintf('--------------------------------------\n');
fprintf('对齐平移量 (Piston Shift): %.4f rad\n', offset);
fprintf('离焦校正量: %.4e, 倾斜校正: %.4e, %.4e\n', coeffs(4), coeffs(2), coeffs(3));
fprintf('--------------------------------------\n');

%% 7. 终极优化：严格圆对称拟合（残差域解包裹自洽版）
% 核心思想：
% 1) Step 5 的 U_final 已经是去除了 piston/tilt/defocus 后的结果
% 2) 此处直接构造“相对理想相位残差”：
%       Phase_Target = angle(U_final .* conj(exp(1i*phi_init)))
% 3) 对残差进行二维解包裹，打破 2pi 的物理截断限制，还原真实偏离趋势
% 4) Step 8 再用 Ideal + Residual 的方式恢复到绝对相位域绘图

fprintf('\n>>> Step 7: 启动严格圆对称优化（残差域解包裹版）...\n');

% -------------------------------------------------------------------------
% 1. 构造残差相位目标：Δphi(x,y) = measured - ideal
% -------------------------------------------------------------------------
U_ideal_phase = exp(1i * phi_init);
Phase_Target_Wrapped = angle(U_final .* conj(U_ideal_phase));   % wrapped residual in [-pi, pi]

% [核心修正 1]：对残差进行二维解包裹，打破 [-pi, pi] 的硬性吸附与跳变限制
Phase_Target = unwrap(unwrap(Phase_Target_Wrapped, [], 1), [], 2); 

% 有效掩膜：口径内 + 振幅足够高
Amp = abs(U_final);
AmpMask  = Amp > 0.20 * max(Amp(:));
Fit_Mask = (Rho_sim <= 0.98 * radius_metalens) & AmpMask;

% 提取数据点（单位：mm）
X_dat = double(X_sim(Fit_Mask)) * 1e3;
Y_dat = double(Y_sim(Fit_Mask)) * 1e3;
Z_dat = Phase_Target(Fit_Mask);

% -------------------------------------------------------------------------
% 2. 自动搜索最佳光轴中心
%    这里固定拟合半径为真实物理半径 radius_metalens
% -------------------------------------------------------------------------
R_fixed_mm = radius_metalens * 1e3;
start_params = [0, 0];

CostFunc = @(p) solve_symmetric_rmse_stable( ...
    p, X_dat, Y_dat, Z_dat, length(A_coeffs_ideal), R_fixed_mm);

options = optimset('Display', 'iter', 'TolX', 1e-4, 'MaxFunEvals', 200);

fprintf('  正在自动搜索最佳光轴中心...\n');
[best_shift, min_rmse] = fminsearch(CostFunc, start_params, options);

fprintf('  [优化完成] 最佳中心偏移: dx = %.6f mm, dy = %.6f mm\n', best_shift(1), best_shift(2));
fprintf('  [优化完成] 最低残差 RMSE: %.6f rad\n', min_rmse);

% -------------------------------------------------------------------------
% 3. 获取最终拟合系数
%    Coeffs_Stable(1) = piston residual (将自动吃掉 unwrap 带来的常数偏移)
%    Coeffs_Stable(2:end) = 残差多项式系数（以 R_fixed_mm 归一化）
% -------------------------------------------------------------------------
[~, Coeffs_Stable, ~, R_max_used] = solve_symmetric_rmse_stable( ...
    best_shift, X_dat, Y_dat, Z_dat, length(A_coeffs_ideal), R_fixed_mm);

Piston_Final = Coeffs_Stable(1);
B_Final      = Coeffs_Stable(2:end);

% 换算回 Zemax 标准归一化半径 radius_norm_mm
A_Final_Zemax = zeros(size(B_Final));
ratio = R_max_used / radius_norm_mm;
for k = 1:length(B_Final)
    A_Final_Zemax(k) = B_Final(k) / (ratio^(2*k));
end

% -------------------------------------------------------------------------
% 4. Step 7 结果输出
% -------------------------------------------------------------------------
fprintf('\n==================================================================\n');
fprintf('  最终残差拟合系数结果 (严格圆对称, R_norm = %.3f mm)\n', radius_norm_mm);
fprintf('==================================================================\n');
fprintf('  %-6s | %-16s | %-16s | %-16s\n', 'Term', 'Design (Ideal)', 'Residual Fit', 'Abs(Final=Ideal+Res)');
fprintf('------------------------------------------------------------------\n');

labels = {'r^2', 'r^4', 'r^6', 'r^8', 'r^10', 'r^12', 'r^14', 'r^16', 'r^18', 'r^20'};
for k = 1:length(A_Final_Zemax)
    if k <= length(labels)
        lb = labels{k};
    else
        lb = sprintf('r^%d', 2*k);
    end

    val_ideal = A_coeffs_ideal(k);
    val_res   = A_Final_Zemax(k);
    val_final = val_ideal + val_res;

    fprintf('  %-6s | %16.6e | %16.6e | %16.6e\n', lb, val_ideal, val_res, val_final);
end
fprintf('------------------------------------------------------------------\n');
fprintf('  Residual piston offset: %.6f rad\n', Piston_Final);
fprintf('==================================================================\n');


%% 8. 全口径三合一可视化（绝对相位域自洽版）
fprintf('\n>>> Step 8: 生成全口径对比图（绝对相位域自洽版）...\n');

% ================= 参数控制 =================
Num_Bins           = 500;
USE_OUTLIER_REJECT = false;   % Zemax 仿真验证通常建议 false
Outlier_Thresh     = 2.5;
USE_SMOOTH         = false;   % Zemax 仿真验证通常建议 false
Scatter_Step       = 10;      % 散点稀疏显示
% ===========================================

% -------------------------------------------------------------------------
% 1. 生成理想绝对相位曲线 Phi_Ideal(r)
% -------------------------------------------------------------------------
r_half = linspace(0, radius_metalens * 1e3, 1000);
r_norm_half = r_half ./ radius_norm_mm;

Phi_Ideal_Half = zeros(size(r_half));
for k = 1:length(A_coeffs_ideal)
    Phi_Ideal_Half = Phi_Ideal_Half + A_coeffs_ideal(k) .* (r_norm_half .^ (2*k));
end

% -------------------------------------------------------------------------
% 2. 生成残差拟合曲线 DeltaPhi_Fit(r)
%    以及绝对拟合曲线 Phi_FitAbs(r) = Phi_Ideal + DeltaPhi_Fit
% -------------------------------------------------------------------------
Phi_ResFit_Half = zeros(size(r_half));
for k = 1:length(A_Final_Zemax)
    Phi_ResFit_Half = Phi_ResFit_Half + A_Final_Zemax(k) .* (r_norm_half .^ (2*k));
end

Phi_FitAbs_Half = Phi_Ideal_Half + Phi_ResFit_Half;

% -------------------------------------------------------------------------
% 3. 散点数据：这里 Z_dat 就是解包裹后的残差 Δphi(x,y)
% -------------------------------------------------------------------------
R_scatter_raw = sqrt((X_dat - best_shift(1)).^2 + (Y_dat - best_shift(2)).^2); % mm
Z_scatter_raw = Z_dat - Piston_Final;   % 去掉残差常数项，让残差以 0 为基准波动

if USE_OUTLIER_REJECT
    Z_ref_vals = interp1(r_half, Phi_ResFit_Half, R_scatter_raw, 'linear', 'extrap');
    residuals  = abs(Z_scatter_raw - Z_ref_vals);
    valid_mask_scatter = residuals < Outlier_Thresh;

    R_scatter = R_scatter_raw(valid_mask_scatter);
    Z_scatter = Z_scatter_raw(valid_mask_scatter);
else
    R_scatter = R_scatter_raw;
    Z_scatter = Z_scatter_raw;
end

% -------------------------------------------------------------------------
% 4. 径向平均，得到残差均值 DeltaPhi_Mean(r)
% -------------------------------------------------------------------------
edges     = linspace(0, max(r_half), Num_Bins + 1);
z_mean_res = zeros(Num_Bins, 1);
bin_count  = zeros(Num_Bins, 1);

for j = 1:Num_Bins
    mask_bin = (R_scatter >= edges(j)) & (R_scatter < edges(j+1));
    bin_count(j) = sum(mask_bin);

    if bin_count(j) > 5
        z_mean_res(j) = mean(Z_scatter(mask_bin));
    else
        if j > 1
            z_mean_res(j) = z_mean_res(j-1);
        else
            z_mean_res(j) = 0;
        end
    end
end

if USE_SMOOTH
    z_mean_res = smoothdata(z_mean_res, 'gaussian', 5);
end

r_centers = (edges(1:end-1) + edges(2:end)) / 2;   % mm

% -------------------------------------------------------------------------
% 5. 由 Ideal + ResidualMean 构造“实测平均绝对相位曲线”
% -------------------------------------------------------------------------
Phi_Ideal_On_Centers = interp1(r_half, Phi_Ideal_Half, r_centers, 'linear', 'extrap');
Phi_MeasAbs_Centers  = Phi_Ideal_On_Centers(:) + z_mean_res(:);

% -------------------------------------------------------------------------
% 6. 统计量计算
% -------------------------------------------------------------------------
% [核心修正 2] 绝对相位域 R^2：比较 实测绝对相位 与 拟合绝对相位
% 将绝对拟合曲线插值到相同的 r_centers 上
Phi_Fit_On_Centers = interp1(r_half, Phi_FitAbs_Half, r_centers, 'linear', 'extrap')';

SSR = sum((Phi_MeasAbs_Centers(:) - Phi_Fit_On_Centers(:)).^2);
SST = sum((Phi_MeasAbs_Centers(:) - mean(Phi_MeasAbs_Centers(:))).^2);

if SST > 0
    R_Square = 1 - SSR / SST;
else
    R_Square = 1;
end

% (2) 1D 残差 RMSE（沿径向平均）
valid_w = bin_count > 0;
if any(valid_w)
    RMSE_phi_1D = sqrt(sum(bin_count(valid_w) .* (z_mean_res(valid_w).^2)) / sum(bin_count(valid_w)));
else
    RMSE_phi_1D = NaN;
end

% (3) 2D 像素级 RMSE（与 Figure 1 口径定义保持一致）
IdealPhase2D = exp(1i * phi_init);
phase_err_2d = angle(U_final .* conj(IdealPhase2D));
Mask2D = (Rho_sim <= 0.98 * radius_metalens);

piston2d = angle(mean(exp(1i * phase_err_2d(Mask2D))));
phase_err_2d0 = angle(exp(1i * (phase_err_2d - piston2d)));
RMSE_phi_2D = sqrt(mean(phase_err_2d0(Mask2D).^2));

% -------------------------------------------------------------------------
% 7. 镜像展开
% -------------------------------------------------------------------------
r_full       = [-fliplr(r_half), r_half];
phi_id_full  = [fliplr(Phi_Ideal_Half), Phi_Ideal_Half];
phi_fit_full = [fliplr(Phi_FitAbs_Half), Phi_FitAbs_Half];

r_mean_full        = [-fliplr(r_centers), r_centers];
phi_meas_mean_full = [flipud(Phi_MeasAbs_Centers(:)); Phi_MeasAbs_Centers(:)]';

r_scat_full = [-R_scatter; R_scatter];
z_scat_full = [ Z_scatter; Z_scatter];

% 散点转回绝对相位域
r_scat_full_row = r_scat_full(:)';
z_scat_full_row = z_scat_full(:)';
phi_id_at_scat  = interp1(r_half, Phi_Ideal_Half, abs(r_scat_full_row), 'linear', 'extrap');
phi_scat_abs    = phi_id_at_scat + z_scat_full_row;

% -------------------------------------------------------------------------
% 8. 绘图
% -------------------------------------------------------------------------
figure('Name', 'Final Polished Profile', 'Color', 'w', 'Position', [100, 100, 1200, 700]);

% [层1] 原始散点云（绝对相位域）
plot(r_scat_full_row(1:Scatter_Step:end), phi_scat_abs(1:Scatter_Step:end), '.', ...
    'Color', [0.90 0.92 0.95], 'MarkerSize', 1, 'DisplayName', 'Raw Data Cloud');
hold on;

% [层2] 实测平均曲线（绝对相位域）
plot(r_mean_full, phi_meas_mean_full, '-', ...
    'Color', [0.00 0.45 0.85], 'LineWidth', 3.0, 'DisplayName', 'Measured Mean Profile');

% [层3] 理想设计
plot(r_full, phi_id_full, '-', ...
    'Color', [0.20 0.20 0.25], 'LineWidth', 2.0, 'DisplayName', 'Ideal Design');

% [层4] 拟合模型（绝对相位域）
plot(r_full, phi_fit_full, '--', ...
    'Color', [1.00 0.50 0.00], 'LineWidth', 2.5, 'DisplayName', 'Fitted Model');

% 样式
y_limits = [min(phi_id_full), max(phi_id_full)];
y_range  = y_limits(2) - y_limits(1);
ylim([y_limits(1) - 0.15*y_range, y_limits(2) + 0.15*y_range]);
xlim([-radius_metalens*1e3, radius_metalens*1e3]);

set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickDir', 'out', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.4);
grid on; box on;

xlabel('Lateral Position $r$ (mm)', 'Interpreter', 'latex', ...
    'FontSize', 18, 'FontWeight', 'bold');
ylabel('Phase $\phi$ (rad)', 'Interpreter', 'latex', ...
    'FontSize', 18, 'FontWeight', 'bold');
title('Full Aperture Phase Profile Comparison', 'Interpreter', 'latex', ...
    'FontSize', 20);

% 标注框
str_stats = {
    ['\bf Fitting Quality:']
    [' $R^2 = ' num2str(R_Square, '%.4f') '$']
    ['']
    ['\bf Deviation from Ideal (Residual Phase):']
    [' $RMSE_{\mathrm{1D}} = ' num2str(RMSE_phi_1D, '%.4f') ' \ \mathrm{rad}$']
    [' $RMSE_{\mathrm{2D}} = ' num2str(RMSE_phi_2D, '%.4f') ' \ \mathrm{rad}$']
};
annotation('textbox', [0.15 0.15 0.28 0.18], 'String', str_stats, ...
    'Interpreter', 'latex', 'FontSize', 14, ...
    'BackgroundColor', 'w', 'FaceAlpha', 0.8, ...
    'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1, ...
    'FitBoxToText', 'on');

hLeg = legend('Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 13);
set(hLeg, 'Interpreter', 'latex', 'Box', 'off');

fprintf('Step 8 完成：\n');
fprintf('  - 横坐标范围 = [%.3f, %.3f] mm\n', -radius_metalens*1e3, radius_metalens*1e3);
fprintf('  - 1D Residual RMSE = %.6f rad\n', RMSE_phi_1D);
fprintf('  - 2D Residual RMSE = %.6f rad\n', RMSE_phi_2D);
fprintf('  - R^2 = %.6f\n', R_Square);

% %% 9. 数据保存 (Data Saving) - 修正维度匹配版
% % 功能：自动以波长命名，将曲线数据保存为 .csv，将系数保存为 .mat
% % 修复：通过插值解决实测数据与拟合数据长度不一致导致的 table 报错
% 
% fprintf('\n>>> Step 9: 保存拟合结果到本地...\n');
% 
% % --- 1. 设置保存路径 ---
% save_dir = fullfile(dataFolder, 'Fit_Results'); 
% if ~exist(save_dir, 'dir')
%     mkdir(save_dir); 
% end
% 
% % --- 2. 构造文件名 ---
% file_name_base = sprintf('%gum_phi_fit', lambda_um);
% 
% % --- 3. 数据维度对齐 (关键修复) ---
% % r_full 是高分辨率坐标 (2000点)，z_hybrid_full 是分箱坐标 (1000点)
% % 我们需要将实测数据插值到 r_full 的坐标上，以便存入同一个表格
% z_meas_interp = interp1(r_mean_full, z_hybrid_full, r_full, 'linear', 'extrap');
% 
% % --- 4. 保存曲线数据 (.csv) ---
% % 使用 (:) 强制转换为列向量，确保方向一致
% Table_Data = table(...
%     r_full(:), ...
%     phi_id_full(:), ...
%     phi_fit_full(:), ...
%     z_meas_interp(:), ...
%     'VariableNames', {'Radius_mm', 'Phase_Ideal_rad', 'Phase_Fitted_rad', 'Phase_Measured_rad'});
% 
% csv_path = fullfile(save_dir, [file_name_base, '.csv']);
% writetable(Table_Data, csv_path);
% fprintf('  [曲线数据] 已保存至: %s\n', csv_path);
% 
% % --- 5. 保存系数与指标 (.mat) ---
% mat_path = fullfile(save_dir, [file_name_base, '_coeffs.mat']);
% save(mat_path, 'A_Final_Zemax', 'Piston_Final', 'R_Square', 'RMSE_Val', 'lambda_um');
% fprintf('  [系数数据] 已保存至: %s\n', mat_path);
% 
% fprintf('  保存完成！\n');

%% --- 辅助函数 ---
function [P_int, t] = intersect_sphere(x0, y0, z0, dx, dy, dz, C, R)
    ox = x0 - C(1); oy = y0 - C(2); oz = z0 - C(3);
    a = dx.^2 + dy.^2 + dz.^2; b = 2*(ox.*dx + oy.*dy + oz.*dz); c = ox.^2 + oy.^2 + oz.^2 - R^2;
    delta = b.^2 - 4*a.*c; t = (-b - sqrt(delta)) ./ (2*a);
    P_int = [x0 + t.*dx, y0 + t.*dy, z0 + t.*dz];
end
function [nx, ny, nz] = refract_sphere(P, dx, dy, dz, C, R, n1, n2, sign_n)
    Nx = (P(:,1)-C(1))/R*sign_n; Ny = (P(:,2)-C(2))/R*sign_n; Nz = (P(:,3)-C(3))/R*sign_n;
    mu = n1/n2; dot_v_n = dx.*Nx + dy.*Ny + dz.*Nz;
    term2 = sqrt(1 - mu.^2 * (1 - dot_v_n.^2)) - mu * dot_v_n;
    nx = mu*dx + term2.*Nx; ny = mu*dy + term2.*Ny; nz = mu*dz + term2.*Nz;
    nv = sqrt(nx.^2 + ny.^2 + nz.^2); nx=nx./nv; ny=ny./nv; nz=nz./nv;
end
function [P_int, t] = intersect_plane(x0, y0, z0, dx, dy, dz, z_plane)
    t = (z_plane - z0) ./ dz; P_int = [x0 + t.*dx, y0 + t.*dy, z0 + t.*dz];
end
function [nx, ny, nz] = refract_plane(dx, dy, dz, n1, n2)
    mu = n1/n2; nx = mu*dx; ny = mu*dy; nz = sqrt(1 - (nx.^2 + ny.^2));
end
function U_out = ASM_Propagate(U_in, z, lambda, dx)
    [Ny, Nx] = size(U_in); df = 1 / (Nx * dx); f = ((1:Nx) - Nx/2 - 1) * df;
    [FX, FY] = meshgrid(f, f); term = 1 - lambda^2 * (FX.^2 + FY.^2); term(term < 0) = 0; 
    H = exp(1i * 2 * pi * z / lambda * sqrt(term));
    U_f = fftshift(fft2(ifftshift(U_in))); U_out = fftshift(ifft2(ifftshift(U_f .* H)));
end

function [rmse, Coeffs_Stable, Z_fit, R_max_used] = solve_symmetric_rmse_stable( ...
    p, X_raw, Y_raw, Z_raw, n_terms, R_fixed_mm)
% 严格圆对称残差拟合（与新版 Step 7 / Step 8 配套）
%
% 输入：
%   p          = [dx, dy]，待优化的中心偏移（单位 mm）
%   X_raw,Y_raw= 原始采样点坐标（单位 mm）
%   Z_raw      = 对应相位数据（单位 rad），这里应为残差相位
%   n_terms    = 拟合项数，例如 length(A_coeffs_ideal)
%   R_fixed_mm = 固定归一化半径（单位 mm），通常取真实物理半径
%
% 输出：
%   rmse         = 当前中心下的拟合残差 RMSE
%   Coeffs_Stable= [piston, a2, a4, a6, ...]'
%   Z_fit        = 每个输入散点对应的拟合值
%   R_max_used   = 实际使用的归一化半径（这里固定为 R_fixed_mm）

    % -----------------------------
    % 1. 中心平移
    % -----------------------------
    dx = p(1);
    dy = p(2);

    Xc = X_raw - dx;
    Yc = Y_raw - dy;
    R  = sqrt(Xc.^2 + Yc.^2);

    % 固定使用物理口径作为归一化半径
    R_max_used = R_fixed_mm;

    % 只使用口径内数据
    valid = (R <= R_max_used) & isfinite(Z_raw);

    % 数据太少时直接惩罚，避免 fminsearch 跑飞
    if nnz(valid) < max(50, n_terms + 5)
        rmse = 1e9;
        Coeffs_Stable = [0; zeros(n_terms,1)];
        Z_fit = zeros(size(Z_raw));
        return;
    end

    Rv = R(valid);
    Zv = Z_raw(valid);

    % -----------------------------
    % 2. 构造圆对称偶次幂多项式
    %    Z(r) = c0 + c1*(rho^2) + c2*(rho^4) + ...
    %    rho = r / R_max_used
    % -----------------------------
    rho = Rv ./ R_max_used;
    rho2 = rho.^2;

    A = zeros(numel(rho), n_terms + 1);
    A(:,1) = 1;   % piston

    for k = 1:n_terms
        A(:,k+1) = rho2.^k;
    end

    % -----------------------------
    % 3. 最小二乘拟合
    % -----------------------------
    coeffs = A \ Zv;

    % 拟合值
    Zv_fit = A * coeffs;

    % 对全体点也生成拟合值，便于后续画图
    Z_fit = zeros(size(Z_raw));
    rho_all = R ./ R_max_used;
    rho2_all = rho_all.^2;

    A_all = zeros(numel(R), n_terms + 1);
    A_all(:,1) = 1;
    for k = 1:n_terms
        A_all(:,k+1) = rho2_all.^k;
    end
    Z_fit(:) = A_all * coeffs;

    % -----------------------------
    % 4. 计算 RMSE
    %    这里用 wrap 后误差更稳健，避免 2pi 跳变影响
    % -----------------------------
    err = angle(exp(1i * (Zv - Zv_fit)));
    rmse = sqrt(mean(err.^2));

    % 输出列向量
    Coeffs_Stable = coeffs(:);
end