function result = run_mdpr_validation(system_params, sampling_params)
% run_mdpr_validation.m
% 参数化的多距离相位反演验证函数
% 对给定的光学系统和离焦采样方案运行 MDPR，返回恢复质量指标
%
% 输入:
%   system_params: 结构体，包含光学系统参数
%     .lambda_um, .n_lens, .n_substrate, .R1_mm, .Tc_lens_mm,
%     .d_air_gap_mm, .d_substrate_mm, .BFL_mm,
%     .D_stop1_mm, .D_metalens_mm,
%     .Grid_Window_mm, .N_sim, .N_det, .det_pitch_um,
%     .A_coeffs_ideal, .radius_norm_mm,
%     .dataFolder  (指向包含 PSF xlsx 文件的文件夹)
%
%   sampling_params: 结构体，包含采样方案参数
%     .delta_z_mm: 离焦步长 [mm]（用于子采样，必须 <= 原始步长）
%     .z_range_mm: [z_min, z_max] 扫描范围 [mm]
%     .max_iter: 最大迭代次数
%     .alpha_start, .alpha_end, .beta_amp: 算法参数
%
% 输出:
%   result: 结构体，包含恢复质量指标
%     .RMSE_2D, .RMSE_1D, .R_square, .A_final, .coeff_error, ...

%% 1. 解包参数
% --- 物理参数 ---
lambda    = system_params.lambda_um * 1e-6;
n_lens    = system_params.n_lens;
n_sub     = system_params.n_substrate;
R1_mm     = system_params.R1_mm;
Tc_mm     = system_params.Tc_lens_mm;
d_air_mm  = system_params.d_air_gap_mm;
d_sub_mm  = system_params.d_substrate_mm;
BFL_mm    = system_params.BFL_mm;
BFL       = BFL_mm * 1e-3;

% --- 孔径 ---
D_stop1_mm = system_params.D_stop1_mm;
D_ml_mm    = system_params.D_metalens_mm;
r_stop1    = (D_stop1_mm * 1e-3) / 2;
r_ml       = (D_ml_mm * 1e-3) / 2;

% --- 仿真网格 ---
N_sim    = system_params.N_sim;
Grid_Win = system_params.Grid_Window_mm * 1e-3;
dx_sim   = Grid_Win / N_sim;

% --- 探测器 ---
N_det      = system_params.N_det;
det_pitch  = system_params.det_pitch_um * 1e-6;
L_det      = N_det * det_pitch;

% --- Zemax 理想系数 ---
A_ideal  = system_params.A_coeffs_ideal;
r_norm   = system_params.radius_norm_mm;

% --- 算法参数 ---
max_iter    = sampling_params.max_iter;
alpha_start = sampling_params.alpha_start;
alpha_end   = sampling_params.alpha_end;
beta_amp    = sampling_params.beta_amp;

fprintf('\n=============================================================\n');
fprintf('  多距离相位反演验证\n');
fprintf('  lambda = %.1f um, D_ml = %.2f mm, BFL = %.1f mm\n', ...
    system_params.lambda_um, D_ml_mm, BFL_mm);
fprintf('  Delta_z = %.3f mm, z_range = [%.2f, %.2f] mm\n', ...
    sampling_params.delta_z_mm, sampling_params.z_range_mm);
fprintf('=============================================================\n');

%% 2. 正向光线追迹 (构建物理入射场)
fprintf('>>> 光线追迹...\n');

r_launch = r_stop1 * 1.05;
N_rays_linear = 5000;
ln = linspace(-r_launch, r_launch, N_rays_linear);
[Ray_X, Ray_Y] = meshgrid(ln, ln);
valid_S1 = (Ray_X.^2 + Ray_Y.^2) <= r_stop1^2;
rx = Ray_X(valid_S1); ry = Ray_Y(valid_S1); rz = zeros(size(rx));
dir_x = zeros(size(rx)); dir_y = zeros(size(rx)); dir_z = ones(size(rx));
OPL = zeros(size(rx));

% Surface 1: 球面
R1 = R1_mm * 1e-3; Center1 = [0, 0, R1];
[P1, t1] = intersect_sphere(rx, ry, rz, dir_x, dir_y, dir_z, Center1, R1);
OPL = OPL + t1 * 1.0;
[dir_x, dir_y, dir_z] = refract_sphere(P1, dir_x, dir_y, dir_z, Center1, R1, 1.0, n_lens, -1);

% Surface 2: 平面
Tc = Tc_mm * 1e-3;
[P2, t2] = intersect_plane(P1(:,1), P1(:,2), P1(:,3), dir_x, dir_y, dir_z, Tc);
OPL = OPL + t2 * n_lens;
[dir_x, dir_y, dir_z] = refract_plane(dir_x, dir_y, dir_z, n_lens, 1.0);

% Surface 3: 基底前表面
z3 = Tc + d_air_mm * 1e-3;
[P3, t3] = intersect_plane(P2(:,1), P2(:,2), P2(:,3), dir_x, dir_y, dir_z, z3);
OPL = OPL + t3 * 1.0;
[dir_x, dir_y, dir_z] = refract_plane(dir_x, dir_y, dir_z, 1.0, n_sub);

% Surface 4: 超透镜面
z4 = z3 + d_sub_mm * 1e-3;
[P4, t4] = intersect_plane(P3(:,1), P3(:,2), P3(:,3), dir_x, dir_y, dir_z, z4);
OPL = OPL + t4 * n_sub;

% 按超透镜口径筛选
R4 = sqrt(P4(:,1).^2 + P4(:,2).^2);
valid_ml = (R4 <= r_ml);
P4_valid  = P4(valid_ml, :);
OPL_valid = OPL(valid_ml);

fprintf(' - 超透镜口径内光线数: %d\n', nnz(valid_ml));

% 在第四面重建波前
x_sim = ((1:N_sim) - N_sim/2 - 1) * dx_sim;
[X_sim, Y_sim] = meshgrid(x_sim, x_sim);
Rho_sim = sqrt(X_sim.^2 + Y_sim.^2);

OPL_rel = OPL_valid - min(OPL_valid);
F_interp = scatteredInterpolant(P4_valid(:,1), P4_valid(:,2), OPL_rel, 'natural', 'linear');
OPD_grid = F_interp(X_sim, Y_sim);
OPD_grid(~isfinite(OPD_grid)) = 0;

Phi_Inc = (2 * pi / lambda) * OPD_grid;

Soft_Mask = exp(-(Rho_sim / r_ml).^60);
Soft_Mask(Rho_sim > r_ml * 1.01) = 0;

U_inc = Soft_Mask .* exp(1i * Phi_Inc);
U_inc = U_inc / sqrt(sum(abs(U_inc(:)).^2));

%% 3. 按采样方案选择离焦面
fprintf('>>> 按照采样方案选择离焦面...\n');

dataFolder = system_params.dataFolder;
files = dir(fullfile(dataFolder, '*.xlsx'));
if isempty(files), error('未找到 xlsx 文件: %s', dataFolder); end

% 提取所有可用的离焦位置
all_dz = [];
for i = 1:numel(files)
    tok = regexp(files(i).name, 'PSF_([-+0-9.eE]+)mm', 'tokens', 'once');
    if ~isempty(tok), all_dz(end+1) = str2double(tok{1}); end
end
all_dz = sort(all_dz);

% 选择符合采样方案的离焦面
z_min = sampling_params.z_range_mm(1);
z_max = sampling_params.z_range_mm(2);
dz_target = sampling_params.delta_z_mm;

% 生成目标离焦位置
target_dz = (z_min : dz_target : z_max)';
% 对每个目标位置，找最近的可用离焦面
selected_dz = zeros(size(target_dz));
for i = 1:length(target_dz)
    [~, idx] = min(abs(all_dz - target_dz(i)));
    selected_dz(i) = all_dz(idx);
end
selected_dz = unique(selected_dz);  % 去重

fprintf(' - 目标离焦范围: [%.2f, %.2f] mm, 步长: %.3f mm\n', z_min, z_max, dz_target);
fprintf(' - 实际选择 %d 个离焦面\n', length(selected_dz));
fprintf(' - 离焦位置: ');
for i = 1:min(10, length(selected_dz))
    fprintf('%.3f ', selected_dz(i));
end
if length(selected_dz) > 10, fprintf('...'); end
fprintf('\n');

%% 4. 读取所选离焦面的数据
% 物理尺寸映射
N_roi = round(L_det / dx_sim);
if mod(N_roi, 2) ~= 0, N_roi = N_roi + 1; end
sim_center = floor(N_sim/2) + 1;

stack_amp = zeros(N_sim, N_sim, length(selected_dz), 'single');
prop_dists = zeros(length(selected_dz), 1);
valid_count = 0;

for i = 1:numel(files)
    fname = files(i).name;
    tok = regexp(fname, 'PSF_([-+0-9.eE]+)mm', 'tokens', 'once');
    if isempty(tok), continue; end
    dz_val = str2double(tok{1});

    % 检查是否在选中的离焦面列表中
    if ~ismember(dz_val, selected_dz), continue; end

    try
        I_raw = readmatrix(fullfile(dataFolder, fname));
    catch
        continue;
    end

    % 尺寸调整
    [rows, cols] = size(I_raw);
    if rows ~= N_det || cols ~= N_det
        if rows > N_det && cols > N_det
            cr = floor(rows/2)+1; cc = floor(cols/2)+1;
            I_raw = I_raw(cr-256:cr+255, cc-256:cc+255);
        else
            I_raw = imresize(I_raw, [N_det, N_det]);
        end
    end

    % 数据清洗
    I_raw = double(I_raw);
    I_raw(isnan(I_raw)) = 0; I_raw(I_raw < 0) = 0;
    bg_thresh = 0.001 * max(I_raw(:));
    I_raw(I_raw < bg_thresh) = 0;

    % 质心对齐
    sum_I = sum(I_raw(:));
    if sum_I > 0
        [rs, cs] = ndgrid(1:size(I_raw,1), 1:size(I_raw,2));
        cy = sum(rs(:).*I_raw(:)) / sum_I;
        cx = sum(cs(:).*I_raw(:)) / sum_I;
        sy = (size(I_raw,1)/2 + 0.5) - cy;
        sx = (size(I_raw,2)/2 + 0.5) - cx;
        I_raw = imtranslate(I_raw, [sx, sy], 'bilinear');
    end

    % 映射到仿真网格
    I_resized = imresize(I_raw, [N_roi, N_roi], 'bilinear');
    I_resized(I_resized < 0) = 0;

    Full_Grid = zeros(N_sim, N_sim, 'single');
    half = floor(N_roi/2);
    idx = (sim_center - half) : (sim_center - half + N_roi - 1);
    Full_Grid(idx, idx) = single(I_resized);

    valid_count = valid_count + 1;
    stack_amp(:,:,valid_count) = sqrt(Full_Grid);
    prop_dists(valid_count) = BFL + dz_val * 1e-3;
end

% 排序
stack_amp = stack_amp(:,:,1:valid_count);
prop_dists = prop_dists(1:valid_count);
[prop_dists, idx] = sort(prop_dists);
stack_amp = stack_amp(:,:,idx);
nPlanes = valid_count;

fprintf(' - 成功加载 %d 个离焦面\n', nPlanes);
if nPlanes < 3
    result.RMSE_2D = NaN;
    result.RMSE_1D = NaN;
    result.R_square = NaN;
    result.coeff_err_low = [NaN, NaN];
    result.coeff_err_high = NaN;
    result.coeff_err_mean_high = NaN;
    result.nPlanes = nPlanes;
    result.selected_dz = selected_dz;
    result.A_final = zeros(size(A_ideal));
    result.A_ideal = A_ideal;
    result.best_shift = [0, 0];
    fprintf(' *** 离焦面数不足，跳过迭代\n');
    return;
end

%% 5. 相位反演迭代
fprintf('>>> 开始迭代 (max_iter=%d)...\n', max_iter);

% 理想相位初始猜测
r_norm_z = double(Rho_sim) ./ (r_norm * 1e-3);
phi_init = zeros(N_sim, N_sim);
r_sq_z = r_norm_z.^2;
for m = 1:length(A_ideal)
    phi_init = phi_init + A_ideal(m) * (r_sq_z.^m);
end

U_total_est = U_inc .* exp(1i * phi_init);

for iter = 1:max_iter
    alpha_w = alpha_start - (alpha_start - alpha_end) * (iter / max_iter);
    U_accum = zeros(N_sim, N_sim);
    E_source = sum(abs(U_inc(:)).^2);

    for k = 1:nPlanes
        dist = prop_dists(k);
        meas_amp_raw = double(stack_amp(:,:,k));

        U_prop = ASM_Propagate(U_total_est, dist, lambda, dx_sim);

        E_calc = sum(abs(U_prop(:)).^2);
        E_meas = sum(meas_amp_raw(:).^2);
        if E_meas > 0
            Scale_Factor = sqrt(E_calc / E_meas);
        else
            Scale_Factor = 0;
        end
        meas_amp_scaled = meas_amp_raw * Scale_Factor;

        Amp_updated = alpha_w * meas_amp_scaled + (1 - alpha_w) * abs(U_prop);
        U_corrected = Amp_updated .* exp(1i * angle(U_prop));
        U_back = ASM_Propagate(U_corrected, -dist, lambda, dx_sim);
        U_accum = U_accum + U_back;
    end

    U_total_temp = U_accum / nPlanes;

    Amp_back = abs(U_total_temp);
    Mask_Binary = double(Rho_sim <= r_ml);
    Amp_Updated = beta_amp * abs(U_inc) + (1 - beta_amp) * Amp_back;
    Amp_Final = Amp_Updated .* Mask_Binary;
    U_total_new = Amp_Final .* exp(1i * angle(U_total_temp));
    U_total_new = U_total_new / sqrt(sum(abs(U_total_new(:)).^2)) * sqrt(E_source);

    err = norm(U_total_new(:) - U_total_est(:), 'fro') / norm(U_total_est(:), 'fro');
    U_total_est = U_total_new;

    if mod(iter, 10) == 0
        fprintf('  Iter %d: Err = %.4f\n', iter, err);
    end
end

%% 6. 结果验证: 提取并评估恢复质量
fprintf('>>> 评估恢复质量...\n');

U_ml_extracted = U_total_est .* conj(U_inc);
H_smooth = fspecial('gaussian', [5 5], 1.0);
U_ml_extracted = imfilter(U_ml_extracted, H_smooth, 'replicate');

U_ideal = Soft_Mask .* exp(1i * phi_init);

% Piston/Tilt/Defocus 校正
[Xg, Yg] = meshgrid(x_sim, x_sim);
U_diff = U_ml_extracted .* conj(U_ideal);
phase_diff_map = angle(U_diff);

Valid_Mask = (Rho_sim <= 0.98 * r_ml);
Norm_Radius = r_ml;
x_val = Xg(Valid_Mask); y_val = Yg(Valid_Mask);
z_val = phase_diff_map(Valid_Mask);
x_n = x_val / Norm_Radius; y_n = y_val / Norm_Radius;
r2_n = x_n.^2 + y_n.^2;

A_fit = [ones(size(x_n)), x_n, y_n, r2_n];
coeffs = A_fit \ z_val;

X_norm = Xg / Norm_Radius; Y_norm = Yg / Norm_Radius;
R2_norm = X_norm.^2 + Y_norm.^2;
Correction_Phase = coeffs(1) + coeffs(2)*X_norm + coeffs(3)*Y_norm + coeffs(4)*R2_norm;
U_final = U_ml_extracted .* exp(-1i * Correction_Phase);

% 2D RMSE (整体)
IdealPhase2D = exp(1i * phi_init);
phase_err_2d = angle(U_final .* conj(IdealPhase2D));
Mask2D = (Rho_sim <= 0.98 * r_ml);
piston2d = angle(mean(exp(1i * phase_err_2d(Mask2D))));
phase_err_2d0 = angle(exp(1i * (phase_err_2d - piston2d)));
RMSE_2D = sqrt(mean(phase_err_2d0(Mask2D).^2));

% 高频系数误差 (高阶项 r^6 及以上)
Phase_Target_Wrapped = angle(U_final .* conj(exp(1i * phi_init)));
Phase_Target = unwrap(unwrap(Phase_Target_Wrapped, [], 1), [], 2);

Amp = abs(U_final);
AmpMask = Amp > 0.20 * max(Amp(:));
Fit_Mask = (Rho_sim <= 0.98 * r_ml) & AmpMask;

X_dat_full = double(X_sim(Fit_Mask)) * 1e3;
Y_dat_full = double(Y_sim(Fit_Mask)) * 1e3;
Z_dat_full = Phase_Target(Fit_Mask);

% 子采样加速fminsearch (每10个点取1个, 仍足够拟合低阶多项式)
subsample_step = 5;
X_dat = X_dat_full(1:subsample_step:end);
Y_dat = Y_dat_full(1:subsample_step:end);
Z_dat = Z_dat_full(1:subsample_step:end);

R_fixed_mm = r_ml * 1e3;

CostFunc = @(p) solve_symmetric_rmse_stable(p, X_dat, Y_dat, Z_dat, length(A_ideal), R_fixed_mm);
options = optimset('Display', 'off', 'TolX', 5e-3, 'MaxFunEvals', 100);
[best_shift, ~] = fminsearch(CostFunc, [0, 0], options);

% 用全数据+最优中心计算最终系数
[~, Coeffs_Stable, ~, R_max_used] = solve_symmetric_rmse_stable(...
    best_shift, X_dat_full, Y_dat_full, Z_dat_full, length(A_ideal), R_fixed_mm);

[~, Coeffs_Stable, ~, R_max_used] = solve_symmetric_rmse_stable(...
    best_shift, X_dat, Y_dat, Z_dat, length(A_ideal), R_fixed_mm);

B_Final = Coeffs_Stable(2:end);
A_Final_Zemax = zeros(size(B_Final));
ratio = R_max_used / r_norm;
for k = 1:length(B_Final)
    A_Final_Zemax(k) = B_Final(k) / (ratio^(2*k));
end

% 分离低阶和高阶系数误差
n_total = length(A_ideal);
n_low = min(2, n_total);  % r^2, r^4 为低阶
coeff_err_low = zeros(1, n_low);
coeff_err_high = zeros(1, n_total - n_low);
for k = 1:n_total
    ideal_val = A_ideal(k);
    fitted_val = ideal_val + A_Final_Zemax(k);
    rel_err = abs(fitted_val - ideal_val) / max(abs(ideal_val), 1e-10);
    if k <= n_low
        coeff_err_low(k) = rel_err;
    else
        coeff_err_high(k - n_low) = rel_err;
    end
end

% 汇总结果
result.RMSE_2D        = RMSE_2D;
result.coeff_err_low  = coeff_err_low;
result.coeff_err_high = coeff_err_high;
result.coeff_err_mean_high = mean(coeff_err_high);
result.nPlanes        = nPlanes;
result.selected_dz    = selected_dz;
result.A_final        = A_Final_Zemax;
result.A_ideal        = A_ideal;
result.best_shift     = best_shift;

fprintf(' - RMSE_2D = %.4f rad\n', RMSE_2D);
fprintf(' - 高阶系数平均相对误差 = %.4f\n', result.coeff_err_mean_high);

%% 辅助函数 (内联, 与主代码保持一致)
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
function [rmse, Coeffs_Stable, Z_fit, R_max_used] = solve_symmetric_rmse_stable(...
    p, X_raw, Y_raw, Z_raw, n_terms, R_fixed_mm)
    dx_p = p(1); dy_p = p(2);
    Xc = X_raw - dx_p; Yc = Y_raw - dy_p;
    R = sqrt(Xc.^2 + Yc.^2);
    R_max_used = R_fixed_mm;
    valid = (R <= R_max_used) & isfinite(Z_raw);
    if nnz(valid) < max(50, n_terms + 5)
        rmse = 1e9; Coeffs_Stable = [0; zeros(n_terms,1)]; Z_fit = zeros(size(Z_raw)); return;
    end
    Rv = R(valid); Zv = Z_raw(valid);
    rho = Rv ./ R_max_used; rho2 = rho.^2;
    A_mat = zeros(numel(rho), n_terms + 1); A_mat(:,1) = 1;
    for k = 1:n_terms, A_mat(:,k+1) = rho2.^k; end
    coeffs = A_mat \ Zv; Zv_fit = A_mat * coeffs;
    Z_fit = zeros(size(Z_raw));
    rho_all = R ./ R_max_used; rho2_all = rho_all.^2;
    A_all = zeros(numel(R), n_terms + 1); A_all(:,1) = 1;
    for k = 1:n_terms, A_all(:,k+1) = rho2_all.^k; end
    Z_fit(:) = A_all * coeffs;
    err = angle(exp(1i * (Zv - Zv_fit)));
    rmse = sqrt(mean(err.^2));
    Coeffs_Stable = coeffs(:);
end
end
