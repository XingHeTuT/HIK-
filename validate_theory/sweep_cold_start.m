% sweep_cold_start.m
% 去热启动 alpha 扫描: 零相位初始猜解, 验证 alpha_opt = pi/2
%
% 与热启动版本的区别:
%   phi_init = zeros(N_sim, N_sim)  而不是从 Zemax 理想系数生成
%   迭代增至 50 次 (零初始需要更多迭代)
%   此实验真正测量"数据里的信息量"而非"初始猜解的质量"

clc; clear all; %#ok<CLALL>

%% ===== 系统A参数 =====
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

% 物理常数
lambda    = sys.lambda_um * 1e-6;
n_lens    = sys.n_lens;
n_sub     = sys.n_substrate;
R1_mm     = sys.R1_mm;
Tc_mm     = sys.Tc_lens_mm;
d_air_mm  = sys.d_air_gap_mm;
d_sub_mm  = sys.d_substrate_mm;
BFL_mm    = sys.BFL_mm;
BFL       = BFL_mm * 1e-3;
D_stop1_mm = sys.D_stop1_mm;
D_ml_mm    = sys.D_metalens_mm;
r_stop1    = (D_stop1_mm * 1e-3) / 2;
r_ml       = (D_ml_mm * 1e-3) / 2;
N_sim    = sys.N_sim;
Grid_Win = sys.Grid_Window_mm * 1e-3;
dx_sim   = Grid_Win / N_sim;
N_det      = sys.N_det;
det_pitch  = sys.det_pitch_um * 1e-6;
L_det      = N_det * det_pitch;
A_ideal  = sys.A_coeffs_ideal;
r_norm   = sys.radius_norm_mm;

z_c_mm = lambda * 1e3 / (pi * ((D_ml_mm/2)/BFL_mm)^2);
a_opt  = pi/2;
fprintf('System A: z_c=%.3fmm, alpha_opt=%.3f (dz_opt=%.3fmm)\n', z_c_mm, a_opt, a_opt*z_c_mm);

%% ===== alpha 配置 =====
export_root = ['D:\LiLinhan_2026\CC-多距离强度相位复原\HIK相位测量\', ...
               'validate_theory\zemax_export'];
alpha_list     = [0.10, 0.20, 0.50, 0.80, 1.20, pi/2, 2.00, 2.50, 3.00];
dz_list_mm     = [0.036, 0.073, 0.181, 0.290, 0.435, 0.570, 0.725, 0.906, 1.088];
output_subdirs = {'alpha_0.10_dz_0.036mm', 'alpha_0.20_dz_0.073mm', ...
                  'alpha_0.50_dz_0.181mm', 'alpha_0.80_dz_0.290mm', ...
                  'alpha_1.20_dz_0.435mm', 'alpha_1.57_dz_0.570mm', ...
                  'alpha_2.00_dz_0.725mm', 'alpha_2.50_dz_0.906mm', ...
                  'alpha_3.00_dz_1.088mm'};
K_expected     = [277, 138, 57, 36, 24, 19, 15, 13, 11];

% 算法参数 (冷启动: 零初始相位, 更多迭代)
max_iter_cold = 50;
alpha_w_start = 1.0;
alpha_w_end   = 0.8;
beta_amp      = 0.8;
z_range_mm = [-5.0, 5.0];

%% ===== 正向光线追迹 (一次, 所有alpha共用) =====
fprintf('\n>>> 光线追迹 (构建入射场)...\n');
r_launch = r_stop1 * 1.05;
N_rays_lin = 5000;
ln = linspace(-r_launch, r_launch, N_rays_lin);
[Ray_X, Ray_Y] = meshgrid(ln, ln);
valid_S1 = (Ray_X.^2 + Ray_Y.^2) <= r_stop1^2;
rx = Ray_X(valid_S1); ry = Ray_Y(valid_S1); rz = zeros(size(rx));
dir_x = zeros(size(rx)); dir_y = zeros(size(rx)); dir_z = ones(size(rx));
OPL = zeros(size(rx));

R1 = R1_mm * 1e-3; Center1 = [0, 0, R1];
[P1, t1] = intersect_sphere(rx, ry, rz, dir_x, dir_y, dir_z, Center1, R1);
OPL = OPL + t1 * 1.0;
[dir_x, dir_y, dir_z] = refract_sphere(P1, dir_x, dir_y, dir_z, Center1, R1, 1.0, n_lens, -1);

Tc = Tc_mm * 1e-3;
[P2, t2] = intersect_plane(P1(:,1), P1(:,2), P1(:,3), dir_x, dir_y, dir_z, Tc);
OPL = OPL + t2 * n_lens;
[dir_x, dir_y, dir_z] = refract_plane(dir_x, dir_y, dir_z, n_lens, 1.0);

z3 = Tc + d_air_mm * 1e-3;
[P3, t3] = intersect_plane(P2(:,1), P2(:,2), P2(:,3), dir_x, dir_y, dir_z, z3);
OPL = OPL + t3 * 1.0;
[dir_x, dir_y, dir_z] = refract_plane(dir_x, dir_y, dir_z, 1.0, n_sub);

z4 = z3 + d_sub_mm * 1e-3;
[P4, t4] = intersect_plane(P3(:,1), P3(:,2), P3(:,3), dir_x, dir_y, dir_z, z4);
OPL = OPL + t4 * n_sub;

R4 = sqrt(P4(:,1).^2 + P4(:,2).^2);
valid_ml = (R4 <= r_ml);
P4_valid = P4(valid_ml, :);
OPL_valid = OPL(valid_ml);
fprintf('  口径内光线: %d\n', nnz(valid_ml));

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

%% ===== 冷启动 MDPR 扫描 =====
fprintf('\n===== COLD START: phi_init = 0 =====\n');
fprintf('max_iter=%d, 9 alpha configurations\n\n', max_iter_cold);

results_cold = cell(length(alpha_list), 1);
rmse_cold = zeros(length(alpha_list), 1);
err_high_cold = zeros(length(alpha_list), 1);
k_actual_cold = zeros(length(alpha_list), 1);
t_total = tic;

for cfg_idx = 1:length(alpha_list)
    dataFolder = fullfile(export_root, output_subdirs{cfg_idx}, 'Zemax_xlsx');

    if ~exist(dataFolder, 'dir')
        fprintf('[%d/9] SKIP: xlsx not found\n', cfg_idx);
        results_cold{cfg_idx} = struct('RMSE_2D', NaN);
        rmse_cold(cfg_idx) = NaN; err_high_cold(cfg_idx) = NaN;
        continue;
    end

    fprintf('\n[%d/9] alpha=%.3f, dz=%.3fmm\n', cfg_idx, alpha_list(cfg_idx), dz_list_mm(cfg_idx));
    t_cfg = tic;

    % ---- 数据加载 ----
    files = dir(fullfile(dataFolder, '*.xlsx'));
    N_roi = round(L_det / dx_sim);
    if mod(N_roi, 2) ~= 0, N_roi = N_roi + 1; end
    sim_center = floor(N_sim/2) + 1;

    stack_amp = zeros(N_sim, N_sim, numel(files), 'single');
    prop_dists = zeros(numel(files), 1);
    valid_count = 0;

    for j = 1:numel(files)
        fname = files(j).name;
        tok = regexp(fname, 'PSF_([-+0-9.eE]+)mm', 'tokens', 'once');
        if isempty(tok), continue; end
        dz_val = str2double(tok{1});

        try
            I_raw = readmatrix(fullfile(dataFolder, fname));
        catch, continue; end

        [rows, cols] = size(I_raw);
        if rows ~= N_det || cols ~= N_det
            if rows > N_det && cols > N_det
                cr = floor(rows/2)+1; cc = floor(cols/2)+1;
                I_raw = I_raw(cr-256:cr+255, cc-256:cc+255);
            else
                I_raw = imresize(I_raw, [N_det, N_det]);
            end
        end

        I_raw = double(I_raw);
        I_raw(isnan(I_raw)) = 0; I_raw(I_raw < 0) = 0;
        bg = 0.001 * max(I_raw(:)); I_raw(I_raw < bg) = 0;

        sum_I = sum(I_raw(:));
        if sum_I > 0
            [rs, cs] = ndgrid(1:size(I_raw,1), 1:size(I_raw,2));
            cy = sum(rs(:).*I_raw(:)) / sum_I;
            cx = sum(cs(:).*I_raw(:)) / sum_I;
            sy = (size(I_raw,1)/2 + 0.5) - cy;
            sx = (size(I_raw,2)/2 + 0.5) - cx;
            I_raw = imtranslate(I_raw, [sx, sy], 'bilinear');
        end

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

    stack_amp = stack_amp(:,:,1:valid_count);
    prop_dists = prop_dists(1:valid_count);
    [prop_dists, idx_sort] = sort(prop_dists);
    stack_amp = stack_amp(:,:,idx_sort);
    nPlanes = valid_count;
    fprintf('  K=%d planes loaded\n', nPlanes);

    % ===== [关键] 冷启动: phi_init = 0 =====
    phi_init = zeros(N_sim, N_sim);
    U_total_est = U_inc .* exp(1i * phi_init);

    % MDPR 迭代
    for iter = 1:max_iter_cold
        alpha_w = alpha_w_start - (alpha_w_start - alpha_w_end) * (iter / max_iter_cold);
        U_accum = zeros(N_sim, N_sim);
        E_source = sum(abs(U_inc(:)).^2);

        for k = 1:nPlanes
            dist = prop_dists(k);
            meas_amp_raw = double(stack_amp(:,:,k));
            U_prop = ASM_Propagate(U_total_est, dist, lambda, dx_sim);

            E_calc = sum(abs(U_prop(:)).^2);
            E_meas = sum(meas_amp_raw(:).^2);
            Scale_Factor = sqrt(E_calc / max(E_meas, 1e-30));
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
        U_total_new = Amp_Updated .* Mask_Binary .* exp(1i * angle(U_total_temp));
        U_total_new = U_total_new / sqrt(sum(abs(U_total_new(:)).^2)) * sqrt(E_source);

        err = norm(U_total_new(:) - U_total_est(:), 'fro') / norm(U_total_est(:), 'fro');
        U_total_est = U_total_new;

        if mod(iter, 25) == 0
            fprintf('    Iter %d: Err=%.4f\n', iter, err);
        end
    end

    % ---- 评估 ----
    U_ml_ext = U_total_est .* conj(U_inc);
    H_smooth = fspecial('gaussian', [5 5], 1.0);
    U_ml_ext = imfilter(U_ml_ext, H_smooth, 'replicate');

    U_ideal = Soft_Mask .* exp(1i * phi_init);  % phi_init=0 here
    % 但实际上理想相位还是从 Zemax 系数算 (用于误差计算)
    phi_ideal_temp = zeros(N_sim, N_sim);
    r_norm_z = double(Rho_sim) ./ (r_norm * 1e-3);
    r_sq_z = r_norm_z.^2;
    for m = 1:length(A_ideal)
        phi_ideal_temp = phi_ideal_temp + A_ideal(m) * (r_sq_z.^m);
    end
    U_ideal_true = Soft_Mask .* exp(1i * phi_ideal_temp);

    % Piston/Tilt/Defocus 校正
    U_diff = U_ml_ext .* conj(U_ideal_true);
    phase_diff_map = angle(U_diff);
    Valid_Mask = (Rho_sim <= 0.98 * r_ml);
    Norm_Radius = r_ml;
    x_val = X_sim(Valid_Mask); y_val = Y_sim(Valid_Mask);
    z_val = phase_diff_map(Valid_Mask);
    x_n = x_val / Norm_Radius; y_n = y_val / Norm_Radius;
    r2_n = x_n.^2 + y_n.^2;
    A_fit = [ones(size(x_n)), x_n, y_n, r2_n];
    coeffs = A_fit \ z_val;
    X_norm = X_sim / Norm_Radius; Y_norm = Y_sim / Norm_Radius;
    R2_norm = X_norm.^2 + Y_norm.^2;
    Correction_Phase = coeffs(1) + coeffs(2)*X_norm + coeffs(3)*Y_norm + coeffs(4)*R2_norm;
    U_final = U_ml_ext .* exp(-1i * Correction_Phase);

    % RMSE
    phase_err_2d = angle(U_final .* conj(U_ideal_true));
    Mask2D = (Rho_sim <= 0.98 * r_ml);
    piston2d = angle(mean(exp(1i * phase_err_2d(Mask2D))));
    phase_err_2d0 = angle(exp(1i * (phase_err_2d - piston2d)));
    RMSE_val = sqrt(mean(phase_err_2d0(Mask2D).^2));

    % 系数拟合
    Phase_Target_Wrapped = angle(U_final .* conj(exp(1i * phi_ideal_temp)));
    Phase_Target = unwrap(unwrap(Phase_Target_Wrapped, [], 1), [], 2);
    Amp = abs(U_final);
    AmpMask = Amp > 0.20 * max(Amp(:));
    Fit_Mask = (Rho_sim <= 0.98 * r_ml) & AmpMask;
    X_dat = double(X_sim(Fit_Mask)) * 1e3;
    Y_dat = double(Y_sim(Fit_Mask)) * 1e3;
    Z_dat = Phase_Target(Fit_Mask);
    R_fixed_mm = r_ml * 1e3;

    X_dat_sub = X_dat(1:5:end); Y_dat_sub = Y_dat(1:5:end); Z_dat_sub = Z_dat(1:5:end);
    CostFunc = @(p) solve_sym(p, X_dat_sub, Y_dat_sub, Z_dat_sub, length(A_ideal), R_fixed_mm);
    options = optimset('Display', 'off', 'TolX', 5e-3, 'MaxFunEvals', 100);
    [best_shift, ~] = fminsearch(CostFunc, [0, 0], options);

    [~, Coeffs_Stable, ~, R_max_used] = solve_sym(best_shift, X_dat, Y_dat, Z_dat, length(A_ideal), R_fixed_mm);
    B_Final = Coeffs_Stable(2:end);
    A_Final_Zemax = zeros(size(B_Final));
    ratio = R_max_used / r_norm;
    for k = 1:length(B_Final)
        A_Final_Zemax(k) = B_Final(k) / (ratio^(2*k));
    end

    % 只对非零理想值项计算相对误差
    err_high_val = 0; n_high = 0;
    err_low_val = 0; n_low = 0;
    for k = 1:length(A_ideal)
        if abs(A_ideal(k)) < 1e-12, continue; end
        rel = abs(A_Final_Zemax(k)) / abs(A_ideal(k));
        if k <= 2
            err_low_val = err_low_val + rel; n_low = n_low + 1;
        else
            err_high_val = err_high_val + rel; n_high = n_high + 1;
        end
    end
    if n_low > 0, err_low_val = err_low_val / n_low; end
    if n_high > 0, err_high_val = err_high_val / n_high; end

    rmse_cold(cfg_idx) = RMSE_val;
    err_high_cold(cfg_idx) = err_high_val;
    k_actual_cold(cfg_idx) = nPlanes;

    results_cold{cfg_idx} = struct('RMSE_2D', RMSE_val, 'err_high', err_high_val, ...
        'err_low', err_low_val, 'nPlanes', nPlanes, 'A_final', A_Final_Zemax, 'A_ideal', A_ideal);

    fprintf('  RMSE=%.4f  Err_high=%.4f  Err_low=%.4f  time=%.0fs\n', ...
        RMSE_val, err_high_val, err_low_val, toc(t_cfg));

    % 中间保存
    save_dir = fileparts(mfilename('fullpath'));
    save(fullfile(save_dir, 'cold_start_intermediate.mat'), ...
        'alpha_list', 'rmse_cold', 'err_high_cold', 'k_actual_cold', ...
        'results_cold', 'cfg_idx');
end

fprintf('\nTotal: %.1f hours\n', toc(t_total)/3600);

%% ===== 结果 =====
fprintf('\n===== COLD START RESULTS =====\n');
fprintf('%-10s | %-10s | %-6s | %-12s | %-12s\n', 'alpha', 'dz[mm]', 'K', 'RMSE_2D', 'Err_high');
fprintf('%s\n', repmat('-', 1, 58));
for i = 1:length(alpha_list)
    fprintf('%-10.3f | %-10.3f | %-6d | %-12.4f | %-12.4f\n', ...
        alpha_list(i), dz_list_mm(i), k_actual_cold(i), rmse_cold(i), err_high_cold(i));
end

valid = ~isnan(rmse_cold);
valid_idx = find(valid);
[~, idx_best] = min(err_high_cold(valid));
fprintf('\nBest (Err_high): alpha=%.3f (dz=%.3fmm), theory: alpha=%.3f\n', ...
    alpha_list(valid_idx(idx_best)), dz_list_mm(valid_idx(idx_best)), a_opt);

% 保存最终结果
save(fullfile(fileparts(mfilename('fullpath')), 'cold_start_final.mat'), ...
    'alpha_list', 'dz_list_mm', 'rmse_cold', 'err_high_cold', 'k_actual_cold', ...
    'results_cold', 'z_c_mm', 'a_opt');
fprintf('\nResults saved.\n');

%% ===== 辅助函数 =====
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
function [rmse, Coeffs, Z_fit, R_max_used] = solve_sym(p, X_raw, Y_raw, Z_raw, n_terms, R_fixed_mm)
    Xc = X_raw - p(1); Yc = Y_raw - p(2); R = sqrt(Xc.^2 + Yc.^2);
    R_max_used = R_fixed_mm;
    valid = (R <= R_max_used) & isfinite(Z_raw);
    if nnz(valid) < max(50, n_terms + 5)
        rmse = 1e9; Coeffs = [0; zeros(n_terms,1)]; Z_fit = zeros(size(Z_raw)); return;
    end
    Rv = R(valid); Zv = Z_raw(valid);
    rho = Rv ./ R_max_used; rho2 = rho.^2;
    A_mat = zeros(numel(rho), n_terms + 1); A_mat(:,1) = 1;
    for k = 1:n_terms, A_mat(:,k+1) = rho2.^k; end
    coeffs = A_mat \ Zv;
    Z_fit = zeros(size(Z_raw));
    rho_all = R ./ R_max_used; rho2_all = rho_all.^2;
    A_all = zeros(numel(R), n_terms + 1); A_all(:,1) = 1;
    for k = 1:n_terms, A_all(:,k+1) = rho2_all.^k; end
    Z_fit(:) = A_all * coeffs;
    err = angle(exp(1i * (Zv - A_mat * coeffs)));
    rmse = sqrt(mean(err.^2));
    Coeffs = coeffs(:);
end
