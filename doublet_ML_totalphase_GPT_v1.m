clc; clear; close all;

%% ============================================================
%  双面超透镜合成相位（严格复振幅传播版）
%  说明：
%  1) 前表面：t1 = exp(i*phi1)
%  2) 衬底内传播：ASM，传播距离 = t_sub
%  3) 后表面：t2 = exp(i*phi2)
%  4) 输出：后表面紧后方的复振幅与相位
%
%  注意：
%  - 本代码采用标量衍射模型（ASM），已比“相位直接相加”严格得多
%  - 若两面微结构周期接近波长、矢量效应很强，则需RCWA/FDTD，不属于本代码范围
%  - 你必须确认 Binary 2 系数的物理单位：'rad' 或 'waves'
%% ============================================================

%% 1. 基本参数
lambda0_um = 10.6;                 % 真空波长 [um]
lambda0    = lambda0_um * 1e-6;    % [m]

n_sub      = 3.4177615;            % 衬底折射率（例如 Si）
t_sub_mm   = 0.725;                % 衬底厚度 [mm]
t_sub      = t_sub_mm * 1e-3;      % [m]

D_clear_mm = 44.64;                % 通光口径 [mm]
R_clear    = D_clear_mm * 1e-3 / 2;% 通光半径 [m]

% ---------- 网格参数 ----------
N          = 2048;                 % 采样点数（建议 1024 / 2048）
L_mm       = 50;                   % 计算窗口边长 [mm]，应略大于口径
L          = L_mm * 1e-3;          % [m]
dx         = L / N;                % 采样间隔 [m]

% ---------- 入射场 ----------
input_type = 'plane';              % 'plane' 或 'gaussian'
gauss_w_mm = 50;                   % 若用高斯入射，1/e振幅半径 [mm]

% ---------- 系数单位 ----------
% 'rad'   : 系数直接生成相位(单位 rad)
% 'waves' : 系数生成的是“波数”，最终乘 2*pi 变成 rad
coef_unit = 'rad';

% ---------- 归一化半径 ----------
% Zemax/设计里若相位多项式使用 rho = r / R_norm
Rnorm1_mm = 1;                 % 前表面相位归一化半径 [mm]
Rnorm2_mm = 1;                 % 后表面相位归一化半径 [mm]
Rnorm1    = Rnorm1_mm * 1e-3;      % [m]
Rnorm2    = Rnorm2_mm * 1e-3;      % [m]

%% 2. 双面相位系数（你把这里替换成自己的）
% 多项式形式：
% phi(r) = c1*rho^2 + c2*rho^4 + c3*rho^6 + ...
% 其中 rho = r / Rnorm
%
% 若 coef_unit='rad'   ：phi 直接就是 rad
% 若 coef_unit='waves' ：phi = 2*pi*( c1*rho^2 + c2*rho^4 + ... )

% ===== 示例：请替换为你的前表面 Binary 2 系数 =====
coef_front = [ ...
    -4.758394571634000E+000, ...
    2.194016497211000E-002, ...
    -5.340807755416000E-005, ...
    3.323255468294000E-008, ...
    0 ...
];

% ===== 示例：请替换为你的后表面 Binary 2 系数 =====
coef_back = [ ...
    5.318079102338000E+000, ...
    -2.351626465626000E-002, ...
    5.795669256364000E-005, ...
    -3.663428239962000E-008, ...
    0 ...
];

%% 3. 坐标与孔径
x = ((1:N) - N/2 - 1) * dx;
[X, Y] = meshgrid(x, x);
R = sqrt(X.^2 + Y.^2);

% 软边光阑，减少频域振铃
soft_edge_order = 60;
Aperture_soft = exp(-(R / R_clear).^soft_edge_order);
Aperture_soft(R > 1.02*R_clear) = 0;

Aperture_hard = double(R <= R_clear);

%% 4. 构造入射场
switch lower(input_type)
    case 'plane'
        U_in = Aperture_soft .* exp(1i * 0);

    case 'gaussian'
        w0 = gauss_w_mm * 1e-3;
        A0 = exp(-(R.^2) / (w0^2));
        U_in = A0 .* Aperture_soft;

    otherwise
        error('input_type 只能是 ''plane'' 或 ''gaussian''');
end

% 能量归一化
U_in = U_in / sqrt(sum(abs(U_in(:)).^2));

%% 5. 计算前后两面的相位分布
phi_front = build_phase_even_poly(R, coef_front, Rnorm1, coef_unit);
phi_back  = build_phase_even_poly(R, coef_back,  Rnorm2, coef_unit);

% 若你怀疑背面坐标需要镜像，可尝试以下几种替换之一：
% phi_back = build_phase_even_poly(sqrt((-X).^2 + Y.^2), coef_back, Rnorm2, coef_unit);
% 对于纯径向对称相位，镜像与否没有区别

t_front = Aperture_soft .* exp(1i * phi_front);
t_back  = Aperture_soft .* exp(1i * phi_back);

%% 6. 严格复振幅传播：前面 -> 衬底传播 -> 后面
% 前表面后方的场
U_after_front = U_in .* t_front;

% 在衬底内传播
lambda_sub = lambda0 / n_sub;
U_before_back = ASM_propagate(U_after_front, t_sub, lambda_sub, dx);

% 后表面调制
U_after_back = U_before_back .* t_back;

%% 7. 计算“简单相加相位”的对照结果
phi_sum_naive = angle(exp(1i * (phi_front + phi_back)));
U_naive = abs(U_in) .* Aperture_soft .* exp(1i * (phi_front + phi_back));

%% 8. 提取严格传播后的合成相位
phi_strict_wrap = angle(U_after_back);

% 仅用于显示的解包裹（二维显示通常还是看 wrap 相位更稳）
phi_strict_unwrap_x = unwrap(phi_strict_wrap(round(N/2)+1, :));
phi_naive_unwrap_x  = unwrap(phi_sum_naive(round(N/2)+1, :));

% 中心归零，便于比较
mid = round(N/2)+1;
phi_strict_unwrap_x = phi_strict_unwrap_x - phi_strict_unwrap_x(mid);
phi_naive_unwrap_x  = phi_naive_unwrap_x  - phi_naive_unwrap_x(mid);

%% 9. 计算严格传播相位与简单相加相位的差异
phase_diff = angle(U_after_back .* conj(exp(1i*(phi_front + phi_back))));

mask_valid = (R <= 0.98*R_clear);
rmse_diff = sqrt(mean(phase_diff(mask_valid).^2));

fprintf('============================================================\n');
fprintf('严格传播相位 vs 简单相加相位 的差异 RMSE = %.6f rad\n', rmse_diff);
fprintf('若该值很小，说明衬底内传播对合成相位影响不大。\n');
fprintf('============================================================\n');

%% 10. 画图
x_mm = x * 1e3;

figure('Color','w','Position',[100,80,1600,900]);

subplot(2,3,1);
imagesc(x_mm, x_mm, phi_front);
axis image; colorbar; colormap jet;
title('前表面相位 \phi_{front} (rad)');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(2,3,2);
imagesc(x_mm, x_mm, phi_back);
axis image; colorbar; colormap jet;
title('后表面相位 \phi_{back} (rad)');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(2,3,3);
imagesc(x_mm, x_mm, phi_strict_wrap);
axis image; colorbar; colormap jet;
title('严格传播后的合成相位 angle(U_{after back})');
xlabel('x (mm)'); ylabel('y (mm)');

subplot(2,3,4);
plot(x_mm, phi_naive_unwrap_x, 'k-', 'LineWidth', 1.8); hold on;
plot(x_mm, phi_strict_unwrap_x, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('x (mm)');
ylabel('Phase (rad)');
title('中心剖面：简单相加 vs 严格传播');
legend('简单相加', '严格传播', 'Location', 'best');

subplot(2,3,5);
imagesc(x_mm, x_mm, phase_diff);
axis image; colorbar; colormap jet;
title(sprintf('相位差：strict - naive，RMSE = %.4g rad', rmse_diff));
xlabel('x (mm)'); ylabel('y (mm)');

subplot(2,3,6);
imagesc(x_mm, x_mm, abs(U_after_back).^2);
axis image; colorbar; colormap hot;
title('后表面紧后方强度 |U|^2');
xlabel('x (mm)'); ylabel('y (mm)');

sgtitle('双面超透镜：严格复振幅传播合成相位计算','FontWeight','bold');

%% 11. 输出结果到工作区
Result = struct();
Result.lambda0_um       = lambda0_um;
Result.n_sub            = n_sub;
Result.t_sub_mm         = t_sub_mm;
Result.D_clear_mm       = D_clear_mm;
Result.x_m              = x;
Result.X_m              = X;
Result.Y_m              = Y;
Result.R_m              = R;
Result.U_in             = U_in;
Result.phi_front        = phi_front;
Result.phi_back         = phi_back;
Result.U_after_front    = U_after_front;
Result.U_before_back    = U_before_back;
Result.U_after_back     = U_after_back;
Result.phi_strict_wrap  = phi_strict_wrap;
Result.phi_sum_naive    = phi_sum_naive;
Result.phase_diff       = phase_diff;
Result.rmse_diff        = rmse_diff;

assignin('base', 'Result_double_metalens', Result);
fprintf('结果已输出到工作区变量：Result_double_metalens\n');

%% ========================= 子函数 ==============================

function phi = build_phase_even_poly(R, coeffs, Rnorm, coef_unit)
% 构造偶次幂旋转对称相位：
% phi = c1*(r/Rnorm)^2 + c2*(r/Rnorm)^4 + ...
%
% coeffs: [c1 c2 c3 ...]
% coef_unit:
%   'rad'   -> 直接得到 rad
%   'waves' -> 最终乘 2*pi

    rho = R / Rnorm;
    rho2 = rho.^2;

    phi_base = zeros(size(R));
    for k = 1:length(coeffs)
        phi_base = phi_base + coeffs(k) * (rho2.^k);
    end

    switch lower(coef_unit)
        case 'rad'
            phi = phi_base;
        case 'waves'
            phi = 2*pi*phi_base;
        otherwise
            error('coef_unit 只能是 ''rad'' 或 ''waves''');
    end
end

function U2 = ASM_propagate(U1, z, lambda, dx)
% 角谱法传播（均匀介质中）
% U1     : 输入复振幅
% z      : 传播距离 [m]
% lambda : 介质内波长 [m]
% dx     : 采样间隔 [m]

    [Ny, Nx] = size(U1);
    if Nx ~= Ny
        error('当前 ASM_propagate 仅支持方阵输入');
    end
    N = Nx;

    k = 2*pi / lambda;

    fx = ((1:N) - N/2 - 1) / (N*dx);
    [FX, FY] = meshgrid(fx, fx);

    KX = 2*pi*FX;
    KY = 2*pi*FY;

    KZ_sq = k^2 - KX.^2 - KY.^2;

    % 处理传播波与倏逝波
    KZ = zeros(size(KZ_sq));
    idx_prop = (KZ_sq >= 0);
    idx_evan = (KZ_sq < 0);

    KZ(idx_prop) = sqrt(KZ_sq(idx_prop));
    KZ(idx_evan) = 1i * sqrt(-KZ_sq(idx_evan));

    H = exp(1i * KZ * z);

    U1_f = fftshift(fft2(ifftshift(U1)));
    U2_f = U1_f .* H;
    U2   = fftshift(ifft2(ifftshift(U2_f)));
end