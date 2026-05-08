function MasterProgramScientific()

clear; clc; close all;

if ~exist('output', 'dir')
    mkdir('output');
end
if ~exist(fullfile('output','frames'), 'dir')
    mkdir(fullfile('output','frames'));
end

%% ==========================================================
% SSI MASTER PROGRAM - GAZETAS / NIST VERSION
% Function-file version for Octave compatibility
%% ==========================================================

%% 1. INPUT DATA
E = 30e9;              % Pa
D = 1.0;               % m
h = 3.0;               % m
T_target = 0.10;       % s
F = 1.0e3;             % N

% Soil properties
G = 80e6;              % Pa
nu = 0.35;
rho_s = 1800;          % kg/m^3
Vs = sqrt(G/rho_s);

% Foundation geometry
B = 2.0;               % m
L = 3.0;               % m
H = 6.0;               % m

% Damping parameters
beta_i  = 0.05;
beta_s  = 0.05;
beta_x  = 0.02;
beta_yy = 0.02;
n_s = 2;
n_x = 2;
n_yy = 2;

% Foundation selector:
% 1 = rectangular surface foundation
% 2 = square foundation
% 3 = strip foundation
foundation_type = 1;

foundation_labels = {'Rectangular', 'Square', 'Strip'};
cases = {'Fixed base', 'Mild SSI', 'Moderate SSI', 'Strong SSI'};
severity = [1.00, 0.85, 0.70, 0.55];

%% 2. COLUMN PROPERTIES
I = pi * D^4 / 64;
k_struct = 3 * E * I / h^3;
m = (T_target / (2*pi))^2 * k_struct;
T_fixed = 2*pi*sqrt(m/k_struct);
Delta_fixed = F / k_struct;
omega_fixed = 2*pi / T_fixed;

%% 3. AUXILIARY GEOMETRY
chi = A_ratio_rect(B, L);
a0 = omega_fixed * B / Vs;

Ibx = foundation_Ibx(B, L);
Iby = foundation_Iby(B, L);
Jb  = foundation_Jb(B, L);

%% 4. STATIC + DYNAMIC STIFFNESS
[Kx_stat, Ky_stat, Kz_stat, Krx_stat, Kry_stat, Kt_stat, ...
 kx_dyn, ky_dyn, kz_dyn, krx_dyn, kry_dyn, kt_dyn] = ...
    gazetas_foundation_stiffness(foundation_type, B, L, H, G, nu, chi, a0, Ibx, Iby, Jb);

Kx = Kx_stat * kx_dyn;
Ky = Ky_stat * ky_dyn;
Kz = Kz_stat * kz_dyn;
Krx = Krx_stat * krx_dyn;
Kry = Kry_stat * kry_dyn;
Kt = Kt_stat * kt_dyn;

kx_ref  = Kx;
kyy_ref = Kry;

%% 5. REFERENCE FLEXIBLE-BASE RESPONSE
T_ratio_ref = sqrt(1 + k_struct/kx_ref + (k_struct*h^2)/kyy_ref);
T_flex_ref = T_fixed * T_ratio_ref;

u_f_ref = F / kx_ref;
theta_ref = (F*h) / kyy_ref;
Delta_flex_ref = F/k_struct + u_f_ref + theta_ref*h;

Tx_ref  = 2*pi*sqrt(m/kx_ref);
Tyy_ref = 2*pi*sqrt((m*h^2)/kyy_ref);

beta_f_ref = (((T_flex_ref/T_fixed)^n_s - 1) / (T_flex_ref/T_fixed)^n_s) * beta_s + ...
             (1 / (T_flex_ref/Tx_ref)^n_x) * beta_x + ...
             (1 / (T_flex_ref/Tyy_ref)^n_yy) * beta_yy;

beta_0_ref = beta_f_ref + (1 / (T_flex_ref/T_fixed)^n_s) * beta_i;

%% 6. SCENARIO ANALYSIS
numCases = length(cases);

T_fixed_cases = T_fixed * ones(1, numCases);
T_flex_cases = zeros(1, numCases);
ratio_T_cases = zeros(1, numCases);
ratio_T2_cases = zeros(1, numCases);
delta_cases = zeros(1, numCases);
damping_cases = zeros(1, numCases);
kx_cases = zeros(1, numCases);
kyy_cases = zeros(1, numCases);

for i = 1:numCases
    kx_i = kx_ref * severity(i);
    kyy_i = kyy_ref * severity(i);

    kx_cases(i) = kx_i;
    kyy_cases(i) = kyy_i;

    ratio_T_cases(i) = sqrt(1 + k_struct/kx_i + (k_struct*h^2)/kyy_i);
    ratio_T2_cases(i) = ratio_T_cases(i)^2;
    T_flex_cases(i) = T_fixed * ratio_T_cases(i);

    uf_i = F / kx_i;
    theta_i = (F*h) / kyy_i;
    delta_cases(i) = F/k_struct + uf_i + theta_i*h;

    Tx_i = 2*pi*sqrt(m/kx_i);
    Tyy_i = 2*pi*sqrt((m*h^2)/kyy_i);

    beta_f_i = (((ratio_T_cases(i))^n_s - 1) / ((ratio_T_cases(i))^n_s)) * beta_s + ...
               (1 / (T_flex_cases(i)/Tx_i)^n_x) * beta_x + ...
               (1 / (T_flex_cases(i)/Tyy_i)^n_yy) * beta_yy;

    damping_cases(i) = beta_f_i + (1 / (ratio_T_cases(i)^n_s)) * beta_i;
end

%% 7. PRINT RESULTS
disp('============================================================');
disp('SSI MASTER PROGRAM - GAZETAS / NIST VERSION');
disp('============================================================');
fprintf('Foundation type = %s\n', foundation_labels{foundation_type});
fprintf('T_fixed = %.4f s\n', T_fixed);
fprintf('T_flex_ref = %.4f s\n', T_flex_ref);
fprintf('Delta_fixed = %.6e m\n', Delta_fixed);
fprintf('Delta_flex_ref = %.6e m\n', Delta_flex_ref);
fprintf('beta_0_ref = %.4f\n', beta_0_ref);

%% 8. CSV EXPORT
main_results = [E, D, h, I, k_struct, m, T_fixed, Delta_fixed, ...
                B, L, H, G, nu, rho_s, Vs, a0, ...
                Kx_stat, Ky_stat, Kz_stat, Krx_stat, Kry_stat, Kt_stat, ...
                kx_dyn, ky_dyn, kz_dyn, krx_dyn, kry_dyn, kt_dyn, ...
                Kx, Ky, Kz, Krx, Kry, Kt, ...
                T_flex_ref, T_ratio_ref, Delta_flex_ref, Tx_ref, Tyy_ref, beta_f_ref, beta_0_ref];

csvwrite(fullfile('output','ssi_gazetas_nist_main_results.csv'), main_results);

%% 9. PLOTS
x = 1:numCases;
bw = 0.38;

fig1 = figure('Color','w');
bar(x - bw/2, T_fixed_cases, bw, 'FaceColor', [0.20 0.72 0.84]);
hold on;
bar(x + bw/2, T_flex_cases, bw, 'FaceColor', [0.85 0.25 0.25]);
set(gca, 'XTick', x, 'XTickLabel', cases);
ylabel('Period (s)');
xlabel('Case');
title('Fixed vs flexible period (Gazetas/NIST SSI)');
legend({'T fixed', 'T flex'}, 'Location', 'northwest');
grid on; box on;
print(fig1, fullfile('output','ssi_gazetas_nist_period_compare.png'), '-dpng', '-r200');

fig2 = figure('Color','w');
plot(x, ratio_T2_cases, '-o', 'LineWidth', 1.8, 'Color', [0.20 0.72 0.84], ...
    'MarkerFaceColor', [0.20 0.72 0.84]);
hold on;
plot(x, damping_cases, '-o', 'LineWidth', 1.8, 'Color', [0.85 0.25 0.25], ...
    'MarkerFaceColor', [0.85 0.25 0.25]);
set(gca, 'XTick', x, 'XTickLabel', cases);
ylabel('Ratio');
xlabel('Case');
title('SSI ratio trends (Gazetas/NIST SSI)');
legend({'T~^2/T^2', 'Damping'}, 'Location', 'northwest');
grid on; box on;
print(fig2, fullfile('output','ssi_gazetas_nist_ratio_trends.png'), '-dpng', '-r200');

fig3 = figure('Color','w');
bar(x, delta_cases, 0.5, 'FaceColor', [0.45 0.60 0.90]);
set(gca, 'XTick', x, 'XTickLabel', cases);
ylabel('Top displacement (m)');
xlabel('Case');
title('Flexible-base displacement by scenario');
grid on; box on;
print(fig3, fullfile('output','ssi_gazetas_nist_displacement_trends.png'), '-dpng', '-r200');

fig4 = figure('Color','w');
plot(x, kx_cases, '-s', 'LineWidth', 1.8, 'Color', [0.15 0.55 0.20], ...
    'MarkerFaceColor', [0.15 0.55 0.20]);
hold on;
plot(x, kyy_cases, '-d', 'LineWidth', 1.8, 'Color', [0.55 0.20 0.60], ...
    'MarkerFaceColor', [0.55 0.20 0.60]);
set(gca, 'XTick', x, 'XTickLabel', cases);
ylabel('Stiffness');
xlabel('Case');
title('Foundation stiffness trends (Gazetas/NIST SSI)');
legend({'k_x', 'k_{yy}'}, 'Location', 'northwest');
grid on; box on;
print(fig4, fullfile('output','ssi_gazetas_nist_stiffness_trends.png'), '-dpng', '-r200');

disp('Done. Files saved in output/.');

end

%% ==========================================================
%% SUBFUNCTIONS
%% ==========================================================

function chi = A_ratio_rect(B, L)
Ab = 4 * B * L;
chi = Ab / (4 * L^2);
end

function Ibx = foundation_Ibx(B, L)
width = 2*B;
length_ = 2*L;
Ibx = width * length_^3 / 12;
end

function Iby = foundation_Iby(B, L)
width = 2*B;
length_ = 2*L;
Iby = length_ * width^3 / 12;
end

function Jb = foundation_Jb(B, L)
width = 2*B;
length_ = 2*L;
Jb = (width * length_) * (width^2 + length_^2) / 12;
end

function [Kx_stat, Ky_stat, Kz_stat, Krx_stat, Kry_stat, Kt_stat, ...
          kx_dyn, ky_dyn, kz_dyn, krx_dyn, kry_dyn, kt_dyn] = ...
          gazetas_foundation_stiffness(ftype, B, L, H, G, nu, chi, a0, Ibx, Iby, Jb)

switch ftype

    case 1
        Kz_stat = (2*G*L/(1-nu)) * (0.73 + 1.54 * chi^0.75);
        Ky_stat = (2*G*L/(2-nu)) * (2 + 2.5 * chi^0.85);
        Kx_stat = Ky_stat - (0.2/(0.75 - nu)) * G * L * (1 - B/L);
        Krx_stat = (G/(1-nu)) * Ibx^0.75 * (L/B)^0.25 * (2.4 + 0.5 * B/L);
        Kry_stat = (3*G/(1-nu)) * Iby^0.75 * (L/B)^0.15;
        Kt_stat = G * Jb^0.75 * (4 + 11 * (1 - B/L)^10);

        kz_dyn = 1.0;
        ky_dyn = 1.0;
        kx_dyn = 1.0;
        krx_dyn = max(0.05, 1 - 0.20 * a0);

        if nu < 0.45
            kry_dyn = max(0.05, 1 - 0.30 * a0);
        else
            kry_dyn = max(0.05, 1 - 0.25 * a0 * (L/B)^0.30);
        end

        kt_dyn = max(0.05, 1 - 0.14 * a0);

    case 2
        Kz_stat = 4.54 * G * B / (1 - nu);
        Ky_stat = 9 * G * B / (2 - nu);
        Kx_stat = Ky_stat;
        Krx_stat = 3.6 * G * B^3 / (1 - nu);
        Kry_stat = Krx_stat;
        Kt_stat = 8.3 * G * B^3;

        kz_dyn = 1.0;
        ky_dyn = 1.0;
        kx_dyn = 1.0;
        krx_dyn = max(0.05, 1 - 0.20 * a0);

        if nu < 0.45
            kry_dyn = max(0.05, 1 - 0.30 * a0);
        else
            kry_dyn = max(0.05, 1 - 0.25 * a0 * (L/B)^0.30);
        end

        kt_dyn = max(0.05, 1 - 0.14 * a0);

    case 3
        Kz_stat = 2*L * (0.73*G/(1-nu)) * (1 + 3.5 * B/H);
        Ky_stat = 2*L * (2*G/(2-nu)) * (1 + 2 * B/H);
        Kx_stat = Ky_stat;
        Krx_stat = 2*L * (pi * G * B^2 / (2*(1-nu))) * (1 + 0.2 * B/H);
        Kry_stat = Krx_stat;
        Kt_stat = Krx_stat;

        kz_dyn = 1.0;
        ky_dyn = 1.0;
        kx_dyn = 1.0;
        krx_dyn = max(0.05, 1 - 0.20 * a0);
        kry_dyn = max(0.05, 1 - 0.20 * a0);
        kt_dyn = max(0.05, 1 - 0.14 * a0);

    otherwise
        error('Invalid foundation type. Use 1, 2, or 3.');
end

end