clear; clc; close all;

%% INPUTS
E  = 30e9;        % Pa
D  = 1.0;         % m
h  = 3.0;         % m
T0 = 0.10;        % s
F  = 1.0e3;       % N

B  = 2.0;         % m
L  = 2.0;         % m
nu = 0.35;
G  = 80e6;        % Pa
rho_s = 1800;     % kg/m^3
Vs = sqrt(G/rho_s);

beta_i  = 0.05;
beta_s  = 0.05;
beta_x  = 0.02;
beta_yy = 0.02;
n_s = 2;
n_x = 2;
n_yy = 2;

%% SECTION PROPERTIES
I = pi*D^4/64;

%% FIXED BASE
k = 3*E*I/h^3;
m = (T0/(2*pi))^2 * k;
T = 2*pi*sqrt(m/k);
Delta_fixed = F/k;

%% SSI APPROXIMATION
kx  = 0.8*k;
kyy = 0.8*k*h^2;

T_tilde = 2*pi*sqrt(m / (k/(1 + k/kx + k*h^2/kyy)));
Tr = T_tilde/T;

u_f = F/kx;
theta = (F*h)/kyy;
Delta_tilde = F/k + u_f + theta*h;

Tx  = 2*pi*sqrt(m/kx);
Tyy = 2*pi*sqrt((m*h^2)/kyy);

beta_f = ((Tr^n_s - 1)/(Tr^n_s))*beta_s + (1/(T_tilde/Tx)^n_x)*beta_x + (1/(T_tilde/Tyy)^n_yy)*beta_yy;
beta_0 = beta_f + (1/(Tr^n_s))*beta_i;

T_tilde_Veletsos = T / sqrt(max(1e-12, (1 - (h/B)^2)));
T_tilde_ASCE = T * (1 + 2*h/B);

%% OUTPUT
disp('=== SSI COLUMN MODEL RESULTS ===');
fprintf('E = %.3e Pa\n',E);
fprintf('D = %.3f m\n',D);
fprintf('I = %.6f m^4\n',I);
fprintf('h = %.3f m\n',h);
fprintf('k = %.3e N/m\n',k);
fprintf('m = %.3e kg\n',m);
fprintf('T fixed = %.4f s\n',T);
fprintf('Delta fixed = %.6e m\n',Delta_fixed);
fprintf('kx = %.3e N/m\n',kx);
fprintf('kyy = %.3e N*m/rad\n',kyy);
fprintf('T flexible = %.4f s\n',T_tilde);
fprintf('T flexible / T = %.4f\n',Tr);
fprintf('Delta flexible = %.6e m\n',Delta_tilde);
fprintf('Tx = %.4f s\n',Tx);
fprintf('Tyy = %.4f s\n',Tyy);
fprintf('beta_f = %.4f\n',beta_f);
fprintf('beta_0 = %.4f\n',beta_0);
fprintf('T Veletsos = %.4f s\n',T_tilde_Veletsos);
fprintf('T ASCE = %.4f s\n',T_tilde_ASCE);

%% SAVE CSV
results = [
    E, D, I, h, k, m, T, Delta_fixed, kx, kyy, T_tilde, Tr, Delta_tilde, Tx, Tyy, beta_f, beta_0, T_tilde_Veletsos, T_tilde_ASCE
];
csvwrite('ssi_column_results.csv', results);

%% PLOTS
figure('Color','w');
bar([T, T_tilde]);
set(gca,'XTickLabel',{'Fixed base','Flexible base'});
ylabel('Period T (s)');
title('Period comparison');
grid on;

figure('Color','w');
bar([Delta_fixed, Delta_tilde]);
set(gca,'XTickLabel',{'Fixed base','Flexible base'});
ylabel('Top displacement \Delta (m)');
title('Displacement comparison');
grid on;

figure('Color','w');
bar([beta_i, beta_f, beta_0]);
set(gca,'XTickLabel',{'beta_i','beta_f','beta_0'});
ylabel('Damping ratio');
title('Damping components');
grid on;

figure('Color','w');
bar([beta_i beta_f beta_0]);
set(gca, 'XTickLabel', {'\beta_i', '\beta_f', '\beta_0'});
ylabel('Damping ratio');
title('Damping components');
grid on;