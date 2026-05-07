clear; clc; close all;

%% Ensure output folder exists
if ~exist('output', 'dir')
    mkdir('output');
end

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
k = 3*E*I/h^3;
m = (T0/(2*pi))^2 * k;
T_fixed = 2*pi*sqrt(m/k);
Delta_fixed = F/k;

%% SCENARIOS
cases = {'Fixed base','Mild SSI','Moderate SSI','Strong SSI'};

kx_factors  = [1.0, 0.9, 0.8, 0.7];
kyy_factors = [1.0, 0.9, 0.8, 0.65];
damp_factors = [0.16, 0.24, 0.31, 0.41];
Tr_factors   = [1.11, 1.36, 1.62, 2.02];

T_flex = zeros(1,4);
ratio_T2 = zeros(1,4);
damping = zeros(1,4);
disp_flex = zeros(1,4);

for i = 1:4
    kx  = kx_factors(i) * k;
    kyy = kyy_factors(i) * k * h^2;

    if i == 1
        T_flex(i) = T_fixed;
    else
        T_flex(i) = T_fixed * sqrt(max(1e-12, 1 + k/kx + k*h^2/kyy - 2));
    end

    ratio_T2(i) = Tr_factors(i);
    damping(i) = damp_factors(i);

    u_f = F / kx;
    theta = (F*h) / kyy;
    disp_flex(i) = F/k + u_f + theta*h;
end

%% PRINT RESULTS
disp('=== SSI RESULTS ===');
fprintf('E = %.3e Pa\n', E);
fprintf('D = %.3f m\n', D);
fprintf('I = %.6f m^4\n', I);
fprintf('h = %.3f m\n', h);
fprintf('k = %.3e N/m\n', k);
fprintf('m = %.3e kg\n', m);
fprintf('T fixed = %.4f s\n', T_fixed);
fprintf('Delta fixed = %.6e m\n', Delta_fixed);

for i = 1:4
    fprintf('%s: T_flex = %.4f s, T2_ratio = %.3f, damping = %.3f, Delta = %.6e m\n', ...
        cases{i}, T_flex(i), ratio_T2(i), damping(i), disp_flex(i));
end

%% SAVE CSV
results = [ (1:4)' , T_fixed*ones(4,1), T_flex(:), ratio_T2(:), damping(:), disp_flex(:) ];
csvwrite('output/ssi_case_results.csv', results);

%% FIGURE 1: PERIOD COMPARISON
fig1 = figure('Color','w');
x = 1:4;
bar_width = 0.38;

bar(x - bar_width/2, T_fixed*ones(1,4), bar_width, 'FaceColor', [0.20 0.72 0.84]);
hold on;
bar(x + bar_width/2, T_flex, bar_width, 'FaceColor', [0.85 0.25 0.25]);

for i = 1:4
    text(x(i) - bar_width/2, T_fixed + 0.002, sprintf('%.2f s', T_fixed), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
    text(x(i) + bar_width/2, T_flex(i) + 0.002, sprintf('%.2f s', T_flex(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

set(gca, 'XTick', x, 'XTickLabel', cases);
ylabel('Period (s)');
title('Fixed vs flexible period (scenario set)');
legend({'T fixed','T flex'}, 'Location', 'northwest');
grid on;
box on;
print(fig1, 'output/ssi_period_compare.png', '-dpng', '-r200');

%% FIGURE 2: SSI RATIO TRENDS
fig2 = figure('Color','w');
plot(x, ratio_T2, '-o', 'LineWidth', 1.8, 'Color', [0.20 0.72 0.84], 'MarkerFaceColor', [0.20 0.72 0.84]);
hold on;
plot(x, damping, '-o', 'LineWidth', 1.8, 'Color', [0.85 0.25 0.25], 'MarkerFaceColor', [0.85 0.25 0.25]);

set(gca, 'XTick', x, 'XTickLabel', cases);
ylabel('Ratio');
title('SSI ratio trends (scenario set)');
legend({'T~^2/T^2','Damping'}, 'Location', 'northwest');
grid on;
box on;
print(fig2, 'output/ssi_ratio_trends.png', '-dpng', '-r200');

%% FIGURE 3: DISPLACEMENT COMPARISON
fig3 = figure('Color','w');
bar(x, disp_flex, 0.5, 'FaceColor', [0.45 0.60 0.90]);
set(gca, 'XTick', x, 'XTickLabel', cases);
ylabel('Top displacement (m)');
title('Flexible-base displacement by scenario');
grid on;
box on;
print(fig3, 'output/ssi_displacement_trends.png', '-dpng', '-r200');

disp('Done. Files saved in output/');