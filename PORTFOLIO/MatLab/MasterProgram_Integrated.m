function MasterProgram_Integrated()
    clear; clc; close all;

    %% ==========================================================
    %% OUTPUT FOLDER
    %% ==========================================================
    if ~exist('output', 'dir')
        mkdir('output');
    end

    %% ==========================================================
    %% 1. STRUCTURAL INPUT
    %% ==========================================================
    fprintf('=============================================================\n');
    fprintf(' SSI MASTER PROGRAM - WINKLER ESTIMATE & RIGOROUS NIST METHOD\n');
    fprintf('=============================================================\n\n');

    fprintf('STRUCTURAL INPUT\n');
    fprintf('----------------\n');
    E = ask_user('Enter column modulus E (Pa)', 30e9);
    D_col = ask_user('Enter circular column diameter D_col (m)', 1.0);
    h = ask_user('Enter effective height h (m)', 3.0);
    T_target = ask_user('Enter fixed-base target period T_fixed (s)', 0.10);
    F = ask_user('Enter lateral force F (N)', 1000);

    if E <= 0 || D_col <= 0 || h <= 0 || T_target <= 0 || F <= 0
        error('Invalid structural input.');
    end

    %% ==========================================================
    %% 2. FOUNDATION INPUT
    %% ==========================================================
    fprintf('\nFOUNDATION INPUT\n');
    fprintf('----------------\n');
    Bp = 0.305;     % reference plate width (m)
    B = ask_user('Enter footing width B (m)', 2.0);
    L = ask_user('Enter footing length L (m)', 3.0);
    Df = ask_user('Enter embedment depth Df (m)', 1.0);

    if B <= 0 || L <= 0 || Df < 0
        error('Invalid foundation geometry input.');
    end

    if B > L
        temp = B; B = L; L = temp;
        fprintf('\nNote: dimensions were swapped so that B <= L.\n');
    end

    %% ==========================================================
    %% 3. SOIL TYPE
    %% ==========================================================
    fprintf('\nSOIL TYPE SELECTION\n');
    fprintf('-------------------\n');
    fprintf('  [1] Sand\n');
    fprintf('  [2] Clay\n');
    soil_type = ask_user('Enter soil type number [1 or 2]', 1);

    soil_type_name = '';
    soil_class_name = '';
    Ks_table_MN = 0;
    nu = 0;
    rho_t = 0;      % t/m^3
    G_kPa = 0;      % kPa

    if soil_type == 1
        soil_type_name = 'Sand';
        fprintf('\nSelect sand class:\n');
        fprintf('  [1] Loose sand\n');
        fprintf('  [2] Medium sand\n');
        fprintf('  [3] Dense sand\n');
        soil_class = ask_user('Enter sand class number', 2);

        switch soil_class
            case 1
                soil_class_name = 'Loose sand'; Ks_table_MN = 12.9; nu = 0.275; rho_t = 1.8; G_kPa = 20000;
            case 2
                soil_class_name = 'Medium sand'; Ks_table_MN = 41.7; nu = 0.325; rho_t = 1.9; G_kPa = 55000;
            case 3
                soil_class_name = 'Dense sand'; Ks_table_MN = 161.0; nu = 0.350; rho_t = 2.0; G_kPa = 140000;
            otherwise
                error('Invalid sand class.');
        end

    elseif soil_type == 2
        soil_type_name = 'Clay';
        fprintf('\nSelect clay class:\n');
        fprintf('  [1] Stiff clay\n');
        fprintf('  [2] Very stiff clay\n');
        fprintf('  [3] Hard clay\n');
        soil_class = ask_user('Enter clay class number', 2);

        switch soil_class
            case 1
                soil_class_name = 'Stiff clay'; Ks_table_MN = 24.1; nu = 0.40; rho_t = 1.75; G_kPa = 25000;
            case 2
                soil_class_name = 'Very stiff clay'; Ks_table_MN = 48.2; nu = 0.44; rho_t = 1.85; G_kPa = 60000;
            case 3
                soil_class_name = 'Hard clay'; Ks_table_MN = 96.4; nu = 0.47; rho_t = 1.90; G_kPa = 115000;
            otherwise
                error('Invalid clay class.');
        end
    else
        error('Invalid soil type.');
    end

    %% ==========================================================
    %% 4. Ks INPUT MODE & DAMPING
    %% ==========================================================
    fprintf('\nKs INPUT MODE\n');
    fprintf('-------------\n');
    fprintf('  [1] Use Ks from soil table\n');
    fprintf('  [2] Enter Ks manually\n');
    Ks_mode = ask_user('Choose Ks input mode', 1);

    if Ks_mode == 1
        Ks_MN = Ks_table_MN;
        Ks_source = 'Soil table';
    elseif Ks_mode == 2
        Ks_MN = ask_user('Enter manual Ks value (MN/m^3)', 40);
        Ks_source = 'Manual input';
    end

    fprintf('\nDAMPING PARAMETERS\n');
    fprintf('------------------\n');
    beta_i = ask_user('Enter structural damping beta_i (e.g. 0.05)', 0.05);
    beta_s = ask_user('Enter soil hysteretic damping beta_s (e.g. 0.05)', 0.05);

    %% ==========================================================
    %% 5. UNIT CONVERSIONS & FIXED-BASE SDOF
    %% ==========================================================
    Ks    = Ks_MN * 1000;       % kN/m^3
    rho_s = rho_t * 1000;       % kg/m^3
    G     = G_kPa * 1000;       % Pa
    Vs    = sqrt(G / rho_s);    % m/s

    I = pi * D_col^4 / 64;
    k_struct = 3 * E * I / h^3;                     % N/m
    m = (T_target / (2*pi))^2 * k_struct;           % kg
    T_fixed = 2*pi*sqrt(m/k_struct);                % s
    Delta_fixed = F / k_struct;                     % m
    omega_fixed = 2*pi / T_fixed;                   % rad/s
    a0 = omega_fixed * B / (2 * Vs);

    %% ==========================================================
    %% 6. PART A: WINKLER SPRINGS ESTIMATE
    %% ==========================================================
    if soil_type == 1
        n_delta = ((B + Bp) / (2 * B))^2;
        n_beta = 1 + 2 * (Df / B);
    else
        n_delta = Bp / B;
        n_beta = 1.0;
    end
    n_sigmax = (1/3) * (2 + B / L);
    Ks_tel = n_delta * n_sigmax * n_beta * Ks;      % kN/m^3

    Kz_kN  = Ks_tel * L * B;                        % kN/m
    Krl_kN = Ks_tel * L * B^3 / 12;                 % kN*m/rad
    Krb_kN = Ks_tel * L^3 * B / 12;                 % kN*m/rad
    Kx_kN  = 0.30 * Kz_kN;                          % kN/m (rough estimate)

    kz_dyn = 1.00;
    kx_dyn = 1.00;
    krl_dyn = max(0.25, 1 - 0.15 * a0);
    if nu < 0.45
        krb_dyn = max(0.25, 1 - 0.20 * a0);
    else
        krb_dyn = max(0.25, 1 - 0.15 * a0 * (L / B)^0.30);
    end

    Kx_dynamic_kN  = kx_dyn  * Kx_kN;
    Kz_dynamic_kN  = kz_dyn  * Kz_kN;
    Krl_dynamic_kN = krl_dyn * Krl_kN;
    Krb_dynamic_kN = krb_dyn * Krb_kN;

    Kx_dynamic  = Kx_dynamic_kN  * 1000;          % N/m
    Krb_dynamic = Krb_dynamic_kN * 1000;          % N*m/rad
    kyy_ref = Krb_dynamic;

    T_ratio = sqrt(1 + k_struct/Kx_dynamic + (k_struct*h^2)/kyy_ref);
    T_flex_est = T_fixed * T_ratio;

    u_f = F / Kx_dynamic;
    theta = (F*h) / kyy_ref;
    Delta_flex_est = F/k_struct + u_f + theta*h;

    Tx_est  = 2*pi*sqrt(m/Kx_dynamic);
    Tyy_est = 2*pi*sqrt((m*h^2)/kyy_ref);

    %% ==========================================================
    %% 7. PART B: RIGOROUS GAZETAS / NIST METHOD & SCENARIOS
    %% ==========================================================
    [K_dyn, beta_rad] = gazetas_impedance_full(B, L, Df, G, nu, Vs, omega_fixed);
    Kx_rig  = K_dyn.x;
    Kyy_rig = K_dyn.ry;   
    beta_x  = beta_rad.x;
    beta_yy = beta_rad.ry; 

    cases = {'Fixed base', 'Mild SSI', 'Moderate SSI', 'Strong SSI'};
    severity = [1.00, 0.85, 0.70, 0.55]; 
    numCases = length(cases);

    T_fixed_cases     = T_fixed * ones(1, numCases);
    
    % Arrays for Winkler
    T_flex_wink_cases = zeros(1, numCases);
    delta_wink_cases  = zeros(1, numCases);

    % Arrays for Gazetas / NIST
    T_flex_cases      = zeros(1, numCases);
    ratio_T_cases     = zeros(1, numCases);
    damping_cases     = zeros(1, numCases);
    delta_cases       = zeros(1, numCases);
    kx_cases          = zeros(1, numCases);
    kyy_cases         = zeros(1, numCases);
    
    n_s = 2; n_x = 2; n_yy = 2; % Wolf / Givens exponents

    for i = 1:numCases
        % ================= WINKLER SCENARIOS =================
        kx_wink_i = Kx_dynamic * severity(i);
        kyy_wink_i = Krb_dynamic * severity(i);
        
        ratio_T_wink = sqrt(1 + k_struct/kx_wink_i + (k_struct*h^2)/kyy_wink_i);
        T_flex_wink_cases(i) = T_fixed * ratio_T_wink;
        delta_wink_cases(i)  = F/k_struct + (F / kx_wink_i) + ((F*h) / kyy_wink_i)*h;

        % ================= GAZETAS/NIST SCENARIOS =================
        kx_i = Kx_rig * severity(i);
        kyy_i = Kyy_rig * severity(i);
        kx_cases(i) = kx_i; kyy_cases(i) = kyy_i;

        ratio_T = sqrt(1 + k_struct/kx_i + (k_struct*h^2)/kyy_i);
        T_flex = T_fixed * ratio_T;
        T_flex_cases(i) = T_flex;
        ratio_T_cases(i) = ratio_T;

        Tx = 2*pi*sqrt(m/kx_i);
        Tyy = 2*pi*sqrt((m*h^2)/kyy_i);

        % Wolf / Roesset / Givens Damping Rule
        term_hysteretic = beta_s * ( (ratio_T^n_s - 1) / ratio_T^n_s );
        term_rad_x      = beta_x * (1 / (T_flex/Tx)^n_x);
        term_rad_yy     = beta_yy * (1 / (T_flex/Tyy)^n_yy);
        
        beta_f = term_hysteretic + term_rad_x + term_rad_yy;
        damping_cases(i) = beta_f + beta_i * (1 / ratio_T^n_s);

        delta_cases(i) = F/k_struct + (F / kx_i) + ((F*h) / kyy_i)*h;
    end

    %% ==========================================================
    %% 8. EXACT PRINTOUT OF RESULTS
    %% ==========================================================
    fprintf('\n=============================================================\n');
    fprintf('RESULTS (PART A: WINKLER METHOD)\n');
    fprintf('=============================================================\n');
    fprintf('Soil type            = %s\n', soil_type_name);
    fprintf('Soil class           = %s\n', soil_class_name);
    fprintf('Ks source            = %s\n', Ks_source);
    fprintf('B                    = %.3f m\n', B);
    fprintf('L                    = %.3f m\n', L);
    fprintf('Df                   = %.3f m\n', Df);
    fprintf('Bp                   = %.3f m\n', Bp);

    fprintf('\n--- Soil properties ---\n');
    fprintf('Ks                   = %.3f MN/m^3\n', Ks_MN);
    fprintf('Ks                   = %.3f kN/m^3\n', Ks);
    fprintf('Ks_tel               = %.3f kN/m^3\n', Ks_tel);
    fprintf('nu                   = %.4f\n', nu);
    fprintf('rho                  = %.1f kg/m^3\n', rho_s);
    fprintf('G                    = %.3f kPa\n', G_kPa);
    fprintf('Vs                   = %.3f m/s\n', Vs);
    fprintf('a0                   = %.4f\n', a0);

    fprintf('\n--- Correction factors ---\n');
    fprintf('n_delta              = %.4f\n', n_delta);
    fprintf('n_sigmax             = %.4f\n', n_sigmax);
    fprintf('n_beta               = %.4f\n', n_beta);

    fprintf('\n--- Structural properties ---\n');
    fprintf('I                    = %.6e m^4\n', I);
    fprintf('k_struct             = %.6e N/m\n', k_struct);
    fprintf('m                    = %.6f kg\n', m);
    fprintf('T_fixed              = %.4f s\n', T_fixed);
    fprintf('omega_fixed          = %.4f rad/s\n', omega_fixed);
    fprintf('Delta_fixed          = %.6e m\n', Delta_fixed);

    fprintf('\n--- Foundation spring constants ---\n');
    fprintf('Kx_dynamic           = %.3f kN/m\n', Kx_dynamic_kN);
    fprintf('Kz_dynamic           = %.3f kN/m\n', Kz_dynamic_kN);
    fprintf('Krl_dynamic          = %.3f kN*m/rad\n', Krl_dynamic_kN);
    fprintf('Krb_dynamic          = %.3f kN*m/rad\n', Krb_dynamic_kN);

    fprintf('\n--- SSI period pre-estimate ---\n');
    fprintf('Tx_est               = %.4f s\n', Tx_est);
    fprintf('Tyy_est              = %.4f s\n', Tyy_est);
    fprintf('T_ratio              = %.4f\n', T_ratio);
    fprintf('T_flex_est           = %.4f s\n', T_flex_est);
    fprintf('Delta_flex_est       = %.6e m\n', Delta_flex_est);

    fprintf('\n=============================================================\n');
    fprintf('RESULTS (PART B: RIGOROUS NIST / GAZETAS METHOD)\n');
    fprintf('=============================================================\n');
    fprintf('Kx_rigorous          = %.2e N/m\n', Kx_rig);
    fprintf('Kyy_rigorous         = %.2e N*m/rad\n', Kyy_rig);
    fprintf('T_flex_NIST          = %.4f s\n', T_flex_cases(1));
    fprintf('System Damping B_0   = %.4f\n', damping_cases(1));
    fprintf('=============================================================\n');

    %% ==========================================================
    %% 9. GENERATE PLOTS (Octave Compatible)
    %% ==========================================================
    disp('Generating comparative plots...');
    x = 1:numCases;

    % Γράφημα 1: Σύγκριση Περιόδων (Τρεις μπάρες: Fixed, Winkler, NIST)
    fig1 = figure('Color','w');
    y_data = [T_fixed_cases; T_flex_wink_cases; T_flex_cases]';
    b = bar(x, y_data, 0.85); 
    b(1).FaceColor = [0.20 0.72 0.84]; % Fixed
    b(2).FaceColor = [0.85 0.55 0.20]; % Winkler
    b(3).FaceColor = [0.85 0.25 0.25]; % NIST
    
    set(gca, 'XTick', x, 'XTickLabel', cases, 'FontSize', 10);
    ylabel('Περίοδος T (s)', 'FontWeight', 'bold');
    title('Σύγκριση Ιδιοπεριόδου: Winkler vs. NIST/Gazetas', 'FontSize', 12);
    legend({'T_{fixed} (Πακτωμένη)', 'T_{flex} (Winkler)', 'T_{flex} (NIST/Gazetas)'}, 'Location', 'northwest');
    grid on; box on;
    print(fig1, fullfile('output','01_period_comparison_combined.png'), '-dpng', '-r300');

    % Γράφημα 2: SSI Ratio & Απόσβεση (Μόνο για NIST/Gazetas)
    fig2 = figure('Color','w');
    [ax, h1, h2] = plotyy(x, ratio_T_cases, x, damping_cases);
    set(h1, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2, 'Color', [0.20 0.40 0.84], 'MarkerFaceColor', [0.20 0.40 0.84]);
    set(h2, 'LineStyle', '-', 'Marker', 's', 'LineWidth', 2, 'Color', [0.85 0.25 0.25], 'MarkerFaceColor', [0.85 0.25 0.25]);
    ylabel(ax(1), 'Λόγος T_{flex} / T_{fixed}', 'FontWeight', 'bold');
    ylabel(ax(2), 'Συνολική Απόσβεση Συστήματος \beta_0', 'FontWeight', 'bold');
    set(ax(1), 'XTick', x, 'XTickLabel', cases, 'FontSize', 10);
    set(ax(2), 'XTick', x, 'XTickLabel', []);
    title('Τάσεις Επιμήκυνσης και Απόσβεσης (NIST/Wolf)', 'FontSize', 12);
    grid(ax(1), 'on'); 
    print(fig2, fullfile('output','02_ratio_and_damping.png'), '-dpng', '-r300');

    % Γράφημα 3: Σύγκριση Μετακινήσεων Κορυφής (Δύο μπάρες: Winkler vs NIST)
    fig3 = figure('Color','w');
    y_disp = [delta_wink_cases; delta_cases]' * 1000; % Σε mm
    b_disp = bar(x, y_disp, 0.7);
    b_disp(1).FaceColor = [0.85 0.55 0.20]; % Winkler
    b_disp(2).FaceColor = [0.45 0.60 0.90]; % NIST

    set(gca, 'XTick', x, 'XTickLabel', cases, 'FontSize', 10);
    ylabel('Συνολική Μετακίνηση Κορυφής (mm)', 'FontWeight', 'bold');
    title('Σύγκριση Μετακίνησης \Delta_{total}: Winkler vs NIST', 'FontSize', 12);
    legend({'\Delta_{total} (Winkler)', '\Delta_{total} (NIST)'}, 'Location', 'northwest');
    grid on; box on;
    print(fig3, fullfile('output','03_displacement_combined.png'), '-dpng', '-r300');

    % Γράφημα 4: Τάσεις Δυσκαμψίας (Μόνο για NIST/Gazetas)
    fig4 = figure('Color','w');
    plot(x, kx_cases, '-s', 'LineWidth', 2, 'Color', [0.15 0.55 0.20], 'MarkerSize', 8, 'MarkerFaceColor', [0.15 0.55 0.20]);
    hold on;
    plot(x, kyy_cases, '-d', 'LineWidth', 2, 'Color', [0.55 0.20 0.60], 'MarkerSize', 8, 'MarkerFaceColor', [0.55 0.20 0.60]);
    set(gca, 'XTick', x, 'XTickLabel', cases, 'FontSize', 10);
    ylabel('Δυσκαμψία Θεμελίου (N/m ή N\cdotm/rad)', 'FontWeight', 'bold');
    title('Μείωση Δυσκαμψίας Εδάφους-Θεμελίωσης ανά Σενάριο', 'FontSize', 12);
    legend({'Οριζόντια Δυσκαμψία (K_x)', 'Στροφική Δυσκαμψία (K_{yy})'}, 'Location', 'northeast');
    grid on; box on;
    print(fig4, fullfile('output','04_stiffness_trends.png'), '-dpng', '-r300');

    disp('ΟΛΟΚΛΗΡΩΘΗΚΕ! Τα συγκριτικά γραφήματα αποθηκεύτηκαν στον φάκελο "output".');
end

%% ==========================================================
%% SUBFUNCTIONS 
%% ==========================================================
function val = ask_user(prompt_text, default_val)
    str = sprintf('%s [%g]: ', prompt_text, default_val);
    user_input = input(str);
    if isempty(user_input)
        val = default_val;
    else
        val = user_input;
    end
end

function [K_dyn, beta_rad] = gazetas_impedance_full(B_total, L_total, D_emb, G, nu, Vs, omega)
    b = B_total / 2; 
    l = L_total / 2;
    if l < b
        temp = l; l = b; b = temp;
    end
    
    a0 = omega * b / Vs;
    psi = min(sqrt(2*(1-nu)/max(1e-5,1-2*nu)), 2.5);

    chi = (4*b*l) / (4*l^2);
    Ibx = (2*l) * (2*b)^3 / 12; 
    Iby = (2*b) * (2*l)^3 / 12; 
    
    Ky_sur = (2*G*l / (2-nu)) * (2 + 2.5 * chi^0.85);
    Kx_sur = Ky_sur - (0.2 / (0.75 - nu)) * G * l * (1 - b/l);
    Kry_sur = (3*G / (1-nu)) * Iby^0.75 * (l/b)^0.15;

    eta_y = 1.0 + (0.33 + 1.34/(1 + l/b)) * (D_emb/b)^0.8;
    eta_x = eta_y; 
    eta_ry = 1.0 + D_emb/b + 1.6 / (0.35 + (l/b)^4) * (D_emb/b)^2; 

    Kx_emb = Kx_sur * eta_x;
    Kry_emb = Kry_sur * eta_ry;

    alpha_x = 1.0;
    if nu < 0.45
        alpha_ry = max(0.01, 1 - 0.30 * a0);
    else
        alpha_ry = max(0.01, 1 - 0.25 * a0 * (l/b)^0.30);
    end

    K_dyn.x  = Kx_emb * alpha_x;
    K_dyn.ry = Kry_emb * alpha_ry;

    beta_rad.x = (4 * (l/b + (D_emb/b)*(psi + l/b)) / (Kx_emb / (G*b))) * (a0 / (2*alpha_x));

    term1_ry = (4/3) * ( (l/b)^3 * (D_emb/b) + psi*(D_emb/b)*(l/b) + (D_emb/b)^3 + ...
               psi*(l/b)^3 + 3*(D_emb/b)*(l/b)^2 + psi*(l/b) ) * a0^2;
    term2_ry = (Kry_emb / (G*b^3)) * ( 1.8 / (1 + 1.75*(l/b - 1)) + a0^2 );
    term3_ry = (4/3) * (l/b + psi) * (D_emb/b)^3 / (Kry_emb / (G*b^3));
    beta_rad.ry = (term1_ry / term2_ry + term3_ry) * (a0 / (2*alpha_ry));
end