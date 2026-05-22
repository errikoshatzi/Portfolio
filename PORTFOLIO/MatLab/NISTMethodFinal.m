function MasterProgramScientific()
    clear; clc; close all;

    % Δημιουργία φακέλου για τα αποτελέσματα
    if ~exist('output', 'dir')
        mkdir('output');
    end

    disp('============================================================');
    disp(' SSI MASTER PROGRAM - ΔΙΑΔΡΑΣΤΙΚΗ ΕΚΔΟΣΗ (WOLF/GIVENS RULE) ');
    disp('============================================================');
    disp('Πατήστε ENTER για να κρατήσετε την προεπιλεγμένη τιμή [σε αγκύλες].');
    disp('------------------------------------------------------------');

    %% 1. ΔΙΑΔΡΑΣΤΙΚΗ ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ (USER INPUT)
    E = ask_user('Μέτρο Ελαστικότητας Κατασκευής E (Pa)', 30e9);
    D_col = ask_user('Διάμετρος/Διάσταση Υποστυλώματος D_col (m)', 1.0);
    h = ask_user('Ενεργό ύψος κατασκευής h (m)', 3.0);
    T_target = ask_user('Στοχευόμενη Ιδιοπερίοδος Κατασκευής T_target (s)', 0.10);
    F = ask_user('Οριζόντια Δύναμη Σχεδιασμού F (N)', 1000);

    disp('--- Ιδιότητες Εδάφους & Θεμελίωσης ---');
    G = ask_user('Μέτρο Διάτμησης Εδάφους G (Pa)', 80e6);
    nu = ask_user('Λόγος Poisson Εδάφους ν', 0.35);
    rho_s = ask_user('Πυκνότητα Εδάφους ρ_s (kg/m^3)', 1800);
    
    B = ask_user('Συνολικό Πλάτος Θεμελίου B (m)', 2.0);
    L = ask_user('Συνολικό Μήκος Θεμελίου L (m)', 3.0);
    D_emb = ask_user('Βάθος Εγκιβωτισμού Θεμελίου D_emb (m)', 1.0);

    disp('--- Παράμετροι Απόσβεσης ---');
    beta_i = ask_user('Δομική Απόσβεση β_i (ως δεκαδικό, π.χ. 0.05)', 0.05);
    beta_s = ask_user('Υστερητική Απόσβεση Εδάφους β_s (π.χ. 0.05)', 0.05);
    
    % Εκθέτες σύμφωνα με Givens (2013) και Wolf
    n_s = 2; n_x = 2; n_yy = 2;

    %% 2. ΥΠΟΛΟΓΙΣΜΟΙ ΑΝΩΔΟΜΗΣ (FIXED BASE)
    Vs = sqrt(G/rho_s);
    I = pi * D_col^4 / 64;
    k_struct = 3 * E * I / h^3;
    m = (T_target / (2*pi))^2 * k_struct;
    T_fixed = 2*pi*sqrt(m/k_struct);
    omega_fixed = 2*pi / T_fixed;
    Delta_fixed = F / k_struct;

    %% 3. ΑΚΡΙΒΗΣ ΥΠΟΛΟΓΙΣΜΟΣ ΕΜΠΕΔΗΣΗΣ (GAZETAS / PAIS & KAUSEL)
    [K_dyn, beta_rad] = gazetas_impedance_full(B, L, D_emb, G, nu, Vs, omega_fixed);
    Kx  = K_dyn.x;
    Kyy = K_dyn.ry;   
    beta_x  = beta_rad.x;
    beta_yy = beta_rad.ry; 

    %% 4. ΣΕΝΑΡΙΑ ΑΝΑΛΥΣΗΣ SSI
    cases = {'Fixed base', 'Mild SSI', 'Moderate SSI', 'Strong SSI'};
    severity = [1.00, 0.85, 0.70, 0.55]; % Μείωση δυσκαμψίας εδάφους λόγω παραμορφώσεων
    numCases = length(cases);

    T_fixed_cases = T_fixed * ones(1, numCases);
    T_flex_cases  = zeros(1, numCases);
    ratio_T_cases = zeros(1, numCases);
    damping_cases = zeros(1, numCases);
    delta_cases   = zeros(1, numCases);
    kx_cases      = zeros(1, numCases);
    kyy_cases     = zeros(1, numCases);

    for i = 1:numCases
        % Μειωμένη δυσκαμψία ανά σενάριο
        kx_i = Kx * severity(i);
        kyy_i = Kyy * severity(i);
        kx_cases(i) = kx_i;
        kyy_cases(i) = kyy_i;

        % Επιμήκυνση Περιόδου
        ratio_T = sqrt(1 + k_struct/kx_i + (k_struct*h^2)/kyy_i);
        T_flex = T_fixed * ratio_T;
        T_flex_cases(i) = T_flex;
        ratio_T_cases(i) = ratio_T;

        % Πλασματικές περίοδοι Tx και Tyy
        Tx = 2*pi*sqrt(m/kx_i);
        Tyy = 2*pi*sqrt((m*h^2)/kyy_i);

        % Κανόνας Ανάμειξης Απόσβεσης (Wolf / Roesset / Givens)
        term_hysteretic = beta_s * ( (ratio_T^n_s - 1) / ratio_T^n_s );
        term_rad_x      = beta_x * (1 / (T_flex/Tx)^n_x);
        term_rad_yy     = beta_yy * (1 / (T_flex/Tyy)^n_yy);
        
        beta_f = term_hysteretic + term_rad_x + term_rad_yy;
        beta_0 = beta_f + beta_i * (1 / ratio_T^n_s);
        damping_cases(i) = beta_0;

        % Μετακίνηση Κορυφής (Kinematics)
        uf_i = F / kx_i;
        theta_i = (F*h) / kyy_i;
        delta_cases(i) = F/k_struct + uf_i + theta_i*h;
    end

    %% 5. ΠΑΡΑΓΩΓΗ ΚΑΙ ΑΠΟΘΗΚΕΥΣΗ ΓΡΑΦΗΜΑΤΩΝ (4 Γραφήματα)
    disp('--- Παραγωγή Γραφημάτων ---');
    x = 1:numCases;
    bw = 0.38;

    % Γράφημα 1: Σύγκριση Περιόδων (Fixed vs Flexible)
    fig1 = figure('Color','w');
    bar(x - bw/2, T_fixed_cases, bw, 'FaceColor', [0.20 0.72 0.84]); hold on;
    bar(x + bw/2, T_flex_cases, bw, 'FaceColor', [0.85 0.25 0.25]);
    set(gca, 'XTick', x, 'XTickLabel', cases, 'FontSize', 10);
    ylabel('Περίοδος T (s)', 'FontWeight', 'bold');
    title('Επιμήκυνση Ιδιοπεριόδου λόγω SSI', 'FontSize', 12);
    legend({'T_{fixed} (Πακτωμένη)', 'T_{flex} (Ευλύγιστη)'}, 'Location', 'northwest');
    grid on; box on;
    print(fig1, fullfile('output','01_period_comparison.png'), '-dpng', '-r300');

    % Γράφημα 2: SSI Ratio & Απόσβεση (ΔΙΟΡΘΩΜΕΝΟ ΓΙΑ OCTAVE ΜΕ PLOTYY)
    fig2 = figure('Color','w');
    [ax, h1, h2] = plotyy(x, ratio_T_cases, x, damping_cases);
    
    set(h1, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2, 'Color', [0.20 0.40 0.84], 'MarkerFaceColor', [0.20 0.40 0.84]);
    set(h2, 'LineStyle', '-', 'Marker', 's', 'LineWidth', 2, 'Color', [0.85 0.25 0.25], 'MarkerFaceColor', [0.85 0.25 0.25]);
    
    ylabel(ax(1), 'Λόγος T_{flex} / T_{fixed}', 'FontWeight', 'bold');
    ylabel(ax(2), 'Συνολική Απόσβεση Συστήματος \beta_0', 'FontWeight', 'bold');
    
    set(ax(1), 'XTick', x, 'XTickLabel', cases, 'FontSize', 10);
    set(ax(2), 'XTick', x, 'XTickLabel', []); % Για να μην κάνουν επικάλυψη τα x-labels
    
    title('Τάσεις Επιμήκυνσης και Απόσβεσης (Κανόνας Wolf)', 'FontSize', 12);
    grid(ax(1), 'on'); 
    
    print(fig2, fullfile('output','02_ratio_and_damping.png'), '-dpng', '-r300');

    % Γράφημα 3: Μετακινήσεις Κορυφής (Top Displacement)
    fig3 = figure('Color','w');
    bar(x, delta_cases*1000, 0.5, 'FaceColor', [0.45 0.60 0.90]); % Μετατροπή σε mm για καλύτερη οπτικοποίηση
    set(gca, 'XTick', x, 'XTickLabel', cases, 'FontSize', 10);
    ylabel('Συνολική Μετακίνηση Κορυφής (mm)', 'FontWeight', 'bold');
    title('Επίδραση SSI στη Συνολική Μετακίνηση (\Delta_{total})', 'FontSize', 12);
    grid on; box on;
    print(fig3, fullfile('output','03_displacement_trends.png'), '-dpng', '-r300');

    % Γράφημα 4: Τάσεις Δυσκαμψίας Θεμελίου (Stiffness Trends)
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

    disp('============================================================');
    disp('ΟΛΟΚΛΗΡΩΘΗΚΕ! Τα 4 γραφήματα αποθηκεύτηκαν στον φάκελο "output".');
    disp('============================================================');
end

%% ==========================================================
%% SUBFUNCTIONS (ΒΟΗΘΗΤΙΚΕΣ ΣΥΝΑΡΤΗΣΕΙΣ)
%% ==========================================================

% Συνάρτηση για διαδραστική εισαγωγή με προεπιλογή (default)
function val = ask_user(prompt_text, default_val)
    str = sprintf('%s [%g]: ', prompt_text, default_val);
    user_input = input(str);
    if isempty(user_input)
        val = default_val;
    else
        val = user_input;
    end
end

% Συνάρτηση πλήρους υπολογισμού ελατηρίων (Gazetas / Pais & Kausel)
function [K_dyn, beta_rad] = gazetas_impedance_full(B_total, L_total, D_emb, G, nu, Vs, omega)
    b = B_total / 2; 
    l = L_total / 2;
    
    % Προσανατολισμός L >= B
    if l < b
        temp = l; l = b; b = temp;
    end
    
    a0 = omega * b / Vs;
    psi = min(sqrt(2*(1-nu)/max(1e-5,1-2*nu)), 2.5);

    % Γεωμετρικά χαρακτηριστικά
    chi = (4*b*l) / (4*l^2);
    Ibx = (2*l) * (2*b)^3 / 12; 
    Iby = (2*b) * (2*l)^3 / 12; 
    Jb  = Ibx + Iby;

    % Στατική Δυσκαμψία Επιφάνειας (Πίνακας 2.2)
    Ky_sur = (2*G*l / (2-nu)) * (2 + 2.5 * chi^0.85);
    Kx_sur = Ky_sur - (0.2 / (0.75 - nu)) * G * l * (1 - b/l);
    Kry_sur = (3*G / (1-nu)) * Iby^0.75 * (l/b)^0.15;

    % Συντελεστές Εγκιβωτισμού (Πίνακας 2.5)
    eta_y = 1.0 + (0.33 + 1.34/(1 + l/b)) * (D_emb/b)^0.8;
    eta_x = eta_y; 
    eta_ry = 1.0 + D_emb/b + 1.6 / (0.35 + (l/b)^4) * (D_emb/b)^2; 

    % Εγκιβωτισμένη Στατική Δυσκαμψία
    Kx_emb = Kx_sur * eta_x;
    Kry_emb = Kry_sur * eta_ry;

    % Δυναμικοί Συντελεστές (Πίνακας 2.2 / 2.6)
    alpha_x = 1.0;
    if nu < 0.45
        alpha_ry = max(0.01, 1 - 0.30 * a0);
    else
        alpha_ry = max(0.01, 1 - 0.25 * a0 * (l/b)^0.30);
    end

    % Τελική Δυναμική Δυσκαμψία
    K_dyn.x  = Kx_emb * alpha_x;
    K_dyn.ry = Kry_emb * alpha_ry;

    % Απόσβεση Ακτινοβολίας για Εγκιβωτισμένα Πέδιλα (Πίνακας 2.7)
    beta_rad.x = (4 * (l/b + (D_emb/b)*(psi + l/b)) / (Kx_emb / (G*b))) * (a0 / (2*alpha_x));

    term1_ry = (4/3) * ( (l/b)^3 * (D_emb/b) + psi*(D_emb/b)*(l/b) + (D_emb/b)^3 + ...
               psi*(l/b)^3 + 3*(D_emb/b)*(l/b)^2 + psi*(l/b) ) * a0^2;
    term2_ry = (Kry_emb / (G*b^3)) * ( 1.8 / (1 + 1.75*(l/b - 1)) + a0^2 );
    term3_ry = (4/3) * (l/b + psi) * (D_emb/b)^3 / (Kry_emb / (G*b^3));
    beta_rad.ry = (term1_ry / term2_ry + term3_ry) * (a0 / (2*alpha_ry));
end