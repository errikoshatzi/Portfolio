function SSI_2Methods_Estimation_Square()
    clear; clc; close all;

    if ~exist('output', 'dir')
        mkdir('output');
    end

    fprintf('=============================================================\n');
    fprintf(' SSI MASTER PROGRAM - WINKLER ESTIMATE & RIGOROUS NIST METHOD\n');
    fprintf('=============================================================\n\n');

    fprintf('STRUCTURAL INPUT\n');
    fprintf('----------------\n');
    E = ask_user('Enter column modulus E (Pa)', 30e9);
    D_col = ask_user('Enter square column side length a (m)', 1.0);
    h = ask_user('Enter effective height h (m)', 3.0);
    T_target = ask_user('Enter fixed-base target period T_fixed (s)', 0.10);
    F = ask_user('Enter lateral force F (N)', 1000);

    if E <= 0 || D_col <= 0 || h <= 0 || T_target <= 0 || F <= 0
        error('Invalid structural input.');
    end

    fprintf('\nFOUNDATION INPUT\n');
    fprintf('----------------\n');
    Bp = 0.305;
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

    fprintf('\nSOIL TYPE SELECTION\n');
    fprintf('-------------------\n');
    fprintf('  [1] Sand\n');
    fprintf('  [2] Clay\n');
    soil_type = ask_user('Enter soil type number [1 or 2]', 1);

    soil_type_name = '';
    soil_class_name = '';
    Ks_table_MN = 0;
    nu = 0;
    rho_t = 0;
    G_kPa = 0;

    if soil_type == 1
        soil_type_name = 'Sand';
        fprintf('\nSelect sand class:\n');
        fprintf('  [1] Loose sand\n');
        fprintf('  [2] Medium sand\n');
        fprintf('  [3] Dense sand\n');
        soil_class = ask_user('Enter sand class number', 2);

        if soil_class == 1
            soil_class_name = 'Loose sand'; Ks_table_MN = 12.9; nu = 0.275; rho_t = 1.8; G_kPa = 20000;
        elseif soil_class == 2
            soil_class_name = 'Medium sand'; Ks_table_MN = 41.7; nu = 0.325; rho_t = 1.9; G_kPa = 55000;
        elseif soil_class == 3
            soil_class_name = 'Dense sand'; Ks_table_MN = 161.0; nu = 0.350; rho_t = 2.0; G_kPa = 140000;
        else
            error('Invalid sand class.');
        end

    elseif soil_type == 2
        soil_type_name = 'Clay';
        fprintf('\nSelect clay class:\n');
        fprintf('  [1] Stiff clay\n');
        fprintf('  [2] Very stiff clay\n');
        fprintf('  [3] Hard clay\n');
        soil_class = ask_user('Enter clay class number', 2);

        if soil_class == 1
            soil_class_name = 'Stiff clay'; Ks_table_MN = 24.1; nu = 0.40; rho_t = 1.75; G_kPa = 25000;
        elseif soil_class == 2
            soil_class_name = 'Very stiff clay'; Ks_table_MN = 48.2; nu = 0.44; rho_t = 1.85; G_kPa = 60000;
        elseif soil_class == 3
            soil_class_name = 'Hard clay'; Ks_table_MN = 96.4; nu = 0.47; rho_t = 1.90; G_kPa = 115000;
        else
            error('Invalid clay class.');
        end
    else
        error('Invalid soil type.');
    end

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

    Ks    = Ks_MN * 1000;
    rho_s = rho_t * 1000;
    G     = G_kPa * 1000;
    Vs    = sqrt(G / rho_s);

    I = (D_col^4) / 12;
    k_struct = 3 * E * I / h^3;
    m = (T_target / (2*pi))^2 * k_struct;
    T_fixed = 2*pi*sqrt(m/k_struct);
    Delta_fixed = F / k_struct;
    omega_fixed = 2*pi / T_fixed;
    a0 = omega_fixed * B / (2 * Vs);

    if soil_type == 1
        n_delta = ((B + Bp) / (2 * B))^2;
        n_beta = 1 + 2 * (Df / B);
    else
        n_delta = Bp / B;
        n_beta = 1.0;
    end
    
    n_sigmax = (1/3) * (2 + B / L);
    Ks_tel = n_delta * n_sigmax * n_beta * Ks;

    Kz_kN  = Ks_tel * L * B;
    Krl_kN = Ks_tel * L * B^3 / 12;
    Krb_kN = Ks_tel * L^3 * B / 12;
    Kx_kN  = 0.30 * Kz_kN;

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

    Kx_dynamic  = Kx_dynamic_kN  * 1000;
    Krb_dynamic = Krb_dynamic_kN * 1000;
    kyy_ref = Krb_dynamic;

    T_ratio = sqrt(1 + k_struct/Kx_dynamic + (k_struct*h^2)/kyy_ref);
    T_flex_est = T_fixed * T_ratio;

    u_f = F / Kx_dynamic;
    theta = (F*h) / kyy_ref;
    Delta_flex_est = F/k_struct + u_f + theta*h;

    Tx_est  = 2*pi*sqrt(m/Kx_dynamic);
    Tyy_est = 2*pi*sqrt((m*h^2)/kyy_ref);

    [K_dyn, beta_rad] = gazetas_impedance_full(B, L, Df, G, nu, Vs, omega_fixed);
    Kx_rig  = K_dyn.x;
    Kyy_rig = K_dyn.ry;   
    beta_x  = beta_rad.x;
    beta_yy = beta_rad.ry; 

    n_s = 2; n_x = 2; n_yy = 2;
    
    ratio_T_NIST = sqrt(1 + k_struct/Kx_rig + (k_struct*h^2)/Kyy_rig);
    T_flex_NIST = T_fixed * ratio_T_NIST;
    
    Tx_NIST = 2*pi*sqrt(m/Kx_rig);
    Tyy_NIST = 2*pi*sqrt((m*h^2)/Kyy_rig);

    term_hysteretic = beta_s * ( (ratio_T_NIST^n_s - 1) / ratio_T_NIST^n_s );
    term_rad_x      = beta_x * (1 / (T_flex_NIST/Tx_NIST)^n_x);
    term_rad_yy     = beta_yy * (1 / (T_flex_NIST/Tyy_NIST)^n_yy);
    beta_f_NIST = term_hysteretic + term_rad_x + term_rad_yy;
    damping_NIST = beta_f_NIST + beta_i * (1 / ratio_T_NIST^n_s);
    
    delta_NIST = F/k_struct + (F / Kx_rig) + ((F*h) / Kyy_rig)*h;

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
    fprintf('T_ratio_NIST         = %.4f\n', ratio_T_NIST);
    fprintf('T_flex_NIST          = %.4f s\n', T_flex_NIST);
    fprintf('Delta_flex_NIST      = %.6e m\n', delta_NIST);
    fprintf('System Damping B_0   = %.4f\n', damping_NIST);
    fprintf('=============================================================\n');

end

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

    % Υπολογισμός απόσβεσης ακτινοβολίας
    beta_rad.x = (4 * (l/b + (D_emb/b)*(psi + l/b)) / (Kx_emb / (G*b))) * (a0 / (2*alpha_x));

    % Σπάσιμο της εξίσωσης σε μικρά κομμάτια για να μην κοπεί κατά την αντιγραφή
    term1a = (l/b)^3 * (D_emb/b);
    term1b = psi*(D_emb/b)*(l/b);
    term1c = (D_emb/b)^3;
    term1d = psi*(l/b)^3;
    term1e = 3*(D_emb/b)*(l/b)^2;
    term1f = psi*(l/b);
    
    term1_ry = (4/3) * (term1a + term1b + term1c + term1d + term1e + term1f) * (a0^2);
    term2_ry = (Kry_emb / (G*b^3)) * ( 1.8 / (1 + 1.75*(l/b - 1)) + a0^2 );
    term3_ry = (4/3) * (l/b + psi) * (D_emb/b)^3 / (Kry_emb / (G*b^3));
    
    beta_rad.ry = (term1_ry / term2_ry + term3_ry) * (a0 / (2*alpha_ry));
end