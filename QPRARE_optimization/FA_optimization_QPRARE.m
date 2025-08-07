function [flips_opt, txdev_opt, rxdev_opt] = FA_optimization_QPRARE6()
% FA_optimization_QPRARE5
% Robust variable flip-angle and phase optimization for Quadrature-Phase RARE (QP-RARE) sequences.
%
%   Author: Fleur de Haar, 2025
%   Inspired by Sofie Rahbek's SPLICE optimization, with multi-phase robustness extension.
%
%   This script optimizes the refocusing flip angle and phase schedule for a QP-RARE sequence
%   using an SNR-based cost function, following the approach of Rahbek et al. The cost function
%   is evaluated over a set of phase offsets, and the optimization maximizes the worst-case SNR.
%   Results, including filter flatness and SNR, are visualized for both optimized and reference (constant flip) schemes.

%% 1. Tissue and Sequence Parameters ----------------------------------------------------------
clearvars; close all;
% Define relaxation times (ms)
T1brain = 900;
T2brain = 95;

T1choice = T1brain;   % T1 [ms]
T2choice = T2brain;   % T2 [ms]

esp         = 3;      % Echo spacing [ms]
L           = 30;     % Number of k-space lines (phase encoding steps)
Nx          = 55;     % Number of points for PSF calculation
startups    = 7;      % Number of prepulses (3 or 7)
profileorder= 'linear'; % k-space ordering
FOV         = 20;     % Field of View [cm]
ETL         = 2*L;    % Echo train length
N           = round(ETL) + startups; % Total number of refocusing pulses
M           = N - startups;          % Imaging pulses under optimization

eps_val     = 0.01;   % Small epsilon to avoid division by zero

fprintf('Total pulses: %d (imaging = %d)\n', N, M);

%% 2. Build Target k-Space Weighting and PSF -------------------------------------------------
% Modified Hann window for k-space target
beta  = 1; alpha = 1.5;
dx    = FOV/Nx;
dk    = 1/FOV;
kgrid = ([1:Nx] - (Nx+1)/2)*dk;
target = (beta/2)*(1 + cos((2*pi*kgrid*dx)/alpha)); % Modified Hann window

[~, psf_target, ax, axfine] = MTF2PSF(target, 'linear', FOV);

if L < Nx
    target = target(1:L);
    N = ETL + startups;
    M = N - startups;
end

target = target(:); % Column vector
psf_target = psf_target / sum(psf_target);

figure;
subplot(1,2,1); plot(target); ylim([0 1]); title('Target T(k)');
subplot(1,2,2); plot(axfine, abs(psf_target)); title('Resulting PSF');
sgtitle(['Alpha = ' num2str(alpha)]);

% k-space sampling order
if strcmp(profileorder,'lowhigh') || strcmp(profileorder,'lh')
    Nc = round((numel(target)+1)/2);
    order(1)      = Nc;
    order(2:2:end)= Nc-1:-1:1;
    order(3:2:end)= Nc+1:numel(target);
    target_ord    = target(order);
else
    target_ord = target;
end
[~, Nc] = max(target);

%% 3. Optimization Constraints and Initialization ---------------------------------------------
if startups==3
    minFlip = deg2rad(1);
    [E_LR, R_LR] = phase_mod_Leroux3(N);
else
    minFlip = deg2rad(1);
    [E_LR, R_LR] = phase_mod_Leroux_7prepulsescheme(N);
end
maxFlip = pi + eps_val;

num_flip_poly = 6; % Polynomial order for flip angle trajectory
p_tx = 6;          % Polynomial order for Tx phase
p_rx = 6;          % Polynomial order for Rx phase

x_img = linspace(-1,1,M)';

% Indices for parameter vector
flips_prep_idx = 1 : startups;
flips_poly_idx = startups + (1:num_flip_poly);
tx_prep_idx    = startups + num_flip_poly + (1:startups);
tx_poly_idx    = startups + num_flip_poly + startups + (1:p_tx);
rx_prep_idx    = startups + num_flip_poly + startups + p_tx + (1:startups);
rx_poly_idx    = startups + num_flip_poly + 2*startups + p_tx + (1:p_rx);

poly_row = x_img(1).^( (numel(p_tx)-1) : -1 : 0 );

% Parameter vector: [prep_flips; flip_poly; prep_tx; tx_poly; prep_rx; rx_poly]
p0 = [ repmat(pi, startups, 1);                % flips_prep
       pi; zeros(num_flip_poly-1,1);           % flips_poly
       zeros(startups,1);                      % tx_prep
       zeros(p_tx,1);                          % tx_poly
       zeros(startups,1);                      % rx_prep
       zeros(p_rx,1) ];                        % rx_poly

Aeq = zeros(2, numel(p0)); beq = zeros(2, 1);
Aeq(1, tx_poly_idx) = poly_row; 
Aeq(2, rx_poly_idx) = poly_row; 
A = []; b = []; lb = []; ub = [];

% Parameter mapping functions
build_flips = @(p) [p(flips_prep_idx); polyval(p(flips_poly_idx), x_img)];
build_txdev = @(p) [p(tx_prep_idx);   polyval(p(tx_poly_idx), x_img)];
build_rxdev = @(p) [p(rx_prep_idx);   polyval(p(rx_poly_idx), x_img)];
build_tx    = @(p) E_LR(:) + build_txdev(p);
build_rx    = @(p) R_LR(:) + build_rxdev(p);

opts = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
    'MaxFunctionEvaluations',1e5,'OptimalityTolerance',1e-9);

%% 4. Run Optimization -----------------------------------------------------------------------
[p_opt, ~] = fmincon(@snr_obj, p0, A, b, Aeq, beq, lb, ub, @bound_constraints, opts);

%% 5. Post-process and Evaluate at Worst-case Phase Offset -----------------------------------
% Extract optimized parameters
flips_opt = build_flips(p_opt);
txdev_opt = build_tx(p_opt);
rxdev_opt = build_rx(p_opt);

phase_offsets = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4];
nPhase = numel(phase_offsets);

% Evaluate SNR for all phase offsets
SNR_vals_opt = zeros(1, nPhase);
Add_all_opt = cell(1, nPhase);
Sub_all_opt = cell(1, nPhase);
filt_add_norm_opt = cell(1, nPhase);
filt_sub_norm_opt = cell(1, nPhase);

for idx = 1:nPhase
    phi = phase_offsets(idx);
    [S, ~, ~] = epg_QPRARE2_override(flips_opt', N, T1choice, T2choice, esp, startups, phi, txdev_opt, rxdev_opt, false);
    sig = S(startups+1:end);
    S_even = sig(1:2:end);
    S_odd  = sig(2:2:end);
    Addkspace = S_even + S_odd;
    Subkspace = S_even - S_odd;
    Add_all_opt{idx} = Addkspace;
    Sub_all_opt{idx} = Subkspace;
    % SNR weighting
    E1 = sum(abs(Addkspace).^2); E2 = sum(abs(Subkspace).^2); Etot = E1 + E2;
    w1 = E1 / Etot; w2 = E2 / Etot;
    % Correction filters
    filt_add = abs(target_ord(:)) ./ (abs(Addkspace(:)) + eps_val);
    filt_sub = abs(target_ord(:)) ./ (abs(Subkspace(:)) + eps_val);
    filt_add_norm_opt{idx} = filt_add ./ filt_add(Nc);
    filt_sub_norm_opt{idx} = filt_sub ./ filt_sub(Nc);
    % SNR metric (higher = better)
    SNR_vals_opt(idx) = 1 / sqrt(w1 * sum(filt_add.^2) + w2 * sum(filt_sub.^2));
end

% Worst-case phase offset (minimum SNR)
[~, worst_idx] = min(SNR_vals_opt);
phi_worst = phase_offsets(worst_idx);

% Store worst-case optimized results
opt.flipdeg = rad2deg(flips_opt);
opt.txdeg = rad2deg(txdev_opt);
opt.rxdeg = rad2deg(rxdev_opt);
opt.startups = startups;
opt.Add = Add_all_opt{worst_idx};
opt.Sub = Sub_all_opt{worst_idx};
opt.filt_add_norm = filt_add_norm_opt{worst_idx};
opt.filt_sub_norm = filt_sub_norm_opt{worst_idx};
opt.SNR = SNR_vals_opt(worst_idx);

fprintf('\nWorst-case phase offset (optimized): phi = %.2f * pi\n', phi_worst/pi);
fprintf('SNR at worst-case phi (optimized): %.4f\n', opt.SNR);

%% 6. Reference (Constant Flip) Evaluation ---------------------------------------------------
flips_ref = [repmat(pi, startups, 1); repmat(pi, M, 1)];
tx_ref = E_LR(:); rx_ref = R_LR(:);

[Sref, ~, ~] = epg_QPRARE2_override(flips_ref', N, T1choice, T2choice, esp, startups, phi_worst, tx_ref, rx_ref, false);
sig_ref = Sref(startups+1:end);
S_even_ref = sig_ref(1:2:end); S_odd_ref = sig_ref(2:2:end);
Add_ref = S_even_ref + S_odd_ref;
Sub_ref = S_even_ref - S_odd_ref;
E1r = sum(abs(Add_ref).^2); E2r = sum(abs(Sub_ref).^2); Etot1r = E1r + E2r;
w1r = E1r / Etot1r; w2r = E2r / Etot1r;
filt_ref_add = abs(target_ord(:)) ./ (abs(Add_ref(:)) + eps_val);
filt_ref_sub = abs(target_ord(:)) ./ (abs(Sub_ref(:)) + eps_val);
filt_ref_add_norm = filt_ref_add ./ filt_ref_add(Nc);
filt_ref_sub_norm = filt_ref_sub ./ filt_ref_sub(Nc);
SNR_ref = 1 / sqrt(w1r * sum(filt_ref_add.^2) + w2r * sum(filt_ref_sub.^2));

ref.flipdeg = rad2deg(flips_ref);
ref.txdeg = rad2deg(tx_ref);
ref.rxdeg = rad2deg(rx_ref);
ref.startups = startups;
ref.Add = Add_ref;
ref.Sub = Sub_ref;
ref.filt_add_norm = filt_ref_add_norm;
ref.filt_sub_norm = filt_ref_sub_norm;
ref.SNR = SNR_ref;

fprintf('SNR at worst-case phi (reference): %.4f\n', ref.SNR);

%% 7. Results Summary Table -------------------------------------------------------------------
phase_deg = rad2deg(phase_offsets);

fprintf('\n=========================================================\n');
fprintf('Summary Table: SNR (higher = better)\n');
fprintf('---------------------------------------------------------\n');
fprintf(' Scheme      | Phase Offset | SNR (worst-case) | Filter Flatness (In-phase, Quadrature)\n');
fprintf('---------------------------------------------------------------------------------------\n');
fprintf('Optimized    | %6.1f deg   |     %.4f      |   %.4f, %.4f\n', ...
    phase_deg(worst_idx), opt.SNR, std(opt.filt_add_norm), std(opt.filt_sub_norm));
fprintf('Reference    | %6.1f deg   |     %.4f      |   %.4f, %.4f\n', ...
    phase_deg(worst_idx), ref.SNR, std(ref.filt_add_norm), std(ref.filt_sub_norm));
fprintf('---------------------------------------------------------------------------------------\n');

fprintf('\nSNR for all phase offsets (Optimized):\n');
for i=1:nPhase, fprintf('%.1f deg: %.4f  ', phase_deg(i), SNR_vals_opt(i)); end
fprintf('\nSNR for all phase offsets (Reference):\n');
% Calculate SNR for reference at all phase offsets
SNR_vals_ref = zeros(1, nPhase);
for idx = 1:nPhase
    phi = phase_offsets(idx);
    [Sref_tmp, ~, ~] = epg_QPRARE2_override(flips_ref', N, T1choice, T2choice, esp, startups, phi, tx_ref, rx_ref, false);
    sig_ref_tmp = Sref_tmp(startups+1:end);
    S_even_ref_tmp = sig_ref_tmp(1:2:end); S_odd_ref_tmp = sig_ref_tmp(2:2:end);
    Add_ref_tmp = S_even_ref_tmp + S_odd_ref_tmp;
    Sub_ref_tmp = S_even_ref_tmp - S_odd_ref_tmp;
    E1r_tmp = sum(abs(Add_ref_tmp).^2); E2r_tmp = sum(abs(Sub_ref_tmp).^2); Etot1r_tmp = E1r_tmp + E2r_tmp;
    w1r_tmp = E1r_tmp / Etot1r_tmp; w2r_tmp = E2r_tmp / Etot1r_tmp;
    filt_ref_add_tmp = abs(target_ord(:)) ./ (abs(Add_ref_tmp(:)) + eps_val);
    filt_ref_sub_tmp = abs(target_ord(:)) ./ (abs(Sub_ref_tmp(:)) + eps_val);
    SNR_vals_ref(idx) = 1 / sqrt(w1r_tmp * sum(filt_ref_add_tmp.^2) + w2r_tmp * sum(filt_ref_sub_tmp.^2));
    fprintf('%.1f deg: %.4f  ', phase_deg(idx), SNR_vals_ref(idx));
end
fprintf('\n');

%% 8. Plot Results ----------------------------------------------------------------------------
plot_qprare_comparison(opt, ref, target_ord);

%% ------------------- Nested Functions --------------------------

    function c = snr_obj(p)
        % Cost function: maximize worst-case SNR across phase offsets
        flips = build_flips(p);
        txdev = build_tx(p);
        rxdev = build_rx(p);
        phase_offsets = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4];

        SNR_vals = zeros(1, length(phase_offsets));
        for idx = 1:length(phase_offsets)
            phi = phase_offsets(idx);
            [S, ~, ~] = epg_QPRARE2_override(flips', N, T1choice, T2choice, esp, startups, phi, txdev, rxdev, false);
            sig = S(startups+1:end);
            S_even = sig(1:2:end);
            S_odd  = sig(2:2:end);

            Addkspace = S_even + S_odd;
            Subtractkspace = S_even - S_odd;

            E1 = sum(abs(Addkspace).^2);
            E2 = sum(abs(Subtractkspace).^2);
            Etot = E1 + E2;
            w1 = E1 / Etot; w2 = E2 / Etot;

            filt_mx = target ./ (abs(Addkspace(:)) + eps_val);
            filt_my = target ./ (abs(Subtractkspace(:)) + eps_val);

            SNR_vals(idx) = 1 / sqrt(w1 * sum(filt_mx.^2) + w2 * sum(filt_my.^2));
        end

        c = -min(SNR_vals); % Minimize negative worst-case SNR (maximizes SNR)
    end

    function [c, ceq] = bound_constraints(p)
        flips = build_flips(p);
        c_min = minFlip - flips;        % imaging α ≥ minFlip
        c_max = flips - maxFlip;        % all α ≤ maxFlip
        c = [c_min; c_max];
        ceq = [];
    end


end

    % =================== Supporting Plot Function ======================
    function plot_qprare_comparison(opt, ref, target_ord)
        % Visual comparison between optimized and reference QP-RARE results
    
        figure('Position', [100 100 1200 800]);
    
        % --- Row 1: Flip angle and phase schemes
        subplot(3,3,1);
        plot(opt.flipdeg,'-o','LineWidth',1.2); xline(opt.startups,'k--');
        title('Optimized flip angles','FontWeight','bold'); ylabel('Degrees');
    
        subplot(3,3,2);
        plot(opt.txdeg,'-','LineWidth',1.2); title('Excitation phase deviation (opt, °)','FontWeight','bold');
    
        subplot(3,3,3);
        plot(opt.rxdeg,'-','LineWidth',1.2); title('Receiver phase deviation (opt, °)','FontWeight','bold');
    
        % --- Row 2: k-space signal profiles
        subplot(3,3,4); hold on;
        plot(target_ord / max(target_ord), 'k--','LineWidth',1.5);
        plot(abs(opt.Add), 'ro-');
        plot(abs(opt.Sub), 'bx-');
        title('Optimized k-space weighting','FontWeight','bold');
        legend('Target','In-phase','Quadrature');
    
        subplot(3,3,5); hold on;
        plot(target_ord / max(target_ord), 'k--','LineWidth',1.5);
        plot(abs(ref.Add), 'ro-');
        plot(abs(ref.Sub), 'bx-');
        title('Reference k-space weighting','FontWeight','bold');
        legend('Target','In-phase','Quadrature');
    
        % --- Row 3: Normalized filter gains with SNR
        subplot(3,3,7); hold on;
        plot(opt.filt_add_norm,'ro-');
        plot(opt.filt_sub_norm,'bx-');
        ylim([0 max(opt.filt_sub_norm)*1.2])
        title(['\bf Opt filters, SNR ∝ ' num2str(opt.SNR, '%.4f')]);
        xlabel('k-line'); ylabel('Gain');
    
        subplot(3,3,8); hold on;
        plot(ref.filt_add_norm,'ro-');
        plot(ref.filt_sub_norm,'bx-');
        ylim([0 max(ref.filt_sub_norm)*1.2])
        title(['\bf Ref filters, SNR ∝ ' num2str(ref.SNR, '%.4f')]);
        xlabel('k-line'); ylabel('Gain');
    end 
