%% QPRARE Robustness Analysis with Sigma of Derivatives
clear all; close all; clc

%% Settings
T1 = 10; T2 = 1; esp = 0.05; 
blib = false;
prepulse_schemes = [3,7];
flip_angles_deg = [180 170 160 150 140 130 120 110 100];
etl_lengths = [21 41 61 81 101 121];
Nmc = 100;  % Monte Carlo samples for phase offset
phi0_samples = 2*pi*rand(1,Nmc);

%% === 1. Robustness vs φ₀ (Monte Carlo over phase offset, std of derivative) ===
flip0 = 160; etl0 = 60;
fprintf('\n=== Robustness vs φ₀ (FA=%d, ETL=%d): ===\n', flip0, etl0);
%fprintf('Prep | Mean σ(dS1) | σ(σ(dS1)) | Mean σ(d|S2|) | σ(σ(d|S2|)) | Mean σ(dI) | σ(σ(dI)) | Mean DecayRatio | σ(DecayRatio)\n');
for p = 1:length(prepulse_schemes)
    prep = prepulse_schemes(p);
    fa_train = deg2rad(flip0 * ones(1, etl0));
    dS1 = zeros(1, Nmc); dS2 = zeros(1, Nmc); dI = zeros(1, Nmc); decay_ratio = zeros(1, Nmc);
    for k = 1:Nmc
        phi0 = phi0_samples(k);
        [~, S1, S2, ~, ~] = epg_QPRARE2(fa_train, etl0, T1, T2, esp, prep, blib, phi0, false);
        idx = (prep+1):etl0;
        I = sqrt(S1(idx).^2 + S2(idx).^2);
        dS1(k) = std(diff(S1(idx)));
        dS2(k) = std(diff(abs(S2(idx))));
        dI(k)  = std(diff(I));
        decay_ratio(k) = I(end)/I(1);
    end
    % Print summary for this prepulse scheme:
    fprintf('Summary for all φ₀ for %d prepulses:\n', prep);
    fprintf('  dS1: mean=%.4g, sigma=%.4g\n', mean(dS1), std(dS1));
    fprintf(' d|S2|: mean=%.4g, sigma=%.4g\n', mean(dS2), std(dS2));
    fprintf('   dI: mean=%.4g, sigma=%.4g\n', mean(dI), std(dI));
    fprintf('Decay ratio: mean=%.4g, sigma=%.4g\n', mean(decay_ratio), std(decay_ratio));
end


%% === 2. Robustness vs Flip Angle (no phase MC, φ₀ = 0, std of derivative) ===
etl0 = 60;  % use the same ETL for flip angle comparison
fprintf('\n=== Robustness vs Flip Angle (ETL=%d, φ₀=0): ===\n', etl0);
%fprintf('Prep | FA | σ(dS1) | σ(d|S2|) | σ(dI) | DecayRatio\n');
for p = 1:length(prepulse_schemes)
    prep = prepulse_schemes(p);
    dS1_all = zeros(1,length(flip_angles_deg));
    dS2_all = zeros(1,length(flip_angles_deg));
    dI_all  = zeros(1,length(flip_angles_deg));
    decay_all = zeros(1,length(flip_angles_deg));
    for i = 1:length(flip_angles_deg)
        fa_deg = flip_angles_deg(i);
        fa_train = deg2rad(fa_deg * ones(1, etl0));
        phi0 = 0;
        [~, S1, S2, ~, ~] = epg_QPRARE2(fa_train, etl0, T1, T2, esp, prep, blib, phi0, false);
        idx = (prep+1):etl0;
        In = sqrt(S1(idx).^2 + S2(idx).^2);
        dS1_all(i)    = std(diff(S1(idx)));
        dS2_all(i)    = std(diff(abs(S2(idx))));
        dI_all(i)     = std(diff(In));
        decay_all(i) = In(end)/In(1);
        %fprintf(' %3d | %3d | %.4f | %.4f | %.4f | %.4f\n', prep, fa_deg, dS1_all(i), dS2_all(i), dI_all(i), decay_all(i));
    end
    % Print summary
    fprintf('Summary over all FAs for %d prepulses:\n', prep);
    fprintf('  dS1: mean=%.4g, sigma=%.4g\n', mean(dS1_all), std(dS1_all));
    fprintf(' d|S2|: mean=%.4g, sigma=%.4g\n', mean(dS2_all), std(dS2_all));
    fprintf('   dI: mean=%.4g, sigma=%.4g\n', mean(dI_all), std(dI_all));
    fprintf('Decay ratio: mean=%.4g, sigma=%.4g\n', mean(decay_all), std(decay_all));
end

%% === 3. Robustness vs ETL (no phase MC, φ₀ = 0, std of derivative) ===
flip0 = 160;
fprintf('\n=== Robustness vs ETL (FA=%d, φ₀=0): ===\n', flip0);
%fprintf('Prep | ETL | σ(dS1) | σ(d|S2|) | σ(dI) | DecayRatio\n');
for p = 1:length(prepulse_schemes)
    prep = prepulse_schemes(p);
    dS1_all = zeros(1,length(etl_lengths));
    dS2_all = zeros(1,length(etl_lengths));
    dI_all  = zeros(1,length(etl_lengths));
    decay_all = zeros(1,length(etl_lengths));
    for i = 1:length(etl_lengths)
        etl = etl_lengths(i);
        fa_train = deg2rad(flip0 * ones(1, etl));
        phi0 = 0;
        [~, S1, S2, ~, ~] = epg_QPRARE2(fa_train, etl, T1, T2, esp, prep, blib, phi0, false);
        idx = (prep+1):etl;
        In = sqrt(S1(idx).^2 + S2(idx).^2);
        dS1_all(i)    = std(diff(S1(idx)));
        dS2_all(i)    = std(diff(abs(S2(idx))));
        dI_all(i)     = std(diff(In));
        decay_all(i) = In(end)/In(1);
        %fprintf(' %3d | %3d | %.4f | %.4f | %.4f | %.4f\n', prep, etl, dS1_all(i), dS2_all(i), dI_all(i), decay_all(i));
    end
    % Print summary
    fprintf('Summary over all ETLs for %d prepulses:\n', prep);
    fprintf('  dS1: mean=%.4g, sigma=%.4g\n', mean(dS1_all), std(dS1_all));
    fprintf(' d|S2|: mean=%.4g, sigma=%.4g\n', mean(dS2_all), std(dS2_all));
    fprintf('   dI: mean=%.4g, sigma=%.4g\n', mean(dI_all), std(dI_all));
    fprintf('Decay ratio: mean=%.4g, sigma=%.4g\n', mean(decay_all), std(decay_all));
end
