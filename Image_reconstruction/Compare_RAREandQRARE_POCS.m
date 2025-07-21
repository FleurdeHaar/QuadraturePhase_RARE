%% compare_quad_POCS_vs_golden.m
% Vergelijk meerdere QP-RARE POCS-reconstructies met één 'gouden standaard' RARE-image

clear all; close all;

%% Bestanden
% Cell-array met je quadrature-RARE datasets
quad_files = { ...
    'QP_RARE_180train_3prepuls.mat', ...
    'QP_RARE_160train_3prepuls_lowerdens.mat', ...
    'QP_RARE_160train_3prepuls.mat',  ...
    'QP_RARE_160train_3prepuls_blib.mat'...
};

% Gouden standaard RARE-bestand
ref_file = 'RARE_normal.mat';

%% Parameters
Nx = 32;    % matrixgrootte in x

%% Laad en reconstrueer referentiebeeld
ref_data = load(ref_file);
Mx_ref     = ref_data.M(:,1);
My_ref     = ref_data.M(:,2);
I_ref      = complex(Mx_ref, My_ref);
Ny_ref     = length(I_ref)/Nx;
k_ref      = reshape(I_ref, [Nx, Ny_ref]);
ref_image  = abs(ifftshift(ifft2(ifftshift(k_ref))));

%% Voorbereiding resultaat-tabel
nQ = numel(quad_files);
results = table('Size',[nQ 5], ...
    'VariableTypes', {'string','double','double','double','double'}, ...
    'VariableNames', {'File','MAE_POCS','MSE_POCS','Leakage_POCS','BetterThanThresh'} );

%% Loop over quadrature-bestanden
for i = 1:nQ
    fname = quad_files{i};
    fprintf('\n=== Processing %s ===\n', fname);
    
    % Run POCS-reconstructie (gebruik je eigen functie)
    recon = reconstruct_magnitude_image_POCS(fname);
    
    % Zorg dat recon dezelfde grootte heeft als ref_image
    if ~isequal(size(recon), size(ref_image))
        error('Grootte mismatch tussen referentie en reconstructie voor %s', fname);
    end
    
    % MAE & MSE
    abs_diff       = abs(ref_image - recon);
    mae_pocs       = mean(abs_diff(:));
    mse_pocs       = mean((abs_diff(:)).^2);
    
    % Leakage (achtergrond)
    thresh         = 0.01 * max(ref_image(:));
    mask_fg        = ref_image >= thresh;
    mask_bg        = ~mask_fg;
    leakage_pocs   = mean(recon(mask_bg));
    
    % Bepaal of recon beter is dan een simpele threshold-reconstructie
    better_flag    = mae_pocs < (mean(ref_image(:))*0.05);  % vb: 5% error
    
    % Vul resultaat-tabel
    results.File(i)            = fname;
    results.MAE_POCS(i)        = mae_pocs;
    results.MSE_POCS(i)        = mse_pocs;
    results.Leakage_POCS(i)    = leakage_pocs;
    results.BetterThanThresh(i)= better_flag;
    
    % Console output
    fprintf('  MAE (POCS): %.4e\n', mae_pocs);
    fprintf('  MSE (POCS): %.4e\n', mse_pocs);
    fprintf('  Leakage:    %.4e\n', leakage_pocs);
    fprintf('  <5%% MAE?    %d\n', better_flag);
end

%% Toon samenvatting
disp('=== Samenvatting POCS vergelijking ===');
disp(results);
