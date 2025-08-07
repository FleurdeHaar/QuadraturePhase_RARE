%% compare_quad_POCS_vs_golden.m
% Compare multiple QP-RARE POCS reconstructions and a standard RARE reconstruction (with blip)
% to a 'golden standard' RARE image (fully sampled, no blip).
% This script evaluates reconstruction performance by calculating MAE, MSE, and background leakage.
% Author: Fleur De Haar, 2025

clear all; close all;

%% Files
% List with your QP-RARE datasets (to be reconstructed with POCS)
quad_files = { ...
    'QP_RARE_180train_7prepuls_blib.mat', ...
    'QP_RARE_VF_7prepuls_blib.mat',  ...
};

% Standard RARE with blip (to be reconstructed with its own normal method)
rare_blib_file = 'RARE_normal_blibx.mat';

% Golden standard (fully sampled standard RARE)
ref_file = 'RARE_normal.mat';

%% Parameters
Nx = 32;    % Image matrix size in x-direction

%% Load and reconstruct the reference image
ref_image = reconstruct_magnitude_image_standard(ref_file, Nx);

%% Prepare results table: 2 QP-RARE + 1 RARE_blib
all_files = [quad_files, {rare_blib_file}];
nF = numel(all_files);
results = table('Size',[nF 5], ...
    'VariableTypes', {'string','double','double','double','double'}, ...
    'VariableNames', {'File','MAE','MSE','Leakage','BetterThanThresh'} );

%% Loop over all files and compute metrics
for i = 1:nF
    fname = all_files{i};
    fprintf('\n=== Processing %s ===\n', fname);
    
    % Decide on reconstruction method
    if ismember(fname, quad_files)
        % QP-RARE: use POCS
        recon = reconstruct_magnitude_image_POCS(fname);
    elseif strcmp(fname, rare_blib_file)
        % Standard RARE with blip: use standard reconstruction
        recon = reconstruct_magnitude_image_standard(fname, Nx);
    else
        error('Unknown file type!');
    end
    
    % Ensure reconstruction and reference have the same size
    if ~isequal(size(recon), size(ref_image))
        error('Size mismatch between reference and reconstruction for %s', fname);
    end
    
    % Calculate metrics
    abs_diff   = abs(ref_image - recon);
    % Background leakage: Calculate mean signal outside the foreground mask
    thresh     = 0.1 * max(ref_image(:));
    mask_fg    = ref_image >= thresh;
    mask_bg    = ~mask_fg;
    leakage    = mean(recon(mask_bg));
    
    mae_fg = mean(abs_diff(mask_fg));
    mse_fg = mean(abs_diff(mask_fg).^2);
  
    
    % Fill results table
    results.File(i)           = fname;
    results.MAE(i)            = mae_fg;
    results.MSE(i)            = mse_fg;
    results.Leakage(i)        = leakage;
    
    % Console output
    fprintf('  MAE:      %.4e\n', mae_fg);
    fprintf('  MSE:      %.4e\n', mse_fg);
    fprintf('  Leakage:  %.4e\n', leakage);

    figure;
    imshow(recon,[]);
    hold on;
    contour(mask_fg, [0.5 0.5], 'r'); % Shows mask outline
    title('Reference image with foreground mask');
end

%% Show summary
disp('=== Summary of reconstructions ===');
disp(results);

% Helper function: standard magnitude image reconstruction from .mat file
function image = reconstruct_magnitude_image_standard(mat_file, Nx)
    data = load(mat_file);
    if isfield(data, 'M')
        Mx = data.M(:,1);
        My = data.M(:,2);
        I  = complex(Mx, My);
    else
        error('Field M not found in %s', mat_file);
    end
    Ny = length(I)/Nx;
    kspace = reshape(I, [Nx, Ny]);
    image = abs(ifftshift(ifft2(ifftshift(kspace))));
end



