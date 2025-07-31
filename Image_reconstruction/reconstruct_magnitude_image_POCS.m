function [magnitude_image , data] = reconstruct_magnitude_image_POCS(file, plotstate)
% RUN_HOMODYNE_RECONSTRUCTION Perform homodyne reconstruction on MRI data.
%
%   magnitude_image = RUN_HOMODYNE_RECONSTRUCTION(file)
%
%   INPUTS:
%       file  - filename of the .mat file containing data.M (Nx*Ny_total x 2 + overscanlines)
%
%   OUTPUT:
%       magnitude_image - reconstructed combined magnitude image
    close all 


    if nargin<2
         plotstate= true;
    end
    %% Load the data
    data = load(file);
    Nx = 32;
    Ny = 32;
    
    % Extract magnetization
    Mx = data.M(:,1);
    My = data.M(:,2);
    I = complex(Mx, My);

    % Calculate basic image parameters
    total_samples = length(I);
    Ny_total = total_samples / Nx;
    overscan_lines = Ny_total - Ny;

    % Reshape: [Nx, Ny_total]
    k = reshape(I, [Nx, Ny_total]);

    %% Step 1: Apply 1D IFFT along kx + sum/diff to obtain S1 and S2
    Ny_half = Ny_total / 2;
    S1_line = zeros(Nx, Ny_half);
    S2_line = zeros(Nx, Ny_half);
    
    % Summing of signals is done in image space 
    for i = 1:Ny_half
        Sn  = k(:, 2*i - 1);
        Sn1 = k(:, 2*i);
        img_line_n  = ifftshift(ifft(ifftshift(Sn)));
        img_line_n1 = ifftshift(ifft(ifftshift(Sn1)));
        S1_line(:, i) = img_line_n1 + img_line_n;
        S2_line(:, i) = img_line_n1 - img_line_n;
    end

    %% Step 2: Transform S1 and S2 back to (half) k-space 
    S1_kspace = zeros(Nx, Ny_half);
    S2_kspace = zeros(Nx, Ny_half);

    for i = 1:Ny_half
        S1_kspace(:, i) = fftshift(fft(ifftshift(S1_line(:, i))));
        S2_kspace(:, i) = fftshift(fft(ifftshift(S2_line(:, i))));
    end

    S_12 = {S1_kspace, S2_kspace};

    %% Step 3: Homodyne reconstruction for both S1 and S2
    partial_fraction = Ny_half / Ny;
    hnover = (1 - partial_fraction) * Ny;
    images = cell(1, 2);

    for idx = 1:2
        S = S_12{idx};
        % set treshhold for iteration
        treshold_pocs = 0.001;
        
        % Expand to full k-space size
        data_full = zeros(Nx, Ny);
        data_full(:, 1:Ny_half) = S; % fill measured lower ky half

        data_pk = data_full;

        % Set part of data not sampled to zero, zero padding for initial guess
        data_pk(:,1 + Ny - hnover:end) = 0; % zero high ky not measured

        % Image of full zero padded k space data 
        im_init = fftshift(ifftn(fftshift(data_pk)));

        % Determine symmetric/centric data with low frequencies
        data_center = data_pk;
        data_center(:,1:hnover-1) = 0;

        % Set the low frequency data to image space to determine the phase
        im_ph = fftshift(ifftn(fftshift(data_center)));

        % Take only magnitude and do the phase correction on the initial guess 
        im_init = abs(im_init) .* exp(1i * angle(im_ph));

        % First guess in k space 
        tmp_k = fftshift(fftn(fftshift(im_init)));
        diff_im = treshold_pocs + 1;
        
        % Iterative to reconstruct part no
        while (abs(diff_im) > treshold_pocs)

            % Replace measured k space data
            tmp_k(:,1:Ny - hnover) = data_pk(:,1:Ny - hnover);

            % Inverse DFT, to image space 
            tmp_im = fftshift(ifftn(fftshift(tmp_k)));

            % take only magnitude and apply the phase term
            tmp_im = abs(tmp_im) .* exp(1i * angle(im_ph));

            % Back to k space
            tmp_k = fftshift(fftn(fftshift(tmp_im)));

            % Compare the reconstructed image 
            diff_im = sum(abs(tmp_im(:) - im_init(:)).^2);
           % fprintf('S%d: Difference is %f\n', idx, diff_im);

            diff_im = max(abs(tmp_im(:)-im_init(:)));
           %fprintf('max diff = %g\n', diff_im);
            im_init = tmp_im;
        end
        im_pocs = tmp_im;
    
        % Save result
        images{idx} = im_pocs;

        % Zero padded image
        %im_zeropadded = fftshift(ifftn(fftshift(data_pk)));

        % Original image for comparison
        im_original = fftshift(ifftn(fftshift(data_full)));

    end

    %% Step 4: Combine components
    % %  Apply a phase correction for the global phase 
    % % Compute global phase difference
    % phase_diff = angle(sum(sum(images{1} .* conj(images{2}))));
    % 
    % %Align S2 phase
    % images{2} = images{2} * exp(-1i * phase_diff);
    % NO NEED NOW


    % Get in and out of phase images
    m_x = images{1};
    m_y = images{2};
    
    
    % Combine the in and out of phase image by SOS
    magnitude_image = sqrt(abs(m_x).^2 + abs(m_y).^2);

    in_phase = imag(m_x); % In phase image, real part
    out_phase = imag(m_y); % Out of phase image, real part 

    %% Step 5: Plot combined results
    if plotstate
       % for idx = 1:2

        show_phase(in_phase, 'In-of-Phase Image S1 \phi');
        show_phase(out_phase, 'Out-Phase Image S2 \phi');
        show_mri(magnitude_image, 'Combined |S1|^2 + |S2|^2');
    end 
 
%         figure;
%         subplot(1,2,1); imagesc(abs(ifftshift(ifftn(fftshift(S_12{idx}))))); axis image off; title(sprintf('Orig. beeld S%d',idx));
%         %subplot(2,1,2); imagesc(abs(fftshift(ifftn(fftshift(data_full(:,1:Ny_half)))))); axis image off; title('Zero-padded');
%         subplot(1,2,2); imagesc(abs(images{idx})); axis image off; title(sprintf('POCS resultaat', images{idx}));
%         %subplot(2,2,4); imagesc(angle(fftshift(ifftn(fftshift(data_full(:,1:hnover)))))); axis image off; title('Phase low frequency');
%         colormap gray; colorbar;
    %end
    %figure; imagesc(magnitude_image); axis image off; title('Combined magnitude', file); colormap gray; colorbar;
end 
 function show_mri(img, ttl)
    % take magnitude
    im = abs(img);

    % threshold at 5th percentile → zero background
    bg_th = prctile(im(:), 5);
    im(im < bg_th) = 0;

    % display
    figure;
    imagesc(im);
    axis image off;
    colormap gray;

    % stretch full dynamic range [0 → max]
    caxis([0 max(im(:))]);

    % simple gray‐scale colorbar
    cb = colorbar;
    cb = colorbar('Ticks',[]);  % no tick marks or numbers
    cb.TickDirection = 'out';

    title(ttl, 'Interpreter','none');
 end

 function show_phase(phi, ttl)
    % 1) Zero out the noise floor (bottom 5%)
    th = prctile(abs(phi(:)), 5);
    im = phi;
    im(abs(im) < th) = 0;

    % 2) Stretch the remaining range to [0…1]
    mn = min(im(:));
    mx = max(im(:));
    im = (im - mn) / (mx - mn);

    % 3) Display as a pure gray‐scale image
    figure;
    imagesc(im);
    axis image off;
    colormap gray;

    % 4) No colorbar (just blank background + full gray map)
    % (If you really want a legend bar, you can add a patch legend manually.)

    % 5) Title the panel
    title(ttl, 'Interpreter','none');
end
 