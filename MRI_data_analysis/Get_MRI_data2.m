clear all;
clc

cd('C:\Users\fleur\OneDrive - TU Eindhoven\Documents\Tue\Master\Stage\Onderzoek\code')
filename = 'P13312.7'; %RARE with imaging 30 Ny lines, first read out used as callibration 


% Reconstruct 2D GRE magnitude image from data obtained with gre2d.tar
% load data from P−file (set CV2 = 1 to save P−file)
d = toppe.utils.loadpfile(filename);

[nfid,ncoil,nslice,necho,nview] = size(d); % necho = 1 by construction

d = permute(d, [1 5 3 2 4]); % [nfid nview nslice ncoil]

d = squeeze(d); % [nfid nview ncoils] = [256 (pislquant+256) ncoils]

% Option: Remove first scan since callibration
d = d(:,2:end,:); 

size(d)
for ic = 1:ncoil
    ims(:,:,ic) = fftshift(ifftn(fftshift(d(:,:,ic))));
end

%% Only first image 
I_1 = ims(:,:,1);

figure;
imagesc(abs(flipdim(flipdim(I1',1),2)))

%% Coil combination of images 
I = sqrt(sum(abs(ims).^2, ndims(d))); % root sum of squares coil combination
figure;
imagesc(flipdim(flipdim(I',1),2))
axis image off; 
colormap gray;
title('MRI RARE sequence ', 'FontSize', 14) % white title text, size 14