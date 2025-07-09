% Example of optimizing the flipangle scheme for a Quadrature Phase RARE sequence 
clear, close all;

% Tissue and seq parameters:
% Several tissues:
% T1csf = 2000;
% T2csf = 250;
% T1gm = 1000;
% T2gm = 100;
% T1wm = 800;
% T2wm = 90;
% T1fat = 300;
% T2fat = 85;
T1phantom = 1000; % [ms]
T2phantom = 100; % [ms]
% T1brain = 900;
% T2brain = 95;

% Choose a tissue: 
T1choice = T1phantom;
T2choice = T2phantom;

% Acquisition settings
esp =  16;  % [ms]
L = 18; % Ammount of k lines sampled (when only part dont fill in full)
startups = 7;           % First prepulses to prepare eigenstate (Le Roux)
profileorder = 'linear';

partkspace = 18 / 32;   % Partial Fourier factor (if used)
FOV = 20;               % Field of view, only for defining target resolution
% adjust L number of K lines if partkspace
%L = round(L*partkspace);
ETL = 2*L; % echo train length, witouth prepulses, sampling each line twice     

% Total echoes for full k-space plus preparation pulses 
N = (ETL/partkspace)+startups; 
TotalRefocussingPulsesForFullGrid = N;

% save ID: 
ID = 'Phantom_LH_st7_PK_2';

%% Define target function (modified Hanning window)
% You want to design your window on an oversampled grid and then crop to
% get the right shape over your actual lines

% Full nx to build the filter on
Nx = (N - startups)/2; %Overshoot grid to build MTF
beta = 1;
alpha = 1.5;
dx = FOV / Nx;
dk = 1 / FOV;
Ksamples = ([1:Nx] - (Nx + 1) / 2) * dk; 
target = beta / 2 * (1 + cos((2 * pi * Ksamples * dx) / alpha)); % Hanning function


% Visualize target (MTF), and the resulting PSF: 
figure('position', [400, 400, 700, 250]);
subplot(1, 2, 1);
plot(target); ylim([0 1]);
title('Target MTF');

[PSF, PSFtarget, ax, axfine] = MTF2PSF(target, 'linear', FOV);
subplot(1, 2, 2);
plot(axfine, abs(PSFtarget));
title('Resulting PSF from Full Grid');
sgtitle(['Alpha = ' num2str(alpha)]);

[~,center] = max(target);

% Adjust for half-scan if needed

%–– Crop the target so it has exactly acquiredLines entries:
if L < Nx
   % Here’s the simplest “linear” cropping:
   target = target(1: L );

   % New ammount of pulses 
   N = ETL+startups; 
  
end

figure;
plot(target); ylim([0 1]);
title('Target Used During Acquisition');


%% Optimization using fmincon:
% Settings:
% lower and upper bounds for the flip angles:
if startups == 3
    minFlip = 140 * pi/180+ eps;    % 140° in radians, lowerbound for Quadrature Phase can notbe too low 
elseif startups == 7
    minFlip = 110 * pi/180+ eps;    % 110° in radians, lowerbound for Quadrature Phase can notbe too low 
end 

% USER OPTION: Fix or optimize prepulse flip angles
fix_prepulses = true;  % <<< Set to false if you want to optimize prepulses too

% Fixed angle for prepulses:
fixed_prep_angle = 180 * pi / 180;

[A,b,Aeq,beq]=deal([]); %Set extra constraints 

% Define optimizer options, maximum iterations, function evaluation max, tolarance for convergence
opt = optimset('MaxIter',1000,'MaxFunEvals',10^6,'TolFun',10^(-9),'AlwaysHonorConstraints','none'); %,'PlotFcns','optimplotfval'

if fix_prepulses
    % --------------------------------------------------------
    % Fix prepulses, optimize only remaining flip angles
    % --------------------------------------------------------

    % Variable part = ETL pulses (N - startups)
    x0_var = 100 * pi/180 * ones(1, N - startups);       % initial guess
    lb_var = minFlip * ones(1, N - startups);            % lower bound
    ub_var = pi * ones(1, N - startups);                 % upper bound

    % Compose full flip angle array with prepulses fixed
    x0 = [repmat(fixed_prep_angle, 1, startups), x0_var];
    lb = [repmat(fixed_prep_angle, 1, startups), lb_var];
    ub = [repmat(fixed_prep_angle, 1, startups), ub_var];

    % Cost function: inserts fixed prepulses automatically
    costfun = @(x_var) my_QPRARE2([repmat(fixed_prep_angle,1,startups), x_var], ...
                                  T1choice,T2choice,esp,N,target./sum(target),startups,profileorder,partkspace);

    % Run optimization only on the ETL part
    [x1_var, fval, exitflag, output]  = fmincon(costfun, x0_var, [], [], [], [], lb_var, ub_var, [], opt);

    % Final result: full flip angle vector including fixed prepulses
    x1 = [repmat(fixed_prep_angle, 1, startups), x1_var];

    disp('>> Prepulse flip angles fixed to 160° and not optimized.');

else
    % --------------------------------------------------------
    % Optimize all flip angles (prepulses + ETL)
    % --------------------------------------------------------

    x0 = 100 * pi/180 * ones(1, N);        % initial guess
    lb = minFlip * ones(1, N);             % lower bound
    ub = pi * ones(1, N);                  % upper bound

    % Cost function uses full x directly
    costfun = @(x) my_QPRARE2(x, T1choice, T2choice, esp, N, target./sum(target), startups, profileorder, partkspace);

    % Run full optimization
    [x1, fval, exitflag, output] = fmincon(costfun, x0, [], [], [], [], lb, ub, [], opt);

    disp('>> All flip angles (including prepulses) are optimized.');
end


%% Look at result
flipdeg = x1*180/pi; 

% The flip angles:
figure('position',[100 100 1100 300])
subplot(1,3,1)
hold on
plot(flipdeg)
xline(startups,'k--','Prep pulses','LabelHorizontalAlignment','left','DisplayName','Preperation pulses')
hold off
ylim([0 180])
xlim([0 N])
title('FAs'), xlabel('Echo # (plus prepulses)'), ylabel('Degrees')

% The MTF for S before filtering
[S,phasediag,P] = epg_QPRARE2(x1,N,T1choice,T2choice,esp, startups);

Signal     = S(startups+1:end);      % length = 2·L
e1  = Signal(1:2:end);           % echoes 1,3,5…  (L points)
e2  = Signal(2:2:end);           % echoes 2,4,6…  (L points)

% Compute the “sum” and “diff” images
Mx  = e2 + e1;                % in‐phase image
My  = e2 - e1;                % quadrature image

% Unfiltered MTF (magnitude‐sum): No start ups!
MTF_unfilt = abs(Mx).^2 + abs(My).^2;

% compute no‐filter SNR
SNR_nofilt = mean( abs(e1).^2 + abs(e2).^2 );  

subplot(1,3,2)
hold on
plot(1:L, abs(Mx), '-o','DisplayName','echo 1')
plot(1:L, abs(My), '-x','DisplayName','echo 2')
plot(1:L, MTF_unfilt, '--s','LineWidth',1.5,'DisplayName','combined')
hold off
xlabel('k‐line index'); ylabel('Magnitude (a.u.)')
title(sprintf('Raw echoes & MTF (SNR_{nofilt}=%.2f)', SNR_nofilt))
legend('Location','best')
ylim([0 max(MTF_unfilt)*1.1])

% Ordering of samples, if center-out (low-high) readout:
if strcmp(profileorder,'lowhigh') || strcmp(profileorder,'lh')
    order = zeros(numel(target),1);
    
    Nc = round((numel(target)+1)/2);
    order(1) = Nc;
    order(2:2:numel(target)) = [Nc-1:-1:1];
    order(3:2:numel(target)) = [Nc+1:numel(target)];
   
    
    if ~isempty(hsfactor) && hsfactor<1
        order = zeros(numel(target),1);
        [~,Nc] = max(target);
        
        order(1:2:Nc*2-1) = [Nc:-1:1];
        order(2:2:Nc*2-1)=[Nc+1:Nc*2-1];
        order(Nc*2:end)=[Nc*2:numel(target)];
    end
    
    target_order = target(order);
    [~,reverse_order] = sort(order); 
else
    target_order = target;
    reverse_order=[1:numel(target)]; 
end 

filt     = abs(target_order) ./ MTF_unfilt;
[~, Nc] = max(target);
filt_apply = filt(reverse_order); % the filter to apply to recorded data - ss the sampling order then is irrelevant.
filt_norm = filt_apply./filt_apply(Nc); % the normalized filter to apply to recorded data. 


% Correction filter & filtered MTF
% (assumes you've already computed & normalized `target_order`)
SNR_filt = 1 / sqrt( sum( filt.^2 ) );  

subplot(1,3,3)
plot(1:numel(filt), filt_norm,'-d','LineWidth',1.5)
xlabel('k‐line index'); ylabel('Filter gain')
title(sprintf('Filter (SNR_{filt}=%.2f)', SNR_filt))
ylim([0 max(filt_norm)*1.2])
xticks([1, round(L/2), L])
xticklabels({'start','k_{center}','end'})

% check-up: just to check that target function is obtained after filtering
%(compare with target_order): No startups since no k sampling
MTF_filt = MTF_unfilt .* filt;  
figure,
plot(1:L, target_order,'-','LineWidth',1.5,'DisplayName','target')
hold on
plot(1:L, MTF_filt,'--','LineWidth',1.5,'DisplayName','filtered MTF')
hold off
xlabel('k-line index'); ylabel('Normalized amplitude')
legend('Location','best')
title('Filtered MTF vs. Target')
ylim([0 1])
%% save result:
close all

% save all parameters: 
%save(['FAopt_' ID '.mat']) % saving all parameters used for optimization

% OR save filter and flips: 
filtername = 'QuadraturephaseRAREPK'; 
% saving filter and relevant parameters:
save(['FAfilt_' ID '.mat'],'filt','filt_apply','filt_norm','flipdeg','startups','partkspace','esp','target','target_order') 

% write flip angles to txt-file: 
fileID = fopen('FILE_FOR_FLIPS_QUADPHASE_PK_2.txt','w');
fprintf(fileID,'%f\r\n',flipdeg) % use '%d\n' if it should save/print integers. 
fclose(fileID)



