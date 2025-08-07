% === EPG vs JEMRIS signal comparison ===

clear all
close all

%% Settings
% NOTE: In JEMRIS these values are in ms
T1 = 100;        % [s]
T2 = 10;        % [s]
esp = 0.5;   % [s]
prepulses = 3; % Number of prepulses (e.g., 3 or 7)
imaging_pulses = 18*2;
etl = prepulses + imaging_pulses; % total pulses
blip = false; % Use blip true or false 
Phi0 = pi;

%% --- Flip angle selection ---
use_variable_flip = false;  % If variable flip angle set this to true 
constant_flipangle = 180; % Set the constant flip angle 

if use_variable_flip
    % === Variable flip angle scheme ===
    % Prepulses: 160 deg
    prep_flips = 180 * ones(1, prepulses);

    N      = etl - prepulses;         % total number of refocusing pulses
    i      = (0:N-1)/(N-1);           % normalized echo index 0→1
        
    % a U-shaped parabola: 4*(i–0.5)^2 goes 1→0→1 as i goes 0→0.5→1
    w      = 4*(i - 0.5).^2;          
    
    % now map that onto [180°→constant_flipangle→180°]
    flips_v = constant_flipangle + (180 - constant_flipangle) * w;

    flip_angles_rad = deg2rad([prep_flips, flips_v]);
    flip_type = 'Variable';
else
    % === Constant flip angle scheme ===
    flip_const_deg = constant_flipangle; % or whatever constant angle you want
    flip_angles_deg = flip_const_deg * ones(1, etl);
    flip_angles_rad = deg2rad(flip_angles_deg);
    flip_type = 'Constant';
end

%% --- EPG simulation ---
doPlot = true;
[S_e, S1_e, S2_e, phasediag, ~] = epg_QPRARE2(flip_angles_rad, etl, T1, T2, esp, prepulses, Phi0,blip,doPlot);

S_epg = (S_e);
inphase_EPG = S1_e;      % Real/in-phase
outphase_EPG = S2_e;     % Imag/quadrature

%% --- Load JEMRIS signal corresponding with EPG simulation ---
file_JEMRIS = 'QP_RARE_160train_7prepuls_K0_final.mat';
data_JEMRIS = load(file_JEMRIS);
Mx = data_JEMRIS.M(:,1);
My = data_JEMRIS.M(:,2);
S_JEMRIS = (Mx - 1i*My);
S1_JEMRIS = Mx;
S2_JEMRIS = My;


%% --- Make all signals same length if not ---
N = min(length(S1_JEMRIS), length(inphase_EPG));
idx_postprep = (prepulses+1):N;

% Clip ALL signals to length N
S1_JEMRIS   = S1_JEMRIS(1:N);
inphase_EPG = inphase_EPG(1:N);
S2_JEMRIS   = S2_JEMRIS(1:N);
outphase_EPG= outphase_EPG(1:N);
S_JEMRIS    = S_JEMRIS(1:N);
S_epg       = S_epg(1:N);

% For x-axis
echo_axis = 1:N;
x_prep = prepulses; % Unchanged

%% === Total signal plot ===
figure;
subplot(2,2,1);
hold on;
plot(echo_axis, abs(S_epg), 'o-', 'DisplayName', 'EPG');
plot(echo_axis, abs(S_JEMRIS), 'x--', 'DisplayName', 'JEMRIS');
ylims = ylim;  % Get current y-axis limits
plot([x_prep x_prep], ylims, 'k--', 'DisplayName', 'Start up');
text(x_prep + 0.5, ylims(2)*.98, sprintf('%d prepulses', prepulses), ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Echo number');
ylabel('|Signal|');
title(['Magnitude of Total signal (S1 - 1i*S2)' ]);
legend show;
grid on;
hold off;

% === In-phase signal plot (S1) ===
subplot(2,2,2);
hold on;
plot(echo_axis, (inphase_EPG), 'o-', 'DisplayName', 'EPG');
plot(echo_axis, (S1_JEMRIS), 'x--', 'DisplayName', 'JEMRIS');
ylims = ylim;
plot([x_prep x_prep], ylims, 'k--', 'DisplayName', 'Start up');
text(x_prep + 0.5, ylims(2)*.98, sprintf('%d prepulses', prepulses), ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Echo number');
ylabel('S1');
title(['In-Phase Component (S1) ']);
legend show;
grid on;
% Annotate: mean/std after prepulses
mean_S1 = mean(abs(inphase_EPG(idx_postprep)));
std_S1 = std(abs(inphase_EPG(idx_postprep)));

ylims = ylim;
xlims = xlim;

x_pos = xlims(1) + 0.95 * (xlims(2) - xlims(1));
y_pos = ylims(1) + 0.05 * (ylims(2) - ylims(1));

%text(x_pos,y_pos, sprintf('Mean = %.3g\nStd = %.3g', mean_S1, std_S1), ...
%    'BackgroundColor','w','EdgeColor','k','FontSize',10,'HorizontalAlignment','right', 'VerticalAlignment','bottom');
hold off;

% === Quadrature-phase signal plot (S2) ===
subplot(2,2,3);
hold on;
plot(echo_axis, (outphase_EPG), 'o-', 'DisplayName', 'EPG');
plot(echo_axis, (S2_JEMRIS), 'x--', 'DisplayName', 'JEMRIS');
ylims = ylim;
plot([x_prep x_prep], ylims, 'k--', 'DisplayName', 'Start up');
text(x_prep + 0.5, ylims(2)*.98, sprintf('%d prepulses', prepulses), ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Echo number');
ylabel('S2');
title(['Out-phase Component (S2)'  ]);
legend show;
grid on;
% Annotate: mean/std after prepulses
mean_S2 = mean(abs(outphase_EPG(idx_postprep)));
std_S2 = std(abs(outphase_EPG(idx_postprep)));

ylims = ylim;
xlims = xlim;

x_pos = xlims(1) + 0.95 * (xlims(2) - xlims(1));
y_pos = ylims(1) + 0.05 * (ylims(2) - ylims(1));

%text(x_pos, y_pos, sprintf('Mean = %.3g\nStd = %.3g', mean_S2, std_S2), ...
  %  'BackgroundColor','w','EdgeColor','k','FontSize',10,'HorizontalAlignment','right', 'VerticalAlignment','bottom');

hold off;



% === Ratio plot |S1| / |S2| ===
subplot(2,2,4);
hold on;
ratio_epg = abs(inphase_EPG) ./ abs(outphase_EPG) ;  % Avoid division by zero
ratio_jemris = abs(S1_JEMRIS) ./ abs(S2_JEMRIS);
plot(echo_axis, ratio_epg, 'o-', 'DisplayName', 'EPG');
plot(echo_axis, ratio_jemris, 'x--', 'DisplayName', 'JEMRIS');
ylims = ylim;
plot([x_prep x_prep], ylims, 'k--', 'DisplayName', 'Start up');
text(x_prep + 0.5, ylims(2)*.98, sprintf('%d prepulses', prepulses), ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Echo number');
ylabel('|S1| / |S2|');
title(['Ratio of In-Phase to Out-phase Signal ']);
legend show;
grid on;
hold off;

%% === Compute SD of ratio after prepulses ===
idx_postprep = (prepulses+1):N;

% Standard deviations
std_ratio_epg = std(ratio_epg(idx_postprep));
std_ratio_jemris = std(ratio_jemris(idx_postprep));

% Print to console
fprintf('   Std. Dev. of |S1/S2| ratio (after %d prepulses):\n', prepulses);
fprintf('   EPG    : %.5f\n', std_ratio_epg);
fprintf('   JEMRIS : %.5f\n', std_ratio_jemris);

%% Compute standard deviation for S2 after prepulses 
% You want to suppress this signal, so should be as low as possible 
% Standard deviations
std_ratio_epg_S2 = std(abs(outphase_EPG(idx_postprep)));
std_ratio_jemris_S2 = std(abs(S2_JEMRIS(idx_postprep)));

% Print to console
fprintf('   Std. Dev. of |S2| (after %d prepulses):\n', prepulses);
fprintf('   EPG    : %.5f\n', std_ratio_epg_S2);
fprintf('   JEMRIS : %.5f\n', std_ratio_jemris_S2);

