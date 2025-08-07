% === EPG vs JEMRIS signal comparison ===

clear all
close all

%% Settings
% NOTE: In JEMRIS these values are in ms
T1 = 1000;        % [ms]
T2 = 100;        % [ms]
esp = 3;   % [ms]
prepulses = 7; % Number of prepulses (e.g., 3 or 7)
imaging_pulses = 90;
etl = prepulses + imaging_pulses; % total pulses
blip = false; % Use blip true or false 
Phi0 = 0;
Phi90 = pi/2;

flip_const_deg = 140; % or whatever constant angle you want
flip_angles_deg = flip_const_deg * ones(1, etl);
flip_angles_rad = deg2rad(flip_angles_deg);
    
%% --- EPG simulation ---
doPlot = false;
[S_e1, S1_e1, S2_e1, ~, ~] = epg_QPRARE2(flip_angles_rad, etl, T1, T2, esp, prepulses, Phi0,blip,doPlot);

S_epg1 = (S_e1);
S1_EPG1 = S1_e1;      % Real/in-phase
S2_EPG1 = S2_e1;     % Imag/quadrature % or ohter way around when other offset

[S_e2, S1_e2, S2_e2, ~, ~] = epg_QPRARE2(flip_angles_rad, etl, T1, T2, esp, prepulses, Phi90,blip,doPlot);

S_epg2 = (S_e2);
S1_EPG2 = S1_e2;      % Real/in-phase
S2_EPG2 = S2_e2;     % Imag/quadrature % or ohter way around when other offset


N = length(S2_EPG1);

% For x-axis
echo_axis = 1:N;
x_prep = prepulses;

figure;
subplot(1,2,1);
hold on;
plot(echo_axis,S1_EPG1, 'x-', 'DisplayName', 'S1 EPG k space data');
plot(echo_axis, S2_EPG1, 'o-', 'DisplayName', 'S2 EPG K space data');
ylims = ylim;  % Get current y-axis limits
plot([x_prep x_prep], ylims, 'k--', 'DisplayName', 'Start up');
text(x_prep + 0.5, ylims(2)*.98, sprintf('%d prepulses', prepulses), ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Echo number');
ylabel('Signal');
title(['In and out of phase signal (all signal in out phase)' ]);
legend show;
grid on;
hold off;
subplot(1,2,2)
hold on;
plot(echo_axis,S1_EPG2, 'x-', 'DisplayName', 'S1 EPG k space data');
plot(echo_axis, S2_EPG2, 'o-', 'DisplayName', 'S2 EPG K space data');
ylims = ylim;  % Get current y-axis limits
plot([x_prep x_prep], ylims, 'k--', 'DisplayName', 'Start up');
text(x_prep + 0.5, ylims(2)*.98, sprintf('%d prepulses', prepulses), ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Echo number');
ylabel('Signal');
title(['In and out of phase signal (all signal in in-phase)' ]);
legend show;
grid on;
hold off;

%% Plot total signal
figure;
hold on;
plot(echo_axis, abs(S_epg1), 'DisplayName','Total Signal Sodd k space')
plot(echo_axis, abs(S_epg2), 'DisplayName','Total Signal Sodd k space')
ylims = ylim;  % Get current y-axis limits
plot([x_prep x_prep], ylims, 'k--', 'DisplayName', 'Start up');
text(x_prep + 0.5, ylims(2)*.98, sprintf('%d prepulses', prepulses), ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('Echo number');
ylabel('Signal');
title(['Total k space signal' ]);
legend show;
grid on;
hold off;

%% Seperate odd and even k space data
sig1   = S_epg1(prepulses+1:end);
S_even1 = abs(sig1(1:2:end));
S_odd1 = abs(sig1(2:2:end));

sig2   = S_epg2(prepulses+1:end);
S_even2 = abs(sig2(1:2:end));
S_odd2 = abs(sig2(2:2:end));


N = length(S_even1);
echo_axis = 1:N;

figure;
subplot(1,2,1);
hold on;
plot(echo_axis,S_even1, 'x-', 'DisplayName', 'Even k space data phi0');
plot(echo_axis, S_even2, 'o-', 'DisplayName', 'Even k space data phi90');
ylims = ylim;  % Get current y-axis limits
xlabel('Echo number');
ylabel('|Signal|');
title(['In and out of absolute phase signal (all signal in out phase)' ]);
legend show;
grid on;
hold off;
subplot(1,2,2)
hold on;
plot(echo_axis,S_odd1, 'x-', 'DisplayName', 'Odd k space data phi0');
plot(echo_axis, S_odd2, 'o-', 'DisplayName', 'Odd k space data phi90');
ylims = ylim;  % Get current y-axis limits
xlabel('Echo number');
ylabel('|Signal|');
title(['In and out of absolute phase signal (all signal in in-phase)' ]);
legend show;
grid on;
hold off;

%% Add and subtract for image reconstruction
S_even1 = (sig1(1:2:end));
S_odd1 = (sig1(2:2:end));

S_even2 = (sig2(1:2:end));
S_odd2 = (sig2(2:2:end));

Addkspace1 = S_even1 + S_odd1;
Subtractkspace1 = S_even1 - S_odd1;

Addkspace2 = S_even2 + S_odd2;
Subtractkspace2 = S_even2 - S_odd2;

N = length(S_even2);
echo_axis = 1:N;

figure;
subplot(1,2,1);
hold on;
plot(echo_axis, abs(Addkspace1), 'x-', 'DisplayName', 'Added k space data phi0');
plot(echo_axis, abs(Subtractkspace1), 'o-', 'DisplayName', 'Subtract k space data phi0');
ylims = ylim;  % Get current y-axis limits
xlabel('Echo number');
ylabel('Signal');
title(['K space data phi0' ]);
legend show;
grid on;
hold off;
subplot(1,2,2)
hold on;
plot(echo_axis,abs(Addkspace2), 'x-', 'DisplayName', 'Added k space data phi90');
plot(echo_axis,abs(Subtractkspace2), 'o-', 'DisplayName', 'Subtract k space data phi90');
ylims = ylim;  % Get current y-axis limits
xlabel('Echo number');
ylabel('Signal');
title(['K space data phi90' ]);
legend show;
grid on;
hold off;

%% Calculate the total signal power
power_I = sum(abs(Addkspace1).^2);
power_O = sum(abs(Subtractkspace1).^2);
power_total = power_I + power_O;

percent_I = 100 * power_I / power_total;
percent_O = 100 * power_O / power_total;

power_I2 = sum(abs(Addkspace2).^2);
power_O2 = sum(abs(Subtractkspace2).^2);
power_total2 = power_I2 + power_O2;

percent_I2 = 100 * power_I2 / power_total2;
percent_O2 = 100 * power_O2 / power_total2;

fprintf('In-phase (I) signal power phi0:  %.2f%%\n', percent_I);
fprintf('Out-of-phase (O) signal power phi0: %.2f%%\n', percent_O);
fprintf('In-phase (I) signal power phi0:  %.2f%%\n', percent_I2);
fprintf('Out-of-phase (O) signal power phi0: %.2f%%\n', percent_O2);



