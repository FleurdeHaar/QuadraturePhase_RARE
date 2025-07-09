% Makes XML file for JEMRIS with QP RARE sequence following original
% formula of le Roux with self defined formula for higher pulses
% Can add you own flip angle scheme
% No phase and frequency encoding --> sampling k = 0
% Crushers are in the sequence
% 

clear all 
close all
% Define echo train length
n_img = 18*2; %Two overscan lines 
n_prep = 7; % Ammount of prepulses defines the scheme 
etl = n_prep + n_img;
use_variable_flip = false; % <-- set to false for constant flip angle
constant_flipangle = 160;
blib = false;  % Set this to true if you want to add a blib to show diffusion effects
%%
% Set initial phase for excitation pulse
excitation_initialphase = 0;

% Choose the phase calculation function based on prepulses
if n_prep == 3
    [E_rad, R_rad] = phase_mod_Leroux3(etl);           % Your 3-prepulse Le Roux function
else 
    [E_rad, R_rad] = phase_mod_Leroux_7prepulsescheme(etl);      % Your 7-prepulse version 
end

E_deg = rad2deg(E_rad); 
R_deg = rad2deg(R_rad);

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

    lips_v = deg2rad([prep_flips, flips_v]);
else
    % --- Constant flip angles ---
    lips_v = constant_flipangle * ones(1, etl); % or set any other constant angle here
end


%% Make XML file 
% Open file for writing
fid = fopen('RARE_quadphase_schemefunction_160_3pre.xml', 'w'); % Name the file 
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<Parameters FOVx="132" FOVy="132" FOVz="1" Name="P" GradMaxAmpl="0.08" GradSlewRate="0.2" Nx="32" Ny="32" Nz="1" TE="50">\n');
fprintf(fid, '   <ConcatSequence Name="QRARE">\n');
fprintf(fid, '      <ConcatSequence Name="O">\n');


% Excitation block 
fprintf(fid, '         <ATOMICSEQUENCE Name="A1">\n');
fprintf(fid, '            <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="90" InitialPhase="%.10g" Name="P1"/>\n', excitation_initialphase);
fprintf(fid, '         </ATOMICSEQUENCE>\n');
if blib
    % GY blip gradient with area = DKy
    fprintf(fid, '         <ATOMICSEQUENCE Name="GYBLIB">\n');
    fprintf(fid, '            <TRAPGRADPULSE Area="DKy" Axis="GY" Name="PBLIB" Observe="DKy=P.DKy"/>\n');
    fprintf(fid, '         </ATOMICSEQUENCE>\n');

    % Zero flip angle hard RF pulse
    fprintf(fid, '         <ATOMICSEQUENCE Name="BLIBZERO">\n');
    fprintf(fid, '            <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="0" Name="PBLIBRF"/>\n');
    fprintf(fid, '         </ATOMICSEQUENCE>\n');
end
fprintf(fid, '         <ATOMICSEQUENCE Name="A2">\n');
fprintf(fid, '            <TRAPGRADPULSE Area="0" Axis="GX" Name="P2"/>\n');
fprintf(fid, '         </ATOMICSEQUENCE>\n');
fprintf(fid, '         <DELAYATOMICSEQUENCE Delay="TE/2" DelayType="B2E" Name="D" Observe="TE=P.TE" StartSeq="A1"/>\n');

% Prepulse blocks 
for i = 1:n_prep
    fprintf(fid, '         <ConcatSequence Name="PREP%d">\n', i);
    fprintf(fid, '            <ATOMICSEQUENCE Name="Aprep1%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE ADCFlag="1" Area="0.27" Axis="GX" Name="Pprep1%d"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="RFPprep%d">\n', i);
    fprintf(fid, '               <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="%.10g" InitialPhase="%.10g" Name="PprepRF%d" Refocusing="1"/>\n', ...
        lips_v(i), E_deg(i), i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="Aprep2%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE ADCFlag="1" Area="0.27" Axis="GX" Name="Pprep2%d"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AprepPE%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Area="0*E" Axis="GY" Name="PEprep%d" Observe="KMy=P.KMAXy, DKy=P.DKy, Ny=P.Ny"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AprepRO%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE ADCFlag="2" ADCs="1" Axis="GX" InitialDelay="TE" InitialPhase="%.10g" Name="ROprep%d" Observe="TE=P.TE"/>\n', ...
        R_deg(i), i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AprepREPH%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Area="0*E" Axis="GY" Name="REPHprep%d" Observe="A=PEprep%d.Area,E=P.TD"/>\n', i, i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '         </ConcatSequence>\n');
end


% Imaging echo train, unique names, no loops or underscores
for i = n_prep+1:etl
    fprintf(fid, '         <ConcatSequence Name="ECHO%d">\n', i);
    fprintf(fid, '            <ATOMICSEQUENCE Name="Aimg1%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE ADCFlag="1" Area="0.27" Axis="GX" Name="Pimg1%d"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="RFPimg%d">\n', i);
    fprintf(fid, '               <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="%.10g" InitialPhase="%.10g" Name="PimgRF%d" Refocusing="1"/>\n', ...
        lips_v(i), E_deg(i), i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="Aimg2%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE ADCFlag="1" Area="0.27" Axis="GX" Name="Pimg2%d"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AimgPE%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Area="0*E" Axis="GY" Name="PEimg%d" Observe="KMy=P.KMAXy, DKy=P.DKy, Ny=P.Ny"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AimgRO%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE ADCFlag="2" ADCs="1" Axis="GX" InitialDelay="TE" InitialPhase="%.10g" Name="ROimg%d" Observe="TE=P.TE"/>\n', ...
        R_deg(i), i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AimgREPH%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Area="0*E" Axis="GY" Name="REPHimg%d" Observe="A=PEimg%d.Area,E=P.TD"/>\n', i, i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '         </ConcatSequence>\n');
end

% Close everything
fprintf(fid, '      </ConcatSequence>\n');
fprintf(fid, '   </ConcatSequence>\n');
fprintf(fid, '</Parameters>\n');
fclose(fid);

disp('✔ XML with variable flip angles, Le Roux phases, no underscores, and linear echo train generated!');

