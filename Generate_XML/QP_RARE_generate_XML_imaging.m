clear all 
close all

% Settings
name = 'RARE_quadphase_schemefunction_180_7pre_scanner_blib_fullkspace.xml' ;% Set name for file saving
n_img = 18*2; % Imaging echoes (bijv. 36), two overscan lines 
n_prep = 7;   % Prepulses (3 of 7)
etl = n_prep + n_img;
use_variable_flip = false;
constant_flipangle = 180;
blib = true;  % Set this to true if you want to add a blib to show diffusion effects
excitation_initialphase = 0; % Set initial phase of 90 pulse

% Le Roux phases
if n_prep == 3
    [E_rad, R_rad] = phase_mod_Leroux3(etl);
else 
    [E_rad, R_rad] = phase_mod_Leroux_7prepulsescheme(etl);
end
E_deg = rad2deg(E_rad);
R_deg = rad2deg(R_rad);

if use_variable_flip
    % === Variable flip angle scheme ===
    % Prepulses: 160 deg
    prep_flips = 180 * ones(1, n_prep);

    N      = etl - n_prep;         % total number of refocusing pulses
    i      = (0:N-1)/(N-1);           % normalized echo index 0→1
        
    % a U-shaped parabola: 4*(i–0.5)^2 goes 1→0→1 as i goes 0→0.5→1
    w      = 4*(i - 0.5).^2;          
    
    % now map that onto [180°→constant_flipangle→180°]
    flips_v = constant_flipangle + (180 - constant_flipangle) * w;

    lips_v = deg2rad([prep_flips, flips_v]);
else
    lips_v = constant_flipangle * ones(1, etl);
end

% Open XML
fid = fopen(name, 'w');
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<Parameters FOVx="132" FOVy="132" FOVz="1" Name="P" GradMaxAmpl="0.08" GradSlewRate="0.2" Nx="32" Ny="32" Nz="1" TE="50">\n');
fprintf(fid, '   <ConcatSequence Name="QRARE">\n');
fprintf(fid, '      <ConcatSequence Name="O">\n');

% Excitation block
fprintf(fid, '         <ATOMICSEQUENCE Name="A1">\n');
fprintf(fid, '            <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="90" InitialPhase="%.10g" Name="P1"/>\n', excitation_initialphase);
fprintf(fid, '         </ATOMICSEQUENCE>\n');
if blib
    % GY blip gradient with area = DKy, phase variatie van 2pi en
    % rephasing
    fprintf(fid, '         <ATOMICSEQUENCE Name="GYBLIB">\n');
    fprintf(fid, '            <TRAPGRADPULSE Area="DKy" Axis="GY" Name="PBLIB" Observe="DKy=P.DKy"/>\n');
    fprintf(fid, '         </ATOMICSEQUENCE>\n');

    % Zero flip angle hard RF pulse
    fprintf(fid, '         <ATOMICSEQUENCE Name="BLIBZERO">\n');
    fprintf(fid, '            <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="0" Name="PBLIBRF"/>\n');
    fprintf(fid, '         </ATOMICSEQUENCE>\n');

end
fprintf(fid, '         <ATOMICSEQUENCE Name="A2">\n');
fprintf(fid, '            <TRAPGRADPULSE Area="0.5*A" Axis="GX" Name="P2" Observe="A=ROimg%d.Area"/>\n', 1+n_prep);
fprintf(fid, '         </ATOMICSEQUENCE>\n');
fprintf(fid, '         <DELAYATOMICSEQUENCE Delay="TE/2" DelayType="B2E" Name="D" Observe="TE=P.TE" StartSeq="A1"/>\n');

% Prepulses
for i = 1:n_prep
    fprintf(fid, '         <ConcatSequence Name="PREP%d">\n', i);
    % Crusher vóór RF
    fprintf(fid, '            <ATOMICSEQUENCE Name="Aprep1%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Area="0.27" Axis="GX" Name="Pprep1%d"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    % RF
    fprintf(fid, '            <ATOMICSEQUENCE Name="RFPprep%d">\n', i);
    fprintf(fid, '               <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="%.10g" InitialPhase="%.10g" Name="PprepRF%d" Refocusing="1"/>\n', ...
        lips_v(i), E_deg(i), i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    % Crusher na RF
    fprintf(fid, '            <ATOMICSEQUENCE Name="Aprep2%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Area="0.27" Axis="GX" Name="Pprep2%d"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    % Dummy PE & RO
    fprintf(fid, '            <ATOMICSEQUENCE Name="AprepPE%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Axis="GY" Name="PEprep%d" Observe="KMy=P.KMAXy, DKy=P.DKy, Ny=P.Ny"/>\n', i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AprepRO%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE ADCFlag="0" Axis="GX" InitialPhase="%.10g" Name="ROprep%d" FlatTopArea="2*KMx" FlatTopTime="20" Observe="KMx=P.KMAXx, Nx=P.Nx, TE=P.TE"/>\n', ...
        R_deg(i), i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <ATOMICSEQUENCE Name="AprepREPH%d">\n', i);
    fprintf(fid, '               <TRAPGRADPULSE Axis="GY" Name="REPHprep%d" Observe="A=PEprep%d.Area"/>\n', i, i);
    fprintf(fid, '            </ATOMICSEQUENCE>\n');
    fprintf(fid, '            <DELAYATOMICSEQUENCE Delay="TE" DelayType="B2E" Name="D%d" Observe="TE=P.TE" StartSeq="RFPprep%d"/>\n', i,i);
    fprintf(fid, '         </ConcatSequence>\n');
end

% Imaging echo train: per PE-lijn twee echoes (linear), met crushers etc.
Ny = 32;
Nx = 32;     

echo_nr = n_prep + 1; % echo-teller

for pe = 1:Ny%((Ny/2)+((n_img-Ny)/2))
    for rep = 1:2   % Twee keer per PE-lijn
        i = echo_nr;
        rf_phase = E_deg(i);
        adc_phase = R_deg(i);
        flip_angle = lips_v(i);
        % ---- Begin echo ----
        fprintf(fid, '         <ConcatSequence Name="ECHO%d">\n', i);
        
        % Crusher vóór RF
        fprintf(fid, '            <ATOMICSEQUENCE Name="CRUSH_PRE%d">\n', i);
        fprintf(fid, '               <TRAPGRADPULSE Area="0.27" Axis="GX" Name="CRUSHER_PRE%d"/>\n', i);
        fprintf(fid, '            </ATOMICSEQUENCE>\n');
        % RF
        fprintf(fid, '            <ATOMICSEQUENCE Name="RFPimg%d">\n', i);
        fprintf(fid, '               <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="%.10g" InitialPhase="%.10g" Name="PimgRF%d" Refocusing="1"/>\n', ...
            flip_angle, rf_phase, i);
        fprintf(fid, '            </ATOMICSEQUENCE>\n');
        % Crusher na RF
        fprintf(fid, '            <ATOMICSEQUENCE Name="CRUSH_POST%d">\n', i);
        fprintf(fid, '               <TRAPGRADPULSE Area="0.27" Axis="GX" Name="CRUSHER_POST%d"/>\n', i);
        fprintf(fid, '            </ATOMICSEQUENCE>\n');
        % PE-gradient: Area="(-KMy + DKy * (pe))"
        fprintf(fid, '            <ATOMICSEQUENCE Name="AimgPE%d">\n', i);
        fprintf(fid, '               <TRAPGRADPULSE Area="(-KMy + DKy*%d)" Axis="GY" Name="PEimg%d" Observe="KMy=P.KMAXy, DKy=P.DKy, Ny=P.Ny"/>\n', ...
            (pe), i);
        fprintf(fid, '            </ATOMICSEQUENCE>\n');
        % RO-gradient & ADC
        fprintf(fid, '            <ATOMICSEQUENCE Name="AimgRO%d">\n', i);
        fprintf(fid, '               <TRAPGRADPULSE ADCFlag="2" ADCs="%d" Axis="GX" FlatTopArea="2*KMx" FlatTopTime="20" InitialPhase="%.10g" Name="ROimg%d" Observe="KMx=P.KMAXx, Nx=P.Nx"/>\n', ...
            Nx, adc_phase, i);
        fprintf(fid, '            </ATOMICSEQUENCE>\n');
        % Rephaser (neutr. PE)
        fprintf(fid, '            <ATOMICSEQUENCE Name="AimgREPH%d">\n', i);
        fprintf(fid, '               <TRAPGRADPULSE Area="-A" Axis="GY" Name="REPHimg%d" Observe="A=PEimg%d.Area"/>\n', i, i);
        fprintf(fid, '            </ATOMICSEQUENCE>\n');
        % Delay
        fprintf(fid, '            <DELAYATOMICSEQUENCE Delay="TE" DelayType="B2E" Name="D%dE" Observe="TE=P.TE" StartSeq="RFPimg%d"/>\n', i, i);
        
        fprintf(fid, '         </ConcatSequence>\n');
        echo_nr = echo_nr + 1;
    end
end

% Close XML
fprintf(fid, '      </ConcatSequence>\n');
fprintf(fid, '   </ConcatSequence>\n');
fprintf(fid, '</Parameters>\n');
fclose(fid);

disp('✔ XML met volledige imaging encodings, crushers en Le Roux fases gegenereerd!');
