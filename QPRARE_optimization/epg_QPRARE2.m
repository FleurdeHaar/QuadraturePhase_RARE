function [S, S1, S2, phasediag, P] = epg_QPRARE2(    flipangle, etl, T1, T2, esp, prepulses, phi0, applyBlip, doPlot)
% THIS ONE WORKS!
%   INPUTS:
%     flipangle  —  either
%                   • a scalar (in radians), meaning “use the same 180° flip for 
%                     all refocusing pulses,” 
%                   • or a 1×ETL vector of radian‐values, giving each RF pulse’s 
%                     flip angle individually.  
%
%     etl        —  total number of RF pulses, including both “prepulses” 
%                   (n_prep = 3)  and  the 2×(number_of_k‐lines) readouts.  
%                   By convention here we choose n_prep=3, so:
%                       number_of_k_lines  =  (etl – 3)/2,  which must be integer.  
%
%     T1, T2     —  relaxation times (in arbritrairy units). 
%     esp        —  echo spacing (in arbritrairy units).  
%
%   OUTPUTS:
%     S        — a 1x(ETL) vector containing both S1 and S2
%     S1       —  a 1×(number_of_k_lines) real vector containing the "in phase echoes" 
%                 (CPMG‐family) signals after phase rotation by R.  
%     S2       —  a 1×(number_of_k_lines) real vector containing the quadrature phase component 
%                 (non-CPMG family) signals after phase rotation by R.  
%
%     phasediag  —  a real‐valued matrix of size [4·ETL  ×  ETL], obtained by stacking 
%                   the EPG states  [F+(−N..+N),  F−(−N..+N),  Z(−N..+N)]  before 
%                   relaxation.  The first 2·(2·ETL) rows are [flipud(F+); F−], the next 
%                   (2·ETL) rows are Z.  This matches your plotting format.  
%
%     P  —  the final 3×(2·ETL) EPG matrix after all ETL pulses, if you wish to inspect 
%            it.  Typically unused.  
%
%   --------------------------------------------------------------------------
%   USAGE EXAMPLE:
%       % 180° constant pulses, ETL=39, brain‐like T1/T2, 5 ms echo spacing:
%       [S1, S2, phasediag, P] = epg_qprare_leRoux(pi, 39, 0.9, 0.095, 0.005)
%
%   --------------------------------------------------------------------------
%   REFERENCE:
%     “Non‐CPMG Fast Spin‐Echo with Full Signal”  Patrick Le Roux,  
%     Journal of Magnetic Resonance 155 (2002) 278–292.  
%     Equations (22–25) define E(n) & R(n) for n = 1..ETL. 
%   --------------------------------------------------------------------------
%
%   Author:  Fleur de Haar, adapted from Sofie Rahbekls work and Le Roux (2002).  
%   Date:    June 2025.  
%

% Inputs and defaults
if nargin<1 || isempty(flipangle)
    flipangle = pi * ones(1, etl);   % “180°” for all refocusing pulses by default
end
if (nargin < 2 || length(etl)==0) etl = length(flipangle); end;

% If user gave a vector smaller than ETL, replicate the last value:
if (etl > length(flipangle)) flipangle(end+1:etl) = flipangle(end); end;

if nargin < 3 || isempty(T1),  T1=1.0;   end   % seconds
if nargin < 4 || isempty(T2),  T2=0.1;   end   % seconds
if nargin < 5 || isempty(esp), esp=0.01; end   % seconds
if nargin < 7 || isempty(phi0) , phi0 = 0 ; end 
if nargin<8, applyBlip = false; end
if nargin<9, doPlot    = false; end
 
% Choose the phase calculation function based on prepulses
if prepulses == 3
    [E, R] = phase_mod_Leroux3(etl);           % Your 3-prepulse Le Roux function
elseif prepulses == 7  
    [E, R] = phase_mod_Leroux_7prepulsescheme(etl);      % Your 7-prepulse version 
end

% Prepare EPG
P      = zeros(3,2*etl); P(3,1) = 1;
Pstore = zeros(4*etl, etl);
Zstore = zeros(2*etl, etl);
S      = zeros(1,etl);
S1     = zeros(1,etl);
S2     = zeros(1,etl);
offset = 0;

% Compute gradient‐blip profile if requested
if applyBlip
    Nlines    = (etl - prepulses)/2;
    totalSweep = 2*pi;            % 2π across FOV
    deltaPhi   = totalSweep / Nlines;
    Nstates    = size(P,2);
    idx0       = etl + 1;
    orders     = (1:Nstates) - idx0;
    phi_blip   = orders * deltaPhi;
end

% Excitation pulse (90°) about axis phi0
% Magnetisatie na puls langs +y (phi0 = 0) 
P    = epg_rf(P, pi/2, phi0); 

% apply gradient blip on coherence orders
if applyBlip
    P(1,:) = P(1,:) .* exp(-1i * phi_blip);
    P(2,:) = P(2,:) .* exp( 1i * phi_blip);
end

% === Echo train ===
for ech = 1:etl

    % == Refocussing pulses == 
    % Use last three values to give Diffusion effects to the relaxation
    % If diffusion effects to be wanted, then D should be put to 1.

    P = epg_grelax(P,T1,T2,esp/2,1,0,1,1);     % relaxation with crusher
    P = epg_rf(P, abs(flipangle(ech)), angle(flipangle(ech))+E(ech)); % Refocus pulse with extra phase E
    P = epg_grelax(P,T1,T2,esp/2,1,0,1,1);      % Relaxation with crusher
     
    rawF0 = P(1,1);
    rawF01 = rawF0 * exp(-1i * (R(ech)+offset));  % apply R(ech), quadratic receiver phase

    S(ech)  = rawF01;        % store complex-valued echo= F0 state, the echo
    S1(ech) = real(rawF01);  % real component, in phase echo
    S2(ech) = imag(rawF01);  % imaginary component, out/quadrature phase component
    
    Pstore(2*etl:4*etl-1,ech) = P(2,:).';	% Put in negative states
    Pstore(1:2*etl,ech) = flipud(P(1,:).');  % Put in positive, overwrite center.
    Zstore(:,ech) = P(3,:).';

end

plotstate = cat(1, Pstore, Zstore);

% Only plot if requested
if doPlot
    set(0,'defaultAxesFontSize',14);
    set(0,'DefaultLineLineWidth',2);

    % subplot 1: total signal
    subplot(2,2,1);
    hold on;
    plot([1:etl]*esp, abs(S));
    x_prep = prepulses*esp;
    ylims = ylim;
    plot([x_prep x_prep], ylims, 'k--', 'LineWidth', 1);
    text(x_prep, ylims(2)*.98, 'Prepulses', 'HorizontalAlignment','left','VerticalAlignment','top');
    xlabel('Echo Time'); ylabel('|Total signal| per echo');
    title('Quadrature Phase RARE: total signal');
    grid on;

    % subplot 2: S1 and S2
    subplot(2,2,2);
    hold on;
    plot([1:etl]*esp, abs(S1), '-o');
    plot([1:etl]*esp, abs(S2), '-x');
    plot([x_prep x_prep], ylim, 'k--');
    text(x_prep, ylims(2)*.98, 'Prepulses', 'HorizontalAlignment','left','VerticalAlignment','top');
    hold off;
    xlabel('Echo Time');
    ylabel('Signal S1 (real) and S2 (imag)');
    title('Quadrature Phase RARE: Signal Components');
    legend('S1','S2','Location','best');
    grid on;

    % subplot 3: state evolution
    subplot(2,2,3);
    dispim(plotstate);
    xlabel('Echo'); ylabel('F (top) and Z (bottom) States');
    title('F and Z states vs time');

    % subplot 4: E and R phases
    subplot(2,2,4);
    hold on;
    plot(1:etl, rad2deg(E), 'o-');
    plot(1:etl, rad2deg(R), 'x-');
    plot([prepulses prepulses ], ylim, 'k--');
    text(x_prep, ylims(2)*.98, 'Prepulses', 'HorizontalAlignment','left','VerticalAlignment','top');
    hold off;
    xlabel('RF Pulse # (n)');
    ylabel('Phase (°)');
    title('Le Roux Quadratic‐Phase Schedule');
    legend('E(n) Tx','R(n) Rx','Location','best');
    grid on;

    % Get magnitude
    S_abs  = abs(S);
    S1_abs = abs(S1);
    S2_abs = abs(S2);
    
    % Fluctuating Signal Index (relative to amplitude)
    FSI_S  = std(S_abs)  / mean(S_abs);
    FSI_S1 = std(S1_abs) / mean(S1_abs);
    FSI_S2 = std(S2_abs) / mean(S2_abs);
    
    % Optionally show
    fprintf('FSI Total Signal (|S|) : %.4f\n', FSI_S);
    fprintf('FSI S1 (|S1|, CPMG-like): %.4f\n', FSI_S1);
    fprintf('FSI S2 (|S2|, QP-like)  : %.4f\n', FSI_S2);

end

phasediag = plotstate;	


