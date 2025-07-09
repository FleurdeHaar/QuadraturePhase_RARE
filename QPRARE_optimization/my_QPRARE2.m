function [c, ceq, GC, GCeq] = my_QPRARE2(x,T1,T2,esp,ETL,target,sup,profileorder,hs)
% function for optimizing the flip angles in a QUADRATURE PHASE RARE sequence. 
% INPUT:
% x: the optimization variable (vector of flip angles)
% T1,T2,esp,ETL: tissue and sequence parameters relevant in the optim.
% target: The desired (normalized) k-space signal weighting
% sup: number of start-up-echoes
% profileorder: determining the k-space sampling order. Either 'linear' (default) or 'lowhigh' (same as center-out). 
% hs: the halfscan factor if halfscan, i.e. partial fourier imaging, is
% used. If left empty, halfscan is NOT used. 
% OBS! when using halfscan together with 'lowhigh', it is assumed that the 
% highest value of the target function (max(target)) applies to the center of k-space.  
% 
% Note: a single filter for the sum of the two echo-families is calculated.
% I.e. not taking oscillations in the individual families into account, as these in reality are not present to
% the extend they appaer when simulating with EPG. (Because of imperfect
% sliceprofile effect). 
% 
% By Sofie Rahbek, June 2022
% QUADRATURE PHASES calculated from paper Le Roux (2002)
% Adjusted for Quadrature Phase RARE by Fleur de Haar, 2025
% Adjusted with different signal for optimization (only CPMG component)
% -------------------------------------------------------------------------
% default values, if no input: 
if nargin > 9
    hs = [];
    if nargin < 8
        profileorder = 'linear'; %default profile order is linear ;
        if nargin < 7
            sup=3; 
        end
    end
end



% Normalizing target (if user input is not normalized):
target=target./sum(target);

% reordering target according to sampling order. E.g. if center-put (low-high) sampling
% is used, the center of the target fct is sampled first. 
N = ETL-sup;
if strcmp(profileorder,'lowhigh') || strcmp(profileorder,'lh')
    order = zeros(numel(target),1);
    Nc = round((N+1)/2);
    order(1) = Nc;
    order(2:2:N) = [Nc-1:-1:1];
    order(3:2:N) = [Nc+1:N];
    
    
    if ~isempty(hs) && hs<1
        order = zeros(numel(target),1);
        [~,Nc] = max(target);
        order(1:2:Nc*2-1) = [Nc:-1:1];
        order(2:2:Nc*2-1)=[Nc+1:Nc*2-1];
        order(Nc*2:end)=[Nc*2:numel(target)];
    end
    
    target = target(order);

end

% With S total signal
[S,~,~] = epg_QPRARE2(x,ETL,T1,T2,esp, sup);


% Remove unstable start up echeos
Signal = S(sup+1:end);

% Summarize two echoes in a row to find combined signal for one k line
% echo_sig = [s1 s2 s3 … s64]; Keep complex information, second echo 90 out
% of phase
e1 = Signal(1:2:end);    % echo i
e2 = Signal(2:2:end);    % echo i+1

Mx = e2 + e1;            % in‐phase combined echo
My = e2 - e1;            % quadrature combined echo

% === MTF filter on combined magnitude image ===
I_k = abs(Mx).^2 +abs(My).^2; 
% Using SOS since signals do not contain same ammount of signal 
filt = target ./ I_k;
SNR_mx = 1 / sum(filt.^2);  % High if Mx matches target

% Penalize fluctuation in the signal
std_Mx = std(abs(Mx));      % low if Mx is stable
ratioS2S1 = abs(My)./abs(Mx); 
std_ratio = std(ratioS2S1); % low if S2 behaves like S1

alpha = 0.5;
beta = 1.5;
gamma = 0.5;

% Cost: maximize SNR, minimize instabilities
c = -(alpha * SNR_mx - beta * std_Mx - gamma * std_ratio);
%c = -SNR_mx;
ceq = [];

%self-defined gradients: (could be implemented to improve computation time)
GC = []; 
GCeq = [];

end