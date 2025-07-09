function [E_rad, R_rad] = phase_mod_Leroux_7prepulsescheme(N)
% Le Roux 7-prepulse quadratic phase scheme, with correct sweep velocities
% and small deltas per the paper.
% N = number of refocusing pulses
% Outputs:
%   E(i): excitation phases, radians
%   R(i): receiver phases, radians
%   delta: small angle delta_i sequence (radians)
%   Delta: sweep velocities (radians)
% References:
% @article{LeRoux2002Non-CPMGSignal,
%     title = {{Non-CPMG Fast Spin Echo with full signal}},
%     year = {2002},
%     journal = {Journal of Magnetic Resonance},
%     author = {Le Roux, Patrick},
%     number = {2},
%     pages = {278--292},
%     volume = {155},
%     publisher = {Academic Press Inc.},
%     doi = {10.1006/jmre.2002.2523},
%     issn = {10907807},
%     pmid = {12036339},
%     keywords = {Carr, Purcell, Meiboom, Gill, Quadratic phase modulation, Spin echo}
% }



% Table 1 sweep velocities (as fractions of 2*pi)
Delta_table = [0.191438, 0.192650, 0.225601, 0.197626, ...
               0.129640, 0.197671, 0.282091];
Delta_const = 957/4999; % = 0.191438

% Fill up Delta_i for N pulses
Delta = zeros(1, N);
if N <= 7
    Delta(1:N) = Delta_table(1:N);
else
    Delta(1:7) = Delta_table;
    Delta(8:N) = Delta_const;
end
Delta = Delta * 2 * pi; % Convert to radians

% Now build phase increments (delta_i) from sweep velocities
delta = zeros(1, N);
delta(1) = Delta(1);              % δ₁ = Δ₁
for i = 2:N
    delta(i) = delta(i-1) + Delta(i); % δᵢ = δᵢ₋₁ + Δᵢ
end

% Initialize phases
E = zeros(1, N); % Rx
R = zeros(1, N);   % Tx

% initialize Receiver 0
R_0 = 0;
% Compute E(i) and R(i)
for i = 1:N
     if i == 1
        E(i) = R_0 + delta(i); 
        R(i) = E(i) + delta(i);

     else 
        E(i) = R(i-1) + delta(i); 
        R(i) = E(i) + delta(i);
      

     end 
end

R = R(1:N);

% Return unwrapped (raw) radians
E_rad = E;
R_rad = R;

end