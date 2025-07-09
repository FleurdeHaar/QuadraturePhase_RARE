function [E_rad, R_rad] = phase_mod_Leroux3(numberpulses)
% phase_mod_Leroux: Computes Le Roux quadratic phase modulation series
%   and prints excitation/receiver phases in radians with high precision.
%
% INPUT:
%   numberpulses - number of pulses (ETL or total refocusing pulses)
%
% OUTPUT:
%   E_rad - excitation phases E(i) in radians
%   R_rad - receiver phases R(i) in radians
% References:
% @article{Bastin2002OnBrain,
%     title = {{On the application of a non‐CPMG single‐shot fast spin‐echo sequence to diffusion tensor MRI of the human brain}},
%     year = {2002},
%     journal = {Magnetic Resonance in Medicine},
%     author = {Bastin, Mark E. and Le Roux, Patrick},
%     number = {1},
%     month = {7},
%     pages = {6--14},
%     volume = {48},
%     doi = {10.1002/mrm.10214},
%     issn = {0740-3194}
% }

    % Force long display precision
    format long

    % Define constant F
    F = (20/49) * pi; % 1.28

    % Preallocate arrays
    E = zeros(1, numberpulses);
    R = zeros(1, numberpulses);
    delta_array = zeros(1, numberpulses + 1);
    
   
    % Precompute delta values, phase increments
    for i = 1:numberpulses + 1
        if i == 1
            delta_array(i) = 2 * F - 0.25*pi;
            
        elseif i == 2
            delta_array(i) = -3 * F - 0.25*pi;
        elseif i == 3
             delta_array(i) = 4 * F - 0.25*pi;
        else
            delta_array(i) = F * (i-1);
        end
    end
    
    % Compute E(i) and R(i)
    for i = 1:numberpulses
         if i == 1
          % assume --> R_0 = 0 
          % First excitation and receiver axis 
              E(i) =  delta_array(i);
              R(i) =  E(i)+ delta_array(i);
             
         else
         % From the next pulse
              E(i) = E(i-1) + delta_array(i-1)+ delta_array(i);
              R(i) = E(i) + delta_array(i);
            
        end
    end

    % Return unwrapped (raw) radians
    E_rad = E;
    R_rad = R;
    
    % Print each value with 12 decimal places
%     fprintf('\nHigh-precision Le Roux phases:\n');
%     fprintf('  i    \tE_rad (rad)        \tR_rad (rad)\n');
%     for i = 1:numberpulses
%         fprintf('  %2d   \t%.12f\t%.12f\n', i, E_rad(i), R_rad(i));
%     end
%     fprintf('\n');

end
