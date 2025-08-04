function [S, S1, S2, phasediag, P] = epg_QPRARE2_override( ...
    flipangle, etl, T1, T2, esp, prepulses, ...
    phi0, Tx_override, Rx_override, applyBlip, doPlot)
% EPG_QPRARE2_OVERRIDE
%   Like epg_QPRARE2 but lets you pass full-length Tx- and Rx-phase vectors.
%   If Tx_override or Rx_override are empty, uses zero deviations.

  %% 1) Parse inputs & defaults
  if nargin<1 || isempty(flipangle),    flipangle   = pi*ones(1,etl); end
  if nargin<2 || isempty(etl),           etl         = numel(flipangle); end
  if etl>numel(flipangle)
      flipangle(end+1:etl) = flipangle(end);
  end
  if nargin<3 || isempty(T1),            T1           = 1.0; end
  if nargin<4 || isempty(T2),            T2           = 0.1; end
  if nargin<5 || isempty(esp),           esp          = 0.01; end
  if nargin<6,                           prepulses    = 3;    end
  if nargin<7 || isempty(phi0),          phi0         = 0;    end
  if nargin<8 || isempty(Tx_override),   Tx_override  = zeros(etl,1); end
  if numel(Tx_override)~=etl, error('Tx_override must have length etl'); end
  if nargin<9 || isempty(Rx_override),   Rx_override  = zeros(etl,1); end
  if numel(Rx_override)~=etl, error('Rx_override must have length etl'); end
  if nargin<10,                          applyBlip    = false; end
  if nargin<11,                          doPlot       = false; end

  %% 2) Le Roux E(n), R(n) baseline
  if prepulses==3
      [E_base, R_base] = phase_mod_Leroux3(etl);
  elseif prepulses==7
      [E_base, R_base] = phase_mod_Leroux_7prepulsescheme(etl);
  else
      error('prepulses must be 3 or 7');
  end

  %% 3) Prepare storage
  P      = zeros(3,2*etl);   P(3,1)=1;
  Pstore = zeros(4*etl, etl);
  Zstore = zeros(2*etl, etl);
  S      = zeros(1, etl);
  S1     = zeros(1, etl);
  S2     = zeros(1, etl);

  if applyBlip
      Nlines    = (etl - prepulses)/2;
      totalSweep= pi;
      deltaPhi  = totalSweep / Nlines;
      Nstates   = size(P,2);
      idx0      = etl + 1;
      orders    = (1:Nstates) - idx0;
      phi_blip  = orders * deltaPhi;
  end

  %% 4) 90° excitation with Tx override
  tx0 = phi0 + Tx_override(1);
  P = epg_rf(P, pi/2, pi + tx0);
  if applyBlip
      P(1,:) = P(1,:).*exp(-1i*phi_blip);
      P(2,:) = P(2,:).*exp( 1i*phi_blip);
  end

  %% 5) Echo train
  for ech = 1:etl
      % pre‐RF relaxation
      P = epg_grelax(P, T1, T2, esp/2, 1,0,1,1);

      % RF with both E(n) and Tx_override
      txph = angle(flipangle(ech)) + E_base(ech) + Tx_override(ech);
      P    = epg_rf(P, abs(flipangle(ech)), txph);

      % post‐RF relaxation
      P = epg_grelax(P, T1, T2, esp/2, 1,0,1,1);

      % record echo with Rx override
      rawF0  = P(1,1);
      rph    = R_base(ech) + Rx_override(ech);
      rawF01 = rawF0 * exp(-1i * rph);
      S(ech) = rawF01;
      S1(ech)= real(rawF01);
      S2(ech)= imag(rawF01);

      % store EPG states
      Pstore(1:2*etl,ech)       = flipud(P(1,:).');
      Pstore(2*etl:4*etl-1,ech) = P(2,:).';
      Zstore(:,ech)             = P(3,:).';
  end

  phasediag = [Pstore; Zstore];

  %% 6) Optional plotting
  if doPlot
      figure;
      subplot(2,1,1);
      plot(1:etl, real(S), '-o', 1:etl, imag(S), '-x');
      legend('S1','S2'); title('Echo signals');
      subplot(2,1,2);
      plot(1:etl, unwrap(angle(S)));
      title('Echo phase evolution');
  end
end
