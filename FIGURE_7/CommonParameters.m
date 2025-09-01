% ========================================================================
%      MODE SELECTION: Define Buffers and Simulation Prefix
%      Note: maintain compatibility with both MATLAB and CalC
% ========================================================================

if mode == 1        % ---- Low EGTA (200 µM = 0.2 mM) ----
        Prefix      = 'EGTA0p2';
        BufferStr   = '0.2 mM EGTA';
        EGTA.total  = 200;                 % in µM
        BAPTA.total = 0;                   % in µM
end

if mode == 2        % ---- High EGTA (10 mM) ----
        Prefix      = 'EGTA10';
        BufferStr   = '10 mM EGTA';
        EGTA.total  = 10000;               % in µM
        BAPTA.total = 0;                   % in µM
end

if mode == 3        % ---- Moderate BAPTA (2 mM) ----
        Prefix      = 'BAPTA2';
        BufferStr   = '2 mM BAPTA';
        EGTA.total  = 0;                   % in µM
        BAPTA.total = 2000;                % in µM
end

% ========================================================================
%                 GEOMETRIC PARAMETERS OF RIBBON SYNAPSE
% ========================================================================

Pulse      = 10;                   % Pulse durations (ms)
TotalTime  = 200;                  % Total simulation times (ms)
postPulse  = TotalTime - Pulse;    % Post-pulse duration (ms)

ICa        = 1.0;                  % Max Ca2+ current (pA)
KD         = 2;                    % Endogenous buffer affinity
Diff       = 0;                    % Buffer diffusivity (immobile)
Kplus      = 0.10;                 % Buffer forward binding rate (µM⁻¹ ms⁻¹)

% Ribbon ellipsoid dimensions

rStalk   = 0.015;                  % Stalk radius (µm)
zAxis    = 0.190;                  % Ellipsoid semi-axis along z (µm)
xAxis    = zAxis;                  % Semi-axis along x (µm)
yAxis    = zAxis * 0.375;          % Semi-axis along y (µm)

% Vertical displacement of ribbon
zElevate = 0.03;                   % Shift above base (µm)
z0       =   zAxis + zElevate;     % Ribbon center height
zTop     = 2*zAxis + zElevate;     % Top surface of ribbon

% Simulation box dimensions (µm)
Lx = 0.64;
Ly = Lx;
Lz = 1.1;
