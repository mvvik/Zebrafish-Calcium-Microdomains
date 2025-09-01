
% ========================================================================
%           BUFFERING AND GEOMETRIC PARAMETERS OF RIBBON SYNAPSE
% ========================================================================

Pulse      = 10;                   % Pulse durations (ms)
TotalTime  = 200;                  % Total simulation times (ms)
postPulse  = TotalTime - Pulse;    % Post-pulse duration (ms)

ICA        = 1.0;                  % Max Ca2+ current (pA)
KD         = 2;                    % Endogenous buffer affinity
Diff       = 0.05;                 % Buffer diffusivity (immobile)
Kplus      = 0.10;                 % Buffer forward binding rate (µM⁻¹ ms⁻¹)

% --- Ribbon ellipsoid dimensions

rStalk   = 0.015;                  % Stalk radius (µm)
zAxis    = 0.190;                  % Ellipsoid semi-axis along z (µm)
xAxis    = zAxis;                  % Semi-axis along x (µm)
yAxis    = zAxis * 0.375;          % Semi-axis along y (µm)

% --- Vertical displacement of ribbon
zElevate = 0.03;                   % Shift above base (µm)
z0       =   zAxis + zElevate;     % Ribbon center height
zTop     = 2*zAxis + zElevate;     % Top surface of ribbon

% --- Simulation box dimensions (µm)
Lx = 0.64;
Ly = Lx;
Lz = 1.1;

% --- No exogenous buffers in Supp Fig 4
EGTA.total  = 0;    
BAPTA.total = 0; 
