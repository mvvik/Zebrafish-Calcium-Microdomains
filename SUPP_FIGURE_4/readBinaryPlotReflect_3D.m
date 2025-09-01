% =========================================================================
%     Read and plot 3D (or 4D) field dump produced by CalC binary plot
% =========================================================================

logMode = true;       % take log10 of Ca values
nSigDgts = 0;         % number of significant digits for color bar

% -----------------------------
%   Read geometry and grid
% -----------------------------
G  = fread(f, 1, 'int');   % geometry (lowest 2 bits = # of dimensions)
N1 = fread(f, 1, 'int');   % x-nodes
N2 = fread(f, 1, 'int');   % y-nodes
N3 = fread(f, 1, 'int');   % z-nodes

X1 = fread(f, N1, 'double'); % x-axis coordinates
X2 = fread(f, N2, 'double'); % y-axis coordinates
X3 = fread(f, N3, 'double'); % z-axis coordinates

% Symmetric reflection in x-axis
X1 = [-X1(end:-1:1), X1];
[X, Z] = meshgrid(X1, X3);

% -----------------------------
%   Choose z-slice
% -----------------------------
yNode  = 1;           % index of y-node for z-slice
yCoord = X2(yNode);   % y-coordinate of slice

N     = N1 * N2 * N3; % total number of nodes
state = false;        % loop control

% -----------------------------
%   Read and plot loop
% -----------------------------
while ~state
    t = fread(f, 1, 'double');   % read time
    if feof(f)
        break;
    end

    A = fread(f, N, 'double');   % read data
    if numel(A) < N
        warning('Incomplete data at time %g. Skipping.', t);
        break;
    end

    B = reshape(A, N1, N2, N3);           % 3D array
    Ca0 = squeeze(B(yNode, :, :))';       % z-slice with transpose
    Ca = zeros(N3, 2*N1);                 % doubled for symmetry

    for k = 1:N3
        Ca(k, :) = [Ca0(k, end:-1:1), Ca0(k, :)]; % mirrored slice
    end

    % Skip t=0 if needed
    if t == 0
        continue;
    end

    % -----------------------------
    %   Logarithmic transform (if enabled)
    % -----------------------------
    if logMode
        Ca(Ca <= 0) = 1e-3;
        Ca = log10(Ca);
    end

    % -----------------------------
    %   Determine color levels
    % -----------------------------
    cornerCa = [ Ca(1,1), Ca(1,end), Ca(end,1), Ca(end,end) ];
    minZ = min(cornerCa);
    maxZ = max(Ca(:));
    ZZ   = linspace(minZ, maxZ, 10);

    CaLevel = floor(10.^ZZ);

    % ensure unique levels
    for shift = 1:2
        if any(diff(CaLevel) == 0)
            CaLevel = floor(10.^(ZZ + shift)) / 10^shift;
        end
    end

    % round large values
    inds = CaLevel > 1.8;
    CaLevel(inds) = round(CaLevel(inds));
    CaLevelLog = log10(CaLevel);

    % -----------------------------
    %   Plot contour
    % -----------------------------
    contourf(X, Z, Ca, CaLevelLog);  % plot z-slice
    caxis([min(CaLevelLog), max(CaLevelLog)]);
    title(titleStr);
    xlabel('X (\mum)');

    h = colorbar;
    shading flat;
    set(h, 'Ticks', CaLevelLog);

    % -----------------------------
    %   Set colorbar labels
    % -----------------------------
    labelStr = arrayfun(@(v) num2str(v), CaLevel, 'UniformOutput', false);
    labelStr{end} = [labelStr{end}, ' \muM'];  % append units
    set(h, 'TickLabels', labelStr);

    % -----------------------------
    %   Exit if end of file reached
    % -----------------------------
    state = feof(f);
end

fclose(f);
