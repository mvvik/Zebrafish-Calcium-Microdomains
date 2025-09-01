%
%   Main MATLAB Script for Figure 7 in the eLife manuscript:
%
%  "Nanophysiology Approach Reveals Diversity in Ca2+ Microdomains ..."
%   Rameshkumar, Shrestha, Boff, Hoon, Matveev, Zenisek, Vaithianathan
%   https://elifesciences.org/reviewed-preprints/105875#s2
%                 Code: Victor Matveev, Aug 31, 2025
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Recompute or just plot?
recomputeFlag = 1;

% --- Special Locations (X, Z coordinates in microns)
Locations = [ 0.02, 0.01; ...
              0.15, 0.01; ...
              0.15, 0.25; ...
              0.00, 0.45; ...
              0.00, 0.80 ];

% --- Colors for plotting special locations (RGB triplets)
Clrs = [ 1,   0,   0;  ...  % red
         0,   0,   1;  ...  % blue
         1, 0.7,   0;  ...  % orange
         1,   0,   1;  ...  % magenta
         0,   0,   0 ];     % black

% ---  Parameters for two simulation runs

ArrayKappa =  [      100,       720,      100,      720,      100,     720 ];   % Buffering capacities
ArrayMode  =  [        1,         1,        2,        2,        3,       3 ];   % EGTA/BAPTA modes
PrefixArray = {'EGTA0p2', 'EGTA0p2', 'EGTA10', 'EGTA10', 'BAPTA2', 'BAPTA2'};   % Identifier file prefix

% --- Panel indices for 2D contour plots
Panels2D = [2  3  8  9 14 15 20 21 26 27; ...
            5  6 11 12 17 18 23 24 29 30];

% --- Plot labels
symbol = 'ABC';
nLocs  = size(Locations, 1);   % number of probe locations

% =========================================================================
%           Compute [Ca2+] across all conditions (if requested)
% =========================================================================

if recomputeFlag
    % Start parallel pool with 6 workers (safeguarded)
    if isempty(gcp('nocreate')), parpool(6); end
    
	% --- Simulation call setup -------------------------------------------
	CalC_version = '6109';
	Script = ' RibbonModel3D.par';

	if contains(computer, 'WIN')
		 prog = ['..\CALC\cwin', CalC_version, 'x64.exe '];
	elseif contains(computer('arch'), 'maca64')
		 prog = ['../CALC/cmac', CalC_version, 'xM1 '];
		 system(['chmod +x ', prog]);
	elseif contains(computer('arch'), 'maci64')
		 prog = ['../CALC/cmac', CalC_version, 'x86 '];
		 system(['chmod +x ', prog]);
	else
		 fprintf('Do not recognize this architecture: %s\n', computer('arch'));
		 fprintf('Compile CalC first: github.com/mvvik/CalC-simple-buffer');
		 return;
    end

    parfor n = 1 : 6
        Kappa  = ArrayKappa(n);  % Select endogenous buffer capacity
        mode   =  ArrayMode(n);  % Select EGTA / BAPTA condition

        % Compute [Ca2+] for the 10 ms pulse
        cmd1 = sprintf('%s %s 1 %d %g', prog, Script, mode, Kappa);
        [~, Msg1] = system(cmd1);
        fprintf('Pulse done: mode = %d, Kappa = %g\n', mode, Kappa);
        if numel(Msg1) ~= 0, warning('Error: %s', Msg1); end

        % Compute [Ca2+] for the 190 ms post-pulse
        cmd2 = sprintf('%s %s 2 %d %g ', prog, Script, mode, Kappa);
        [~, Msg2] = system(cmd2);
        fprintf('Post-pulse done: mode = %d, Kappa = %g\n', mode, Kappa);
        if numel(Msg2) ~= 0, warning('Error: %s', Msg2); end
    end
end

% =========================================================================
% --------               Plot traces & contour maps
% =========================================================================
for n = 1 : 6
    Kappa  = ArrayKappa(n);
    mode   = ArrayMode(n);
    nn     = 1 - bitand(n,1);
    CommonParameters;
    figure(mode);
    
    % --- Loop over probe locations
    for k = 1:nLocs
        subplot(nLocs, 6, 1 + 3*nn + 6*(k-1)); 
        hold on;

        % --- Pulse trace
        fname1 = sprintf('DATA/%s_Kappa%d_KD%d_iter1_Ca%d', Prefix, Kappa, KD, k);
        [E1, XX, Y] = processAndPlot(fname1, 2, 0, 1, 2, Clrs(k,:));
        minY = Y(1); maxY = Y(end);

        % --- Post-pulse trace (concatenated with offset)
        fname2 = sprintf('DATA/%s_Kappa%d_KD%d_iter2_Ca%d', Prefix, Kappa, KD, k);
        [E2, XX, Y] = processAndPlot(fname2, 2, XX(end), 1, 2, Clrs(k,:));
        maxY = max([max(Y), maxY]);

        % --- Tick placement logic (more robust)
        dY   = maxY - minY;
        dgts = min(round(-log10(dY)) + 2, 3);
        Y1   = round((minY + 0.5*dY) * 10^dgts) / 10^dgts;
        Y2   = round(maxY * 10^dgts) / 10^dgts;
        maxY = max([Y2, maxY]);

        if Y1 == Y2
            yticks([minY, Y2]);
        elseif minY == Y1
            yticks([Y1, Y2]);
        else
            yticks([minY, Y1, Y2]);
        end
        
        % --- Legend & labels
        TitleStr = sprintf('X = %.0fnm, Z = %.0fnm', 1e3*Locations(k,1:2) );
        lgnd = legend(TitleStr, 'Location', 'northeast');
        lgnd.ItemTokenSize = [3 10];
        set(gca, 'FontSize', 9);
        axis([0, 200, minY, minY + dY*1.01]);

        if k == 1
            ylabel('\muM', 'FontSize', 9);
            title([symbol(mode), num2str(nn+1)], 'FontSize', 16);
        elseif k == nLocs
            xlabel('Time (ms)', 'FontSize', 9);
        end
    end

    % --- 2D Contour plot (simulation dump)
    subplot(nLocs, 6, Panels2D(nn+1,:));
    fnameBin = sprintf('DATA/CaPulse%s_Pulse%dms_Kappa%d_KD2.dat', Prefix, Pulse, Kappa);
    f  = fopen(fnameBin, 'rb');
    if f < 0
        warning('Could not open %s', fnameBin);
    else
        hold on;
        titleStr = sprintf('%s, \\kappa=%d, I_{Ca}=%.3gpA, Time=%dms', BufferStr, Kappa, ICa, Pulse);
        readBinaryPlotReflect_3D;  % generates contour plot
        axis([-0.6 0.6001 0 0.9]);
    end

    % --- Mark probe locations
    for m = 1:nLocs
        plot(Locations(m,1), Locations(m,2), 'o', ...
             'MarkerFaceColor', Clrs(m,:), ...
             'MarkerEdgeColor', Clrs(m,:), ...
             'MarkerSize', 5);
    end
end
