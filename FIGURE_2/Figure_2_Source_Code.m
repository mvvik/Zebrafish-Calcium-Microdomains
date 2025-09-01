%
%   MATLAB Script for kinetic data fit in Fig. 2 in the eLife manuscript:
%
%  "Nanophysiology Approach Reveals Diversity in Ca2+ Microdomains ..."
%   Rameshkumar, Shrestha, Boff, Hoon, Matveev, Zenisek, Vaithianathan
%   https://elifesciences.org/reviewed-preprints/105875#s2
%               Code: Victor Matveev, Sep 1, 2025
% =========================================================================

nTrials = 800;     % --- Number of random starting points for optimization
T0      = 24;      % --- Stimulus onset time (ms)
tfs     = 12;      % --- Font size for figure titles

% --------  Columns to analyze (2 = proximal, 4 = distal)
COLS  = 2:2:4;
nCols = numel(COLS);

% --------  Parameter search ranges (log-spaced)
minP = log([0.004  1    50   1e-5   1e-3 ]);
maxP = log([0.4   20   200   10     100  ]);
nPars = numel(minP);

% --------  Utility for formatting numbers with ~3 significant digits
sigDgts = @(x, n) (x > 1) * max([0, n - floor(log10(x)) - 1]) + ...
                  (x < 1) * (n + floor(abs(log10(x - floor(x)))) );
DS      = @(x) num2str(sigDgts(x, 3));

% --------  Colors and labels
Str = {'Proximal', 'Distal'};
Clr = [0.9 0   0;   % red
       0    0 1;   % blue
       0   0.7 0.15]; % green (not used here, but available)

figure;

% =================== Loop over experimental conditions ===================
for mode = 1:2
    switch mode
        case 1   % --- High-Affinity Indicator Dye
            dataStr  = 'Fig 2A: Cal520HA';
            dataFile = 'Figure_2A_Source_Data.txt'; 
        case 2   % --- Low-Affinity Indicator Dye
            dataStr  = 'Fig 2B: Cal520LA';  
            dataFile = 'Figure_2B_Source_Data.txt';
    end

    % ================  Load data (time + measurements)  ==================
    x = processNcolumns(dataFile, 4);
    T = x(1, :);   % First row is time

    % ===========  Baseline correction: subtract mean before T0  ==========
    ind2 = find(T > T0, 1) - 1;
    if isempty(ind2) || ind2 < 1
        error('T0 (%.2f) is outside time range.', T0);
    end
    for jjj = 2:2:4
        x(jjj, :) = x(jjj, :) - mean(x(jjj, 4:ind2));
    end

    % ============  Oversampled time axis for smoother fitting ============
    MM = round(5 * max(T));
    TT = linspace(min(T), max(T), MM);

    % Optimization settings
    opt = optimset('TolX', 1e-5, 'TolFun', 1e-5, 'Display', 'off');

    % ==================  Loop over proximal/distal signals ===============
    for jjj = 1:nCols
        COL = COLS(jjj);
        subplot(2, 2, (mode - 1)*2 + jjj);

        % ===================== Interpolated data trace  ==================
        YY  = interp1(T, x(COL,:), TT, 'linear');

        % ==========================  Find peak ===========================
        % Given S/N level, finding max by data fitting wouldn't be advisable

        [AMP, indMax] = max(YY);   % Direct peak detection / extraction
        Tmax = TT(indMax);
        DDT  = Tmax - T0;

        % ====================  Define model components  ==================

        sigma = @(p) AMP * (tanh(p(5)*(TT - T0 - p(4))) + tanh(p(5)*p(4))) ./ ...
                           (tanh(p(5)*(DDT - p(4)))   + tanh(p(5)*p(4)));

        C1    = @(p) (TT >= Tmax) .* abs(p(1))       .* exp(-(TT-Tmax)./p(2));
        C2    = @(p) (TT >= Tmax) .* abs(AMP-p(1))   .* exp(-(TT-Tmax)./p(3));

        Y0    = @(p,flag) sigma(p).*(TT>=T0 & TT<Tmax) + C1(p) + ...
                          flag*C2(p) + 2*(1-flag)*abs(AMP-abs(p(1)));

        Error = @(p,flag) sum(abs(Y0(p,flag) - YY).^2);

        % ========== Parameter optimization with multiple starts ==========
        ResultsOut = zeros(nTrials, nPars+1);
        flag = 1;  % Two-component fit (set 0 for single exp)

        parfor ind = 1:nTrials
            P0  = exp(minP + rand(1,nPars).*(maxP-minP)); % log-uniform init
            P1  = fminsearch(@(pp) Error(pp, flag), P0, opt);
            ERR = Error(P1, flag);
            ResultsOut(ind,:) = [P1, ERR];
        end

        % =======================  Best fit  ==============================
        [ERR, I] = min(ResultsOut(:,end));
        P = abs(ResultsOut(I,1:nPars));

        % ==================  Extract rise time  ==========================
        ttt = linspace(T0, Tmax, 500);
        sigFun = (tanh(P(5)*(ttt - T0 - P(4))) + tanh(P(5)*P(4))) ./ ...
                 (tanh(P(5)*(DDT - P(4)))     + tanh(P(5)*P(4)));
        iRise = find(sigFun > 0.5, 1);
        tRise = ttt(iRise) - T0;

        % ======================== Reporting ==============================
        if flag
            if P(2) > P(3), P = [1-P(1), P(3), P(2)]; end
            fStr1 = sprintf('%.3f exp(-t/%.3f) + %.3f exp(-t/%.3f)', ...
                            P(1)/AMP, P(2), 1-P(1)/AMP, P(3));
        else
            fStr1 = sprintf('exp(-t/%.3f)', P(2));
        end

        fStr2 = sprintf('%s-%s: t_{rise}=%.3f ms', dataStr, Str{jjj}, tRise);
        fprintf('\n Best Err = %g, tRise = %.3f ms\n', ERR, tRise);

        % ========================== Plot =================================
        hold off;
        plot(TT, Y0(P,flag), 'k-', 'LineWidth', 2); hold on;
        plot(TT, YY, '-', 'LineWidth', 1, 'Color', Clr(jjj,:));
        title({fStr2, fStr1}, 'FontSize', tfs);
        axis tight;
        if mode == 2, xlabel('Time (ms)'); end
    end
end

