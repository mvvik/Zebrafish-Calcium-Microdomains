%
%   MATLAB Script for kinetic data fit in Fig. 3 in the eLife manuscript:
%
%  "Nanophysiology Approach Reveals Diversity in Ca2+ Microdomains ..."
%   Rameshkumar, Shrestha, Boff, Hoon, Matveev, Zenisek, Vaithianathan
%   https://elifesciences.org/reviewed-preprints/105875#s2
%               Code: Victor Matveev, Sep 1, 2025
% =========================================================================

clear;
minT    = 1;     % --- Beginning of data read time window
maxT    = 295;   % --- End of data read time window
nTrials = 400;   % --- Number of parameter optimization trials
T0      = 24;    % --- Stimulus onset time (ms)

%--- Helper functions for digit formatting --------------------------------

sigDgts = @(x, n) (x > 1) * max([0, n - floor(log10(x)) - 1]) + ...
                  (x < 1) * (n + floor(abs(log10(x - floor(x)))) );
DS      = @(x) num2str(sigDgts(x, 3));

%--------------------------------------------------------------------------

COLS    = 2:2:4;                             % Columns to analyze 
nCols   = numel(COLS);                       % Total number of columns

minP    = log([0.004  1  50  1e-5  1e-3 ]);  % Lower param bounds
maxP    = log([0.4   20 200    10  100 ]);   % Upper param bounds 
nPars   = numel(minP);

figure;
tfs = 12;                                   % Font size for titles

for mode = 1 : 4
    %------ Select dataset and title string suffix ------------------------

    switch mode
        case 1
            dataStr = 'Fig 3C: 3Cal520HA'; 
            fname   = 'Figure_3C_Source_Data.txt';
        case 2
            dataStr = 'Fig 3D: Cal520LA';   
            fname   = 'Figure_3D_Source_Data.txt';
        case 3
            dataStr = 'Fig 3E: Cal520LA';  
            fname   = 'Figure_3E_Source_Data.txt';
        case 4
            dataStr = 'Fig 3F: Cal520LA';  
            fname   = 'Figure_3F_Source_Data.txt';
    end

    %--- Load and preprocess data -----------------------------------------

    x    = processNcolumns(fname, 4);
    Str  = {'RBP-Prox', 'RBP-Dist',  'RBP-Prox', 'RBP-Dist', ...
            'RBP-Prox', 'Free-Prox', 'RBP-Dist', 'Free-Dist'};

    Clr  = [0.9 0 0;   0 0 1;  0 0.7 0.15;  0.9 0 0;  0 0 1];
    T    = x(1, :);                         % First column = time (ms)

    ind2 = find(x(1,:) > T0, 1) - 1;
    for jjj = 2:2:4
        x(jjj, :) = x(jjj, :) - mean(x(jjj, 4:ind2));
    end

    %--- Oversampling and optimization setup ------------------------------

    MM  = round(5 * max(T));                % # interpolation points
    TT  = linspace(min(T), max(T), MM);
    opt = optimset('TolX', 1e-5, 'TolFun', 1e-5, 'Display', 'off');

    %--- Loop over selected columns ---------------------------------------

    for jjj = 1 : nCols
        flag = 1;
        COL  = COLS(jjj);
        subplot(4, 2, (mode - 1)*2 + jjj);

        fname = ['Data/DataFit_Results_TwoExp_COL_', ...
                  num2str(COL), '_NEW.mat'];

        %--- Interpolate trace --------------------------------------------

        YY  = interp1(T, x(COL, :), TT, 'linear');
        [AMP, indMax] = max(YY);
        Tmax  = TT(indMax);
        DDT   = Tmax - T0;

        %--- Model components ---------------------------------------------

        sigma = @(p) AMP * (tanh(p(5)*(TT - T0 - p(4))) + tanh(p(5)*p(4))) ...
                       / (tanh(p(5)*(DDT - p(4))) + tanh(p(5)*p(4)));
        C1    = @(p) (TT >= Tmax) .* abs(p(1))      .* exp(-(TT-Tmax)./p(2));
        C2    = @(p) (TT >= Tmax) .* abs(AMP-p(1))  .* exp(-(TT-Tmax)./p(3));
        Y0    = @(p,flag) sigma(p).*(TT >= T0).*(TT < Tmax) + ...
                          C1(p) + flag*C2(p) + ...
                          2*(1-flag)*abs((AMP-abs(p(1))));
        Error = @(p,flag) sum(abs(Y0(p,flag) - YY).^2);

        %--- Parameter optimization ---------------------------------------

        ResultsOut = zeros(nTrials, nPars+1);
        parfor ind = 1 : nTrials
            P1  = exp(minP + rand(1, nPars) .* (maxP - minP));
            P1  = fminsearch(@(x) Error(x, flag), P1, opt);
            ERR = Error(P1, flag);
            ResultsOut(ind, :) = [P1, ERR];
        end

        [~, I] = min(ResultsOut(:, end));
        P      = abs(ResultsOut(I, 1:nPars));
        ERR    = ResultsOut(I, end);

        %--- Extract timing metrics ---------------------------------------

        [~, ind] = max(Y0(P, flag));
        tMax = TT(ind);
        DT   = tMax - T0;

        ttt   = linspace(T0, Tmax, 500);
        sig   = @(p) (tanh(p(5)*(ttt - T0 - p(4))) + tanh(p(5)*p(4))) ...
                       / (tanh(p(5)*(DDT - p(4))) + tanh(p(5)*p(4)));
        iRise = find(sig(P) > 0.5, 1);
        tRise = ttt(iRise) - T0;

        %--- Format fitting results ---------------------------------------

        if flag
            if P(2) > P(3), P = [1-P(1); P(3); P(2); P(4); P(5)]; end
            formStr1 = [' %.', DS(P(1)/AMP), 'f exp(-t/%.', DS(P(2)), 'f) + ' ...
                        '%.', DS(1-P(1)/AMP), 'f exp(-t/%.', DS(P(3)), 'f)'];
            fStr1    = sprintf(formStr1, P(1)/AMP, P(2), 1-P(1)/AMP, P(3));
        else
            formStr1 = ['exp(-t/%.', DS(P(2)), 'f)'];
            fStr1    = sprintf(formStr1, P(2));
        end

        formStr2 = ['%s-%s: t_{rise}=%.', DS(tRise), 'f ms'];
        fStr2    = sprintf(formStr2, dataStr, Str{jjj+(mode-1)*2}, tRise);
        fprintf('\n Best Err = %g  DT = %g\n', ERR, tRise);

        %--- Plot data and fit --------------------------------------------

        hold off;
        plot(TT, Y0(P, flag), '-', 'LineWidth', 2, 'Color', 'k'); hold on;
        plot(TT, YY,          '-', 'LineWidth', 1, 'Color', Clr(jjj,:));
        title({fStr2, fStr1}, 'FontSize', tfs);
        axis tight; drawnow;
        if mode == 4, xlabel('Time (ms)'); end
    end
end
