function [A, X, Y] = processAndPlot(fileName, N, offset, plotFlag, lineWidth, clr)
%   Load data, process N columns, and optionally plot results.
%
%   [A, X, Y] = process_and_plot(fileName, N, offset, plotFlag, lineWidth, clr)
%
%   Inputs:
%       fileName   - Name of the file to load data from
%       N          - Number of columns per block (used in reshaping)
%       offset     - Horizontal offset applied to X-values when plotting
%       plotFlag   - If true (nonzero), plots the data
%       lineWidth  - Width of plotted lines
%       clr        - Color of plotted lines (RGB triplet or MATLAB color char)
%
%   Outputs:
%       A   - Processed data matrix (N rows, variable number of columns)
%       X   - First row of A (independent variable)
%       Y   - Last dependent variable plotted (or 2nd row if plotFlag = false)

    % -----------------------------
    % Load and reshape data
    % -----------------------------
    f = fopen(fileName, 'r');
    if f < 0
        warning('Could not open file: %s', fileName);
        A = []; X = []; Y = [];
        return;
    end

    % Read raw data as a column vector
    B = fscanf(f, '%f');
    fclose(f);

    % Trim to a multiple of N and reshape into N rows
    K = floor(numel(B) / N); 
    A = reshape(B(1:N*K), N, K);

    % Extract independent (X) and dependent (Y) variables
    X = A(1, :);
    Y = A(2, :);  % Default: second row (overwritten during plotting if plotFlag=1)

    % -----------------------------
    % Plotting (if requested)
    % -----------------------------
    if plotFlag
        holdState = ishold;  % remember current hold state
        hold on;
        for k = 2:size(A, 1)
            Y = A(k, :); % overwrite Y with current row
            plot(X + offset, Y, '-', 'Color', clr, 'LineWidth', lineWidth);
        end
        if ~holdState, hold off; end  % restore original hold state
    end
end
