function A = processNcolumns(fname, P)
%-- PROCESSNCOLUMNS  Read numeric data from file and reshape into P rows
%   A = processNcolumns(fname, P) reads all floating-point numbers
%   from the text file 'fname' and reshapes them into a P×K matrix,
%   truncating excess values if the count is not a multiple of P.

    % Try opening file for reading
    f = fopen(fname,'r');
    if f < 0
        A = [];  % Return empty if file not found
        return;
    end

    % Read all floats into a column vector
    B = fscanf(f, '%f');
    fclose(f);

    % Reshape into P×K (discard remainder if needed)
    K = floor(numel(B) / P);
    A = reshape(B(1:P*K), P, K);
end