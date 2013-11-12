function h = make_plot(M, mis_per_epoch)
% MAKE_PLOT Make convergence plot for one matrix.
%   h = make_plot(M, mis_per_epoch)
%   Plots the convergence for the matrix M and returns h, a handle of the
%   figure. The parameter mis_per_epoch indicates how often, in number of
%   major iterations (sequences of n steps), the residual is printed.
%

    h = 0;
    % thread counts
    T = [];
    % tall matrix of residual norms
    R = [];
    read_data;
    plot_data
    
    function read_data
        files = dir([M '*.txt']);
        for j = 1:length(files)
            X = dlmread(files(j).name);
            T(j) = sscanf(files(j).name, [M '-T%d.txt']);
            R(:, j) = X(:, 2);
        end
    end

    function plot_data
        h = figure;
        X = (0:size(R, 1) - 1) * mis_per_epoch;
        semilogy(X, R);
        xlabel('major iteration');
        ylabel('norm(b - Ax) / norm(b)');
        title(M);
        legend(strtrim(cellstr(num2str(T', 'T = %d'))));
    end
end
