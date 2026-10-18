function R = build_uniform_covariance(N, L)

% =========================================================
% BUILD UNIFORM CHANNEL COVARIANCE MATRIX
%
% Based on Eq. (20) of:
% H.-R. Park,
% "A Low-Complexity Channel Estimation for OFDM Systems
% Based on CIR Length Adaptation"
%
% gamma(k1,k2) =
%
%       1,                                      k1 = k2
%
%       [1 - exp(-j*2*pi*L*(k1-k2)/N)]
%       ---------------------------------        k1 ~= k2
%              j*2*pi*L*(k1-k2)/N
%
% =========================================================

    % -----------------------------------------------------
    % Input validation
    % -----------------------------------------------------

    if ~isscalar(N) || ~isfinite(N) || ...
            N < 1 || N ~= floor(N)

        error('N must be a positive integer.');

    end

    if ~isscalar(L) || ~isfinite(L) || ...
            L <= 0

        error('L must be a positive scalar.');

    end

    if L > N

        error('L cannot be greater than N.');

    end


    % -----------------------------------------------------
    % Create covariance matrix
    % -----------------------------------------------------

    R = zeros(N,N);


    % -----------------------------------------------------
    % Calculate every covariance entry
    % -----------------------------------------------------

    for k1 = 1:N

        for k2 = 1:N

            if k1 == k2

                % Eq. (20), diagonal case
                R(k1,k2) = 1;

            else

                difference = k1 - k2;

                numerator = ...
                    1 - exp(-1j*2*pi*L*difference/N);

                denominator = ...
                    1j*2*pi*L*difference/N;

                R(k1,k2) = numerator / denominator;

            end

        end

    end


    % -----------------------------------------------------
    % Information
    % -----------------------------------------------------

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' Uniform Channel Covariance Matrix\n');
    fprintf('============================================\n');

    fprintf('N                : %d\n', N);
    fprintf('CIR length L     : %.4f\n', L);

    fprintf('Matrix dimensions: %d x %d\n', ...
        size(R,1), size(R,2));

    fprintf('============================================\n');

end