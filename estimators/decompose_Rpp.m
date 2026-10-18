function [U, Lambda] = decompose_Rpp(Rpp)

% =========================================================
% EIGENDECOMPOSITION OF R_pp
%
% R_pp = U * Lambda * U^H
%
% The covariance matrix R_pp should theoretically be
% Hermitian positive semidefinite.
%
% Because MATLAB uses finite-precision floating-point
% arithmetic, eigenvalues that are theoretically zero can
% sometimes appear as tiny negative numbers such as:
%
%       -1e-15
%
% These are numerical round-off errors.
%
% We therefore:
%
%   1. Verify Hermitian structure
%   2. Compute eigendecomposition
%   3. Sort eigenvalues in descending order
%   4. Reject genuinely negative eigenvalues
%   5. Clamp only numerical negative values to zero
%
% =========================================================

    %% -----------------------------------------------------
    % Input validation
    % -----------------------------------------------------

    if nargin < 1

        error('Rpp must be provided.');

    end

    if ~isnumeric(Rpp) || ~ismatrix(Rpp)

        error('Rpp must be a numeric matrix.');

    end

    [Nr,Nc] = size(Rpp);

    if Nr ~= Nc

        error('Rpp must be a square matrix.');

    end

    if Nr < 1

        error('Rpp cannot be empty.');

    end


    %% -----------------------------------------------------
    % Hermitian verification
    % -----------------------------------------------------

    hermitianError = ...
        norm(Rpp-Rpp','fro');

    hermitianScale = ...
        max(1,norm(Rpp,'fro'));

    hermitianTolerance = ...
        1e-10 * hermitianScale;

    if hermitianError > hermitianTolerance

        error( ...
            ['Rpp is not Hermitian. ' ...
             'Hermitian error = %.6e'], ...
            hermitianError);

    end


    %% -----------------------------------------------------
    % Force exact Hermitian symmetry numerically
    %
    % This does NOT change the theoretical covariance model.
    % It only removes tiny floating-point asymmetry.
    % -----------------------------------------------------

    RppHermitian = ...
        (Rpp + Rpp')/2;


    %% -----------------------------------------------------
    % Eigendecomposition
    % -----------------------------------------------------

    [U,D] = eig(RppHermitian);

    eigenvalues = ...
        real(diag(D));


    %% -----------------------------------------------------
    % Numerical PSD validation
    % -----------------------------------------------------

    eigenvalueScale = ...
        max(1,max(abs(eigenvalues)));

    negativeTolerance = ...
        1e-12 * eigenvalueScale;

    minEigenvalue = ...
        min(eigenvalues);


    % Genuine negative eigenvalue
    if minEigenvalue < -negativeTolerance

        error( ...
            ['Rpp has a genuinely negative eigenvalue. ' ...
             'Minimum eigenvalue = %.6e. ' ...
             'This indicates a covariance construction problem.'], ...
            minEigenvalue);

    end


    %% -----------------------------------------------------
    % Clamp numerical negative values
    % -----------------------------------------------------
    %
    % Example:
    %
    %     -2.35e-15 -> 0
    %
    % but:
    %
    %     -0.01 -> ERROR above
    %
    % -----------------------------------------------------

    numClamped = ...
        sum(eigenvalues < 0);

    eigenvalues(eigenvalues < 0) = 0;


    %% -----------------------------------------------------
    % Sort eigenvalues from largest to smallest
    % -----------------------------------------------------

    [eigenvalues,order] = ...
        sort(eigenvalues,'descend');

    U = U(:,order);

    Lambda = ...
        diag(eigenvalues);


    %% -----------------------------------------------------
    % Final verification
    % -----------------------------------------------------

    reconstructionError = ...
        norm( ...
            RppHermitian - ...
            U*Lambda*U', ...
            'fro');


    %% -----------------------------------------------------
    % Information
    % -----------------------------------------------------

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' R_pp Eigendecomposition\n');
    fprintf('============================================\n');

    fprintf('Matrix size       : %d x %d\n',Nr,Nr);

    fprintf('Hermitian error   : %.6e\n', ...
        hermitianError);

    fprintf('Minimum raw eigenvalue : %.6e\n', ...
        minEigenvalue);

    fprintf('Clamped numerical negatives : %d\n', ...
        numClamped);

    fprintf('Reconstruction error : %.6e\n', ...
        reconstructionError);

    fprintf('\nEigenvalues:\n');

    for k = 1:length(eigenvalues)

        fprintf( ...
            'lambda(%d) = %.12e\n', ...
            k,eigenvalues(k));

    end

    fprintf('============================================\n');

end