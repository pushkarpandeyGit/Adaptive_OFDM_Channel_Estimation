function Dp = build_Dp(eigenvalues,pRank,beta,gamma_p)
%BUILD_DP Construct the SVD filtering matrix D_p.
%
% Paper equation:
%
%   delta_k = 1 / (lambda_k + beta/gamma_p),   k = 1,...,p
%   delta_k = 0,                               otherwise
%
% Inputs:
%   eigenvalues : eigenvalues of R_pp
%   pRank       : number of retained SVD modes
%   beta        : regularization parameter
%   gamma_p     : pilot SNR in linear scale
%
% Output:
%   Dp          : NxN diagonal filtering matrix

    % ---------------------------------------------------------
    % STEP 1: Make sure eigenvalues are a vector
    % ---------------------------------------------------------
    %
    % If a diagonal matrix is accidentally supplied, extract
    % its diagonal first.
    %
    if ismatrix(eigenvalues) && size(eigenvalues,1) == size(eigenvalues,2)
        eigenvalues = diag(eigenvalues);
    end

    % Force column-vector form
    eigenvalues = eigenvalues(:);

    % Number of eigenvalues
    N = numel(eigenvalues);

    % ---------------------------------------------------------
    % STEP 2: Validate inputs
    % ---------------------------------------------------------

    if N == 0
        error('build_Dp: eigenvalues cannot be empty.');
    end

    if ~isscalar(pRank) || ~isnumeric(pRank)
        error('build_Dp: pRank must be a numeric scalar.');
    end

    pRank = round(pRank);

    if pRank < 1 || pRank > N
        error('build_Dp: pRank must satisfy 1 <= pRank <= N.');
    end

    if ~isscalar(beta) || ~isnumeric(beta)
        error('build_Dp: beta must be a numeric scalar.');
    end

    if ~isscalar(gamma_p) || ~isnumeric(gamma_p) || gamma_p <= 0
        error('build_Dp: gamma_p must be a positive numeric scalar.');
    end

    % ---------------------------------------------------------
    % STEP 3: Regularization term
    % ---------------------------------------------------------

    regularization = beta / gamma_p;

    % ---------------------------------------------------------
    % STEP 4: Create NxN diagonal matrix
    % ---------------------------------------------------------

    Dp = zeros(N,N);

    % ---------------------------------------------------------
    % STEP 5: Keep only first pRank eigenmodes
    % ---------------------------------------------------------

    for k = 1:pRank

        denominator = eigenvalues(k) + regularization;

        % Avoid division by an effectively zero denominator
        if abs(denominator) < eps
            error('build_Dp: denominator is too close to zero.');
        end

        Dp(k,k) = 1 / denominator;

    end

    % ---------------------------------------------------------
    % STEP 6: Verification
    % ---------------------------------------------------------

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' SVD Filtering Matrix Dp\n');
    fprintf('============================================\n');
    fprintf('Number of eigenvalues : %d\n',N);
    fprintf('Selected rank         : %d\n',pRank);
    fprintf('beta                  : %.6e\n',beta);
    fprintf('gamma_p               : %.6e\n',gamma_p);
    fprintf('regularization        : %.6e\n',regularization);
    fprintf('Dp dimensions         : %d x %d\n',size(Dp,1),size(Dp,2));
    fprintf('============================================\n');

end