function preamble = preamble_sequences(Pall)
%PREAMBLE_SEQUENCES Generate IEEE 802.16 STS and LTS frequency sequences.
%
% INPUT:
%   Pall : 201-element IEEE 802.16 QPSK sequence
%          corresponding to k = -100,...,100
%
% OUTPUT:
%   preamble.k
%   preamble.Pall
%   preamble.PSTS
%   preamble.PLTS

% ---------------------------------------------------------
% 1. Subcarrier indices
% ---------------------------------------------------------
k = -100:100;

% ---------------------------------------------------------
% 2. Check the input
% ---------------------------------------------------------
if length(Pall) ~= 201
    error('Pall must contain exactly 201 elements.');
end

Pall = Pall(:).';     % Make sure it is a row vector

% ---------------------------------------------------------
% 3. Generate STS frequency-domain sequence
%
% PSTS,k = 2 Pall,k       if mod(k,4) = 0
%           0             otherwise
% ---------------------------------------------------------
PSTS = zeros(1,201);

stsmask = mod(k,4) == 0;

PSTS(stsmask) = 2 * Pall(stsmask);

% ---------------------------------------------------------
% 4. Generate LTS frequency-domain sequence
%
% PLTS,k = sqrt(2) Pall,k  if mod(k,2) = 0
%           0              otherwise
% ---------------------------------------------------------
PLTS = zeros(1,201);

ltsmask = mod(k,2) == 0;

PLTS(ltsmask) = sqrt(2) * Pall(ltsmask);

% ---------------------------------------------------------
% 5. Store results
% ---------------------------------------------------------
preamble.k    = k;
preamble.Pall = Pall;
preamble.PSTS = PSTS;
preamble.PLTS = PLTS;

% ---------------------------------------------------------
% 6. Verification
% ---------------------------------------------------------
fprintf('\n============================================\n');
fprintf(' IEEE 802.16 Preamble Sequences\n');
fprintf('============================================\n');

fprintf('Pall length : %d\n', length(Pall));

fprintf('Nonzero Pall : %d\n', nnz(Pall));

fprintf('Nonzero PSTS : %d\n', nnz(PSTS));

fprintf('Nonzero PLTS : %d\n', nnz(PLTS));

fprintf('============================================\n');

end