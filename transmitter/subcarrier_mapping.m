function mapping = subcarrier_mapping(p)
%SUBCARRIER_MAPPING Creates the IEEE 802.16 OFDM subcarrier mapping.
%
% Output:
%   mapping.usedIndex   -> 200 used subcarriers
%   mapping.dataIndex   -> 192 data subcarriers
%   mapping.pilotIndex  -> 8 pilot subcarriers
%   mapping.nullIndex   -> 56 null subcarriers

% ---------------------------------------------------------
% 1. All subcarrier positions
% ---------------------------------------------------------
subcarrierIndex = -p.N/2 : p.N/2-1;

% ---------------------------------------------------------
% 2. Pilot locations specified by the paper
% ---------------------------------------------------------
pilotLocations = p.pilotLocations;

% ---------------------------------------------------------
% 3. Used subcarriers
%
% The IEEE 802.16 setup uses 200 used subcarriers.
% Zero/DC is not a used subcarrier.
% ---------------------------------------------------------
usedIndex = subcarrierIndex( ...
    subcarrierIndex >= -100 & ...
    subcarrierIndex <= 100 & ...
    subcarrierIndex ~= 0);

% ---------------------------------------------------------
% 4. Pilot subcarriers
% ---------------------------------------------------------
pilotIndex = find(ismember(usedIndex, pilotLocations));

% ---------------------------------------------------------
% 5. Data subcarriers
%
% Everything used that is NOT a pilot is data.
% ---------------------------------------------------------
dataIndex = find(~ismember(usedIndex, pilotLocations));

% ---------------------------------------------------------
% 6. Null subcarriers
% ---------------------------------------------------------
nullIndex = find(~ismember(subcarrierIndex, usedIndex));

% ---------------------------------------------------------
% 7. Store results
% ---------------------------------------------------------
mapping.subcarrierIndex = subcarrierIndex;
mapping.usedIndex      = usedIndex;
mapping.pilotIndex     = pilotIndex;
mapping.dataIndex      = dataIndex;
mapping.nullIndex      = nullIndex;

% Actual subcarrier locations
mapping.pilotLocations = usedIndex(pilotIndex);
mapping.dataLocations  = usedIndex(dataIndex);
mapping.nullLocations  = subcarrierIndex(nullIndex);

% ---------------------------------------------------------
% 8. Verification
% ---------------------------------------------------------
assert(length(usedIndex) == p.Nu, ...
    'Number of used subcarriers is incorrect.');

assert(length(dataIndex) == p.Nd, ...
    'Number of data subcarriers is incorrect.');

assert(length(pilotIndex) == p.Np, ...
    'Number of pilot subcarriers is incorrect.');

assert(length(nullIndex) == p.N - p.Nu, ...
    'Number of null subcarriers is incorrect.');

fprintf('\n============================================\n');
fprintf(' OFDM Subcarrier Mapping\n');
fprintf('============================================\n');
fprintf('Total subcarriers : %d\n', p.N);
fprintf('Used subcarriers  : %d\n', length(usedIndex));
fprintf('Data subcarriers  : %d\n', length(dataIndex));
fprintf('Pilot subcarriers : %d\n', length(pilotIndex));
fprintf('Null subcarriers  : %d\n', length(nullIndex));
fprintf('============================================\n');
end