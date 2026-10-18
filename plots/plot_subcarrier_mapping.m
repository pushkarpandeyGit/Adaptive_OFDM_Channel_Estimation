function plot_subcarrier_mapping(mapping)
%PLOT_SUBCARRIER_MAPPING Visualize OFDM subcarrier allocation.

figure;

% Start with all subcarriers as null
allocation = zeros(size(mapping.subcarrierIndex));

% Data = 1
allocation(ismember(mapping.subcarrierIndex, ...
    mapping.dataLocations)) = 1;

% Pilot = 2
allocation(ismember(mapping.subcarrierIndex, ...
    mapping.pilotLocations)) = 2;

% Plot
stem(mapping.subcarrierIndex, allocation, 'filled');

grid on;

xlabel('Subcarrier Index');
ylabel('Type');

title('IEEE 802.16 OFDM Subcarrier Mapping');

yticks([0 1 2]);
yticklabels({'Null', 'Data', 'Pilot'});

xlim([-128 127]);
end