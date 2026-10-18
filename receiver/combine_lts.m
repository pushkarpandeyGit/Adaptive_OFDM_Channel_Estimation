function lts = combine_lts(lts1, lts2)

% =========================================================
% Combine the two repeated LTS symbols
% =========================================================

lts1 = lts1(:).';
lts2 = lts2(:).';

% Average the two received LTS repetitions
lts = (lts1 + lts2) / 2;

fprintf('\n');
fprintf('============================================\n');
fprintf(' LTS Combination\n');
fprintf('============================================\n');

fprintf('LTS 1 length : %d samples\n', length(lts1));
fprintf('LTS 2 length : %d samples\n', length(lts2));
fprintf('Combined LTS : %d samples\n', length(lts));

fprintf('============================================\n');

end