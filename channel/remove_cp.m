function x = remove_cp(y, G, N)

y = y(:).';

if nargin < 2 || isempty(G)
    G = 64;
end

if nargin < 3 || isempty(N)
    N = 128;
end

startIndex = G + 1;
endIndex = G + N;

if endIndex > length(y)
    error('Received signal is too short to remove CP and extract N samples.');
end

x = y(startIndex:endIndex);

end