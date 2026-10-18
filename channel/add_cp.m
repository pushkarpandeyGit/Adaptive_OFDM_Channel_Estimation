function x_cp = add_cp(x, G)

x = x(:).';

if nargin < 2 || isempty(G)
    G = 64;
end

if G > length(x)
    error('CP length G cannot exceed the input length.');
end

cp = x(end-G+1:end);

x_cp = [cp, x];

end