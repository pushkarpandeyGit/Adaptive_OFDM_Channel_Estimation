function R = uniform_covariance(N,L)
%UNIFORM_COVARIANCE Paper Eq. (20) uniform delay power profile covariance.
k = (-N/2:N/2-1).';
d = k-k.';
R = ones(N);
idx = (d~=0);
R(idx) = (1-exp(-1j*2*pi*L*d(idx)/N)) ./ (1j*2*pi*L*d(idx)/N);
R = (R+R')/2;
end
