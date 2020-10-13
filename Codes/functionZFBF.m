function wZFBF = functionZFBF(H,D)

%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<2
    D = repmat( eye(N), [1 1 Kr]);
end

%Pre-allocation of MRT beamforming
wZFBF = zeros(size(H'));

%Computation of ZFBF, based on Definition 3.4
for k = 1:Kr
    effectivechannel = (H*D(:,:,k))'; %Effective channels
    channelinversion = effectivechannel/(effectivechannel'*effectivechannel); %Compute zero-forcing based on channel inversion
    wZFBF(:,k) = channelinversion(:,k)/norm(channelinversion(:,k));  %Normalization of zero-forcing direction
end
