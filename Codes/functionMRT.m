function wMRT = functionMRT(H,D)




%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<2
    D = repmat( eye(N), [1 1 Kr]);
end

%Pre-allocation of MRT beamforming
wMRT = zeros(size(H'));

%Computation of MRT, based on Definition 3.2
for k = 1:Kr
    channelvector = (H(k,:)*D(:,:,k))'; %Useful channel
    wMRT(:,k) = channelvector/norm(channelvector); %Normalization of useful channel
end
