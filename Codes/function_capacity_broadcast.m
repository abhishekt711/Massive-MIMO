function [sumCapacity,p] = function_capacity_broadcast(H,Pmax)

%Extract system dimensions
K = size(H,1); %Number of terminals
M = size(H,2); %Number of service antennas

%Solve the capacity-achieving power optimization problem using CVX, by
%solving the convex problem stated in Theorem in "Sum Capacity of the
%Vector Gaussian Broadcast Channel and Uplink?Downlink Duality" by Pramod
%Viswanath and David Tse 
cvx_begin
cvx_quiet(true); % this suppresses screen output from the solver
variable p(K);
minimize det_inv(eye(M)+H'*diag(p)*H);
subject to
    sum(p)<=Pmax
    min(p)>=0
cvx_end

%Compute the sum capacity using the optimized power allocation
sumCapacity = real(log2(det(eye(M)+H'*diag(p)*H)));
