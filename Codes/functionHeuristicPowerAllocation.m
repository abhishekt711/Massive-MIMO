function powerallocation = functionHeuristicPowerAllocation(rhos,q,weights)



Kt = size(rhos,1); %Number of base stations (BSs)
Kr = size(rhos,2); %Number of users (in total)


%Pre-allocation of matrix for power allocation coefficients
powerallocation=size(Kt,Kr);


%Iteration over base stations to perform power allocation
for j = 1:Kt
    indicesOfNonzero = find(rhos(j,:)>0); %Find which users that are served by BS j
    
    %Case 1: Compute waterlevel if all of the users served by BS j are
    %allocated non-zero power.
    nuAllActive = (q(j)+sum(1./rhos(j,indicesOfNonzero)))/sum(weights(indicesOfNonzero));
    
    %Case 2: Compute waterlevel if only a subset of the users served by BS
    %j are allocated non-zero power. The range of the waterlevel is
    %achieved by checking when there is equality in (3.37); that is, when
    %users are activated. The fminbnd-algorithm finds the waterlevel that
    %minimize the difference between the allocated power and available power.
    nuRangeLower = min(1./(rhos(j,indicesOfNonzero)'.*weights(indicesOfNonzero)));
    nuRangeUpper = max(1./(rhos(j,indicesOfNonzero)'.*weights(indicesOfNonzero)));
    nu = fminbnd(@(x) functionAllocDiff(x,q(j),rhos(j,indicesOfNonzero)',weights(indicesOfNonzero)),nuRangeLower,nuRangeUpper);
    
    %Check if the difference between the allocated power and the available
    %power is minimized by allocating power to all users or only subset.
    if functionAllocDiff(nu,q(j),rhos(j,indicesOfNonzero)',weights(indicesOfNonzero)) < functionAllocDiff(nuAllActive,q(j),rhos(j,indicesOfNonzero)',weights(indicesOfNonzero))
        %Compute power allocation with optimal waterlevel (only a subset of users are active)
        powerallocation(j,indicesOfNonzero) = max([weights(indicesOfNonzero)*nu-1./rhos(j,indicesOfNonzero)' zeros(length(indicesOfNonzero),1)],[],2);
    else
        %Compute power allocation with optimal waterlevel (all users are active)
        powerallocation(j,indicesOfNonzero) = max([weights(indicesOfNonzero)*nuAllActive-1./rhos(j,indicesOfNonzero)' zeros(length(indicesOfNonzero),1)],[],2);
    end
    
    %Scale the power allocation to use full power (to improve numerical accuracy)
    powerallocation(j,:) = q(j)*powerallocation(j,:)/sum(powerallocation(j,:));
end


function difference = functionAllocDiff(nu,q,rhos,weights)
%Computes the power allocation of (3.37) for a given waterlevel and returns
%the absolute difference between the total allocated power and the total
%available power.
%
%INPUT:
%nu      = Waterlevel
%q       = Total available power
%rhos    = K x 1 vector with effective channel gains
%weights = K x 1 vector with positive weights for each user
%
%OUTPUT:
%difference = Absolute difference between allocated and available power

difference = abs( sum( max([nu*weights-1./rhos zeros(size(weights))],[],2) ) - q);
