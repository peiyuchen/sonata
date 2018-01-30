function [ DD ] = accumulate_matrix( cost )
% input  : cost matrix
% output : accumulate matrix

    cost = 1 - cost;
    
    N = size(cost,1);
    M = size(cost,2);

    DD = zeros(size(cost));
    DD(1,1) = cost(1,1);

    for n = 2:N
        DD(n,1) = cost(n,1) + DD(n-1,1);
    end
    for m = 2:M
        DD(1,m)=cost(1,m) + DD(1,m-1);
    end
    for n=2:N
        for m=2:M
            DD(n,m)=cost(n,m)+min([DD(n-1,m),DD(n-1,m-1),DD(n,m-1)]);
        end
    end

end

