% Implementation of the IDA algorithm from the paper
function A_bar = simpleIDA(w_tar,B_init,A_conv,precision,niter)
    % I follow the protocol at page 4 of the paper
    % Step 1 & 2
    N = numel(w_tar);
    A_bar = zeros(N);
    for iter = 1:niter
        for i = 1:N
            for j = i:N % Since it's symmetric, I need to compute only half of the terms
                if i == j
                    for m = 1:N
                    A_bar(i,j) = A_bar(i,j) + B_init(i,m)*(w_tar(m)^2)*B_init(m,j)';
                    A_bar(j,i) = A_bar(i,j) + B_init(i,m)*(w_tar(m)^2)*B_init(m,j)';
                    end
                else % Step 2
                    A_bar(i,j) = A_conv(i,j);
                    A_bar(j,i) = A_conv(j,i);
                end
            end
        end
        % Step 3
        [B_inter,D_inter] = eig(A_bar); 
        w_inter = diag(D_inter);
        if abs(w_inter - w_tar) <= precision
            break
        end
        B_init = B_inter; % Set the newly found B matrix as the initial guess for B in the next cycle
    end
    disp('Convergence not reached :(')
end