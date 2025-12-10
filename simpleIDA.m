% Implementation of the IDA algorithm from the paper
function A_bar = simpleIDA(w_tar,B_init,A_conv,precision,niter,verb)
    % I follow the protocol at page 4 of the paper
    % Step 1 & 2
    cN = numel(w_tar); % because w_tar has size = [1 3*N] (since we have 3 coordinates, x,y,z)
    N = cN/3;
    A_bar = zeros(cN);
    for iter = 1:niter
        for i = 1:cN
            for j = i:cN % Since it's symmetric, I need to compute only half of the terms
                if i == j
                    for m = 1:cN
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
        w_conv = repelem([2*pi*1.2e6,2*pi*1e6,2*pi*0.2e6],N);
        score = abs(w_inter(1:N)./w_conv(1) - w_tar(1:N)./w_conv(1)); % score only along x
        if rem(iter,10) == 0 && verb
            disp(['Iteration # ' num2str(iter) ' with score = ']);
            disp(score');
        end
        if score <= precision
            break
        end
        B_init = B_inter; % Set the newly found B matrix as the initial guess for B in the next cycle
    end
    disp('Convergence not reached :(');
end