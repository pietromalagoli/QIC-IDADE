function A_bar = simpleIDA(w_tar,B_init,A_conv,precision,niter,verb)
% simpleIDA is an iterative method for computing a modified matrix A_bar
%
% Inputs:
%   w_tar     - Target frequencies (1 x 3*N vector)
%   B_init    - Initial guess for matrix B (cN x cN matrix)
%   A_conv    - Convergence matrix (cN x cN matrix)
%   precision - Convergence criterion (scalar)
%   niter     - Maximum number of iterations (scalar)
%   verb      - Verbosity flag (boolean), if true, displays iteration info
%
% Outputs:
%   A_bar     - Computed modified matrix (cN x cN matrix)
%
% Description:
%   This function implements a simple iterative algorithm to compute a 
%   modified matrix A_bar based on the input target frequencies w_tar, 
%   an initial guess for matrix B_init, and a convergence matrix A_conv. 
%   The algorithm iteratively updates the matrix A_bar until the score 
%   based on the target frequencies converges to a specified precision 
%   or the maximum number of iterations is reached. The function also 
%   provides optional verbosity to display progress at every 10th iteration.

    % Step 1 & 2
    N = numel(w_tar); 
    N = N/3;
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