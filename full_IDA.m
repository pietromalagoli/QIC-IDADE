function [out] = full_IDA(w_tar,B_ens,A_conv,precision,niter,verb)
    % full_IDA computes the optimal A matrix for a given set of target eigenfrequencies.
    %
    % Inputs:
    %   w_tar     - A vector of target eigenfrequencies along a specified axis (1 x N).
    %   B_ens     - A 3D matrix containing an ensemble of initial guesses for the B matrix (N x N x P).
    %   A_conv    - A matrix used to construct a modified version of the A matrix (N x N).
    %   precision  - A scalar value indicating the convergence threshold for the optimization.
    %   niter     - An integer specifying the maximum number of iterations for the optimization process.
    %   verb      - A boolean flag to control the verbosity of the output during iterations.
    %
    % Outputs:
    %   out       - A structure containing the following fields:
    %       A_opt   - The optimal A matrix found during the optimization (N x N).
    %       score    - The score associated with the optimal A matrix.
    %       A_bars   - The set of A matrices computed during the iterations (if not converged).
    %
    % The function iteratively refines the B matrix and computes the corresponding A matrices
    % until convergence is reached or the maximum number of iterations is exceeded.
    % Convergence is determined based on the minimum score calculated from the difference
    % between the computed eigenfrequencies and the target eigenfrequencies.

    converged = false; % Flag for convergence
    N = numel(w_tar); % because w_tar has size = [1 N] (since we study only along one direction)
    P = size(B_ens,3);
    w_tars = repelem(w_tar',1,1,P); % size = 1 N P (I added the dummy variable (first dimension) to be able to vetorize the computation of A_bars
    A_bars = zeros(N,N,P);
    w_conv = [2*pi*1.2e6,2*pi*1e6,2*pi*0.2e6]; % convential trap frequencies for Yb from figure 4
    score = zeros(1,P); % store abs(w_inter(p)-w_tar(p)) per ogni p
    %w_inters = zeros(size(w_tars));
    w_inters = w_tars;
    for iter = 1:niter
        mask_conv = ~eye(N) * A_conv; % construct a copy of A_conv with zeros on the diagonal
        A_bars = arrayfun(@(p) (B_ens(:,:,p)*diag(w_inters(:,:,p))*B_ens(:,:,p)') .* eye(N) + mask_conv, 1:P, 'UniformOutput', false);
        A_bars = cat(3,A_bars{:});
            
        % Step 3
        for p = 1:P
            [B_inter,D_inter] = eig(A_bars(:,:,p)); 
            w_inters(:,:,p) = diag(D_inter);
            B_ens(:,:,p) = B_inter;
            score(p) = sum(abs(w_inters(:,:,p)./w_conv(3)-w_tars(:,:,p)./w_conv(3)));
        end
        if min(score) <= precision
            converged = true;
            [winning_score,p_star] = min(score);
            out.A_opt = A_bars(:,:,p_star); % if converged, return the winning matrix and its score
            out.score = winning_score;
            break
        end
        if rem(iter,10) == 0 && verb
            disp(['Iteration # ' num2str(iter) ' with score = ']);
            disp(min(score)');
        end
    end
    if ~converged % In matlab, ~ is the symbol for 'not'
        disp(['Convergence not reached for ensemble. Minimum score: ' num2str(min(score))]);
        out.A_bars = A_bars; % if not converged, return the whole set of matrices
        out.score = -1; % used as a flag for not convergence
    else
        disp(['Convergence reached after ' num2str(iter) ...
            ' iterations with score ' num2str(opt.score) ' :).']);
        out.A_bars = []; % if converged, do not return the full set of matrices
    end
end


