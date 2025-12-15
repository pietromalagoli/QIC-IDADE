function [out] = full_IDA(w_tar,B_ens,A_conv,precision,niter,verb)
    % w_tar: target eigenfrequencies ALONG AN AXIS, not all 3.
    % B_ens: ensemble of intial guesses for the B matrix (NxNxP matrix)
    converged = false; % Flag for convergence
    N = numel(w_tar); % because w_tar has size = [1 N] (since we study only along one direction)
    %W2 = diag(w_tar.^2); % diagonal matrix of w_tar^2
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


