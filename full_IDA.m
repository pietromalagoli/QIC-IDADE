function [out] = full_IDA(w_tar,B_ens,A_conv,precision,niter)
    % B_ens: ensemble of intial guesses for the B matrix (NxNxP matrix)
    converged = false; % Flag for convergence
    N = numel(w_tar);
    P = size(B_ens,3);
    w_tars = repelem(w_tar,1,1,P); % size = 1 N P (I added the dummy variable (first dimension) to be able to vetorize the computation of A_bars
    A_bars = zeros(N,N,P);
    for iter = 1:niter
        for i = 1:N
            for j = i:N % Since it's symmetric, I need to compute only half of the terms
                if i == j
                    for m = 1:N
                    B_t = permute(B_ens,[2,1,3]); % Transpose of B_ens
                    A_bars(i,j,:) = A_bars(i,j,:) + B_ens(i,m,:).*(w_tars(:,m,:).^2).*B_t(m,j,:);
                    A_bars(j,i,:) = A_bars(j,i,:) + B_ens(j,m,:).*(w_tars(:,m,:).^2).*B_t(m,i,:);
                    end
                else % Step 2
                    A_bars(i,j,:) = A_conv(i,j);
                    A_bars(j,i,:) = A_conv(j,i);
                end
            end
        end
        % Step 3
        B_ens = zeros(size(B_ens));
        w_inters = zeros(size(w_tars));
        score = zeros(1,P); % store abs(w_inter(p)-w_tar(p)) per ogni p
        for p = 1:P
            [B_inter,D_inter] = eig(A_bars(:,:,p)); 
            w_inters(:,:,p) = diag(D_inter);
            B_ens(:,:,p) = B_inter;
            score(p) = sum(abs(w_inters(:,:,p)-w_tars(:,:,p)));
        end
        if min(score) <= precision
            converged = true;
            [winning_score,p_star] = min(score);
            out.A_opt = A_bars(:,:,p_star); % if converged, return the winning matrix and its score
            out.score = winning_score;
            break
        end
        B_ens = B_inter; % Set the newly found B matrix as the initial guess for B in the next cycle
    end
    if ~converged % In matlab, ~ is the symbol for 'not'
        disp('Convergence not reached for ensemble :(');
        out.A_bars = A_bars; % if not converged, return the whole set of matrices
    else
        disp(['Convergence reached after ' num2str(iter) ...
            ' iterations with score ' num2str(opt.score) ' :).']);
        out.A_bars = []; % if converged, do not return the full set of matrices
    end
end
