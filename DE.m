function A_bars = DE(A_bars,B_ens,w_tar,A_conv,DEparams,verb)
    % This function implements the Differential Evolution (DE) algorithm for optimizing a set of matrices.
    % 
    % Inputs:
    %   A_bars   - A 3D matrix of size (N, N, P) representing P matrices to be optimized.
    %   B_ens    - A 3D matrix of size (N, N, P) representing the ensemble of matrices used in the optimization.
    %   w_tar    - A 1D array of size (1, P) containing the target eigenvalues for each matrix.
    %   A_conv   - A 2D matrix of size (N, N) representing the conventional matrix that influences the optimization.
    %   DEparams  - A vector containing three parameters for the DE algorithm: 
    %               DEparams(1) - xi (bias), 
    %               DEparams(2) - k (randomizer), 
    %               DEparams(3) - eta (mutation probability).
    %   verb     - A boolean flag indicating whether to display verbose output (true) or not (false).
    %
    % Outputs:
    %   A_bars   - The optimized set of matrices after applying the DE algorithm.
    %
    % The function evaluates the eigenvalues of the input matrices, computes a score based on how close 
    % the eigenvalues are to the target values, and iteratively mutates the matrices to minimize this score.
    % The process continues until the matrices converge towards the target eigenvalues.

    [xi, k, eta] = deal(DEparams(1), DEparams(2), DEparams(3)); % mutation parameters: xi=bias, k=randomizer, eta=mutation probability
    P = size(A_bars,3);
    N = numel(w_tar); % because w_tar has size = [1 N] (since we study only along one direction)
    w_tars = repelem(w_tar,1,P);
    w = zeros(size(w_tars));
    w_primes = zeros(size(w_tars));
    x = zeros(N,P); 
    x_m = zeros(size(x)); % mutated x
    x_prime = zeros(size(x)); % x after mutation
    A_primes = zeros(size(A_bars)); % ensemble of mutated matrices
    score = zeros(1,P);
    new_score = zeros(1,P);
    w_conv = [2*pi*1.2e6,2*pi*1e6,2*pi*0.2e6]; % convential trap frequencies for Yb from figure 4
    for p = 1:P 
        % Step 1
        w(:,p) = eig(A_bars(:,:,p)); 
        x(:,p) = diag(A_bars(:,:,p)); % in x we store the diagonal elements of A, because they're the only ones the tweezer can have an effect on 
        score(p) = sum(abs(w(:,p)./w_conv(3)-w_tars(:,p)./w_conv(3))); % score to asses how close one is to the target 
    end
    % Step 2
    [~, p_star] = min(score);
    % Step 3
    for p = 1:P % P must be at least 3
        indices = setdiff(1:P, p); % Exclude the current index p
        selected = randsample(indices, 2); % Randomly sample 2 indices from the remaining
        q = selected(1);
        r = selected(2);
        x_m(:,p) = x(:,p) + xi*(x(:,p_star)-x(:,p)) + k*(x(:,q)-x(:,r));
        mutation_mask = randsample([0,1],N,true,[1-eta,eta]); % mutate each element with probability eta
        x_prime(:,p) = x(:,p); % the elementes that are not mutated are kept 
        x_prime(logical(mutation_mask),p) = x_m(logical(mutation_mask),p); % mutate the extracted elements
        % Step 4
        A_primes(:,:,p) = diag(x_prime(:,p)) + A_conv.*~eye(N); % here i construct A' as having x' on the diagonal and A_conv elsewhere since the tweezer does not affect the off-diagonal elements
        w_primes(:,p) = eig(A_primes(:,:,p)); % compute the mutated eigenvalues
        new_score(p) = sum(abs(w_primes(:,p)./w_conv(3) - w_tars(:,p)./w_conv(3))); % calculate the new score for the mutated matrix
        if new_score(p) < score(p) % if w'(p) are closer to the target, we return A'(p), else A(p)
            A_bars(:,:,p) = A_primes(:,:,p);
            disp(['Evolution converged for element ' num2str(p) ' with score = ']);
            disp(new_score(p));
        end
    end
    if verb
        disp('Old score:' );
        disp(score);
        disp('New score:');
        disp(new_score);
    end
end