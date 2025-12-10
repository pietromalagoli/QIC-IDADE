function [A_conv,B_init,w_init] = AB_conv(N,init,params)
    % "Function for computing the A_conv and B_conv matrices for the given
    % parameters."
    % Arguments:
    % - init: struct with initialization values for charge and mass,
    % i.e. q,M, respectively named q_init,M_init
    % - params: struct with the parameters for the system, such as eps0,Vdc, ...
    % Returns:
    % - A_conv: A matrix w/o the tweezers (for the off-diagonal constraints)
    % - B_init: Initial guess for the B matrix, as the B matrix of the
    % system w/o tweezers.
    % - w_init: Eigenfrequencies w/o tweezers, useful if one wants to focus
    % IDADE's work on a specific direction (coordinate).

    % Extract initialization values and parameters
    q_init = init.q_init;
    M_init = init.M_init;
    eps0 = params.eps0;
    Vdc = params.Vdc;
    Vrf = params.Vrf;
    OmegaRF = params.OmegaRF;
    xi = params.xi;
    psi = params.psi;
    
    syms q M [1 N] real; % electronic charge and mass (respectively) of each ion
    assume(M>0); % set the condition for positive mass
    syms r [3 N] real; % position vector of each ion
    
    % Define the potentials
    phiC = 0; % Coulomb potential
    for i = 1:N-1 % devo fare 1:N-1 perché per i=N dopo j=i+1:N è vuoto e quindi dà errore. e infatti matematicamente il termine della somma per i = N si ski
        for j = i+1:N
            phiC = phiC + q(i) * q(j) / (4*pi*eps0*norm(r(:,i)-r(:,j)));
        end
    end
    
    % Trap frequencies
    syms wTx(M_i,q_i); 
    syms wTy(M_i,q_i);
    syms wTz(M_i,q_i);
    
    % Define the trap frequencies
    wTx(M_i,q_i) = sqrt(-q_i*Vdc*xi(2)/M_i + q_i^2*Vrf^2*psi(2)^2/(2*M_i^2*OmegaRF^2));
    wTy(M_i,q_i) = sqrt(-q_i*Vdc*xi(1)/M_i + q_i^2*Vrf^2*psi(1)^2/(2*M_i^2*OmegaRF^2));
    wTz(M_i,q_i) = sqrt(q_i*Vdc*xi(3)/M_i);

    phiT = 0; % Trap potential (in the paper, phi_conv)
    for i =  1:N % I directly compute the sum of the potential over the i index
        phiT = phiT + 0.5*M(i)*((wTx(M(i),q(i))^2)*r(1,i)^2 + ...
                                 (wTy(M(i),q(i))^2)*r(2,i)^2 + ...
                                 (wTz(M(i),q(i))^2)*r(3,i)^2);
    end
    phiTotal = phiC + phiT; % Total potential (without tweezers)
    
    %% Stability Analysis
    coords = [r(1,:) r(2,:) r(3,:)];
    gradPhi = gradient(phiTotal,coords);
    
    % Substitute the variables' initial values in the gradient expression
    init_list = [q,M];
    param_inits = [q_init,M_init];
    gradPhi = subs(gradPhi,init_list, param_inits); 
    
    % Find the critical points
    crits = vpasolve(gradPhi == 0.,'Random',true); 
    crits_mat = cell2mat(struct2cell(crits));
    crits_mat(abs(crits_mat) < 1e-25) = 0.; % I set a tolerance for the values of crits, setting to zero too low values (This allows to regain some stability)
    
    % Now I compute the Hessian martix
    H = hessian(phiTotal,coords);
    H = subs(H,init_list, param_inits); % Substitute the parameters' values in the Hessian expression
    H_eq = subs(H,symvar(H),crits_mat'); % compute the value of the Hessian matrix in the critical point previously found
    H_eq(isnan(H_eq)) = 0.; % Substitute NaN values with zeros in the Hessian matrix

    % Find the normal modes and frequencies
    M = repelem(cell2mat(sym2cell(M)), 3);  % each mass appears 3 times because there are 3 variables for each ion (x,y,z)
    Minv2 = diag(1 ./ sqrt(M));
    A_conv = Minv2 * H_eq * Minv2; % A matrix as the mass-weighted Hessian (matrix)
    A_conv = subs(A_conv,M,repelem(M_init,3)); % subtitutes the values of the masses in the expression
    
    % Find B and w by diagonalizing A
    [B_init,W] = eig(A_conv); 
    w_init = diag(W); 
end