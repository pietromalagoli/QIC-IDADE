function [A_conv,B_init,w_init] = AB_conv(N,init,params)
    % This function computes the A_conv and B_init matrices, as well as the 
    % initial guess for the eigenfrequencies (w_init) of a system of ions 
    % influenced by electric fields, without the presence of optical tweezers.
    %
    % Inputs:
    %   N      - Integer representing the number of ions in the system.
    %   init   - Struct containing initialization values:
    %             - q_init: Initial charge of the ions (vector of size N).
    %             - M_init: Initial mass of the ions (vector of size N).
    %   params - Struct containing system parameters:
    %             - eps0: Permittivity of free space.
    %             - Vdc: Direct current voltage.
    %             - Vrf: Radio frequency voltage.
    %             - OmegaRF: Angular frequency of the radio frequency.
    %             - xi: Vector of parameters related to the trap.
    %             - psi: Vector of parameters related to the trap.
    %
    % Outputs:
    %   A_conv - 3D matrix of size [N N 3] representing the A matrix 
    %             without the tweezers, with one NxN matrix for each 
    %             direction (x, y, z).
    %   B_init  - 3D matrix of size [N N 3] representing the initial guess 
    %             for the B matrix of the system without tweezers, with one 
    %             NxN matrix for each direction (x, y, z).
    %   w_init  - 2D matrix of size [N 3] containing the eigenfrequencies 
    %             without tweezers, with one column for each direction 
    %             (x, y, z).

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
    for i = 1:N-1 
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
    coords = [r(1,:);r(2,:); r(3,:)];
    gradPhi = gradient(phiTotal,coords(:));
    
    % Substitute the variables' initial values in the gradient expression
    init_list = [q,M];
    param_inits = [q_init,M_init];
    gradPhi = subs(gradPhi,init_list, param_inits); 
    
    % Find the critical points
    crits = vpasolve(gradPhi == 0.,'Random',true); 
    crits_mat = cell2mat(struct2cell(crits))';
    crits_mat(abs(crits_mat) < 1e-25) = 0.; % I set a tolerance for the values of crits, setting to zero too low values (This allows to regain some stability)
    
    % Now I compute the Hessian martix
    dim = size(coords,1);
    Hs = arrayfun(@(alpha) hessian(phiTotal,coords(alpha,:)), 1:dim,'UniformOutput',false); % compute the Hessian for each direction (x,y,z)
    Hs = cat(3,Hs{:});
    Hs = subs(Hs,init_list,param_inits); % substitute the parameters' values in the Hessian matrix
    Hs_eq = subs(Hs,symvar(Hs),crits_mat); % compute the value of the Hessian matrix in the critical point previously found
    Hs_eq(isnan(Hs_eq)) = 0.; % Substitute NaN values with zeros in the Hessian matrix

    % Find the normal modes and frequencies
    Minv2 = diag(1 ./ sqrt(M_init));
    A_conv = arrayfun(@(alpha) Minv2 * Hs_eq(:,:,alpha) * Minv2, 1:dim, 'UniformOutput', false);
    A_conv = cat(3,A_conv{:});
    
    % Find B and w by diagonalizing A
    [B_init,w_init] = arrayfun(@(alpha) eig(A_conv(:,:,alpha)), 1:dim, 'UniformOutput', false); 
    B_init = cat(3,B_init{:});
    w_init = arrayfun(@(alpha) diag(w_init{alpha}), 1:dim, 'UniformOutput', false);
    w_init = cat(2,w_init{:});
end