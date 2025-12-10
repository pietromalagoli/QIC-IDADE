%% Note: all this analysis is done on 171Yb+ ions
clear
%{
digits(4)
sympref('FloatingPointOutput',true);
%}

N = 4; % number of ions in the chain (tpar: trap parameters)
%eps0 = 8.85418782e-12; % permittivity of free space ([e0] = m-3 kg-1 s4 A2)

syms q M [1 N] real; % electronic charge and mass (respectively) of each ion
assume(M>0); % set the condition for positive mass
syms r [3 N] real; % position vector of each ion
syms eps0 Vdc Vrf OmegaRF real positive; % permittivity of free space, trap's DC and rf peak voltages (respectively)
syms xi psi [1 3] real positive; % Xi (DC electrodes distance?) and psy (rf electrodes distance?) trap geometric factors

eps0 = sym('epsilon_0'); % This ensures that eps0 is displayed as the greek letter
Vdc = str2sym('V_DC');
Vrf = str2sym('V_RF');
OmegaRF = str2sym('Omega_RF');

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
%% 
coords = [r(1,:);r(2,:); r(3,:)];
gradPhi = gradient(phiTotal,coords(:));
% Set values for parameters
q_init = repelem(1e-19,N); % 1 unit of electric charge in Coulomb
M_init = repelem(1e-25,N); % here I used the fact that the mass of Yb atoms (whose ions are used in an example in the paper have mass ~173 u)
xi_init = [0.06, 0.06, 0.06]; % Example values for trap geometric factors
psi_init = [0.06, 0.06]; % Example values for trap geometric factors
Vdc_init = 10; % in volts (reference: Table 1 at http://dx.doi.org/10.1088/0143-0807/34/3/787) 
Vrf_init = 550; % in volts (reference: Table 1 at http://dx.doi.org/10.1088/0143-0807/34/3/787)
OmegaRF_init = 1.8e7; % Frequency in Hz (reference: Table 1 at http://dx.doi.org/10.1088/0143-0807/34/3/787)
eps0_init = 8.85e-12; % permittivity of free space ([e0] = m-3 kg-1 s4 A2)

param_list = [M,sym('Omega_RF'),sym("V_DC"),sym("V_RF"),sym("epsilon_0"),psi1,psi2,q,xi1,xi2,xi3];
param_inits = [M_init,OmegaRF_init,Vdc_init,Vrf_init,eps0_init,psi_init,q_init,xi_init];
gradPhi = subs(gradPhi,param_list, param_inits); % Substitute the parameters' values in the gradient expression
%% 
crits = vpasolve(gradPhi == 0.,'Random',true);
% I set a tolerance for the values of crits, setting to zero too low values
% (This allows to regain some stability)
crits_mat = cell2mat(struct2cell(crits))';
crits_mat(abs(crits_mat) < 1e-25) = 0.;

%% Now I compute the Hessian matrix for each direction
H = hessian(phiTotal,coords(:));
dim = size(coords,1);
Hs = arrayfun(@(alpha) hessian(phiTotal,coords(alpha,:)), 1:dim,'UniformOutput',false); % compute the Hessian for each direction (x,y,z)
Hs = cat(3,Hs{:});
Hs = subs(Hs,param_list,param_inits); % substitute the parameters' values in the Hessian matrix
Hs_eq = subs(Hs,symvar(Hs),crits_mat); % compute the value of the Hessian matrix in the critical point previously found
Hs_eq(isnan(Hs_eq)) = 0.; % Substitute NaN values with zeros in the Hessian matrix

%% Study the eigenvalues of H_eq for stability and modes
e = arrayfun(@(alpha) eig(Hs_eq(:,:,alpha)), 1:dim,'UniformOutput',false);
e = cat(2,e{:})';
% Find the normal modes and frequencies
Minv2 = diag(1 ./ sqrt(M_init));
As = arrayfun(@(alpha) Minv2 * Hs_eq(:,:,alpha) * Minv2, 1:dim, 'UniformOutput', false);
As = cat(3,As{:});
% ora poi A va diagonalizzata e trovato eigenvalues and eigenvectors
[Bs,Ws] = arrayfun(@(alpha) eig(As(:,:,alpha)), 1:dim, 'UniformOutput', false); 
Bs = cat(3,Bs{:});
Ws = arrayfun(@(alpha) diag(Ws{alpha}), 1:dim, 'UniformOutput', false);
Ws = cat(2,Ws{:})';

%% Now we add the tweezers to system (and therefore the associated term to the potential)
% Some parameters, taken from Figure 5 and table 1 (first column) in the paper
N_excited = 5; % number of excited state (in the case of 171Yb+ is 5)
c = 2.998e8; % speed of light in m/s
syms omega_ai Gamma [N N_excited] real positive; % atomic transition frequencies and scattering rates (Einstein A coefficients) of the transition from ground state to each excited state a for each ion
syms omega_li P_i lambda_i real positive; % tweezer frequency, power and wavelength on each ion
syms r_star r [3 N] real; % equilibrium positions
syms sigma0 real positive; % beam waist
syms x y z real; % dummy variables

% Define the tweezer intensity profile 
sigma_x = sigma0*sqrt(1+(x/(pi*(sigma0^2)/lambda_i))^2);
I_i = (2*P_i/(pi*sigma0^2))*((sigma0^2)/sigma_x)^2*exp(-2*(y^2+z^2)/sigma_x^2);

% Compute the complete form of the optical potential (tweezer)
phiTz = sym(0);
for i = 1:N % sum over the ions
    x_i = r(1,i) - r_star(1,i);
    y_i = r(2,i) - r_star(2,i);
    z_i = r(3,i) - r_star(3,i);

    I_sub = subs(I_i, ...   % substitute the values in the intensity profile expression
        [x,   y,   z], ...
        [x_i, y_i, z_i]);
    sum_a = 0.;
    for a = 1:N_excited % sum over the excited states of each ion
        sum_a = sum_a + (3*pi*c^2)/(2*omega_ai(i,a)^3)* ...
        (Gamma(i,a)/(omega_ai(i,a)-omega_li)+Gamma(i,a)/(omega_ai(i,a)+omega_li));
    end
    phiTz = phiTz + sum_a * I_sub;
end

%Stability analysis with the added optical potetential term
new_phiTotal = phiC + phiT + phiTz; % Total potential (with tweezers)
coords = [r(1,:); r(2,:); r(3,:)];
new_gradPhi = gradient(new_phiTotal,coords(:));

% Parameters values for the optical potential term (I use only 171Yb+ ions)
% Some parameters, taken from Figure 5 and table 1 (first column) in the paper
beam_waist = 0.458e-6; % in meters
lambda = 3.75e-7; % wavelength of the tweezer on each ion (in meters)
omega_li_init = c/lambda; % laser frequency (as c/wavelength)
Power = 0.0871; % Tweezer power on each ion (in Watt)
states_wavelength = [3.69419e-7,3.47730e-7,3.28937e-7,2.97056e-7,2.89138e-7]; % wavelengths of each excited state of 171Yb+ (taken from Figure 5 from the paper)
omega_ai_init = repelem(c./states_wavelength,N); % frequency of each excited state as c/wavelength for each ion in Hz
Gamma_init = repelem([1.23e8,1.055e7,1.62e8,2.61e7,3.42e7],N); % scattering rate (Einstein A coefficients of each excited state for each ion in Hz

new_param_list = [M,sym('Omega_RF'),sym("V_DC"),sym("V_RF"),sym("epsilon_0"),psi1,psi2,q,xi1,xi2,xi3, ... % Conventional trap parameters
              omega_li, P_i, lambda_i, sigma0, c,omega_ai(:)',Gamma(:)',r_star(:)']; % Optical potential parameters including phiTz

% Allora, per ora uso per r_star i valori di equilibrio trovati per la
% trappola convenzionale. Probabilmente non è corretto, però.
new_param_inits = [M_init,OmegaRF_init,Vdc_init,Vrf_init,eps0_init,psi_init,q_init,xi_init,...  % Conventional trap parameters
    omega_li_init,Power,lambda,beam_waist,c,omega_ai_init,Gamma_init,crits_mat]; % Optical potential parameters including phiTz
new_gradPhi = subs(new_gradPhi,new_param_list, new_param_inits); % Substitute the parameters' values in the gradient expression

new_crits = vpasolve(new_gradPhi == 0.,'Random',true);
% I set a tolerance for the values of crits, setting to zero too low values
% (This allows to regain some stability)
new_crits_mat = cell2mat(struct2cell(new_crits))';
new_crits_mat(abs(new_crits_mat) < 1e-25) = 0.;

%% Now I compute the Hessian matrix for each direction
dim = size(coords,1);
new_Hs = arrayfun(@(alpha) hessian(new_phiTotal,coords(alpha,:)), 1:dim,'UniformOutput',false); % compute the Hessian for each direction (x,y,z)
new_Hs = cat(3,new_Hs{:});
new_Hs = subs(new_Hs,new_param_list,new_param_inits); % substitute the parameters' values in the Hessian matrix
new_Hs_eq = subs(new_Hs,symvar(new_Hs),new_crits_mat); % compute the value of the Hessian matrix in the critical point previously found
new_Hs_eq(isnan(new_Hs_eq)) = 0.; % Substitute NaN values with zeros in the Hessian matrix

%% Study the eigenvalues of H_eq for stability and modes
new_e = arrayfun(@(alpha) eig(new_Hs_eq(:,:,alpha)), 1:dim,'UniformOutput',false);
new_e = cat(2,new_e{:})';
% Find the normal modes and frequencies
new_As = arrayfun(@(alpha) Minv2 * new_Hs_eq(:,:,alpha) * Minv2, 1:dim, 'UniformOutput', false);
new_As = cat(3,new_As{:});
% ora poi A va diagonalizzata e trovato eigenvalues and eigenvectors
[new_Bs,new_Ws] = arrayfun(@(alpha) eig(new_As(:,:,alpha)), 1:dim, 'UniformOutput', false); 
new_Bs = cat(3,new_Bs{:});
new_Ws = arrayfun(@(alpha) diag(new_Ws{alpha}), 1:dim, 'UniformOutput', false);
new_Ws = cat(2,new_Ws{:})';