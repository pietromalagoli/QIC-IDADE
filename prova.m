clear
%% Note: all this analysis is done on 171Yb+ ions
%{
digits(4)
sympref('FloatingPointOutput',true);
%}

N = 3; % number of ions in the chain (tpar: trap parameters)
%eps0 = 8.85418782e-12; % permittivity of free space ([e0] = m-3 kg-1 s4 A2)

syms q M [1 N] real; % electronic charge and mass (respectively) of each ion
assume(M>0); % set the condition for positive mass
%q = sym('q_%d',[1 N]);
%q
%syms q1 q2 q3 q4 q5 real;
%q = [q1,q2,q3,q4,q5];
syms r [3 N] real; % position vector of each ion
syms eps0 Vdc Vrf OmegaRF real positive; % permittivity of free space, trap's DC and rf peak voltages (respectively)
syms xi psi [1 3] real positive; % Xi (DC electrodes distance?) and psy (rf electrodes distance?) trap geometric factors

eps0 = sym('epsilon_0'); % This unsures that eps0 is displayed as the greek letter
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

disp(phiC);

% Trap frequencies
syms wTx(M_i,q_i); 
syms wTy(M_i,q_i);
syms wTz(M_i,q_i);

% Define the trap frequencies
wTx(M_i,q_i) = sqrt(-q_i*Vdc*xi(2)/M_i + q_i^2*Vrf^2*psi(2)^2/(2*M_i^2*OmegaRF^2));
wTy(M_i,q_i) = sqrt(-q_i*Vdc*xi(1)/M_i + q_i^2*Vrf^2*psi(1)^2/(2*M_i^2*OmegaRF^2));
wTz(M_i,q_i) = sqrt(q_i*Vdc*xi(3)/M_i);
disp(wTx);
disp(wTy);
disp(wTz);

phiT = 0; % Trap potential (in the paper, phi_conv)
for i =  1:N % I directly compute the sum of the potential over the i index
    phiT = phiT + 0.5*M(i)*((wTx(M(i),q(i))^2)*r(1,i)^2 + ...
                             (wTy(M(i),q(i))^2)*r(2,i)^2 + ...
                             (wTz(M(i),q(i))^2)*r(3,i)^2);
end
disp(phiT);
phiTotal = phiC + phiT; % Total potential (without tweezers)
disp(phiTotal);
coords = [r(1,:) r(2,:) r(3,:)];
gradPhi = gradient(phiTotal,coords);
disp(gradPhi);
% Set placeholder values for parameters
q_init = [1e-19,1e-19,1e-19]; % 1 unit of electric charge in Coulomb
M_init = [1e-25,1e-25,1e-25]; % here I used the fact that the mass of Yb atoms (whose ions are used in an example in the paper have mass ~173 u)
xi_init = [0.06, 0.06, 0.06]; % Example values for trap geometric factors
psi_init = [0.06, 0.06]; % Example values for trap geometric factors
Vdc_init = 10; % in volts (reference: Table 1 at http://dx.doi.org/10.1088/0143-0807/34/3/787) 
Vrf_init = 550; % in volts (reference: Table 1 at http://dx.doi.org/10.1088/0143-0807/34/3/787)
OmegaRF_init = 1.8e7; % Frequency in Hz (reference: Table 1 at http://dx.doi.org/10.1088/0143-0807/34/3/787)
eps0_init = 8.85e-12; % permittivity of free space ([e0] = m-3 kg-1 s4 A2)

param_list = [M1,M2,M3,sym('Omega_RF'),sym("V_DC"),sym("V_RF"),sym("epsilon_0"),psi1,psi2,q1,q2,q3,xi1,xi2,xi3];
param_inits = [M_init,OmegaRF_init,Vdc_init,Vrf_init,eps0_init,psi_init,q_init,xi_init];
gradPhi = subs(gradPhi,param_list, param_inits); % Substitute the parameters' values in the gradient expression

crits = vpasolve(gradPhi == 0.,'Random',true);
% I set a tolerance for the values of crits, setting to zero too low values
% (This allows to regain some stability)
crits_mat = cell2mat(struct2cell(crits));
crits_mat(abs(crits_mat) < 1e-25) = 0.;
disp(crits_mat');
% Now I compute the Hessian martix
H = hessian(phiTotal,coords);
H = subs(H,param_list, param_inits); % Substitute the parameters' values in the Hessian expression

H_eq = subs(H,symvar(H),crits_mat'); % compute the value of the Hessian matrix in the critical point previously found 
% Study the eigenvalues of H_eq for stability and modes
H_eq(isnan(H_eq)) = 0.; % Substitute NaN values with zeros in the Hessian matrix
e = eig(H_eq);
disp(e');
% Find the normal modes and frequencies
M = repelem(cell2mat(sym2cell(M)), 3);  % each mass appears 3 times because there are 3 variables for each ion (x,y,z)
Minv2 = diag(1 ./ sqrt(M));
A = Minv2 * H_eq * Minv2;
A = subs(A,M,repelem(M_init,3)); % subtitutes the values of the masses in the expression
% ora poi A va diagonalizzata e trovato eigenvalues and eigenvectors
[B,W] = eig(A); 
% I use diag to extract the diagonal elements of the matrix
modes = diag(B); % Normal modes
freqs = diag(W); % Normal modes' frequencies
disp(modes');
disp(freqs');

%% Now we add the tweezers to system (and therefore the associated term to the potential)
% Some parameters, taken from Figure 5 and table 1 (first column) in the paper
N_excited = 5; % number of excited state (in the case of 171Yb+ is 5)
beam_waist = 0.458e-6; % in meters
lambda = 3.75e-7; % wavelength of the tweezer on each ion (in meters)
P = 0.0871; % Tweezer power on each ion (in Watt)
syms omega_ai Gamma_ai [N N_excited] real positive; % atomic transition frequencies and scattering rates of the transition from ground state to each excited state a for each ion
syms omega_li P_i lambda_i real positive; % tweezer frequency, power and wavelength on each ion
syms r_star r [3 N] real; % equilibrium positions
syms sigma0 real positive; % beam waist
syms x y z real; % dummy variables
pi = sym('pi');
c = sym('c');

% Define the tweezer intensity profile 
%{
syms I_i(x,y,z,P_i,lambda_i); % intensity profile
syms sigma(x); % spot size 
sigma(x) = sigma0*sqrt(1+(x/(pi*sigma0^2/lambda_i))^2);
I_i = (2*P_i/(pi*sigma0^2))*((sigma0^2)/sigma(x))^2*exp(-2*(y^2+z^2)/sigma(x)^2);
% Define a symbolic function for the difference between the ion's position
% and its equilibrium position
syms Delta(alpha,alpha_star);
Delta(alpha, alpha_star) = alpha - alpha_star; 
%}
sigma_x = sigma0*sqrt(1+(x/(pi*(sigma0^2)/lambda_i))^2);
I_i = (2*P_i/(pi*sigma0^2))*((sigma0^2)/sigma_x)^2*exp(-2*(y^2+z^2)/sigma_x^2);


phiTz = sym(0);
for i = 1:N % sum over the ions
    for a = 1:N_excited % sum over the excited states of each ion
        x_i = r(1,i) - r_star(1,i);
        y_i = r(2,i) - r_star(2,i);
        z_i = r(3,i) - r_star(3,i);

        I_sub = subs(I_i, ...
            [x,   y,   z], ...
            [x_i, y_i, z_i]);

        phiTz = phiTz + (3*pi*c^2)/(2*omega_ai(i,a)^3)* ...
        (Gamma_ai(i,a)/(omega_ai(i,a)-omega_li)+Gamma_ai(i,a)/(omega_ai(i,a)+omega_li))* ...
        I_sub;
    end
end
disp(phiTz);
