N = 5; % number of ions in the chain (tpar: trap parameters)
%eps0 = 8.85418782e-12; % permittivity of free space ([e0] = m-3 kg-1 s4 A2)

syms q M [1 N] real; % electronic charge and mass (respectively) of each ion
%q = sym('q_%d',[1 N]);
%q
%syms q1 q2 q3 q4 q5 real;
%q = [q1,q2,q3,q4,q5];
syms r [3 N] real; % position vector of each ion
syms eps0 Vdc Vrf OmegaRF real positive; % permittivity of free space, trap's DC and rf peak voltages (respectively)
syms Xi Psi [1 3] real; % Xi_a and phy_a trap geometric factors

%%
% Define the potentials
phiC = 0; % Coulomb potential
for i = 1:N-1 % devo fare 1:N-1 perché per i=N dopo j=i+1:N è vuoto e quindi dà errore. e infatti matematicamente il termine della somma per i = N si ski
    for j = i+1:N
        phiC = phiC + q(i) * q(j) / (4*pi*eps0*norm(r(:,i)-r(:,j)));
    end
end

%{
% Trap frequencies
syms wTx(M,q); 
syms wTy(M,q);
syms wTz(M,q);

% Define the trap frequencies
wTx = sqrt(q*Vdc*Xi(3)/M);
wTy = sqrt(-q*Vdc*Xi(1)/M + q^2*Vrf^2*Psi(1)^2/(2*M^2*OmegaRF^2));
wTz = sqrt(-q*Vdc*Xi(2)/M + q^2*Vrf^2*Psi(2)^2/(2*M^2*OmegaRF^2));
%}

phiT = 0; % Trap potential (in the paper, phi_conv)
for i =  1:N % I directly compute the sum of the potential over the i index
    phiT = phiT + 0.5*M(i)*((sqrt(q(i)*Vdc*Xi(3)/M(i))^2)*r(1,i)^2 + ...
                             (sqrt(-q(i)*Vdc*Xi(1)/M(i) + q(i)^2*Vrf^2*Psi(1)^2/(2*M(i)^2*OmegaRF^2))^2)*r(2,i)^2 + ...
                             (sqrt(-q(i)*Vdc*Xi(2)/M(i) + q(i)^2*Vrf^2*Psi(2)^2/(2*M(i)^2*OmegaRF^2))^2)*r(3,i)^2);
end

phiTotal = phiC + phiT; % Total potential (without tweezers)
disp(phiTotal);

coord = [r(1,:) r(2,:) r(3,:)];
H = hessian(phiTotal,coord);

