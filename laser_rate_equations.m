

%{
The Rate Equations:

**Carriers**
dN/dt = (eta_i*I)/(q*V) - (Rsp+Rnr) - vg*g*Np

**Photons**
dNp/dt = [Gamma*vg*g - 1/tau_p]*Np + Gamma*Rsp'


**Gain**
g = g0N*ln((N+Ns)/(Ntr+Ns))

For InGaAs/GaAs 8nm QW:
    Ntr = 1.8e18
    Ns = -0.4e18
    g0N = 1800

Rsp ~ BN^2, B ~10^-10 cm^3/s

%}


% Constants (from Coldren et al.)
eta_i = 0.8;
q = 1.602e-19;
V = 4e-12;          % cm^3
B = 0.8e-10;         % cm^3/s
vg = (3/4.2)*1e10;    % cm/s
Gamma = 0.032;
tau_p = 2.77e-12;     % s
Rsp_p = 1.02e23;      % cm^-3/s
am = 45.6;      % cm^-1
lambda = 0.98;  % um
h = 6.626176e-34; % J*s
hv = (1.23985/lambda)*q; % eV
Vp = 1.25e-10;  % cm^3

Ntr = 1.8e18;
Ns = -0.4e18;
g0N = 1800;


% For now, assume a single mode laser. I plan to add additional modes
% later.

% Initial Conditions
N0 = 0;
Np0 = 0;

% Time step for the calculations
dt = 0.5e-12;


% Time vector
t_array = 0:dt:10e-9;
total_ns = max(t_array)./(1e-9);

% Threshold Current
Ith = 1.11e-3;  % A

% Input Current
I = 2.*Ith.*ones(1, length(t_array));

% Initialize carrier density, photon density and gain arrays
N = zeros(1, length(t_array));
Np = zeros(1, length(t_array));
g = zeros(1, length(t_array));

% Set the first element of each array to their initial conditions
N(1) = N0;
Np(1) = Np0;
g(1) = 0;


% Calculate number of iterations per nanosecond
i_ns = floor( (1e-9)/dt);
j = 1;

for i = 2:length(t_array)

    %g = g0N.*log((N(i-1)+Ns)./(Ntr+Ns));
    g(i) = real(g0N.*log((N(i-1)+Ns)./(Ntr+Ns)));

    Rsp = B.*N(i-1).^2;
    Rnr = (3.5e-30).*N(i-1).^3;

    dN = ((eta_i*I(i))./(q*V) - (Rsp+Rnr) - vg.*g(i).*Np(i-1)).*dt;
    
    dNp = ((Gamma*vg.*g(i) - 1/tau_p).*Np(i-1) + Gamma*Rsp_p).*dt;

    N(i) = N(i-1) + dN;
    Np(i) = Np(i-1) + dNp;
    
    % Display the progress at each nanosecond
    if i == j*i_ns
        %fprintf('Completed %d out of %d iterations\n', i, length(t_array));
        fprintf('Simulated for %d ns\n', j);
        j = j+1;
    end
end

% Calculate optical power output from photon density Np
Po = vg.*am.*Np.*hv.*Vp;
%{
figure(1)
yyaxis left
plot(t_array.*1e9, N, 'b-');
ylabel('Carrier Density (cm^{-3})');

yyaxis right
plot(t_array.*1e9, Po.*1e3, 'r-');
ylabel('Power Output (mW)');

xlabel('Time (ns)');
title('Dual Y-axis Plot')
%}

% Plot carrier density and output power together. Current on separate
% subplot
figure(1)
subplot(2, 1, 1);  % top
yyaxis left
plot(t_array.*1e9, N, 'b-');
ylabel('Carrier Density (cm^{-3})');

yyaxis right
plot(t_array.*1e9, Po.*1e3, 'r-');
ylabel('Power Output (mW)');

subplot(2, 1, 2);  % bottom
plot(t_array.*1e9, I.*1e3, 'r-');
ylabel('Current (mA)');
