clear all; close all; clc;

P = oneatm; % Pa
mdot_L = 1.0 / 100 ; % Fuel stream, Kg/s
mdot_R = -16.6 / 100; % Air stream, Kg/s
rhoL = 0.716; % Density of CH4, Kg/m^3
rhoR = 1.3947; % Density of Air, Kg/m^3
uL = mdot_L / rhoL;
uR = mdot_R / rhoR;

fuel = Methane();
air = Air();
gas = GRI30('Mix');
setPressure(gas, P);

N = 1001; % Total num of grid points
K = nSpecies(gas); % Total num of species

MW = molecularWeights(gas); % Kg / Kmol
NAME = speciesNames(gas);

zL = 0.0;
zR = 0.1 / 10; % m
L = zR - zL;
z = linspace(zL, zR, N);
dz = z(2)-z(1);

Tmin = 300.0; % K
Tmax = 1500.0; % K
Tmax_pos = relaxation(zL, zR, 0.5);

PREV = 1;
CUR = 2;

rho = zeros(2, N); % Kg / m^3
u = zeros(2, N); % m/s
V = zeros(2, N);
T = zeros(2, N); % K
Y = zeros(2, K, N);
Nbla = zeros(1, 2); % The eigenvalue

mu = zeros(1, N); % Viscosity, Pa * s = Kg / (m * s)
cp = zeros(1, N); % Specific heat, J / (Kg * K)
lambda = zeros(1, N); % Thermal conductivity, W / (m * K)
D = zeros(K, N); % Binary diffusion coefficients, m^2 / s

RS = zeros(N, 1); % Chemical reaction source term, J / (m^3 * s)
RR = zeros(K, N); % Chemical reaction rate, Kg / (m^3 * s)

%% Init
rho(PREV, :) = linspace(rhoL , rhoR, N);
rho(CUR, :) = rho(PREV, :);

u(PREV, :) = linspace(uL, uR, N);
u(CUR, :) = u(PREV, :);

Nbla(PREV) = -0.1;

Y(PREV, speciesIndex(gas, 'CH4'), :) = linspace(1.0, 0.0, N);
Y(PREV, speciesIndex(gas, 'O2'), :) = linspace(0.0, massFraction(air, 'O2'), N);
Y(PREV, speciesIndex(gas, 'N2'), :) = linspace(0.0, massFraction(air, 'N2'), N);

V(PREV, :) = -df(rho(PREV, :) .* u(PREV, :), dz, N) ./ (2 * rho(PREV, :));

for i = 1:N
    if abs(z(i) - Tmax_pos) < 0.1 * L
        T(PREV, i) = Tmax;
    else
        T(PREV, i) = Tmin;
    end
end
T(CUR, :) = T(PREV, :);

%% Loop
err = 1.0;
iter_cnt = 0;
while(err > 1e-3)
    iter_cnt = iter_cnt + 1;
    fprintf("Iteration %d:\n", iter_cnt);
    
    %====================Calc physical properties=====================
    for i = 1:N
        local_T = T(PREV, i);
        
        set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(PREV, :, i)));
        
        mu(i) = viscosity(gas);
        lambda(i) = thermalConductivity(gas);
        cp(i) = cp_mass(gas);
        D(:, i) = mixDiffCoeffs(gas);      
        
        w = netProdRates(gas); % kmol / (m^3 * s)
        h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol                                              
       
        RS(i) = -dot(w, h); % J / (m^3 * s)
        RR(:, i) = w.* MW; % Kg / (m^3 * s)
        
%         if(i == int32(N/2))
%             fprintf('Position: %d, Local T: %f K\n', i, local_T);
%             fprintf('%16s%24s%32s%32s%40s\n', 'Species', 'Y', 'Enthalpy(J * Kmol^-1)', 'RR(kmol * m^-3 * s^-1)', 'h*wdot(J * m^-3 * s^-1)');
%             st = zeros(K, 1);
%             for k = 1:K
%                 st(k) = h(k)*w(k);
%                 fprintf('%16s%24.6f%32.6e%32.6e%40.6e\n', NAME{1,k}, Y(PREV, k, i), h(k), w(k), st(k));
%             end
%             fprintf('Local Energy Source(J * m^-3 * s^-1): %16.6e\n', -sum(st));
%         end        

    end
    
    %===========================Diagnose==============================
    subplot(3, 4, 1)
    plot(z, T(PREV, :))
    title('$$T$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('K')
    
    subplot(3, 4, 2)
    plot(z, rho(PREV, :))
    title('$$\rho$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('Kg\cdotm^{-3}');
    
    subplot(3, 4, 3)
    plot(z, u(PREV, :) );
    title('$$u$$','Interpreter','latex')
    xlabel('z / m');
    ylabel('m/s');
    
    subplot(3, 4, 4)
    plot(z, RS);
    title('$$-\sum{h_k\dot{\omega}_k}$$','Interpreter','latex')
    xlabel('z / m');
    ylabel('J\cdotm^{-3}\cdots^{-1}');
      
    subplot(3, 4, 5)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'CH4'), :)))
    title('Y_{CH4}');
    xlabel('z / m');
    
    subplot(3, 4, 9)
    plot(z, squeeze(RR(speciesIndex(gas, 'CH4'), :)))
    title('$$\dot{\omega}_{CH_4}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')
    
    subplot(3, 4, 6)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'O2'), :)))
    title('Y_{O2}');
    xlabel('z / m');
    
    subplot(3, 4, 10)
    plot(z, squeeze(RR(speciesIndex(gas, 'O2'), :)))
    title('$$\dot{\omega}_{O_2}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')
    
    subplot(3, 4, 7)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'CO2'), :)))
    title('Y_{CO2}')
    xlabel('z / m')
    
    subplot(3, 4, 11)
    plot(z, squeeze(RR(speciesIndex(gas, 'CO2'), :)))
    title('$$\dot{\omega}_{CO_2}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')
    
    subplot(3, 4, 8)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'H2O'), :)))
    title('Y_{H2O}')
    xlabel('z / m')
    
    subplot(3, 4, 12)
    plot(z, squeeze(RR(speciesIndex(gas, 'H2O'), :)))
    title('$$\dot{\omega}_{H_2O}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')
        
    %==============================Solve V============================
    coef = zeros(N, N);
    for i = 2 : N-1
        coef(i, i-1) = -0.5 * rho(PREV, i) * u(PREV, i) / dz - mu(i) / dz^2;
        coef(i, i) = rho(PREV, i) * V(PREV, i) + 2 * mu(i) / dz^2;
        coef(i, i+1) = 0.5 * rho(PREV, i) * u(PREV, i) / dz - mu(i) / dz^2;
    end
    rhs = -Nbla(PREV)*ones(N, 1);
    
    A = coef(2:N-1, 2:N-1);
    b = rhs(2:N-1);
    b(1) = b(1) - coef(2, 1) * V(PREV, 1);
    b(N-2) = b(N-2) - coef(N-1, N) * V(PREV, N);
    x = linsolve(A, b);
    
    V(CUR, 1) = V(PREV, 1);
    V(CUR, 2:N-1) = x;
    V(CUR, N) = V(PREV, N);
    
    %==============================Solve u============================
    u(CUR, 1) = u(PREV, 1);
    for i = 2 : N - 1
        u(CUR, i) = (-2 * rho(PREV, i) * V(CUR, i) * dz + rho(PREV, i-1) * u(PREV, i-1)) / rho(PREV, i);
    end
    u(CUR, N) = u(PREV, N);
    
    %=============================Correct V===========================
    flux = 0.0;
    flux = flux - 2 * rho(PREV, 1) * V(CUR, 1) * (dz/2);
    for i = 2 : N-1
        flux = flux - 2 * rho(PREV, i) * V(CUR, i) * dz;
    end
    flux = flux - 2 * rho(PREV, N) * V(CUR, N) * (dz/2);
    gain_factor = flux / (mdot_R - mdot_L) ;
    V(CUR, :) = V(CUR, :) / gain_factor;
    
    %===========================Correct Nbla==========================
    dVdz = df(V(CUR, :), dz, N);
    
    lhs1 = dot(rho(PREV, :) .* u(CUR, :), dVdz);
    lhs2 = dot(rho(PREV, :) .* V(CUR, :), V(CUR, :));
    rhs2 = sum(df(mu .* dVdz, dz, N));
    
    Nbla(CUR) = (rhs2 - lhs1 - lhs2) / N;
    err = abs(Nbla(CUR) - Nbla(PREV));
    fprintf("\terr = %f\n", err);
    Nbla(CUR) = relaxation(Nbla(PREV), Nbla(CUR), 0.5);
    
    %=============================Solve T============================   
    coef = zeros(N, N);
    for i = 2 : N-1
        coef(i, i-1) = -0.5 * rho(PREV, i) * u(CUR, i) * cp(i) / dz - lambda(i) / dz^2;
        coef(i, i) = 2 * lambda(i) / dz^2;
        coef(i, i+1) = 0.5 * rho(PREV, i) * u(CUR, i) * cp(i) / dz - lambda(i) / dz^2;
    end
    
    A = coef(2:N-1, 2:N-1);
    b = RS(2:N-1);
    b(1) = b(1) - coef(2, 1) * T(PREV, 1);
    b(N-2) = b(N-2) - coef(N-1, N) * T(PREV, N);
    x = linsolve(A, b);
    
    T(CUR, 1) = T(PREV, 1);
    T(CUR, 2:N-1) = x;
    T(CUR, N) = T(PREV, N);
    
    for i = 1:N
        T(CUR, i) = max(300, T(CUR, i));
    end
        
    %=============================Sovle Yk=============================
    for k=1:K
        coef = zeros(N, N);
        for i = 2 : N-1
            coef(i, i-1) = -rho(PREV, i) * u(CUR, i) / (2*dz) - rho(PREV, i) * D(k, i) / dz^2;
            coef(i, i) = 2 * rho(PREV, i) * D(k, i) / dz^2;
            coef(i, i+1) = rho(PREV, i) * u(CUR, i) / (2*dz) - rho(PREV,i) * D(k, i) / dz^2;
        end
        rhs = transpose(RR(k, :));
        
		A = coef(2:N-1, 2:N-1);
		b = rhs(2:N-1);
		b(1) = b(1) - coef(2, 1) * Y(PREV, k, 1);
		b(N-2) = b(N-2) - coef(N-1, N) * Y(PREV, k, N);
		x = linsolve(A, b);
		
		Y(CUR, k, 1) = Y(PREV, k, 1);
        Y(CUR, k, 2:N-1) = x;
		Y(CUR, k, N) = Y(PREV, k, N);
        
        for i = 2 : N-1
            Y(CUR, k, i) = max(Y(CUR, k, i), 0.0);
        end
    end
    
    for i = 2 : N-1
        tmp = sum(Y(CUR, :, i));
        Y(CUR, :, i) = Y(CUR, :, i) / tmp;
    end
    
    %===========================Update density==========================
    rho(CUR, 1) = rho(PREV, 1);
    for i = 2:N-1
        rho(CUR, i) = P / (gasconstant * T(CUR, i) * sum(squeeze(Y(CUR, :, i)) ./ MW'));
    end
    rho(CUR, N) = rho(PREV, N);
      
	%===========================Swap Index==============================
    PREV = 3 - PREV;
    CUR = 3 - CUR;
end

%% Helpers
function ret = df(f, dx, N)
    ret = zeros(1, N);
    
    ret(1) = f(2) - f(1);
    for i = 2 : N-1
        ret(i) = (f(i+1) - f(i-1))/2;
    end
    ret(N) = f(N) - f(N-1);
    
    ret = ret / dx;
end

function ret = ddf(f, dx, N)
    ret = zeros(1, N);
    
    ret(1) = f(3) - 2*f(2) + f(1);
    for i = 2 : N-1
        ret(i) = f(i+1) - 2*f(i) + f(i-1);
    end
    ret(N) = f(N-2) - 2*f(N-1) + f(N);
    
    ret = ret / dx^2;
end

function ret = relaxation(a, b, alpha)
    ret = (1-alpha) * a + alpha * b;
end
