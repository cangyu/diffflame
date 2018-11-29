clear all; close all; clc;

P = oneatm;
mdot_L = 1.0; %Fuel stream, Kg/s
mdot_R = -16.6; %Air stream, Kg/s
rhoL = 0.716; %Density of CH4, Kg/m^3
rhoR = 1.3947; %Density of Air, Kg/m^3
uL = mdot_L / rhoL;
uR = mdot_R / rhoR;

fuel = Methane();
air = Air();
gas = GRI30('Mix');
setPressure(gas, P);

N = 101; % Total num of grid points
K = nSpecies(gas); % Total num of species

MW = molecularWeights(gas); % Kg/Kmol

zL = 0.0;
zR = 0.1; %10cm
L = zR - zL;
z = linspace(zL, zR, N);
dz = z(2)-z(1);

T_L = 300.0; %K
T_R = 300.0; %K
Tmax = 2000.0; %K
Tmax_pos = lin_dist(zL, zR, 0.5);
Tcoef = polyfit([zL, zR, Tmax_pos], [T_L, T_R, Tmax],2);

PREV = 1;
CUR = 2;

rho = zeros(2, N); % Kg/m^3
u = zeros(2, N); % m/s
V = zeros(2, N);
T = zeros(2, N); % K
Y = zeros(2, K, N);
Nbla = zeros(1, 2); % The eigenvalue

mu = zeros(1, N); % Viscosity, Pa * s = Kg / (m * s)
cp = zeros(1, N); % Specific heat, J / (Kg * K)
lambda = zeros(1, N); % Thermal conductivity, W / (m * K)
D = zeros(K, N); % Binary diffusion coefficients, m^2 / s

RR = zeros(N, 1); % Chemical source term, J / (m^3 * s)

%% Init
rho(PREV, :) = linspace(rhoL , rhoR, N);
% rho(PREV, :) = ones(NumOfPnt, 1) * 1.0;
rho(CUR, :) = rho(PREV, :);
u(PREV, :) = linspace(uL, uR, N);
Nbla(PREV) = -0.1;
Y(PREV, speciesIndex(gas, 'CH4'), :) = linspace(1.0, 0.0, N);
Y(PREV, speciesIndex(gas, 'O2'), :) = linspace(0.0, massFraction(air, 'O2'), N);
Y(PREV, speciesIndex(gas, 'N2'), :) = linspace(0.0, massFraction(air, 'N2'), N);
V(PREV, :) = -df(rho(PREV, :) .* u(PREV, :), dz, N) ./ (2 * rho(PREV, :));
for i = 1:N
    T(PREV, i) = polyval(Tcoef, z(i));
end
T(CUR, :) = T(PREV, :);
% plot(z, T(PREV, :))
% plot(z, squeeze(Y(PREV, in2, :)))

%% Loop
err = 1.0;
iter_cnt = 0;
while(err > 1e-6)
    iter_cnt = iter_cnt + 1;
    
    %====================Calc physical properties=====================
    for i = 1:N
        local_T = T(PREV, i);
        
        setTemperature(gas, local_T);
        setMassFractions(gas, Y(PREV, :, i));
        
        mu(i) = viscosity(gas);
        lambda(i) = thermalConductivity(gas);
        cp(i) = cp_mass(gas);
        D(:, i) = mixDiffCoeffs(gas);
        
        w = netProdRates(gas); % kmol / (m^3 * s)
        h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol                                              
        RR(i) = dot(w, h); % J / (m^3 * s)
    end
    
    %==============================Solve V============================
    coef = zeros(N, N);
    coef(1, 1) = - rho(PREV, 1) * u(PREV, 1) / dz - mu(1) / dz^2 + rho(PREV, 1) * V(PREV, 1);
    coef(1, 2) =  rho(PREV, 1) * u(PREV, 1) / dz + 2 * mu(1) / dz^2;
    coef(1, 3) = -mu(1) / dz^2;
    for i = 2 : N-1
        coef(i, i-1) = -0.5 * rho(PREV, i) * u(PREV, i) / dz - mu(i) / dz^2;
        coef(i, i) = rho(PREV, i) * V(PREV, i) + 2 * mu(i) / dz^2;
        coef(i, i+1) = 0.5 * rho(PREV, i) * u(PREV, i) / dz - mu(i) / dz^2;
    end
    coef(N, N-2) = -mu(N) / dz^2;
    coef(N, N-1) = -rho(PREV, N) * u(PREV, N) / dz + 2 * mu(N) / dz^2;
    coef(N, N) = rho(PREV, N) * u(PREV, N) / dz - mu(N) / dz^2 + rho(PREV, N) * V(PREV, N);
    
    rhs = -Nbla(PREV)*ones(N, 1);
    
    V(CUR, :) = linsolve(coef, rhs);
    
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
    ddVddz = ddf(V(CUR, :), dz, N);
    
    lhs1 = dot(rho(PREV, :) .* u(CUR, :), dVdz);
    lhs2 = dot(rho(PREV, :) .* V(CUR, :), V(CUR, :));
    rhs2 = dot(mu, ddVddz);
    
    Nbla(CUR) = (rhs2 - lhs1 - lhs2) / N;
    err = abs(Nbla(CUR) - Nbla(PREV));
    Nbla(CUR) = lin_dist(Nbla(PREV), Nbla(CUR), 0.5);
    
    %=============================Solve T============================
    coef = zeros(N, N);
    coef(1, 1) = - rho(PREV, 1) * u(CUR, 1) * cp(1) / dz - lambda(1) / dz^2;
    coef(1, 2) = rho(PREV, 1) * u(CUR, 1) * cp(1) / dz + 2 * lambda(1) / dz^2;
    coef(1, 3) = -lambda(1) / dz^2;
    for i = 2 : N-1
        coef(i, i-1) = -0.5 * rho(PREV, i) * u(CUR, i) * cp(i) / dz - lambda(i) / dz^2;
        coef(i, i) = 2 * lambda(i) / dz^2;
        coef(i, i+1) = 0.5 * rho(PREV, i) * u(CUR, i) * cp(i) / dz - lambda(i) / dz^2;
    end
    coef(N, N-2) = -lambda(N) / dz^2;
    coef(N, N-1) = -rho(PREV, N) * u(CUR, N) * cp(N) / dz + 2 * lambda(N) / dz^2;
    coef(N, N) = rho(PREV, N) * u(CUR, N) * cp(N) / dz - lambda(N) / dz^2;    
    
    T(CUR, :) = linsolve(coef, -RR);
        
    %=============================Sovle Yk=============================
    ydot = zeros(K, N);
    for i = 1:N
        local_T = T(CUR, i);
        setTemperature(gas, local_T);
        setMassFractions(gas, Y(PREV, :, i));
        ydot(:, i) = netProdRates(gas) .* MW;
    end
    
    for k=1:K
        coef = zeros(N, N);
        coef(1, 1) = -rho(PREV, 1) * u(CUR, 1) / dz + rho(PREV, 1) * D(k, 1) / dz^2;
        coef(1, 2) = rho(PREV, 1) * u(CUR, 1) / dz - 2 * rho(PREV, 1)* D(k, 1) /dz^2;
        coef(1, 3) = rho(PREV, 1) * D(k, 1) / dz^2;
        for i = 2 : N-1
            coef(i, i-1) = -rho(PREV, i) * u(CUR, i) / (2*dz) + rho(PREV, i) * D(k, i) / dz^2;
            coef(i, i) = -2 * rho(PREV, i) * D(k, i) / dz^2;
            coef(i, i+1) = rho(PREV, i) * u(CUR, i) / (2*dz) + rho(PREV,i) * D(k, i) / dz^2;
        end
        coef(N, N-2) = rho(PREV, N) * D(k, N) / dz^2;
        coef(N, N-1) = -rho(PREV, N) * u(CUR, N) / dz - 2 * rho(PREV, N) * D(k, N) / dz^2;
        coef(N, N) = rho(PREV, N) * u(CUR, N) / dz + rho(PREV, N) * D(k, N) / dz^2;
        
        rhs = transpose(ydot(k, :));
        
        Y(CUR,k, :) = linsolve(coef, rhs);
    end
    
    %===========================Update density==========================
    rho(CUR, 1) = rho(PREV, 1);
    for i = 2:N-1
        tmp = 0.0;
        for k = 1:K
            tmp = tmp + Y(CUR, k, i) / MW(k);
        end
        rho(CUR, i) = P / (gasconstant * T(CUR, i) * tmp);
    end
    rho(CUR, N) = rho(PREV, N);
    
    fprintf("Iteration %d: err = %f\n", iter_cnt, err);
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

function ret = lin_dist(a, b, alpha)
    ret = (1-alpha) * a + alpha * b;
end
