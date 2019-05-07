% Solve opposed diffusion flame using the damped Newton method.
clear all; close all; clc;

global N K Le P rtol atol z zL zR mdot_L mdot_R T_L T_R Y_L Y_R C U phi F;
global gas MW NAME iCH4 iH2 iO2 iN2 iAR iH2O iCO iCO2 iNO iNO2;

%% Mechanism
gas = GRI30('Mix'); % Using the GRI 3.0 mechanism
Le = 1.0; % Unity Lewis
MW = molecularWeights(gas); % Molecular weight, Kg/Kmol
NAME = speciesNames(gas); % Name of each species
K = nSpecies(gas); % Total num of species

iCH4 = speciesIndex(gas, 'CH4');
iH2 = speciesIndex(gas, 'H2');
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
iAR = speciesIndex(gas, 'AR');
iH2O = speciesIndex(gas, 'H2O');
iCO = speciesIndex(gas, 'CO');
iCO2 = speciesIndex(gas, 'CO2');
iNO = speciesIndex(gas, 'NO');
iNO2 = speciesIndex(gas, 'NO2');

%% Settings
P = oneatm; % The constant pressure, Pa
mdot_f = 0.1; % Mass flux of fuel, Kg/(m^2 * s).
mdot_o = 0.3; % Mass flux of oxidizer, Kg/(m^2 * s).
L = 0.05; % Length of domain, m.
rtol = 1e-5; % Relative tolerance
atol = sqrt(eps); % Absolute tolerance

%% B.C.
mdot_L = mdot_f ; % Mass flux at left, Kg/(m^2 * s)
mdot_R = -mdot_o; % Mass flux at right, Kg/(m^2 * s)
T_L = 300.0; % Temperature at left, K
T_R = 300.0; % Temperature at right, K
% CH4:0.5, H2:0.5 at left
Y_L = zeros(K, 1);
Y_L(iCH4) = 0.5*MW(iCH4) / (0.5*MW(iCH4) + 0.5*MW(iH2));
Y_L(iH2) = 1.0 - Y_L(iCH4);
% N2:0.78, O2:0.21, AR:0.01at right
Y_R = zeros(K, 1);
Y_R(iAR) = 0.01*MW(iAR) / (0.78*MW(iN2) + 0.21*MW(iO2) + 0.01*MW(iAR));
Y_R(iO2) = 0.21*MW(iO2) / (0.78*MW(iN2) + 0.21*MW(iO2) + 0.01*MW(iAR));
Y_R(iN2) = 1.0 - (Y_R(iO2) + Y_R(iAR));

%% Initialize
[N, raw_data, trans_data] = load_existing_case(mdot_f, mdot_o, L);
zL = raw_data(1, 1); % Position of left endpoint, m
zR = raw_data(N, 1); % Position of right endpoint, m
if abs(zR - zL - L) > 1e-6
    error("Inconsistent domain size");
end
z = raw_data(:, 1); % Coordinates for each point, m

rho0 = trans_data(:, 1);
u0 = raw_data(:, 2);
if abs(rho0(1) * u0(1) - mdot_L) > 1e-6 
    error("Inconsistent mass flux at left");
end
if abs(rho0(N) * u0(N) - mdot_R) > 1e-6
    error("Inconsistent mass flux at right");
end
V0 = raw_data(:, 3);
T0 = raw_data(:, 4);
if abs(T0(1) - T_L) > 1e-6
    error("Inconsistent temperature at left");
end
if abs(T0(N) - T_R) > 1e-6
    error("Inconsistent temperature at right");
end
Nbla0 = raw_data(:, 5);
Y0 = zeros(K, N);
for k = 1:K
    Y0(k, :) = raw_data(:, 5+k);
    if abs(Y0(k, 1) - Y_L(k)) > 1e-6
        error("Inconsistent Y_%s at left", NAME{k});
    end
    if abs(Y0(k, N) - Y_R(k)) > 1e-6
        error("Inconsistent Y_%s at right", NAME{k});
    end
end

CUR = 1;
NEXT = 2;

C = 4+K;  % Num of unknowns per node
U = C*N; % Total num of unknowns
phi = zeros(U, 2); % Solution vector
F = zeros(U, 2); % Residual vector

phi(:, CUR) = construct_solution_vector(u0, V0, T0, Nbla0, Y0);

%% Solve 
J = zeros(U, U); % Jacobian matrix
global_converged = false;
global_iter_cnt = 0;
while(~global_converged)
    global_iter_cnt = global_iter_cnt + 1;
    
    F(:, CUR) = calculate_ss_residual_vector(phi(:, CUR));
    % Calculate the Jacobian by finite difference perturbations
    tic
    for j = 1:U
        delta = perturbation_delta(phi(j, CUR));
        phi(j, CUR) = phi(j, CUR) + delta;
        F(:, NEXT) = calculate_ss_residual_vector(phi(:, CUR));
        phi(j, CUR) = phi(j, CUR) - delta;
        J(:, j) = (F(:, NEXT) - F(:, CUR))/delta;
    end
    toc
    % Solve the Jacobian
    dphi = linsolve(J, F(:, CUR));
    
    % Update
    phi(:, NEXT) = phi(:, CUR) + dphi(:);
    
    CUR = 3-CUR;
    NEXT = 3-NEXT;
end

%% Auxiliary functions
function ret = construct_solution_vector(u, V, T, Nbla, Y)
    global N K U;
    
    ret = zeros(U, 1);
    cnt = 1;
    for i = 1:N
        ret(cnt) = u(i); cnt = cnt + 1;
        ret(cnt) = V(i); cnt = cnt + 1;
        ret(cnt) = T(i); cnt = cnt + 1;
        ret(cnt) = Nbla(i); cnt = cnt + 1;
        for k = 1:K
            ret(cnt) = Y(k, i); cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
    end
end

function [u, V, T, Nbla, Y] = mapback_solution_vector(phi)
    global N K U;
    
    u = zeros(N, 1); % m/s
    V = zeros(N, 1); % 1/s
    T = zeros(N, 1); % K
    Nbla = zeros(N, 1); %The eigenvalue
    Y = zeros(K, N); %Mass fraction
    
    cnt = 1;
    for i = 1:N
        u(i) = phi(cnt); cnt = cnt + 1;
        V(i) = phi(cnt); cnt = cnt + 1;
        T(i) = phi(cnt); cnt = cnt + 1;
        Nbla(i) = phi(cnt); cnt = cnt + 1;
        for k = 1:K
            Y(k, i) = phi(cnt); cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
    end
end

function ret = calculate_ss_residual_vector(phi)
    global N K P U gas MW z mdot_L mdot_R T_L T_R Y_L Y_R;
    
    [u, V, T, Nbla, Y] = mapback_solution_vector(phi);
    rho = zeros(N, 1); % Density, Kg / m^3
    mu = zeros(N, 1); %Viscosity, Pa * s = Kg / (m * s)
    Cp = zeros(N, 1); %Specific heat, J / (Kg * K)
    Cp_R = zeros(K, N); % Non-dimensionalized specific heats of species at constant pressure by gas constant, whose unit is J/(Kmol*K)
    lambda = zeros(N, 1); %Thermal conductivity, W / (m * K)
    D = zeros(K, N); %Binary diffusion coefficients, m^2 / s
    RS = zeros(N, 1); %Energy source due to chemical reaction, J / (m^3 * s)
    RR = zeros(K, N); %Chemical reaction rate, Kg / (m^3 * s)
    j = zeros(K, N); % Diffusion mass flux, Kg / (m^2 * s)
    
    % Properties
    for i = 1:N
        loc_T = T(i); % K
        set(gas, 'T', loc_T, 'P', P, 'Y', squeeze(Y(:, i)));
        rho(i) = density(gas); %Kg / m^3
        mu(i) = viscosity(gas); % Pa * s = Kg / (m * s)
        lambda(i) = thermalConductivity(gas); % W / (m * K)
        Cp(i) = cp_mass(gas); % J / (Kg * K)
        Cp_R(:, i) = gasconstant * cp_R(gas) ./ MW; % J / (Kg * K)
        D(:, i) = lambda(i) / (rho(i) * Cp(i)); % Compute from Unity Lewis, m^2 / s
        w = netProdRates(gas); % Kmol / (m^3 * s)
        h = enthalpies_RT(gas) * loc_T * gasconstant; % J/Kmol
        RS(i) = dot(w, h); % J / (m^3 * s)
        RR(:, i) = w .* MW; % Kg / (m^3 * s)
    end
    
    % Derivatives
    dVdz = df_upwind(V, z, u);
    ddVddz = ddf(V, z);
    dTdz = df_upwind(T, z, u);
    ddTddz = ddf(T, z);
    dYdz = zeros(K, N);
    ddYddz = zeros(K, N);
    for k = 1:K
        dYdz(k, :) = df_upwind(Y(k, :), z, u);
        ddYddz(k, :) = ddf(Y(k, :), z);
    end
     
     % Species diffusion mass flux
     for i = 1:N
         j(:, i) = -rho(i) * D(:, i) .* dYdz(:, i); % Be careful with sign convention
         j(:, i) = j(:, i) - Y(:, i) * sum(j(:, i)); % Ensure the sum of mass flux is zero.
     end
    
    % Residuals
    ret = zeros(U, 1);
    cnt = 1;
    for i = 1:N
        % Continuity equation
        if i == N
            ret(cnt) = rho(i) * u(i) - mdot_R; % B.C. of mdot at right.
        else
            ret(cnt) = (rho(i+1) * u(i+1) - rho(i) * u(i))/(z(i+1) - z(i)) + rho(i) * V(i) + rho(i+1) * V(i+1);
        end
        cnt = cnt +1;
        % Radial Momentum
        if i == 1 || i == N
            ret(cnt) = V(i); % B.C. of V at both left and right.
        else
            ret(cnt) = rho(i)*u(i)*dVdz(i)+rho(i)*V(i)^2 + Nbla(i) - mu(i)*ddVddz(i);
        end
        cnt = cnt + 1;
        % Energy
        if i == 1
            ret(cnt) = T(i) - T_L; % B.C. of T at left.
        elseif i == N
            ret(cnt) = T(i) - T_R; % B.C. of T at right.
        else
            ret(cnt) = rho(i) * Cp(i) * u(i) * dTdz(i) -lambda(i) * ddTddz(i) + RS(i) + dot(j(:, i), Cp_R(:, i)) * dTdz(i);
        end
        cnt = cnt + 1;
        % Eigenvalue
        if i == 1
            ret(cnt) = rho(i) * u(i) - mdot_L; % B.C. of mdot at left.
        else
            ret(cnt) = Nbla(i) - Nbla(i-1);
        end
        cnt = cnt + 1;
        % Species
        for k = 1:K
            if i == 1
                ret(cnt) = rho(i) * u(i) * Y(k, i) + j(k, i) - mdot_L * Y_L(k); % B.C. of Y_k at left.
            elseif i==N
                ret(cnt) = rho(i) * u(i) * Y(k, i) + j(k, i) - mdot_R * Y_R(k); % B.C. of Y_k at right.
            else
                ret(cnt) = rho(i)*u(i)*dYdz(k, i)-rho(i)*D(k,i)*ddYddz(k,i)-RR(k, i);
            end
            cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
    end
end

function check_residual_vector(F, threshold)
    global U C;
    
    for cnt = 1:U
        if F(cnt) > threshold
            fprintf("Pnt%d_Var%d: %g\n", floor((cnt-1)/C), mod(cnt-1, C) + 1, F(cnt));
        end
    end
end

function ret = relaxation(a, b, alpha)
    ret = (1-alpha) * a + alpha * b;
end

function ret = df_upwind(f, x, upwind_var)
%Upwind 1st order difference
    N = length(x);
    ret = zeros(1, N);
    ret(1) = (f(2) - f(1))/(x(2)-x(1));
    for i = 2 : N-1
        if upwind_var(i) > 0
            ret(i) = (f(i) - f(i-1))/(x(i)-x(i-1));
        else
            ret(i) = (f(i+1) - f(i))/(x(i+1)-x(i));
        end
    end
    ret(N) = (f(N) - f(N-1))/(x(N)-x(N-1));
end

function ret = ddf(f, x)
%Central 2nd order difference
    N = length(x);
    ret = zeros(1, N);
    ret(1) = 2.0/(x(2)-x(3))*((f(2)-f(1))/(x(2)-x(1)) - (f(3)-f(1))/(x(3)-x(1)));
    for i = 2 : N-1
        dxl = x(i) - x(i-1);
        dxr = x(i+1) - x(i);
        dfl = f(i-1) - f(i);
        dfr = f(i+1) - f(i);
        ret(i) = 2.0 / (dxl+dxr) * (dfl/dxl+dfr/dxr);
    end
    ret(N) = 2.0/(x(N-2)-x(N-1))*((f(N)-f(N-2))/(x(N)-x(N-2)) - (f(N)-f(N-1))/(x(N)-x(N-1)));
end

function ret = perturbation_delta(x)
    global rtol atol;
    
    % Preserve sign(x)
    if x >= 0
        ret = rtol * x + atol;
    else
        ret = rtol * x - atol;
    end
end

function [lines, raw_data, trans_data] = load_existing_case(mf, mo, domain_len)
    case_str = sprintf("mf=%g_mo=%g_L=%g", mf, mo, domain_len);
    raw_data_path = "../data/" + case_str + "_raw.txt";
    trans_data_path = "../data/" + case_str + "_transformed.txt";
    
    raw_tbl = importdata(raw_data_path);
    raw_data = raw_tbl.data;
    
    trans_tbl = importdata(trans_data_path);
    trans_data = trans_tbl.data;
    
    a = size(raw_data);
    b = size(trans_data);
    if  a(1) ~= b(1)
        error("Inconsistent input data!");
    end
    lines = a(1);
end
