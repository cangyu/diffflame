% Solve opposed diffusion flame using the damped Newton method.
clear all; close all; clc;

global N K  C U z zL zR mdot_L mdot_R T_L T_R Y_L Y_R;
global Le P gas MW NAME iCH4 iH2 iO2 iN2 iAR iH2O iCO iCO2 iNO iNO2;
global rtol atol ss_atol ss_rtol ts_atol ts_rtol phi_prev;

%% Mechanism
gas = GRI30('Mix'); % Using the GRI 3.0 mechanism
Le = 1.0; % Unity Lewis
MW = molecularWeights(gas); % Molecular weight, Kg/Kmol
NAME = speciesNames(gas); % Name of each species
K = nSpecies(gas); % Total num of species
C = 4+K;  % Num of unknowns per node

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
ss_rtol = 1e-4*ones(C, 1);
ts_rtol = 1e-4*ones(C, 1);
ss_atol = 1e-9*ones(C, 1);
ts_atol = 1e-11*ones(C, 1);

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
U = C*N; % Total num of unknowns
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

phi_prev = construct_solution_vector(u0, V0, T0, Nbla0, Y0); 
phi = phi_prev; % Solution vector

%% Solve 
global_converged = false;
global_iter_cnt = 0;
dt = 1e-6;
while(~global_converged)
    global_iter_cnt = global_iter_cnt + 1;
    
    F = calculate_ss_residual_vector(phi);
    ss_norm2(phi, F)
    
    J = calculate_jacobian(0.0, phi, F);
    
    % Solve the Jacobian
    dphi = linsolve(J, -F);
    
    % Update
    phi = phi + dphi;
    phi_prev = phi;
end

%% Functions
function J = calculate_jacobian(rdt, phi, F)
    % Calculate the Jacobian by finite difference perturbations
    global U phi_prev;
    
    J = zeros(U, U);
    tic
    if rdt == 0.0
        fprintf("Calculating Steady-State Jacobian matrix ...\n");
        for j = 1:U
            delta = perturbation_delta(phi(j));
            phi(j) = phi(j) + delta;
            F_perturb = calculate_ss_residual_vector(phi);
            phi(j) = phi(j) - delta;
            J(:, j) = (F_perturb - F)/delta;
        end
    else
        fprintf("Calculating Time-Stepping Jacobian matrix ...\n");
        for j = 1:U
            delta = perturbation_delta(phi(j));
            phi(j) = phi(j) + delta;
            F_perturb = calculate_ts_residual_vector(rdt, phi_prev, phi);
            phi(j) = phi(j) - delta;
            J(:, j) = (F_perturb - F)/delta;
        end
    end
    toc
end

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
            ret(cnt) = ret(cnt) / rho(i);
        end
        cnt = cnt + 1;
        % Energy
        if i == 1
            ret(cnt) = T(i) - T_L; % B.C. of T at left.
        elseif i == N
            ret(cnt) = T(i) - T_R; % B.C. of T at right.
        else
            ret(cnt) = rho(i) * Cp(i) * u(i) * dTdz(i) -lambda(i) * ddTddz(i) + RS(i) + dot(j(:, i), Cp_R(:, i)) * dTdz(i);
            ret(cnt) = ret(cnt) / (rho(i) * Cp(i));
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
                ret(cnt) = ret(cnt) / rho(i);
            end
            cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
    end
end

function ret = calculate_ts_residual_vector(rdt, phi_prev, phi_cur)
    global N K P U gas MW z mdot_L mdot_R T_L T_R Y_L Y_R;

    [~, V_prev, T_prev, ~, Y_prev] = mapback_solution_vector(phi_prev);
    [u, V, T, Nbla, Y] = mapback_solution_vector(phi_cur);
    rho = zeros(N, 1); % Density, Kg / m^3
    mu = zeros(N, 1); % Viscosity, Pa * s = Kg / (m * s)
    Cp = zeros(N, 1); % Specific heat, J / (Kg * K)
    Cp_R = zeros(K, N); % Specific heat of each species, J / (Kg * K)
    lambda = zeros(N, 1); % Thermal conductivity, W / (m * K)
    D = zeros(K, N); % Binary diffusion coefficients, m^2 / s
    RS = zeros(N, 1); % Energy source due to chemical reaction, J / (m^3 * s)
    RR = zeros(K, N); % Chemical reaction rate, Kg / (m^3 * s)
    j = zeros(K, N); % Diffusion mass flux, Kg / (m^2 * s)
    
    % Properties
    for i = 1:N
        loc_T = T(i); % K
        set(gas, 'T', loc_T, 'P', P, 'Y', squeeze(Y(:, i)));
        rho(i) = density(gas); % Kg / m^3
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
         j(:, i) = j(:, i) - Y(:, i) * sum(j(:, i)); % Ensure the sum of mass flux is zero
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
            ret(cnt) = ret(cnt) / rho(i);
            ret(cnt) = ret(cnt) + rdt*(V(i)-V_prev(i)); % Add transient term
        end
        cnt = cnt + 1;
        % Energy
        if i == 1
            ret(cnt) = T(i) - T_L; % B.C. of T at left.
        elseif i == N
            ret(cnt) = T(i) - T_R; % B.C. of T at right.
        else
            ret(cnt) = rho(i)*Cp(i)*u(i)*dTdz(i)-lambda(i)*ddTddz(i)+RS(i)+dot(j(:, i), Cp_R(:, i))*dTdz(i);
            ret(cnt) = ret(cnt) / (rho(i) * Cp(i));
            ret(cnt) = ret(cnt) + rdt*(T(i)-T_prev(i)); % Add transient term
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
                ret(cnt) = ret(cnt) / rho(i);
                ret(cnt) = ret(cnt) + rdt*(Y(k, i)-Y_prev(k, i)); % Add transient term
            end
            cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
    end
end

function diagnose_vector(F, threshold)
    global U C;
    
    for cnt = 1:U
        if F(cnt) > threshold
            fprintf("pnt%d-var%d: %g\n", floor((cnt-1)/C), mod(cnt-1, C) + 1, F(cnt));
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

function ret = norm1(res_vec)
    ret = log10(norm(res_vec, inf));
end

function ret = ss_norm2(sol_vec, res_vec)
    global ss_atol ss_rtol N K C U;
    
    [u, V, T, Nbla, Y] = mapback_solution_vector(sol_vec);
    
    w = zeros(C, 1);
    w(1) = ss_rtol(1) * norm(u, 1) / N + ss_atol(1);
    w(2) = ss_rtol(2) * norm(V, 1) / N + ss_atol(2);
    w(3) = ss_rtol(3) * norm(T, 1) / N + ss_atol(3);
    w(4) = ss_rtol(4) * norm(Nbla, 1) / N + ss_atol(4);
    for k = 1:K
        n = k+4;
        w(n) = ss_rtol(n) * norm(Y(k, :), 1) / N + ss_atol(n);
    end
    
    [res_u, res_V, res_T, res_Nbla, res_Y] = mapback_solution_vector(res_vec);
    
    ret = 0.0;
    ret = ret + sum(res_u .^ 2) / w(1)^2;
    ret = ret + sum(res_V .^ 2) / w(2)^2;
    ret = ret + sum(res_T .^ 2) / w(3)^2;
    ret = ret + sum(res_Nbla .^ 2) / w(4)^2;
    for k = 1:K
        n = k+4;
        ret = ret + sum(res_Y(k, :) .^ 2) / w(n)^2;
    end
    
    ret = sqrt(ret/U);
end

function ret = ts_norm2(sol_vec, res_vec)
    global ts_atol ts_rtol N K C U;
    
    [u, V, T, Nbla, Y] = mapback_solution_vector(sol_vec);
    
    w = zeros(C, 1);
    w(1) = ts_rtol(1) * norm(u, 1) / N + ts_atol(1);
    w(2) = ts_rtol(2) * norm(V, 1) / N + ts_atol(2);
    w(3) = ts_rtol(3) * norm(T, 1) / N + ts_atol(3);
    w(4) = ts_rtol(4) * norm(Nbla, 1) / N + ts_atol(4);
    for k = 1:K
        n = k+4;
        w(n) = ts_rtol(n) * norm(Y(k, :), 1) / N + ts_atol(n);
    end
    
    [res_u, res_V, res_T, res_Nbla, res_Y] = mapback_solution_vector(res_vec);
    
    ret = 0.0;
    ret = ret + sum(res_u .^ 2) / w(1)^2;
    ret = ret + sum(res_V .^ 2) / w(2)^2;
    ret = ret + sum(res_T .^ 2) / w(3)^2;
    ret = ret + sum(res_Nbla .^ 2) / w(4)^2;
    for k = 1:K
        n = k+4;
        ret = ret + sum(res_Y(k, :) .^ 2) / w(n)^2;
    end
    
    ret = sqrt(ret/U);
end

function draw_profile(sol_vec)
    global P N K z iCH4 iH2 iO2 iN2 iH2O iCO iCO2 iNO iNO2 gas MW;
    
    [u, V, T, Nbla, Y] = mapback_solution_vector(sol_vec);
    rho = zeros(N, 1);
    for i = 1:N
        rho(i) = calc_density(P, T(i), Y(:, i));
    end
    RS = zeros(N, 1); % Energy source due to chemical reaction, J / (m^3 * s)
    RR = zeros(K, N); % Chemical reaction rate, Kg / (m^3 * s)
    for i = 1:N
        loc_T = T(i); % K
        set(gas, 'T', loc_T, 'P', P, 'Y', squeeze(Y(:, i)));
        w = netProdRates(gas); % Kmol / (m^3 * s)
        h = enthalpies_RT(gas) * loc_T * gasconstant; % J/Kmol
        RS(i) = dot(w, h); % J / (m^3 * s)
        RR(:, i) = w .* MW; % Kg / (m^3 * s)
    end
    
    
    h = figure(1);
    set(h, 'position', get(0,'ScreenSize'));

    subplot(3, 6, 1)
    plot(z, T)
    title('$$T$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('K')

    subplot(3, 6, 2)
    plot(z, rho)
    title('$$\rho$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('Kg\cdotm^{-3}');

    subplot(3, 6, 3)
    plot(z, u);
    title('$$u$$','Interpreter','latex')
    xlabel('z / m');
    ylabel('m/s');

    subplot(3, 6, 4)
    plot(z, RS);
    title('$$-\sum{h_k\dot{\omega}_k}$$','Interpreter','latex')
    xlabel('z / m');
    ylabel('J\cdotm^{-3}\cdots^{-1}');

    subplot(3, 6, 5)
    plot(z, Y(iNO, :))
    title('Y_{NO}');
    xlabel('z / m');

    subplot(3, 6, 6)
    plot(z, Y(iCO, :))
    title('Y_{CO}');
    xlabel('z / m');

    subplot(3, 6, 7)
    plot(z, Y(iCH4, :))
    title('Y_{CH4}');
    xlabel('z / m');

    subplot(3, 6, 8)
    plot(z, Y(iO2, :))
    title('Y_{O2}');
    xlabel('z / m');

    subplot(3, 6, 9)
    plot(z, Y(iCO2, :))
    title('Y_{CO2}')
    xlabel('z / m')

    subplot(3, 6, 10)
    plot(z, Y(iH2O, :))
    title('Y_{H2O}')
    xlabel('z / m')

    subplot(3, 6, 11)
    plot(z, Y(iH2, :))
    title('Y_{H2}')
    xlabel('z / m')

    subplot(3, 6, 12)
    plot(z, Y(iN2, :))
    title('Y_{N2}')
    xlabel('z / m')

    subplot(3, 6, 13)
    plot(z, RR(iCH4, :))
    title('$$\dot{\omega}_{CH_4}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')

    subplot(3, 6, 14)
    plot(z, RR(iO2, :))
    title('$$\dot{\omega}_{O_2}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')

    subplot(3, 6, 15)
    plot(z, RR(iCO2, :))
    title('$$\dot{\omega}_{CO_2}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')

    subplot(3, 6, 16)
    plot(z, RR(iH2O, :))
    title('$$\dot{\omega}_{H_2O}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')

    subplot(3, 6, 17)
    plot(z, RR(iH2, :))
    title('$$\dot{\omega}_{H_2}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')
end

function show_solution_diff(sol1, sol2)
    % TODO
end

function ret = calc_density(p, t, y)
    global MW K;
    
    tmp = 0;
    for k = 1:K
        tmp = tmp + y(k) / MW(k);
    end
    
    ret = p / (gasconstant * t * tmp);
end
