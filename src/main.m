% Solve opposed diffusion flame using the damped Newton method.
clear; close all; clc;

global N K  C U z zL zR mdot_L mdot_R T_L T_R Y_L Y_R DampFactor MaxDampRound;
global Le P gas MW NAME iCH4 iH2 iO2 iN2 iAR iH2O iCO iCO2 iNO iNO2 filtering_sigma;
global rtol atol ss_atol ss_rtol ts_atol ts_rtol diag_mask ts_mask phi_prev gConverged;
global MW_C MW_H MW_O MW_N Yc_fu Yh_fu Yo_fu Yc_ox Yh_ox Yo_ox gauss_weight;

%% Mechanism
gas = GRI30('Mix'); % Use the GRI 3.0 mechanism.
Le = 1.0; % Unity Lewis
MW = molecularWeights(gas); % Molecular weight, Kg/Kmol
NAME = speciesNames(gas); % Name of each species
K = nSpecies(gas); % Total num of species
C = 4+K;  % Num of unknowns per node, (u V T Nbla Y) in sequence.

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

MW_C = MW(iCH4) - 2*MW(iH2);
MW_H = MW(iH2)/2;
MW_O = MW(iO2)/2;
MW_N = MW(iN2)/2;

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
DampFactor = sqrt(2);
MaxDampRound = 7;
filtering_sigma = 2.52e-4;

%% B.C.
mdot_L = mdot_f ; % Mass flux at left, Kg/(m^2 * s)
mdot_R = -mdot_o; % Mass flux at right, Kg/(m^2 * s)
T_L = 300.0; % Temperature at left, K
T_R = 300.0; % Temperature at right, K
% CH4:0.5, H2:0.5 at left
Y_L = zeros(K, 1);
Y_L(iCH4) = 0.5*MW(iCH4) / (0.5*MW(iCH4) + 0.5*MW(iH2));
Y_L(iH2) = 1.0 - Y_L(iCH4);
Yc_fu = MW_C/(MW_C+4*MW_H)*Y_L(iCH4);
Yh_fu = 1.0-Yc_fu;
Yo_fu = 0.0;
% N2:0.78, O2:0.21, AR:0.01 at right
Y_R = zeros(K, 1);
Y_R(iAR) = 0.01*MW(iAR) / (0.78*MW(iN2) + 0.21*MW(iO2) + 0.01*MW(iAR));
Y_R(iO2) = 0.21*MW(iO2) / (0.78*MW(iN2) + 0.21*MW(iO2) + 0.01*MW(iAR));
Y_R(iN2) = 1.0 - (Y_R(iO2) + Y_R(iAR));
Yc_ox = 0.0;
Yh_ox = 0.0; 
Yo_ox = Y_R(iO2);

%% Initialize
[N, raw_data, trans_data] = load_existing_case(mdot_f, mdot_o, L);
U = C*N; % Total num of unknowns
zL = raw_data(1, 1); % Position of left endpoint, m
zR = raw_data(N, 1); % Position of right endpoint, m
if abs(zR - zL - L) > 1e-6
    error('Inconsistent domain size');
end
z = raw_data(:, 1); % Coordinates for each point, m
rho0 = trans_data(:, 1);
u0 = raw_data(:, 2);
if abs(rho0(1) * u0(1) - mdot_L) > 1e-6 
    error('Inconsistent mass flux at left');
end
if abs(rho0(N) * u0(N) - mdot_R) > 1e-6
    error('Inconsistent mass flux at right');
end
V0 = raw_data(:, 3);
T0 = raw_data(:, 4);
if abs(T0(1) - T_L) > 1e-6
    error('Inconsistent temperature at left');
end
if abs(T0(N) - T_R) > 1e-6
    error('Inconsistent temperature at right');
end
Nbla0 = raw_data(:, 5);
Y0 = zeros(K, N);
for k = 1:K
    Y0(k, :) = raw_data(:, 5+k);
    if abs(Y0(k, 1) - Y_L(k)) > 1e-6
        error('Inconsistent Y_%s at left', NAME{k});
    end
    if abs(Y0(k, N) - Y_R(k)) > 1e-6
        error('Inconsistent Y_%s at right', NAME{k});
    end
end
gauss_weight = zeros(N, N);
for i = 1:N
    for ii = 1:N
        gauss_weight(i, ii) = gauss1d(z(ii)-z(i));
    end
end
ts_mask = transient_jacobian_mask();
diag_mask = full(blktridiag(ones(C), ones(C), ones(C), N));
phi_prev = construct_solution_vector(u0, V0, T0, Nbla0, Y0); 

%% Solve 
gConverged = false;
gIterCnt = 0;
dt = 1e-6;
phi = phi_prev; % Solution vector
while(~gConverged)
    gIterCnt = gIterCnt + 1;
    fprintf('Iter%d:\n', gIterCnt);
    F = calculate_residual_vector(0.0, phi);
    report_solution(F, phi);

    % Check convergence first
    J = calculate_jacobian(0.0, phi, F); 
    dphi = linsolve(J, -F); % The undamped correction vector.
    gConverged = check_convergence(phi, dphi);
    if gConverged
        break;
    end

    % If not coverged, try Steady-State seeking for a convergent solution
    try
        fprintf('\tTry Steady-State...\n');
        phi_next = calculate_damped_solution_vector(J, phi, dphi);
        phi_prev = phi;
        phi = phi_next;
        fprintf('\tSteady-State success.\n');
        SS_OK = true;
    catch
        fprintf('\tSteady-State failure.\n');
        SS_OK = false;
    end
    
    % If Steady-State fails, try Time-stepping to get closer to the convergent solution
    if ~SS_OK
        fprintf('\tTry Time-Stepping...\n');
        loc_ts_cnt = 0;
        
        while loc_ts_cnt < 3
            try
                rdt = 1/dt;
                loc_Jac = J + rdt * ts_mask;
                loc_F = calculate_residual_vector(rdt, phi);
                loc_dphi = linsolve(loc_Jac, -loc_F);
                phi_next = calculate_damped_solution_vector(loc_Jac, phi, loc_dphi);
                phi_prev = phi;
                phi = phi_next;
                fprintf('\t\tTime-Stepping success on dt=%g(rdt=%g), trying larger dt...\n', dt, rdt);
                dt = dt * 1.5;
            catch
                fprintf('\t\tTime-Stepping failure on dt=%g(rdt=%g).\n, trying smallar dt...', dt, rdt);
                dt = dt / 2;
                if dt < 1e-18
                    throw('Too small time-step.');
                end
            end
            loc_ts_cnt = loc_ts_cnt + 1;
        end
        
        fprintf('\tTime-Stepping success.\n');
    end
    
    fprintf('Iter%d done!\n', gIterCnt);
end
fprintf('Converged!\n');

%% Functions
function report_solution(F, phi)
    [~, ~, loc_Temp, ~, ~] = mapback_solution_vector(phi);
    fprintf('\tTmax=%gK\n', max(loc_Temp));
    fprintf('\tlog10(||F||_inf)=%g\n', log10(norm1(F)));
end

function ret=check_convergence(phi, dphi)
    global atol rtol;

    n2 = norm1(dphi);
    n3 = max(atol, rtol*norm1(phi));
    n4 = norm2(0.0, phi, dphi);
    
    fprintf('\t||dphi||_inf=%g, convergence criteria: %g\n', n2, n3);
    fprintf('\t||dphi||_weighted=%g, convergence criteria: %g\n', n4, 1.0);
    
    ret = n2 < n3 || n4 < 1;
end

function diagnose_residual(F, F_cmp)
    global z NAME iCH4 iH2 iO2 iN2 iAR iH2O iCO iCO2 iNO iNO2;

    [res_u0, res_V0, res_T0, res_A0, res_Y0] = mapback_solution_vector(F);
    [res_u1, res_V1, res_T1, res_A1, res_Y1] = mapback_solution_vector(F_cmp);
    dF = F + F_cmp;
    [res_u, res_V, res_T, res_A, res_Y] = mapback_solution_vector(dF);
    
    figure(1);
    subplot(211);
    plot(z, res_u0, 'r-', z, res_u1, 'b-');
    title('Continuity eqn');
    subplot(212);
    plot(z, res_u, '-.');
    
    figure(2);
    subplot(211);
    plot(z, res_V0, 'r-', z, res_V1, 'b-');
    title('Radial Momentum eqn');
    subplot(212);
    plot(z, res_V, '-.');
    
    figure(3);
    subplot(211);
    plot(z, res_T0, 'r-', z, res_T1, 'b-');
    title('Energy eqn');
    subplot(212);
    plot(z, res_T, '-.');
    
    figure(4);
    subplot(211);
    plot(z, res_A0, 'r-', z, res_A1, 'b-');
    title('Eigenvalue eqn');
    subplot(212);
    plot(z, res_A, '-.');
    
    figure(5);
    subplot(211);
    plot(z, res_Y0(iCH4, :), 'r-', z, res_Y1(iCH4, :), 'b-');
    title(NAME{iCH4});
    subplot(212);
    plot(z, res_Y(iCH4, :), '-.');
    
    figure(6);
    subplot(211);
    plot(z, res_Y0(iH2, :), 'r-', z, res_Y1(iH2, :), 'b-');
    title(NAME{iH2});
    subplot(212);
    plot(z, res_Y(iH2, :), '-.');

    figure(7);
    subplot(211);
    plot(z, res_Y0(iO2, :), 'r-', z, res_Y1(iO2, :), 'b-');
    title(NAME{iO2});
    subplot(212);
    plot(z, res_Y(iO2, :), '-.');
    
    figure(8);
    subplot(211);
    plot(z, res_Y0(iN2, :), 'r-', z, res_Y1(iN2, :), 'b-');
    title(NAME{iN2});
    subplot(212);
    plot(z, res_Y(iN2, :), '-.');
    
    figure(9);
    subplot(211);
    plot(z, res_Y0(iAR, :), 'r-', z, res_Y1(iAR, :), 'b-');
    title(NAME{iAR});
    subplot(212);
    plot(z, res_Y(iAR, :), '-.');
    
    figure(10);
    subplot(211);
    plot(z, res_Y0(iH2O, :), 'r-', z, res_Y1(iH2O, :), 'b-');
    title(NAME{iH2O});
    subplot(212);
    plot(z, res_Y(iH2O, :), '-.');
    
    figure(11);
    subplot(211);
    plot(z, res_Y0(iCO, :), 'r-', z, res_Y1(iCO, :), 'b-');
    title(NAME{iCO});
    subplot(212);
    plot(z, res_Y(iCO, :), '-.');
    
    figure(12);
    subplot(211);
    plot(z, res_Y0(iCO2, :), 'r-', z, res_Y1(iCO2, :), 'b-');
    title(NAME{iCO2});
    subplot(212);
    plot(z, res_Y(iCO2, :), '-.');
    
    figure(13);
    subplot(211);
    plot(z, res_Y0(iNO, :), 'r-', z, res_Y1(iNO, :), 'b-');
    title(NAME{iNO});
    subplot(212);
    plot(z, res_Y(iNO, :), '-.');
    
    figure(14);
    subplot(211);
    plot(z, res_Y0(iNO2, :), 'r-', z, res_Y1(iNO2, :), 'b-');
    title(NAME{iNO2});
    subplot(212);
    plot(z, res_Y(iNO2, :), '-.');
end

function ret = transient_jacobian_mask()
    global U N C;
    
    ret = ones(U, 1);
    
    for i = 1:C
        ret(i) = 0;
        ret(U-i+1) = 0;
    end
    
    cnt = C+1;
    for pnt_idx = 2:N-1
        for var_idx = 1:C
            if var_idx == 1 || var_idx == 4
                ret(cnt) = 0;
            end
            cnt = cnt + 1;
        end
    end
    
    ret = diag(ret);
end

function phi_next = calculate_damped_solution_vector(Jac, phi0, dphi0)
    global MaxDampRound DampFactor gConverged;

    tau = 1.0;
    k = 0;
    dampingOK = false;
    s0 = norm2(0.0, phi0, dphi0);
    
    while ~dampingOK
        % Check if has reached the max trial limit.
        if k > MaxDampRound
            break;
        end
        phi_next = phi0 + tau * dphi0;
        
        % Check if the terminating condition is satisfied.
        F_damp = calculate_residual_vector(0.0, phi_next);
        dphi_damp = linsolve(Jac, -F_damp);
        s1 = norm2(0.0, phi_next, dphi_damp);
        dampingOK = norm1(dphi_damp) < norm1(dphi0) || s1 < 1.0 || s1 < s0;
        if dampingOK
            if s1 < 1.0
                gConverged = true;
            end
            break;
        end
        
        % Update
        tau = tau / DampFactor;
        k = k + 1;
    end
    
    if ~dampingOK
        throw('Max damping iter reached.')
    end
end

function J = calculate_jacobian(rdt, phi, F)
    % Calculate the Jacobian by finite difference perturbations
    global U diag_mask;
    J = zeros(U, U);
    
    if rdt == 0.0
        fprintf('\tCalculating Steady-State Jacobian matrix ...\n');
    else
        fprintf('\tCalculating Time-Stepping Jacobian matrix ...\n');
    end
    
    tic
    for j = 1:U
        delta = perturbation_delta(phi(j));
        phi(j) = phi(j) + delta;
        F_perturb = calculate_residual_vector(rdt, phi);
        phi(j) = phi(j) - delta;
        J(:, j) = (F_perturb - F)/delta;
    end
    toc
    
    %J = J .* diag_mask; % Enforce those off tri-diagnoal blocks being 0.
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
        error('Internal error!');
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
        error('Internal error!');
    end
end

function ret = calculate_residual_vector(rdt, phi)
    global N K P U gas MW z mdot_L mdot_R T_L T_R Y_L Y_R phi_prev;
    
    [~, V_prev, T_prev, ~, Y_prev] = mapback_solution_vector(phi_prev);
    [u, V, T, Nbla, Y] = mapback_solution_vector(phi);
    rho = zeros(N, 1); % Density, Kg / m^3
    mu = zeros(N, 1); %Viscosity, Pa * s = Kg / (m * s)
    Cp = zeros(N, 1); %Specific heat, J / (Kg * K)
    Cp_R = zeros(K, N); % Specific heat of each species, J / (Kg * K)
    lambda = zeros(N, 1); %Thermal conductivity, W / (m * K)
    D = zeros(K, N); %Binary diffusion coefficients, m^2 / s
    unfiltered_wdot = zeros(K, N); % Kmol / (m^3 * s)
    enthalpy = zeros(K, N); % J/Kmol
    RS = zeros(N, 1); %Energy source due to chemical reaction, J / (m^3 * s)
    RR = zeros(K, N); %Chemical reaction rate, Kg / (m^3 * s)
    
    % Properties
    for i = 1:N
        loc_T = T(i); % K
        set(gas, 'T', loc_T, 'P', P, 'Y', Y(:, i));
        rho(i) = density(gas); % Kg / m^3
        mu(i) = viscosity(gas); % Pa * s = Kg / (m * s)
        lambda(i) = thermalConductivity(gas); % W / (m * K)
        Cp(i) = cp_mass(gas); % J / (Kg * K)
        Cp_R(:, i) = gasconstant * cp_R(gas) ./ MW; % J / (Kg * K)
        D(:, i) = lambda(i) / (rho(i) * Cp(i)); % Compute from Unity Lewis, m^2 / s
        unfiltered_wdot(:, i) = netProdRates(gas); % Kmol / (m^3 * s)
        enthalpy(:, i) = enthalpies_RT(gas) * loc_T * gasconstant; % J/Kmol
    end
    
    % Filtering
    filtered_wdot = calculate_filtered_wdot(unfiltered_wdot); % Kmol / (m^3 * s)
    %filtered_wdot = unfiltered_wdot;
    for i = 1:N
        RS(i) = sum(enthalpy(:, i) .* filtered_wdot(:, i)); % J / (m^3 * s)
        RR(:, i) = filtered_wdot(:, i) .* MW; % Kg / (m^3 * s)
    end
    
    % Derivatives
    dmudz0 = df_central(mu, z);
    dlambdadz0 = df_central(lambda, z);
    drhodz0 = df_central(rho, z);
    dVdz = df_upwind(V, z, u);
    dVdz0 =  df_central(V, z);
    dTdz = df_upwind(T, z, u);
    dTdz0 = df_central(T, z);
    dYdz = zeros(K, N);
    dYdz0 = zeros(K, N);
    dDdz0 = zeros(K, N);
    for k = 1:K
        dYdz(k, :) = df_upwind(Y(k, :), z, u);
        dYdz0(k, :) = df_central(Y(k, :), z);
        dDdz0(k, :) = df_central(D(k, :), z);
    end
    j = zeros(K, N); % Diffusion mass flux, Kg / (m^2 * s)
    for i = 1:N
         j(:, i) = -rho(i) * D(:, i) .* dYdz0(:, i); % Be careful with sign convention.
         j(:, i) = j(:, i) - Y(:, i) * sum(j(:, i)); % Ensure the sum of mass flux is zero.
    end

    % Divergence
    ddVddz = ddf(V, z);
    divVisc = dmudz0 .* dVdz0 + mu .* ddVddz;
    %divVisc = df_central(mu .* dVdz0, z);
    ddTddz = ddf(T, z);
    divHeat = dlambdadz0 .* dTdz0 + lambda .* ddTddz;  
    %divHeat = df_central(lambda .* dTdz0, z);
    ddYddz = zeros(K, N);
    for k = 1:K
        ddYddz(k, :) = ddf(Y(k, :), z);
    end
    divDiffus = zeros(K, N);
    for k = 1:K
        %divDiffus(k, :) = -(drhodz0' .* D(k, :) .* dYdz0(k, :) + rho' .* dDdz0(k, :) .* dYdz0(k, :) + rho' .* D(k, :) .* ddYddz(k, :));
        divDiffus(k, :) = df_central(j(k, :), z);
    end

    % Residuals
    ret = zeros(U, 1);
    cnt = 1;
    for i = 1:N
        % Continuity
        if i == N
            ret(cnt) = rho(i) * u(i) - mdot_R; % B.C. of mdot at right.
        else
            ret(cnt) = (rho(i+1) * u(i+1) - rho(i) * u(i))/(z(i+1) - z(i)) + rho(i) * V(i) + rho(i+1) * V(i+1);
        end
        cnt = cnt +1;
        % Continuity eqn done.
        
        % Radial momentum
        if i == 1 || i == N
            ret(cnt) = V(i); % B.C. of V at both left and right.
        else
            ret(cnt) = rho(i)*u(i)*dVdz(i)+rho(i)*V(i)^2 + Nbla(i) - divVisc(i);
            ret(cnt) = ret(cnt) / rho(i);
            ret(cnt) = ret(cnt) + rdt*(V(i)-V_prev(i)); % Add transient term.
        end
        cnt = cnt + 1;
        % Radial momentum eqn done.
        
        % Energy
        if i == 1
            ret(cnt) = T(i) - T_L; % B.C. of T at left.
        elseif i == N
            ret(cnt) = T(i) - T_R; % B.C. of T at right.
        else
            ret(cnt) = rho(i)*Cp(i)*u(i)*dTdz(i)-divHeat(i)+RS(i)+dot(j(:, i),Cp_R(:, i))*dTdz(i);
            ret(cnt) = ret(cnt) / (rho(i)*Cp(i));
            ret(cnt) = ret(cnt) + rdt*(T(i)-T_prev(i)); % Add transient term.
        end
        cnt = cnt + 1;
        % Energy eqn done.
        
        % Eigenvalue
        if i == 1
            ret(cnt) = rho(i) * u(i) - mdot_L; % B.C. of mdot at left.
        else
            ret(cnt) = Nbla(i) - Nbla(i-1);
        end
        cnt = cnt + 1;
        % Eigenvalue eqn done.
        
        % Species
        for k = 1:K
            if i == 1
                ret(cnt) = Y(k, i) - Y_L(k); % B.C. of Y_k at left.
            elseif i==N
                ret(cnt) = Y(k, i) - Y_R(k); % B.C. of Y_k at right.
            else
                ret(cnt) = rho(i)*u(i)*dYdz(k, i)+divDiffus(k, i)-RR(k, i);
                ret(cnt) = ret(cnt) / rho(i);
                ret(cnt) = ret(cnt) + rdt*(Y(k, i)-Y_prev(k, i)); % Add transient term.
            end
            cnt = cnt + 1;
        end
        % Species eqn done.
    end
    
    if cnt ~= U+1
        error('Internal error!');
    end
end

function ret = relaxation(a, b, alpha)
    ret = (1-alpha) * a + alpha * b;
end

function ret = df_backward(f, x)
    % 1st order derivative using backward difference.
    N = length(x);
    ret = zeros(N, 1);
    
    ret(1) = (f(2) - f(1))/(x(2)-x(1));
    for i = 2 : N
        ret(i) = (f(i) - f(i-1))/(x(i)-x(i-1));
    end
end

function ret = df_upwind(f, x, upwind_var)
    % 1st order derivative using upwind difference.
    N = length(x);
    ret = zeros(N, 1);
    
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

function ret = df_central(f, x)
    % 1st order derivative using central difference.
    N = length(x);
    ret = zeros(N, 1);
    
    ret(1) = ((x(2) - x(1)) / (x(3) - x(1)) * (f(3) - f(1)) - (x(3) - x(1)) / (x(2) - x(1)) * (f(2) - f(1))) / (x(2) - x(3));
    for i = 2:N-1
        ret(i) = ((x(i) - x(i - 1)) / (x(i + 1) - x(i)) * (f(i + 1) - f(i)) - (x(i + 1) - x(i)) / (x(i) - x(i - 1)) * (f(i - 1) - f(i))) / (x(i + 1) - x(i - 1));
    end
    ret(N) = ((x(N) - x(N-2)) / (x(N) - x(N-1)) * (f(N-1) - f(N)) - (x(N) - x(N-1)) / (x(N) - x(N-2)) * (f(N-2) - f(N))) / (x(N-2) - x(N-1));
end

function ret = ddf(f, x)
    % 2nd order derivative using central difference.
    N = length(x);
    ret = zeros(N, 1);
    
    ret(1) = 2.0/(x(2)-x(3))*((f(2)-f(1))/(x(2)-x(1)) - (f(3)-f(1))/(x(3)-x(1)));
    for i = 2 : N-1
        dxl = x(i) - x(i-1);
        dxr = x(i+1) - x(i);
        dfl = f(i-1) - f(i);
        dfr = f(i+1) - f(i);
        ret(i) = 2.0 * (dfl/dxl+dfr/dxr) / (x(i+1) - x(i-1));
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
    case_str = sprintf('../data/mf=%g_mo=%g_L=%g', mf, mo, domain_len);
    raw_data_path = sprintf('%s_raw.txt', case_str);
    trans_data_path = sprintf('%s_transformed.txt', case_str);
    
    raw_tbl = importdata(raw_data_path);
    raw_data = raw_tbl.data;
    
    trans_tbl = importdata(trans_data_path);
    trans_data = trans_tbl.data;
    
    a = size(raw_data);
    b = size(trans_data);
    if  a(1) ~= b(1)
        error('Inconsistent input data!');
    end
    lines = a(1);
end

function ret = norm1(res_vec)
    ret = norm(res_vec, inf);
end

function ret = norm2(rdt, sol, step)
    global ts_atol ts_rtol ss_atol ss_rtol N K C U;

    if rdt == 0.0
        loc_atol = ss_atol;
        loc_rtol = ss_rtol;
    else
        loc_atol = ts_atol;
        loc_rtol = ts_rtol;
    end
    
    [u, V, T, Nbla, Y] = mapback_solution_vector(sol);
    
    w = zeros(C, 1);
    w(1) = loc_rtol(1) * norm(u, 1) / N + loc_atol(1);
    w(2) = loc_rtol(2) * norm(V, 1) / N + loc_atol(2);
    w(3) = loc_rtol(3) * norm(T, 1) / N + loc_atol(3);
    w(4) = loc_rtol(4) * norm(Nbla, 1) / N + loc_atol(4);
    for k = 1:K
        n = k+4;
        w(n) = loc_rtol(n) * norm(Y(k, :), 1) / N + loc_atol(n);
    end
    
    [res_u, res_V, res_T, res_Nbla, res_Y] = mapback_solution_vector(step);
    
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

function ret = bilger(Yc, Yh, Yo)
    % Compute the mixture fraction using the Bilger formula.
    global MW_C MW_H MW_O Yc_fu Yh_fu Yo_fu Yc_ox Yh_ox Yo_ox;

    a = 2 * (Yc - Yc_ox) / MW_C + (Yh - Yh_ox) / (2*MW_H) - 2 * (Yo - Yo_ox) / MW_O;
    b = 2 * (Yc_fu - Yc_ox) / MW_C + (Yh_fu - Yh_ox) / (2*MW_H) - 2 * (Yo_fu - Yo_ox) / MW_O;
    ret = a / b;
end

function ret = calculate_mixture_fraction(sol_vec)
    global P gas N;
    
    [~, ~, T, ~, Y] = mapback_solution_vector(sol_vec);
    ret = zeros(N, 1);
    for i = 1:N
        set(gas, 'T', T(i), 'P', P, 'Y', Y(:, i));
        yc = elementalMassFraction(gas, 'C');
        yh = elementalMassFraction(gas, 'H');
        yo = elementalMassFraction(gas, 'O');
        ret(i) = bilger(yc, yh, yo);
    end
end

function ret = calculate_density(sol_vec)
    global MW N P;
    
    ret = zeros(N, 1);
    [~, ~, T, ~, Y] = mapback_solution_vector(sol_vec);
    
    for i = 1:N
        ret(i) = P / (gasconstant * T(i) * sum(Y(:, i) ./ MW));
    end
end

function diagnose_vector(F, threshold)
    global U C;
    
    for cnt = 1:U
        if F(cnt) > threshold
            fprintf('pnt%d-var%d: %g\n', floor((cnt-1)/C), mod(cnt-1, C) + 1, F(cnt));
        end
    end
end

function show_solution_profile(sol_vec)
    global P N K z iCH4 iH2 iO2 iN2 iAR iH2O iCO iCO2 iNO iNO2 gas MW;
    
    Z = calculate_mixture_fraction(sol_vec);
    rho = calculate_density(sol_vec);
    [u, V, T, Nbla, Y] = mapback_solution_vector(sol_vec);
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
    
    ax1 = subplot(4, 7, [1 2 3 8 9 10 15 16 17]);
    plot(ax1, z, Y(iCH4, :), z, Y(iH2, :), z, Y(iN2, :), z, Y(iO2, :), z, Y(iAR, :)*10, z, Y(iH2O, :), z, Y(iCO, :)*10, z, Y(iCO2, :)*10, z, Y(iNO, :)*1e3, z, Y(iNO2, :)*1e4);
    legend(ax1, 'Y_{CH_4}','Y_{H_2}','Y_{N_2}','Y_{O_2}','10\cdotY_{AR}','Y_{H_2O}','10\cdotY_{CO}','10\cdotY_{CO_2}','1e3\cdotY_{NO}','1e4\cdotY_{NO_2}');
    ylim(ax1, [0 1]);
    title(ax1, 'Y');
    
    ax2 = subplot(4, 7, [22 23 24]);
    xlabel('z / m');
    yyaxis left
    plot(ax2, z, Z);
    ylabel(ax2, 'Z');
    yyaxis right
    plot(ax2, z, T);
    ylabel(ax2, 'T/K');
    
    ax3 = subplot(4, 7, [4 5 6]);
    plot(ax3, z, RR(iCH4, :), z, RR(iH2, :), z, RR(iO2, :), z, RR(iH2O, :));
    legend(ax3, 'RR_{CH_4}','RR_{H_2}','RR_{O_2}','RR_{H_2O}');
    title(ax3, '$$\dot{\omega}$$','Interpreter','latex');
    ylabel('Kg\cdotm^{-3}\cdots^{-1}');
    
    ax4 = subplot(4, 7, [11 12 13]);
    plot(ax4, z, RR(iCO, :), z, RR(iCO2, :));
    legend(ax4,'RR_{CO}','RR_{CO_2}');
    title(ax4, '$$\dot{\omega}$$','Interpreter','latex');
    ylabel('Kg\cdotm^{-3}\cdots^{-1}');
    
    ax5 = subplot(4, 7, [18 19 20]);
    plot(ax5, z, RR(iNO, :), z, RR(iNO2, :)*1e2);
    legend(ax5, 'RR_{NO}','1e2\cdotRR_{NO_2}');
    title(ax5, '$$\dot{\omega}$$','Interpreter','latex');
    ylabel('Kg\cdotm^{-3}\cdots^{-1}');
    
    ax6 = subplot(4, 7, [25 26 27]);
    plot(ax6, z, -RS);
    ylabel(ax6, 'J\cdotm^{-3}\cdots^{-1}');
    title(ax6, '$$-\sum{h_k\dot{\omega}_k}$$','Interpreter','latex');
    set(ax6,'YAxisLocation','right');
    xlabel('z / m');
    
    ax7 = subplot(4, 7, 7);
    plot(ax7, z, rho);
    ylabel(ax7, 'kg\cdotm^{-3}');
    title(ax7, '$$\rho$$','Interpreter','latex');
    set(ax7,'YAxisLocation','right');
    
    ax8 = subplot(4, 7, 14);
    plot(ax8, z, u);
    ylabel(ax8, 'm\cdots^{-1}');
    title(ax8, '$$u$$','Interpreter','latex');
    set(ax8,'YAxisLocation','right');
    
    ax9 = subplot(4, 7, 21);
    plot(ax9, z, V);
    ylabel(ax9, 's^{-1}');
    title(ax9, '$$V$$','Interpreter','latex');
    set(ax9,'YAxisLocation','right');
    
    ax10 = subplot(4, 7, 28);
    plot(ax10, z, Nbla);
    ylabel(ax10, 'Kg\cdotm^{-3}\cdots^{-2}');
    title(ax10, '$$\Lambda$$','Interpreter','latex');
    set(ax10,'YAxisLocation','right');
    xlabel('z / m');
end

function show_solution_diff(sol0, sol1)
    global z iCH4 iH2 iO2 iH2O iCO iCO2 iNO iNO2;
    
    rho0 = calculate_density(sol0);
    [u0, V0, T0, Nbla0, Y0] = mapback_solution_vector(sol0);
    rho1 = calculate_density(sol1);
    [u1, V1, T1, Nbla1, Y1] = mapback_solution_vector(sol1);
    
    delta_rho = rho1 - rho0;
    delta_u = u1 - u0;
    delta_V = V1 - V0;
    delta_T = T1 - T0;
    delta_Nbla = Nbla1 - Nbla0; 
    delta_Y = Y1 - Y0;
    
    h = figure(2);
    set(h, 'position', get(0,'ScreenSize'));
    
    subplot(3, 5, 1);
    plot(z, delta_rho);
    title('$$\Delta\rho$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('Kg\cdotm^{-3}');

    subplot(3, 5, 2);
    plot(z, delta_u);
    title('$$\Delta u$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('m\cdots^{-1}');
    
    subplot(3, 5, 3);
    plot(z, delta_V);
    title('$$\Delta V$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('s^{-1}');
    
    subplot(3, 5, 4);
    plot(z, delta_T);
    title('$$\Delta T$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('K');
    
    subplot(3, 5, 5);
    plot(z, delta_Nbla);
    title('$$\Delta\Lambda$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('Kg\cdotm^{-3}\cdots^{-2}');
    
    subplot(3, 5, 6);
    plot(z, delta_Y(iCH4, :));
    title('$$\Delta Y_{CH4}$$','Interpreter','latex');
    xlabel('z / m');
    
    subplot(3, 5, 7);
    plot(z, delta_Y(iH2, :));
    title('$$\Delta Y_{H2}$$','Interpreter','latex');
    xlabel('z / m');

    subplot(3, 5, 8);
    plot(z, delta_Y(iO2, :));
    title('$$\Delta Y_{O2}$$','Interpreter','latex');
    xlabel('z / m');

    subplot(3, 5, 9);
    plot(z, delta_Y(iCO2, :));
    title('$$\Delta Y_{CO2}$$','Interpreter','latex');
    xlabel('z / m');

    subplot(3, 5, 10);
    plot(z, delta_Y(iH2O, :));
    title('$$\Delta Y_{H2O}$$','Interpreter','latex');
    xlabel('z / m');

    subplot(3, 5, 11)
    plot(z, delta_Y(iNO, :))
    title('$$\Delta Y_{NO}$$','Interpreter','latex');
    xlabel('z / m')

    subplot(3, 5, 12)
    plot(z, delta_Y(iNO2, :))
    title('$$\Delta Y_{NO2}$$','Interpreter','latex');
    xlabel('z / m')
    
    subplot(3, 5, 13)
    plot(z, delta_Y(iCO, :))
    title('$$\Delta Y_{CO}$$','Interpreter','latex');
    xlabel('z / m');
end

function A = blktridiag(Amd,Asub,Asup,n)
% BLKTRIDIAG: computes a sparse (block) tridiagonal matrix with n blocks
% usage: A = BLKTRIDIAG(Amd,Asub,Asup,n)  % identical blocks
% usage: A = BLKTRIDIAG(Amd,Asub,Asup)    % a list of distinct blocks
%
% BLKTRIDIAG runs in two distinct modes. The first mode
% supplies three blocks, one for the main diagonal, and
% the super and subdiagonal blocks. These blocks will be
% replicated n times down the main diagonals of the matrix,
% and n-1 times down the sub and super diagonals.
% 
% The second mode is to supply a list of distinct blocks
% for each diagonal, as planes of 3d arrays. No replication
% factor is needed in this mode.
%
% arguments: (input mode 1)
%  Amd  - pxq array, forming the main diagonal blocks
%
%  Asub - pxq array, sub diagonal block
%         Asub must be the same size and shape as Amd
%
%  Asup - pxq array, super diagonal block
%         Asup must be the same size and shape as Amd
%
%  n    - scalar integer, defines the number of blocks
%         When n == 1, only a single block will be formed, A == Amd
%
% arguments: (input mode 2)
%  Amd  - pxqxn array, a list of n distinct pxq arrays
%         Each plane of Amd corresponds to a single block
%         on the main diagonal.
%
%  Asub - pxqx(n-1) array, a list of n-1 distinct pxq arrays,
%         Each plane of Asub corresponds to a single block
%         on the sub-diagonal.
%
%  Asup - pxqx(n-1) array, a list of n-1 distinct pxq arrays,
%         Each plane of Asup corresponds to a single block
%         on the super-diagonal.
%
% Note: the sizes of Amd, Asub, and Asup must be consistent
% with each other, or an error will be generated.
% 
% arguments: (output)
%  A    - (n*p by n*q) SPARSE block tridiagonal array
%         If you prefer that A be full, use A=full(A) afterwards.
%
%
% Example 1:
%  Compute the simple 10x10 tridiagonal matrix, with 2 on the
%  diagonal, -1 on the off diagonal.
%
%  A = blktridiag(2,-1,-1,10);
%
%
% Example 2:
%  Compute the 5x5 lower bi-diagonal matrix, with blocks of
%  [1 1;1 1] on the main diagonal, [2 2;2 2] on the sub-diagonal,
%  and blocks of zeros above.
%
%  A = blktridiag(ones(2),2*ones(2),zeros(2),5);
%
% Example 3:
%  Compute the 3x6 tridiagonal matrix, with non-square blocks
%  that vary along the main diagonal, [2 2] on the sub-diagonal,
%  and [1 1] on the super-diagonal. Note that all blocks must have
%  the same shape.
%
%  A = blktridiag(rand(1,2,3),2*ones(1,2,2),ones(1,2,2));
%
%
% See also: blkdiag, spdiags, diag
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 4.0
% Original release date: 4/01/06
% Current release date: 12/14/07
% Which mode of operation are we in?
if nargin==4
  % replicated block mode
  
  % verify the inputs in this mode are 2-d arrays.
  if (length(size(Amd))~=2) || ...
     (length(size(Asub))~=2) || ...
     (length(size(Asup))~=2) 
    error 'Inputs must be 2d arrays if a replication factor is provided'
  end
  
  % get block sizes, check for consistency
  [p,q] = size(Amd);
  if isempty(Amd)
    error 'Blocks must be non-empty arrays or scalars'
  end
  if any(size(Amd)~=size(Asub)) || any(size(Amd)~=size(Asup))
    error 'Amd, Asub, Asup are not identical in size'
  end
  if isempty(n) || (length(n)>1) || (n<1) || (n~=floor(n))
    error 'n must be a positive scalar integer'
  end
  
  % scalar inputs?
  % since p and q are integers...
  if (p*q)==1
    if n==1
      A = Amd;
    else
      % faster as Jos points out
      A = spdiags(repmat([Asub Amd Asup],n,1),-1:1,n,n);
    end
    % no need to go any farther
    return
  end
  
  % use sparse. the main diagonal elements of each array are...
  v = repmat(Amd(:),n,1);
  % then the sub and super diagonal blocks.
  if n>1
    % sub-diagonal
    v=[v;repmat(Asub(:),n-1,1)];
    
    % super-diagonal
    v=[v;repmat(Asup(:),n-1,1)];
  end
  
elseif nargin==3
  % non-replicated blocks, supplied as planes of a 3-d array
  
  % get block sizes, check for consistency
  [p,q,n] = size(Amd);
  if isempty(Amd)
    error 'Blocks must be (non-empty) arrays or scalars'
  end
  
  if (p~=size(Asub,1)) || (q~=size(Asub,2)) || (p~=size(Asup,1)) || (q~=size(Asup,2))
    error 'Amd, Asub, Asup do not have the same size blocks'
  end
  if (n>1) && (((n-1) ~= size(Asub,3)) || ((n-1) ~= size(Asup,3)))
    error 'Asub and Asup must each have one less block than Amd'
  end
  
  % scalar inputs?
  if (p*q)==1
    if n==1
      A = Amd(1);
    else
      % best to just use spdiags
      A = spdiags([[Asub(:);0], Amd(:), [0;Asup(:)]],-1:1,n,n);
    end
    % no need to go any farther
    return
  end
  
  % The main diagonal elements
  v = Amd(:);
  % then the sub and super diagonal blocks.
  if n>1
    % sub-diagonal
    v=[v;Asub(:)];
    % super-diagonal
    v=[v;Asup(:)];
  end
else
  % must have 3 or 4 arguments
  error 'Must have either 3 or 4 arguments to BLKTRIDIAG'
end
% now generate the index arrays. first the main diagonal
[ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-1);
rind = 1+ind1(:)+p*ind3(:);
cind = 1+ind2(:)+q*ind3(:);
% then the sub and super diagonal blocks.
if n>1
  % sub-diagonal
  [ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-2);
  rind = [rind;1+p+ind1(:)+p*ind3(:)];
  cind = [cind;1+ind2(:)+q*ind3(:)];
  % super-diagonal
  rind = [rind;1+ind1(:)+p*ind3(:)];
  cind = [cind;1+q+ind2(:)+q*ind3(:)];
end
% build the final array all in one call to sparse
A = sparse(rind,cind,v,n*p,n*q);
end

function ret = gauss1d(x)
    global filtering_sigma;

    ret = 1.0 / (sqrt(2*pi)*filtering_sigma) * exp(-0.5*(x/filtering_sigma)^2);
end

function ret = calculate_filtered_wdot(wdot)
    global z N K gauss_weight;
    
    ret = zeros(K, N);
    for i = 1:N
        for k = 1:K
            ret(k, i) = trapz(z, gauss_weight(i, :) .* wdot(k, :)); 
        end
    end
end
