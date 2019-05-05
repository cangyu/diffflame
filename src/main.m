% Solve unsteady opposed diffusion flame using the damped Newton method.
clear; close all; clc;

%% Init
global N;
global K;
global z;
global rho;
global u;
global V;
global P;
global T;
global Nbla;
global Y;
global mu;
global Cp;
global lambda;
global D;
global RS;
global RR;
global dVdz;
global dTdz;
global dYdz;
global ddVddz;
global ddTddz;
global ddYddz;
global C;
global U;
global phi;
global F;
global gas;
global MW;
global NAME;
global mdot_L;
global mdot_R;

mdot_f = 0.1; % Mass flux of fuel at left side, Kg/(m^2 * s).
mdot_o = 0.3; % Mass flux of oxidizer at right side, Kg/(m^2 * s).
L = 0.05; % Length of domain, m.

case_str = sprintf("mf=%g_mo=%g_L=%g", mdot_f, mdot_o, L);
raw_data_path = "../data/" + case_str + "_raw.txt";
trans_data_path = "../data/" + case_str + "_transformed.txt";

raw_data = readtable(raw_data_path);
trans_data = readtable(trans_data_path);
if height(raw_data) ~= height(trans_data)
    error("Inconsistent input data!");
end

N = height(raw_data); % Total num of points
zL = raw_data{1, 1}; % Position of left endpoint, m
zR = raw_data{N, 1}; % Position of right endpoint, m
if zR - zL ~= L
    error("Inconsistent domain!");
end
z = raw_data{:, 1}; % Coordinates for each point, m

gas = GRI30('Mix');
MW = molecularWeights(gas); % Molecular weight, Kg/Kmol
NAME = speciesNames(gas); % Name of each species
K = nSpecies(gas); % Total num of species

ch4_idx = speciesIndex(gas, 'CH4');
h2_idx = speciesIndex(gas, 'H2');
o2_idx = speciesIndex(gas, 'O2');
n2_idx = speciesIndex(gas, 'N2');
h2o_idx = speciesIndex(gas, 'H2O');
co_idx = speciesIndex(gas, 'CO');
co2_idx = speciesIndex(gas, 'CO2');
no_idx = speciesIndex(gas, 'NO');
no2_idx = speciesIndex(gas, 'NO2');

CUR = 1;
NEXT = 2;

rho = zeros(2, N); % Kg / m^3
u = zeros(2, N); % m/s
V = zeros(2, N); % 1/s
P = oneatm; % The constant pressure, Pa
T = zeros(2, N); % K
Nbla = zeros(2, N); %The eigenvalue
Y = zeros(2, K, N); %Mass fraction

mu = zeros(N, 1); %Viscosity, Pa * s = Kg / (m * s)
Cp = zeros(N, 1); %Specific heat, J / (Kg * K)
lambda = zeros(N, 1); %Thermal conductivity, W / (m * K)
D = zeros(K, N); %Binary diffusion coefficients, m^2 / s

RS = zeros(N, 1); %Energy source due to chemical reaction, J / (m^3 * s)
RR = zeros(K, N); %Chemical reaction rate, Kg / (m^3 * s)

dVdz = zeros(N, 1);
dTdz = zeros(N, 1);
dYdz = zeros(K, N);
ddVddz = zeros(N, 1);
ddTddz = zeros(N, 1);
ddYddz = zeros(K, N);

C = 4+K;  % Num of unknowns per node
U = C*N; % Total num of unknowns
phi=zeros(2, U); % Solution vector
F= zeros(2, U); % Residual vector
J = zeros(U, U); % Jacobian matrix

Le = 1.0; % Unity Lewis
mdot_L = mdot_f ; % Fuel stream, Kg/(m^2 * s)
mdot_R = -mdot_o; % Air stream, Kg/(m^2 * s)
rtol = 1e-5; % Relative tolerance
atol = sqrt(eps); % Absolute tolerance

rho(CUR, :) = trans_data{:, 1};
u(CUR, :) = raw_data{:, 2};
V(CUR, :) = raw_data{:, 3};
T(CUR, :) = raw_data{:, 4};
Nbla(CUR, :) = raw_data{:, 5};
for k = 1:K
    Y(CUR, k, :) = raw_data{:, 5+k};
end

%%  Solve 
global_converged = false;
global_iter_cnt = 0;
while(~global_converged)
    global_iter_cnt = global_iter_cnt + 1;

    construct_solution_vector(CUR, CUR);

    % Calculate the Jacobian by finite difference perturbations
    for j = 1:U
        % var base init
        rho(NEXT, :) = rho(CUR, :);
        u(NEXT, :) = u(CUR, :);
        V(NEXT, :) = V(CUR, :);
        Nbla(NEXT, :) = Nbla(CUR, :);
        T(NEXT, :) = T(CUR, :);
        Y(NEXT, :, :) = Y(CUR, :, :);
        % add perturbation
        node_idx = ceil(j / C) + 1;
        var_idx = mod(j-1, C);
        if  var_idx == 0
            delta = perturbation_delta(rtol, atol, u(CUR, node_idx));
            u(NEXT, node_idx) = u(CUR, node_idx) + delta;
        elseif var_idx == 1
            delta = perturbation_delta(rtol, atol, V(CUR, node_idx));
            V(NEXT, node_idx) = V(CUR, node_idx) + delta;
        elseif var_idx == 2
            delta = perturbation_delta(rtol, atol, T(CUR, node_idx));
            T(NEXT, node_idx) = T(CUR, node_idx) + delta;
        elseif var_idx == 3
            delta = perturbation_delta(rtol, atol, Nbla(CUR, node_idx));
            Nbla(NEXT, node_idx) = Nbla(CUR, node_idx) + delta;
        else
            spec_idx = var_idx - 3;
            delta = perturbation_delta(rtol, atol, Y(CUR, spec_idx, node_idx));
            Y(NEXT, spec_idx, node_idx) = Y(CUR, spec_idx, node_idx) + delta;
        end
        % update properties
        for i = 1:N
            local_T = T(NEXT, i);
            set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(NEXT, :, i)));
            mu(NEXT,i) = viscosity(gas);
            lambda(NEXT,i) = thermalConductivity(gas);
            Cp(NEXT,i) = cp_mass(gas);
            D(NEXT,:, i) = lambda(NEXT,i) / (rho(NEXT, i) * Cp(NEXT, i) * Le);
            w = netProdRates(gas); % kmol / (m^3 * s)
            h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol
            RS(NEXT,i) = dot(w, h); % J / (m^3 * s)
            RR(NEXT,:, i) = w.* MW; % Kg / (m^3 * s)
        end
        % compute derivatives
        dVdz(NEXT, :) = df_upwind(V(NEXT, :), z,  u(NEXT, :));
        dTdz(NEXT, :) = df_upwind(T(NEXT, :), z, u(NEXT, :));
        for k = 1:K
            dYdz(NEXT, k, :) = df_upwind(Y(NEXT, k, :), z, u(NEXT, :));
        end
        ddVddz(NEXT, :) = ddf(V(NEXT, :), z);
        ddTddz(NEXT, :) = ddf(T(NEXT, :), z);
         for k = 1:K
            ddYddz(NEXT, k, :) = ddf(Y(NEXT, k, :), z);
        end
        % calculate residuals
        cnt = 1;
        for i = 2:N-1
            % u
            F(NEXT, cnt) = (rho(NEXT, i+1) * u(NEXT, i+1) - rho(NEXT, i) * u(NEXT, i))/(z(i+1) - z(i)) + rho(NEXT, i) * V(NEXT, i) + rho(NEXT, i+1) * V(NEXT, i+1);
            cnt = cnt +1;
            %V
            F(NEXT, cnt) = rho(NEXT, i)*u(NEXT, i)*dVdz(NEXT,i)+rho(NEXT,i)*V(NEXT,i)^2 + Nbla(NEXT, i) - mu(NEXT,i)*ddVddz(NEXT,i);
            cnt = cnt + 1;
            %T
            F(NEXT, cnt) = rho(NEXT, i) * Cp(NEXT,i) * u(NEXT, i) * dTdz(NEXT,i) -lambda(NEXT,i) * ddTddz(NEXT,i) + RS(NEXT,i);
            cnt = cnt + 1;
            % Lambda
            if i == 2
                F(NEXT, cnt) = rho(NEXT, i) * u(NEXT, i) - mdot_L;
            else
                F(NEXT, cnt) = Nbla(NEXT, i) - Nbla(NEXT, i-1);
            end
            cnt = cnt + 1;
            %Yk
            for k = 1:K
                F(NEXT, cnt) = rho(NEXT, i)*u(NEXT, i)*dYdz(NEXT,k, i)-D(NEXT,k,i)*ddYddz(NEXT,k,i)-RR(NEXT,k, i);
                cnt = cnt + 1;
            end
        end
        % update column
        J(:, j) = (F(NEXT, :) - F(CUR, :))/delta;
    end

    %Solve the Jacobian
    dphi = solBlkDiagMat(J, F(CUR, :)', C);

end

%% Aux functions
function construct_solution_vector(src_level, dst_level)
    global u;
    global V;
    global T;
    global Nbla;
    global Y;
    global phi;
    global N;
    global K;
    global U;
    
    cnt = 1;
    for i = 1:N
        phi(dst_level, cnt) = u(src_level, i); cnt = cnt + 1;
        phi(dst_level, cnt) = V(src_level, i); cnt = cnt + 1;
        phi(dst_level, cnt) = T(src_level, i); cnt = cnt + 1;
        phi(dst_level, cnt) = Nbla(src_level, i); cnt = cnt + 1;
        for k = 1:K
            phi(dst_level, cnt) = Y(src_level, k, i); cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
    end
end

function mapback_solution_vector(src_level, dst_level)
    global N;
    global K;
    global u;
    global V;
    global T;
    global Nbla;
    global Y;
    global U;
    global phi;
    
    cnt = 1;
    for i = 1:N
        u(dst_level, i) = phi(src_level, cnt); cnt = cnt + 1;
        V(dst_level, i) = phi(src_level, cnt); cnt = cnt + 1;
        T(dst_level, i) = phi(src_level, cnt); cnt = cnt + 1;
        Nbla(dst_level, i) = phi(src_level, cnt); cnt = cnt + 1;
        for k = 1:K
            Y(dst_level, k, i) = phi(src_level, cnt); cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
    end
end

function calculate_residual_vector(src_level, dst_level)
    global N;
    global K;
    global z;
    global rho;
    global u;
    global V;
    global P;
    global T;
    global Nbla;
    global Y;
    global mu;
    global Cp;
    global lambda;
    global D;
    global RS;
    global RR;
    global dVdz;
    global dTdz;
    global dYdz;
    global ddVddz;
    global ddTddz;
    global ddYddz;
    global U;
    global F;
    global gas;
    global MW;
    global mdot_L;
    
    % Update properties
    for i = 1:N
        loc_T = T(src_level, i);
        set(gas, 'T', loc_T, 'P', P, 'Y', squeeze(Y(src_level, :, i)));
        mu(i) = viscosity(gas);
        lambda(i) = thermalConductivity(gas);
        Cp(i) = cp_mass(gas);
        D(:, i) = lambda(i) / (rho(src_level, i) * Cp(i)); % Compute from Unity Lewis
        w = netProdRates(gas); % kmol / (m^3 * s)
        h = enthalpies_RT(gas) * loc_T * gasconstant; % J/Kmol
        RS(i) = dot(w, h); % J / (m^3 * s)
        RR(:, i) = w.* MW; % Kg / (m^3 * s)
    end
    
    % Calcaulate derivatives
    dVdz(:) = df_upwind(V(src_level, :), z, u(src_level, :));
    dTdz(:) = df_upwind(T(src_level, :), z, u(src_level, :));
    for k = 1:K
        dYdz(k, :) = df_upwind(Y(src_level, k, :), z, u(src_level, :));
    end
    ddVddz(:) = ddf(V(src_level, :), z);
    ddTddz(:) = ddf(T(src_level, :), z);
     for k = 1:K
        ddYddz(k, :) = ddf(Y(src_level, k, :), z);
     end
    
    % Calculate residuals
    cnt = 1;
    for i = 1:N
        F(dst_level, cnt) = (rho(src_level, i+1) * u(src_level, i+1) - rho(src_level, i) * u(src_level, i))/(z(i+1) - z(i)) + rho(src_level, i) * V(src_level, i) + rho(src_level, i+1) * V(src_level, i+1);
        cnt = cnt +1;
        F(dst_level, cnt) = rho(src_level, i)*u(src_level, i)*dVdz(i)+rho(src_level,i)*V(src_level,i)^2 + Nbla(src_level, i) - mu(i)*ddVddz(i);
        cnt = cnt + 1;
        F(dst_level, cnt) = rho(src_level, i) * Cp(i) * u(src_level, i) * dTdz(i) -lambda(i) * ddTddz(i) + RS(i);
        cnt = cnt + 1;
        if i == 2
            F(dst_level, cnt) = rho(src_level, i) * u(src_level, i) - mdot_L;
        else
            F(dst_level, cnt) = Nbla(src_level, i) - Nbla(src_level, i-1);
        end
        cnt = cnt + 1;
        for k = 1:K
            F(dst_level, cnt) = rho(src_level, i)*u(src_level, i)*dYdz(k, i)-D(k,i)*ddYddz(k,i)-RR(k, i);
            cnt = cnt + 1;
        end
    end
    
    if cnt ~= U+1
        error("Internal error!");
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

function ret = perturbation_delta(rel_perturb, abs_perturb, x)
    ret = rel_perturb * x + abs_perturb;
end

function x = solBlkDiagMat(B, b, bandwidth)
    x = B\b;
end
