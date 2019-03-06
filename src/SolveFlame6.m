function SolveFlame6(mdot_f, mdot_o, L, N, output_dir)
%%Solve the counter-flow diffusion flame with the damped Newton method.
%  mdot_f: Mass flux of fuel at left side, Unit: Kg/(m^2 * s).
%  mdot_o: Mass flux of oxidizer at right side, Unit: Kg/(m^2 * s).
%  L: Length of domain, Unit: m.
%  N: Total num of grid points distributed uniformly.
%  output_dir: Target directory where transformed data files are stored.

    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    
    %% Setup solution environment
    Le = 1.0;
    P = oneatm; %The constant pressure, Pa
    mdot_L = mdot_f ; %Fuel stream, Kg/(m^2 * s)
    mdot_R = -mdot_o; %Air stream, Kg/(m^2 * s)
    
    r = sqrt(eps);
    a = sqrt(eps);
    
    gas = GRI30('Mix');
    ch4_idx = speciesIndex(gas, 'CH4');
    o2_idx = speciesIndex(gas, 'O2');
    n2_idx = speciesIndex(gas, 'N2');
    co2_idx = speciesIndex(gas, 'CO2');
    h2o_idx = speciesIndex(gas, 'H2O');
    h2_idx = speciesIndex(gas, 'H2');
    co_idx = speciesIndex(gas, 'CO');
    no_idx = speciesIndex(gas, 'NO');
    MW = molecularWeights(gas); %Molecular weight, Kg/Kmol
    NAME = speciesNames(gas); %Name of each species
    K = nSpecies(gas); %Total num of species
    
    zL = -L/2; %Position of left endpoint, m
    zR = L/2; %Position of right endpoint, m
    z = linspace(zL, zR, N); %Coordinates for each point, m
    
    CUR = 1;
    NEXT = 2;
    
    rho = zeros(2, N); % Kg / m^3
    u = zeros(2, N); % m/s
    V = zeros(2, N);
    Nbla = zeros(2, N); %The eigenvalue
    T = zeros(2, N); % K
    Y = zeros(2, K, N); %Mass fraction
    
    mu = zeros(2, N); %Viscosity, Pa * s = Kg / (m * s)
    cp = zeros(2, N); %Specific heat, J / (Kg * K)
    lambda = zeros(2, N); %Thermal conductivity, W / (m * K)
    D = zeros(2, K, N); %Binary diffusion coefficients, m^2 / s
    
    RS = zeros(2, N); %Energy source due to chemical reaction, J / (m^3 * s)
    RR = zeros(2, K, N); %Chemical reaction rate, Kg / (m^3 * s)
    
    dVdz = zeros(2, N);
    dTdz = zeros(2, N);
    dYdz = zeros(2, K, N);
    ddVddz = zeros(2, N);
    ddTddz = zeros(2, N);
    ddYddz = zeros(2, K, N);
    
    C = 4+K;  % Num of unknowns per node
    U = C*(N-2);
    phi=zeros(2, U);
    F= zeros(2, U);
    J = zeros(U, U);
    
    %% Initialize using the Burke-Shumann solution
    [zc, uc, Tc,  Yc] = DiffFlameSim(L, P, 300.0, mdot_L, -mdot_R);
    tpts = z - zL;

    T(CUR, :) = spline(zc, Tc, tpts);
    u(CUR, :) = spline(zc, uc, tpts);
    for k = 1:K
        Y(CUR, k, :) = spline(zc, Yc(k, :), tpts);
    end

    for i = 1:N
        rho(CUR, i) = P / (gasconstant * T(CUR, i) * sum(squeeze(Y(CUR, :, i)) ./ MW'));
    end

    %V set to 0 at boundary as no vertical slip physically
    V(CUR, :) = -df_upwind(rho(CUR, :) .* u(CUR, :), z, u(CUR, :)) ./ (2 * rho(CUR, :));
    V(CUR,1)=0.0;
    V(CUR,N)=0.0;
    
    for i = 1:N
        set(gas, 'T', T(CUR, i), 'P', P, 'Y', squeeze(Y(CUR, :, i)));
        mu(CUR,i) = viscosity(gas);
    end

    %Select initial guess of the eigenvalue
    dVdz(CUR, :) = df_upwind(V(CUR, :), z, u(CUR, :));
    ddVddz(CUR, :) = ddf(V(CUR, :), z);
    lhs1 = dot(rho(CUR, 2:N-1) .* u(CUR, 2:N-1), dVdz(CUR, 2:N-1));
    lhs2 = dot(rho(CUR, 2:N-1) .* V(CUR, 2:N-1), V(CUR, 2:N-1));
    rhs2 = dot(mu(CUR, 2:N-1),  ddVddz(CUR, 2:N-1));
    Nbla(CUR, :) = (rhs2 - lhs1 - lhs2) / (N-2);
    
    %% Solve 
    global_converged = false;
    global_iter_cnt = 0;
    while(~global_converged)
        global_iter_cnt = global_iter_cnt + 1;
        
        %% Update properties
        for i = 1:N
            local_T = T(CUR, i);
            set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(CUR, :, i)));
            mu(CUR,i) = viscosity(gas);
            lambda(CUR,i) = thermalConductivity(gas);
            cp(CUR,i) = cp_mass(gas);
            D(CUR, :, i) = lambda(CUR,i) / (rho(CUR, i) * cp(CUR, i) * Le);
            w = netProdRates(gas); % kmol / (m^3 * s)
            h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol
            RS(CUR,i) = dot(w, h); % J / (m^3 * s)
            RR(CUR,:, i) = w.* MW; % Kg / (m^3 * s)
        end
                
        %% Calcaulate derivatives
        dVdz(CUR, :) = df_upwind(V(CUR, :), z, u(CUR, :));
        dTdz(CUR, :) = df_upwind(T(CUR, :), z, u(CUR, :));
        for k = 1:K
            dYdz(CUR, k, :) = df_upwind(Y(CUR, k, :), z, u(CUR, :));
        end
        ddVddz(CUR, :) = ddf(V(CUR, :), z);
        ddTddz(CUR, :) = ddf(T(CUR, :), z);
         for k = 1:K
            ddYddz(CUR, k, :) = ddf(Y(CUR, k, :), z);
        end
       
        %% Calculate residuals
        cnt = 1;
        for i = 2:N-1
            % u
            F(CUR, cnt) = (rho(CUR, i+1) * u(CUR, i+1) - rho(CUR, i) * u(CUR, i))/(z(i+1) - z(i)) + rho(CUR, i) * V(CUR, i) + rho(CUR, i+1) * V(CUR, i+1);
            cnt = cnt +1;
            % V
            F(CUR, cnt) = rho(CUR, i)*u(CUR, i)*dVdz(CUR,i)+rho(CUR,i)*V(CUR,i)^2 + Nbla(CUR, i) - mu(CUR,i)*ddVddz(CUR,i);
            cnt = cnt + 1;
            % T
            F(CUR, cnt) = rho(CUR, i) * cp(CUR,i) * u(CUR, i) * dTdz(CUR,i) -lambda(CUR,i) * ddTddz(CUR,i) + RS(CUR,i);
            cnt = cnt + 1;
            % Lambda
            if i == 2
                F(CUR, cnt) = rho(CUR, i) * u(CUR, i) - mdot_L;
            else
                F(CUR, cnt) = Nbla(CUR, i) - Nbla(CUR, i-1);
            end
            cnt = cnt + 1;
            % Yk
            for k = 1:K
                F(CUR, cnt) = rho(CUR, i)*u(CUR, i)*dYdz(CUR,k, i)-D(CUR,k,i)*ddYddz(CUR,k,i)-RR(CUR,k, i);
                cnt = cnt + 1;
            end
        end
        
        %% Calculate the Jacobian by finite difference perturbations
        for j = 1:U
            %% var base init
            rho(NEXT, :) = rho(CUR, :);
            u(NEXT, :) = u(CUR, :);
            V(NEXT, :) = V(CUR, :);
            Nbla(NEXT, :) = Nbla(CUR, :);
            T(NEXT, :) = T(CUR, :);
            Y(NEXT, :, :) = Y(CUR, :, :);
            %% add perturbation
            node_idx = ceil(j / C) + 1;
            var_idx = mod(j-1, C);
            if  var_idx == 0
                delta = perturbation_delta(r, a, u(CUR, node_idx));
                u(NEXT, node_idx) = u(CUR, node_idx) + delta;
            elseif var_idx == 1
                delta = perturbation_delta(r, a, V(CUR, node_idx));
                V(NEXT, node_idx) = V(CUR, node_idx) + delta;
            elseif var_idx == 2
                delta = perturbation_delta(r, a, T(CUR, node_idx));
                T(NEXT, node_idx) = T(CUR, node_idx) + delta;
            elseif var_idx == 3
                delta = perturbation_delta(r, a, Nbla(CUR, node_idx));
                Nbla(NEXT, node_idx) = Nbla(CUR, node_idx) + delta;
            else
                spec_idx = var_idx - 3;
                delta = perturbation_delta(r, a, Y(CUR, spec_idx, node_idx));
                Y(NEXT, spec_idx, node_idx) = Y(CUR, spec_idx, node_idx) + delta;
            end
            %% update properties
            for i = 1:N
                local_T = T(NEXT, i);
                set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(NEXT, :, i)));
                mu(NEXT,i) = viscosity(gas);
                lambda(NEXT,i) = thermalConductivity(gas);
                cp(NEXT,i) = cp_mass(gas);
                D(NEXT,:, i) = lambda(NEXT,i) / (rho(NEXT, i) * cp(NEXT, i) * Le);
                w = netProdRates(gas); % kmol / (m^3 * s)
                h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol
                RS(NEXT,i) = dot(w, h); % J / (m^3 * s)
                RR(NEXT,:, i) = w.* MW; % Kg / (m^3 * s)
            end
            %% compute derivatives
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
            %% calculate residuals
            cnt = 1;
            for i = 2:N-1
                % u
                F(NEXT, cnt) = (rho(NEXT, i+1) * u(NEXT, i+1) - rho(NEXT, i) * u(NEXT, i))/(z(i+1) - z(i)) + rho(NEXT, i) * V(NEXT, i) + rho(NEXT, i+1) * V(NEXT, i+1);
                cnt = cnt +1;
                %V
                F(NEXT, cnt) = rho(NEXT, i)*u(NEXT, i)*dVdz(NEXT,i)+rho(NEXT,i)*V(NEXT,i)^2 + Nbla(NEXT, i) - mu(NEXT,i)*ddVddz(NEXT,i);
                cnt = cnt + 1;
                %T
                F(NEXT, cnt) = rho(NEXT, i) * cp(NEXT,i) * u(NEXT, i) * dTdz(NEXT,i) -lambda(NEXT,i) * ddTddz(NEXT,i) + RS(NEXT,i);
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
            %% update column
            J(:, j) = (F(NEXT, :) - F(CUR, :))/delta;
        end
        
        %% Solve the Jacobian
        dphi = solBlkDiagMat(J, F(CUR, :)', C);
        
    end
    
    %% Output        

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

function [z, u, T, y] = DiffFlameSim(domain_length, p, tin, mdot_f, mdot_o)
    runtime = cputime;  % Record the starting time

    initial_grid = domain_length*linspace(0,1, 51);  % Units: m
    tol_ss    = [1e-4 1e-9];        % [rtol atol] for steady-state problem
    tol_ts    = [1e-4 1e-9];        % [rtol atol] for time stepping
    loglevel  = 1;                      % Amount of diagnostic output (0 to 5)
    refine_grid = 1;                    % 1 to enable refinement, 0 to disable

    fuel = GRI30('Mix');
    ox = GRI30('Mix');
    oxcomp = 'O2:0.21, N2:0.78'; % Air composition
    fuelcomp = 'CH4:0.5, H2:0.5'; % Fuel composition

    set(fuel,'T', tin, 'P', p, 'X', fuelcomp);
    set(ox,'T',tin,'P',p,'X', oxcomp);

    f = AxisymmetricFlow(fuel,'flow');
    set(f, 'P', p, 'grid', initial_grid);
    set(f, 'tol', tol_ss, 'tol-time', tol_ts);

    % Set the oxidizer inlet.
    inlet_o = Inlet('air_inlet');
    set(inlet_o, 'T', tin, 'MassFlux', mdot_o, 'X', oxcomp);

    % Set the fuel inlet.
    inlet_f = Inlet('fuel_inlet');
    set(inlet_f, 'T', tin, 'MassFlux', mdot_f, 'X', fuelcomp);

    fl = CounterFlowDiffusionFlame(inlet_f, f, inlet_o, fuel, ox, 'O2');

    solve(fl, loglevel, 0);

    enableEnergy(f);
    setRefineCriteria(fl, 2, 200.0, 0.1, 0.2);
    solve(fl, loglevel, refine_grid);
    
    writeStats(fl);
    elapsed = cputime - runtime;
    fprintf('Elapsed CPU time: %10.4g\n',elapsed);

    z = grid(fl, 'flow'); % Get grid points of flow
    spec = speciesNames(fuel); % Get species names in gas
    u = solution(fl, 'flow', 'u'); 
    T = solution(fl, 'flow', 'T'); % Get temperature solution
    y = zeros(length(spec), length(z));
    for i = 1:length(spec)
        y(i,:) = solution(fl, 'flow', spec{i}); % Get mass fraction of all species from solution
    end
end
