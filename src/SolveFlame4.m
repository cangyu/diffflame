function SolveFlame4(mdot_f, mdot_o, L, N, output_dir)
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
    
    unknown_per_node = 1+K+1;
    unknown_num = unknown_per_node*N;
    phi=zeros(2, unknown_num);
    F= zeros(2, unkown_num);
    J = zeros(unkown_num, unkown_num);
    
    %% Initialize using the Burke-Shumann solution
    
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
            D(CUR,:, i) = mixDiffCoeffs(gas);
            w = netProdRates(gas); % kmol / (m^3 * s)
            h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol
            RS(CUR,i) = -dot(w, h); % J / (m^3 * s)
            RR(CUR,:, i) = w.* MW; % Kg / (m^3 * s)
        end
                
        %% Calcaulate derivatives
        dVdz(CUR, :) = df(V(CUR, :), z, N);
        dTdz(CUR, :) = df(T(CUR, :), z, N);
        for k = 1:K
            dYdz(CUR, k, :) = df(Y(CUR, k, :), z, N);
        end
        ddVddz(CUR, :) = ddf(V(CUR, :), z, N);
        ddTddz(CUR, :) = ddf(T(CUR, :), z, N);
         for k = 1:K
            ddYddz(CUR, k, :) = ddf(Y(CUR, k, :), z, N);
        end
       
        %% Calculate residuals
        cnt = 1;
        for i = 1:N
            F(CUR, cnt) = rho(CUR, i) * cp(CUR,i) * u(CUR, i) * dTdz(CUR,i) -lambda(CUR,i) * ddTddz(CUR,i) - RS(CUR,i);
            cnt = cnt + 1;
            for k = 1:K
                F(CUR, cnt) = rho(CUR, i)*u(CUR, i)*dYdz(CUR,k, i)-D(CUR,k,i)*ddYddz(CUR,k,i)-RR(CUR,k, i);
                cnt = cnt + 1;
            end
            F(CUR, cnt) = rho(CUR, i)*u(CUR, i)*dVdz(CUR,i)+rho(CUR,i)*V(CUR,i)^2 + Nbla(CUR, i) - mu(CUR,i)*ddVddz(CUR,i);
            cnt = cnt + 1;
        end
        
        %% Calculate the Jacobian by finite difference perturbations
        for j = 1:unknown_num
            %% var base init
            rho(NEXT, :) = rho(CUR, :);
            u(NEXT, :) = u(CUR, :);
            V(NEXT, :) = V(CUR, :);
            Nbla(NEXT, :) = Nbla(CUR, :);
            T(NEXT, :) = T(CUR, :);
            Y(NEXT, :, :) = Y(CUR, :, :);
            %% add perturbation
            node_idx = ceil(j / unknown_per_node);
            var_idx = mod(j, unknown_per_node);
            if  var_idx == 1
                delta = perturbation_delta(r, a, rho(CUR, node_idx));
                rho(NEXT, node_idx) = rho(CUR, node_idx) + delta;
            elseif var_idx == 0
                delta = perturbation_delta(r, a, Nbla(CUR, node_idx));
                Nbla(NEXT, node_idx) = Nbla(CUR, node_idx) + delta;
            else
                spec_idx = var_idx - 1;
                delta = perturbation_delta(r, a, Y(CUR, spec_idx, node_idx));
                Y(NEXT, spec_idx, node_idx) = Y(NEXT, spec_idx, node_idx) + delta;
            end
            %% update properties
            for i = 1:N
                local_T = T(NEXT, i);
                set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(NEXT, :, i)));
                mu(NEXT,i) = viscosity(gas);
                lambda(NEXT,i) = thermalConductivity(gas);
                cp(NEXT,i) = cp_mass(gas);
                D(NEXT,:, i) = mixDiffCoeffs(gas);
                w = netProdRates(gas); % kmol / (m^3 * s)
                h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol
                RS(NEXT,i) = -dot(w, h); % J / (m^3 * s)
                RR(NEXT,:, i) = w.* MW; % Kg / (m^3 * s)
            end
            %% compute derivatives
            dVdz(NEXT, :) = df(V(NEXT, :), z, N);
            dTdz(NEXT, :) = df(T(NEXT, :), z, N);
            for k = 1:K
                dYdz(NEXT, k, :) = df(Y(NEXT, k, :), z, N);
            end
            ddVddz(NEXT, :) = ddf(V(NEXT, :), z, N);
            ddTddz(NEXT, :) = ddf(T(NEXT, :), z, N);
             for k = 1:K
                ddYddz(NEXT, k, :) = ddf(Y(NEXT, k, :), z, N);
            end
            %% calculate residuals
            cnt = 1;
            for i = 1:N
                F(NEXT, cnt) = rho(NEXT, i) * cp(NEXT,i) * u(NEXT, i) * dTdz(NEXT,i) -lambda(NEXT,i) * ddTddz(NEXT,i) - RS(NEXT,i);
                cnt = cnt + 1;
                for k = 1:K
                    F(NEXT, cnt) = rho(NEXT, i)*u(NEXT, i)*dYdz(NEXT,k, i)-D(NEXT,k,i)*ddYddz(NEXT,k,i)-RR(NEXT,k, i);
                    cnt = cnt + 1;
                end
                F(NEXT, cnt) = rho(NEXT, i)*u(NEXT, i)*dVdz(NEXT,i)+rho(NEXT,i)*V(NEXT,i)^2 + Nbla(NEXT, i) - mu(NEXT,i)*ddVddz(NEXT,i);
                cnt = cnt + 1;
            end
            %% update column
            for i = 1:unknown_num
                J(i, j) = (F(NEXT, i) - F(CUR, i))/delta;
            end
        end
        
        %% Solve the Jacobian
        dphi = solBlkDiagMat(J, F(CUR, :));
        
        for j = 1:unknown_num
            
        end
        
        
        
        
    end
    
    %% Output        

end

function ret = relaxation(a, b, alpha)
    ret = (1-alpha) * a + alpha * b;
end

function ret = df(f, x, N)
%Upwind 1st order difference
    ret = zeros(1, N);
    ret(1) = (f(2) - f(1))/(x(2)-x(1));
    for i = 2 : N-1
        if f(i) > 0
            ret(i) = (f(i) - f(i-1))/(x(i)-x(i-1));
        else
            ret(i) = (f(i+1) - f(i))/(x(i+1)-x(i));
        end
    end
    ret(N) = (f(N) - f(N-1))/(x(N)-x(N-1));
end

function ret = ddf(f, x, N)
%Central 2nd order difference
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

function x = solBlkDiagMat(B, b)
    x = B\b;
end
