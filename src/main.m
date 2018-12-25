clear all; close all; clc;

fuel = Methane();
oxidizer = Air();
gas = GRI30('Mix');

global K;
global N;
global PREV;
global CUR;
global rho;
global u;
global V;
global P;
global Nbla;
global T;
global Y;

MW = molecularWeights(gas); %Molecular weight, Kg/Kmol
NAME = speciesNames(gas); %Name of each species
K = nSpecies(gas); %Total num of species

P = oneatm; %The constant pressure, Pa
mdot_L = 1.0/100 ; %Fuel stream, Kg/s
mdot_R = -16.6/100; %Air stream, Kg/s
rhoL = 0.716; %Density of CH4, Kg/m^3
rhoR = 1.3947; %Density of Air, Kg/m^3
S = 1.0; %Cross area, m^2
uL = mdot_L / rhoL / S; %Velocity at left entrance, m/s 
uR = mdot_R / rhoR / S; %Velocity at right entrance, m/s

N = 5001; %Total num of grid points
zL = -0.025; %Position of left endpoint, m
zR = 0.025; %Position of right endpoint, m
L = zR - zL; %Length of domain, m
z = linspace(zL, zR, N); %Coordinates for each point, m
dz = z(2)-z(1); %The uniform gap, m
dz2 = dz^2;

PREV = 1;
CUR = 2;

rho = zeros(2, N); % Kg / m^3
u = zeros(2, N); % m/s
V = zeros(2, N);
Nbla = zeros(2, 1); %The eigenvalue
T = zeros(2, N); % K
Y = zeros(2, K, N);

mu = zeros(1, N); %Viscosity, Pa * s = Kg / (m * s)
cp = zeros(1, N); %Specific heat, J / (Kg * K)
lambda = zeros(1, N); %Thermal conductivity, W / (m * K)
D = zeros(K, N); %Binary diffusion coefficients, m^2 / s

RS = zeros(N, 1); %Energy source due to chemical reaction, J / (m^3 * s)
RR = zeros(K, N); %Chemical reaction rate, Kg / (m^3 * s)

CFL = 0.8;
max_dT = 0.5;
max_dY = zeros(K, 1);

cr = 0.0; %Upwind coef for i+1
cm = 0.0; %Upwind coef for i
cl = 0.0; %Upwind coef for i-1

iter_cnt = 0;

%=============================Init========================================
if exist('../data/iter0.txt','file')
    while exist(sprintf('../data/iter%d.txt', iter_cnt),'file')
        iter_cnt = iter_cnt + 1;
    end
    iter_cnt = iter_cnt - 1;
    report(0, 'Loading existing data ...');
    fin = fopen(sprintf('../data/iter%d.txt', iter_cnt), 'r');
    data_set = fscanf(fin, '%e', [6+K N]);
    rho(PREV, :) = data_set(1, :);
    u(PREV, :) = data_set(2, :);
    V(PREV, :) = data_set(3, :);
    P = data_set(4, 1);
    Nbla(PREV) = data_set(5, 1);
    T(PREV, :) = data_set(6, :);
    for k = 1:K
        Y(PREV, k, :) = data_set(6+k, :);
    end
    fclose(fin);
else
    report(0, 'Initializing ...');
    [zc, uc, Tc,  Yc] = DiffFlameSim(L, P, 300.0, mdot_L, -mdot_R);
    tpts = z - zL;
    
    T(PREV, :) = spline(zc, Tc, tpts);
    u(PREV, :) = spline(zc, uc, tpts);
    for k = 1:K
        Y(PREV, k, :) = spline(zc, Yc(k, :), tpts);
    end
    
    rho(PREV, 1) = rhoL;
    for i = 2:N-1
        rho(PREV, i) = P / (gasconstant * T(PREV, i) * sum(squeeze(Y(PREV, :, i)) ./ MW'));
    end
    rho(PREV, N) = rhoR;
    
    V(PREV, :) = -df(rho(PREV, :) .* u(PREV, :), dz, N) ./ (2 * rho(PREV, :));
    Nbla(PREV) = -0.1;
    
    write_data(0, PREV);
end

%==================================Loop=================================
report(0, 'Main program running ...');
err = 1.0;

while(err > 1e-3)
    iter_cnt = iter_cnt + 1;
    report(1, sprintf('Iteration %d:', iter_cnt));
    
    %% Calc physical properties
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
    end
    
    %% Plot
    h = figure(1);
    set(h, 'position', get(0,'ScreenSize'));
    subplot(3, 6, 1)
    plot(z, T(PREV, :))
    title('$$T$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('K')
    
    subplot(3, 6, 2)
    plot(z, rho(PREV, :))
    title('$$\rho$$','Interpreter','latex');
    xlabel('z / m');
    ylabel('Kg\cdotm^{-3}');
    
    subplot(3, 6, 3)
    plot(z, u(PREV, :) );
    title('$$u$$','Interpreter','latex')
    xlabel('z / m');
    ylabel('m/s');
    
    subplot(3, 6, 4)
    plot(z, RS);
    title('$$-\sum{h_k\dot{\omega}_k}$$','Interpreter','latex')
    xlabel('z / m');
    ylabel('J\cdotm^{-3}\cdots^{-1}');
      
    subplot(3, 6, 7)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'CH4'), :)))
    title('Y_{CH4}');
    xlabel('z / m');
    
    subplot(3, 6, 8)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'O2'), :)))
    title('Y_{O2}');
    xlabel('z / m');
    
    subplot(3, 6, 9)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'CO2'), :)))
    title('Y_{CO2}')
    xlabel('z / m')
    
    subplot(3, 6, 10)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'H2O'), :)))
    title('Y_{H2O}')
    xlabel('z / m')
    
    subplot(3, 6, 11)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'N2'), :)))
    title('Y_{N2}')
    xlabel('z / m')
    
    subplot(3, 6, 12)
    plot(z, squeeze(Y(PREV, speciesIndex(gas, 'NO'), :)))
    title('Y_{NO}')
    xlabel('z / m')
    
    subplot(3, 6, 13)
    plot(z, squeeze(RR(speciesIndex(gas, 'CH4'), :)))
    title('$$\dot{\omega}_{CH_4}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')
    
    subplot(3, 6, 14)
    plot(z, squeeze(RR(speciesIndex(gas, 'O2'), :)))
    title('$$\dot{\omega}_{O_2}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')

    subplot(3, 6, 15)
    plot(z, squeeze(RR(speciesIndex(gas, 'CO2'), :)))
    title('$$\dot{\omega}_{CO_2}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')

    subplot(3, 6, 16)
    plot(z, squeeze(RR(speciesIndex(gas, 'H2O'), :)))
    title('$$\dot{\omega}_{H_2O}$$','Interpreter','latex');
    xlabel('z / m')
    ylabel('Kg\cdotm^{-3}\cdots^{-1}')
    
    fpic = sprintf('../pic/iter%d.png', iter_cnt);
    saveas(h, fpic);
    
    %% Solve V
    coef = zeros(N, N);
    for i = 2 : N-1
        % Upwind
        if u(PREV, i) < 0
            cr = 1.0; cm = -1.0; cl = 0.0;
        else
            cr = 0.0; cm = 1.0; cl = -1.0;
        end
        coef(i, i-1) = rho(PREV, i)*u(PREV, i)*cl/dz - mu(i)/dz2;
        coef(i, i) = rho(PREV, i)*u(PREV, i)*cm/dz + rho(PREV, i)*V(PREV, i) + 2*mu(i)/dz2;
        coef(i, i+1) = rho(PREV, i)*u(PREV, i)*cr/dz - mu(i)/dz2;
    end
    rhs = -Nbla(PREV)*ones(N, 1);
    A = coef(2:N-1, 2:N-1);
    b = rhs(2:N-1);
    b(1) = b(1) - coef(2, 1) * V(PREV, 1);
    b(N-2) = b(N-2) - coef(N-1, N) * V(PREV, N);
    x = solTriDiag(A, b);
    V(CUR, 1) = V(PREV, 1);
    V(CUR, 2:N-1) = x;
    V(CUR, N) = V(PREV, N);
    
    %% Solve u
    u(CUR, 1) = u(PREV, 1);
    for i = 2 : N - 1
        u(CUR, i) = (-2 * rho(PREV, i) * V(CUR, i) * dz + rho(PREV, i-1) * u(PREV, i-1)) / rho(PREV, i);
    end
    u(CUR, N) = u(PREV, N);
    
    %% Correct V
    flux = 0.0;
    flux = flux - 2 * rho(PREV, 1) * V(CUR, 1) * (dz/2);
    for i = 2 : N-1
        flux = flux - 2 * rho(PREV, i) * V(CUR, i) * dz;
    end
    flux = flux - 2 * rho(PREV, N) * V(CUR, N) * (dz/2);
    gain_factor = flux / (mdot_R - mdot_L) ;
    V(CUR, :) = V(CUR, :) / gain_factor;
    
    %% Correct Nbla
    dVdz = df(V(CUR, :), dz, N);
    
    lhs1 = dot(rho(PREV, :) .* u(CUR, :), dVdz);
    lhs2 = dot(rho(PREV, :) .* V(CUR, :), V(CUR, :));
    rhs2 = sum(df(mu .* dVdz, dz, N));
    
    Nbla(CUR) = (rhs2 - lhs1 - lhs2) / N;
    err = abs(Nbla(CUR) - Nbla(PREV));
    report(2, sprintf('errNbla = %f', err));
    Nbla(CUR) = relaxation(Nbla(PREV), Nbla(CUR), 0.5);
    
    %% CFL condition
    dt_cfl = CFL * dz / max(abs(u(CUR, :)));
    report(2, sprintf('Time step given by CFL condition: %es', dt_cfl));
        
    %% Solve T
    report(2, 'Solving T equation ...');
    temp_iter_cnt = 0;
    ok = false;
    while(~ok)
        temp_iter_cnt = temp_iter_cnt + 1;
        
        %Compute energy source term
        for i = 2 : N-1
            local_T = T(PREV, i);
            set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(PREV, :, i)));
            lambda(i) = thermalConductivity(gas);
            cp(i) = cp_mass(gas);
            D(:, i) = mixDiffCoeffs(gas);
            w = netProdRates(gas); % kmol / (m^3 * s)
            h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol
            RS(i) = -dot(w, h); % J / (m^3 * s)
        end
        %Filtering
        %TODO
        
        dYdz = zeros(K, N);
        for k = 1:K
            dYdz(k, :) = df(Y(PREV, k, :), dz, N);
        end
        SpecDiffTerm = zeros(N, 1);
        for i = 2:N-1
            SpecDiffTerm(i) = rho(PREV, i) * cp(i) * dot(D(:, i), dYdz(:, i));
        end
        
        %Choose proper time step
        dt = dt_cfl;
        for i = 2:N-1
            dt_chem = rho(PREV, i)*cp(i)*max_dT/(abs(RS(i))+1e-20);
            dt = min(dt, dt_chem);
        end
        
        %Construct the coefficient matrix
        coef = zeros(N, N);
        for i = 2 : N-1
            % Upwind
            if u(PREV, i) < 0
                cr = 1.0; cm = -1.0; cl = 0.0;
            else
                cr = 0.0; cm = 1.0; cl = -1.0;
            end
            coef(i, i-1) = rho(PREV, i)*cp(i)*u(CUR, i)*cl*dt/dz - lambda(i)*dt/dz2 - SpecDiffTerm(i)*cl*dt/dz;
            coef(i, i) = rho(PREV, i)*cp(i)*(1+u(CUR, i)*cm*dt/dz) + 2*lambda(i)*dt/dz2 - SpecDiffTerm(i)*cm*dt/dz;
            coef(i, i+1) = rho(PREV, i)*cp(i)*u(CUR, i)*cr*dt/dz - lambda(i)*dt/dz2 - SpecDiffTerm(i)*cr*dt/dz;
        end
        A = coef(2:N-1, 2:N-1);
        
        %Construct the RHS
        rhs = zeros(N, 1);
        for i = 2 : N-1
            rhs(i) = rho(PREV, i)*cp(i)*T(PREV, i)+dt*RS(i);
        end
        b = rhs(2:N-1);
        b(1) = b(1) - coef(2, 1) * T(PREV, 1);
        b(N-2) = b(N-2) - coef(N-1, N) * T(PREV, N);
        
        %Solve
        x = solTriDiag(A, b);
        
        %Calc error
        change_of_T = abs(x - squeeze(T(PREV, 2:N-1))');
        macT = max(change_of_T);
        relative_change_of_T = zeros(N-2, 1);
        for i = 2:N-1
            idx = i-1;
            relative_change_of_T(idx) = change_of_T(idx) / T(PREV, i);
        end
        mrcT = max(relative_change_of_T);
        
        %Check convergence
        if temp_iter_cnt > 50
            ok = true;
        else
            ok = mrcT < 2e-4;
        end
        report(3, sprintf('Time step=%es, MaxAbsChange=%eK, MaxRelChange=%e', dt, macT, mrcT));
        
        %Check constraint: no less than 300, no greater than 3000
        for i = 2:N-1
            idx = i - 1;
            if x(idx) < 300.0
                x(idx) = 300.0;
            end
            if x(idx) > 3000.0
                x(idx) = 3000.0;
            end
        end
        
        %Next round
        T(PREV, 2:N-1) = x(:);
    end
    
    %Update
    report(3, sprintf('Converges after %d iterations!', temp_iter_cnt));
    T(CUR, :) = T(PREV, :);
    
    %% Sovle Y
    report(2, 'Solving Y equations ...');
    
    %Update diffusion coefficients and RR
    for i = 1:N
        local_T = T(CUR, i);
        set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(PREV, :, i)));
        D(:, i) = mixDiffCoeffs(gas);    
        w = netProdRates(gas); % kmol / (m^3 * s)
        RR(:, i) = w .* MW; % Kg / (m^3 * s)
    end
    
    %Solve each species
    for k=1:K
        max_dY(k) = 1e-2 * max(Y(PREV, k, :));
        if max_dY(k) < 1e-10
            max_dY(k) = 1e-6;
        end
        
        ok = false;
        y_iter_cnt = 0;
        prev_mrcr = 1e100;
        while(~ok)
            y_iter_cnt = y_iter_cnt + 1;
            
            %Choose proper time step
            dt = dt_cfl;
            for i = 2:N-1
                dt_chem = rho(PREV, i)*max_dY(k)/(abs(RR(k, i))+1e-80);
                dt = min(dt, dt_chem);
            end
                        
            %Construct the coefficient matrix
            coef = zeros(N, N);
            for i = 2 : N-1
                % Upwind
                if u(PREV, i) < 0
                    cr = 1.0; cm = -1.0; cl = 0.0;
                else
                    cr = 0.0; cm = 1.0; cl = -1.0;
                end
                coef(i, i-1) = rho(PREV, i)*(u(CUR, i)*cl*dt/dz-D(k, i)*dt/dz2);
                coef(i, i) = rho(PREV, i)*(1+2*D(k, i)*dt/dz2+u(CUR, i)*cm*dt/dz);
                coef(i, i+1) = rho(PREV, i)*(u(CUR, i)*cr*dt/dz-D(k, i)*dt/dz2);
            end
            A = coef(2:N-1, 2:N-1);
            
            %Construct the RHS
            rhs = zeros(N, 1);
            for i = 2 : N-1
                rhs(i) = rho(PREV, i)*Y(PREV, k, i)+dt*RR(k, i);
            end
            b = rhs(2:N-1);
            b(1) = b(1) - coef(2, 1) * Y(PREV, k, 1);
            b(N-2) = b(N-2) - coef(N-1, N) * Y(PREV, k, N);
            
            %Solve
            x = solTriDiag(A, b);
            
            %Check constraints: no less than 0, no greater than 1.0
            for i = 2:N-1
                idx = i-1;
                if x(idx) < 0
                    x(idx) = 0.0;
                end
                if x(idx) > 1.0
                    x(idx) = 1.0;
                end
            end
            
            %Stat err
            errY = max(abs(squeeze(Y(PREV, k, 2:N-1)) - x));
            cur_mrcr = errY/ (max(Y(PREV, k, 2:N-1)) + 1e-60);
            
            %Check convergence
            if y_iter_cnt > 100
                ok = true;
            elseif cur_mrcr >= prev_mrcr
                ok = true;
            else
                ok = cur_mrcr < 5e-3;
            end
            prev_mrcr = cur_mrcr;
            report(3, sprintf('(%d/%d)%s: dY_max=%e, Time step=%es, MaxAbsChange=%e, MaxRelChange=%e',k, K, NAME{1, k}, max_dY(k), dt, errY, cur_mrcr));

            %Next round
            Y(PREV, k, 2:N-1) = x(:);
        end
        
        %Update
        report(4, sprintf('Converges after %d iterations!', y_iter_cnt));
        Y(CUR, k, :) = Y(PREV, k, :);
    end
    
    %Normalization
    for i = 2 : N-1
        Y(CUR, :, i) = Y(CUR, :, i) / sum(Y(CUR, :, i));
    end
    
    %% Update density
    rho(CUR, 1) = rho(PREV, 1);
    for i = 2:N-1
        rho(CUR, i) = P / (gasconstant * T(CUR, i) * sum(squeeze(Y(CUR, :, i)) ./ MW'));
    end
    rho(CUR, N) = rho(PREV, N);
    
    %% Save current iteration
    report(2, 'Writing data ...');
    write_data(iter_cnt, CUR);
    
    %% Swap Index
    PREV = 3 - PREV;
    CUR = 3 - CUR;
end
report(0, sprintf('Main program converges after %d iterations!', iter_cnt));

%=================================Helpers================================
function report(level, msg)
    n = int32(2 * level);
    k = 0;
    while(k < n)
        fprintf(' ');
        k = k+1;
    end
    fprintf('%s\n', msg);
end

function write_data(iter, idx)
    global N
    global K
    global rho
    global u
    global V
    global P
    global Nbla
    global T
    global Y
    
    fout = fopen(sprintf('../data/iter%d.txt', iter), 'w');
    for i = 1:N
        fprintf(fout, '%24.6e', rho(idx, i));
        fprintf(fout, '%24.6e', u(idx, i));
        fprintf(fout, '%24.6e', V(idx, i));
        fprintf(fout, '%24.6e', P);
        fprintf(fout, '%24.6e', Nbla(idx));
        fprintf(fout, '%24.6e', T(idx, i));
        for k = 1:K
            fprintf(fout, '%24.6e', Y(idx, k, i));
        end
        fprintf(fout, '\n');
    end
    fclose(fout);
end
