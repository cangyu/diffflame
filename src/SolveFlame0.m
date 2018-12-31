function SolveFlame0(mdot_f, mdot_o, L, N, ChemTbl_DIR)
%Iterate seperately to solve the counter-flow diffusion flame without filtering
%  mdot_f: Mass flux of fuel at left side, Unit: Kg/(m^2 * s).
%  mdot_o: Mass flux of oxidizer at right side, Unit: Kg/(m^2 * s).
%  L: Length of domain, Unit: m.
%  N: Total num of grid points distributed uniformly.

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

    P = oneatm; %The constant pressure, Pa
    mdot_L = mdot_f ; %Fuel stream, 
    mdot_R = -mdot_o; %Air stream, Kg/(m^2 * s)

    zL = -L/2; %Position of left endpoint, m
    zR = L/2; %Position of right endpoint, m
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
    Y = zeros(2, K, N); %Mass fraction

    mu = zeros(1, N); %Viscosity, Pa * s = Kg / (m * s)
    cp = zeros(1, N); %Specific heat, J / (Kg * K)
    lambda = zeros(1, N); %Thermal conductivity, W / (m * K)
    D = zeros(K, N); %Binary diffusion coefficients, m^2 / s

    RS = zeros(N, 1); %Energy source due to chemical reaction, J / (m^3 * s)
    RR = zeros(K, N); %Chemical reaction rate, Kg / (m^3 * s)

    CFL = 0.8;
    max_dT = 0.5;
    max_dY = zeros(K, 1);
    
    cl = 0.0; %Upwind coef for i-1
    cm = 0.0; %Upwind coef for i
    cr = 0.0; %Upwind coef for i+1

    iter_cnt = 0;
    err = 1.0;

    %=============================Init========================================
    filefolder = fullfile('../data');
    diroutput = dir(fullfile(filefolder, 'iter*.txt'));
    n = length(diroutput);
    if n > 0
        max_idx = 0;
        for i = 1:n
            t = diroutput(i).name;
            a = isstrprop(t, 'digit');
            b = t(a);
            c = str2num(b);
            max_idx = max(max_idx, c);
        end
        iter_cnt = max_idx;

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

        for i = 1:N
            rho(PREV, i) = P / (gasconstant * T(PREV, i) * sum(squeeze(Y(PREV, :, i)) ./ MW'));
        end

        V(PREV, :) = -df(rho(PREV, :) .* u(PREV, :), dz, N) ./ (2 * rho(PREV, :));
        Nbla(PREV) = -0.1;
    end

    rho(CUR, :) = rho(PREV, :);
    u(CUR, :) = u(PREV, :);
    V(CUR, :) = V(PREV, :);
    Nbla(CUR) = Nbla(PREV);
    T(CUR, :) = T(PREV, :);
    Y(CUR, :, :) = Y(PREV, :, :);

    %==================================Loop=================================
    report(0, 'Main program running ...');
    while(err > 1e-3 && max(T(PREV, :)) > 500)
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
                ok = mrcT < 1e-4;
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
            max_dY(k) = 1e-3 * max(Y(PREV, k, :));
            if max_dY(k) < 1e-10
                max_dY(k) = 1e-10;
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
                    ok = cur_mrcr < 1e-3;
                end
                prev_mrcr = cur_mrcr;
                report(3, sprintf('(%d/%d)%s: dY_max=%e, Time step=%es, MaxAbsChange=%e, MaxRelChange=%e',k, K, NAME{1, k}, max_dY(k), dt, errY, cur_mrcr));

                %Next round
                Y(PREV, k, 2:N-1) = x(:);
            end

            %Update
            report(4, sprintf('Converges after %d iterations!', y_iter_cnt));
            Y(CUR, k, :) = relaxation(Y(PREV, k, :), Y(CUR, k, :), 0.5);
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
        fout = fopen(sprintf('../data/iter%d.txt', iter_cnt), 'w');
        for i = 1:N
            fprintf(fout, '%24.6e', rho(CUR, i));
            fprintf(fout, '%24.6e', u(CUR, i));
            fprintf(fout, '%24.6e', V(CUR, i));
            fprintf(fout, '%24.6e', P);
            fprintf(fout, '%24.6e', Nbla(CUR));
            fprintf(fout, '%24.6e', T(CUR, i));
            for k = 1:K
                fprintf(fout, '%24.6e', Y(CUR, k, i));
            end
            fprintf(fout, '\n');
        end
        fclose(fout);

        %% Swap Index
        PREV = 3 - PREV;
        CUR = 3 - CUR;
    end
    report(0, sprintf('Main program converges after %d iterations!', iter_cnt));

    %==================Transform to Z space and output============================
    report(0, 'Transforming to Z space ...');
    MixFrac = zeros(N, 1);
    for i = 1:N
       MixFrac(i) = calcZ(Y(CUR, ch4_idx, i), Y(CUR, o2_idx, i));
    end
    dMixFrac = df(MixFrac, dz, N);

    DiffCoef = zeros(N, 1);
    for i = 1:N
        DiffCoef(i) = 0.33 * D(ch4_idx, i) + 0.14 * D(o2_idx, i) + 0.52 * D(n2_idx, i); 
    end

    kai = zeros(N, 1);
    for i = 1:N
        kai(i) = 2.0 * DiffCoef(i) * dMixFrac(i)^2;
    end

    fout = fopen(sprintf('%s/mf=%.6e_mo=%.6e.txt', ChemTbl_DIR, mdot_f, mdot_o), 'w');
    fprintf(fout, '%18s\t%18s\t%18s\t%18s\t%18s\t%18s\t%18s\t%18s\t%18s\t%18s\t%18s\n', ...
        'D', 'kai', 'Z', 'dZ', 'Y_CH4', 'Y_CO', 'Y_H2', 'Y_H2O', 'Y_CO2', 'T', 'Y_NO');
    for j = 1:N
        fprintf(fout, '%18.8e\t', DiffCoef(j));
        fprintf(fout, '%18.8e\t', kai(j));
        fprintf(fout, '%18.8e\t', MixFrac(j));
        fprintf(fout, '%18.8e\t', dMixFrac(j));
        fprintf(fout, '%18.8e\t', Y(CUR, ch4_idx, j));
        fprintf(fout, '%18.8e\t', Y(CUR, co_idx, j));
        fprintf(fout, '%18.8e\t', Y(CUR, h2_idx, j));
        fprintf(fout, '%18.8e\t', Y(CUR, h2o_idx, j));
        fprintf(fout, '%18.8e\t', Y(CUR, co2_idx, j));
        fprintf(fout, '%18.8e\t', T(CUR, j));
        fprintf(fout, '%18.8e\n', Y(CUR, no_idx, j));
    end
    fclose(fout);
    
    %% Task Done
    report(0, 'Done!');
	delete('*.xml');
end

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

function ret = calcZ(Yf, Yo)
    s = 4;
    Yo_0 = 0.232;% 空气中氧化剂质量分数    
    Yf_0 = 1.0;% 燃料中燃料质量分数
    ret = (s*Yf-Yo+Yo_0)/(s*Yf_0+Yo_0);% 混合分数
    if (ret < 0)
        ret = 0.0;
    end
end

function x = solTriDiag(B, b)
    [n, ~] = size(B);
    A = spdiags(spdiags(B, -1:1), -1:1, n, n);
    x = A\b;
end

function ret = relaxation(a, b, alpha)
    ret = (1-alpha) * a + alpha * b;
end

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

function [z, u, T, y] = DiffFlameSim(domain_length, p, tin, mdot_f, mdot_o)
    runtime = cputime;  % Record the starting time

    initial_grid = domain_length*[0.0 0.2 0.4 0.6 0.8 1.0];  % Units: m
    tol_ss    = [1.0e-5 1.0e-9];        % [rtol atol] for steady-state problem
    tol_ts    = [1.0e-3 1.0e-9];        % [rtol atol] for time stepping
    loglevel  = 1;                      % Amount of diagnostic output (0 to 5)
    refine_grid = 1;                    % 1 to enable refinement, 0 to disable

    fuel = GRI30('Mix');
    ox = GRI30('Mix');
    oxcomp     =  'O2:0.21, N2:0.78';   % Air composition
    fuelcomp   =  'CH4:1';             % Fuel composition

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
    saveSoln(fl,'ch4.xml','energy',['solution with energy equation']);
    
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
