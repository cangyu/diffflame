function SolveFlame3(mdot_f, mdot_o, L, N, ChemTbl_DIR, MAX_ITER)
%Iterate together to solve the counter-flow diffusion flame without filtering.
%  mdot_f: Mass flux of fuel at left side, Unit: Kg/(m^2 * s).
%  mdot_o: Mass flux of oxidizer at right side, Unit: Kg/(m^2 * s).
%  L: Length of domain, Unit: m.
%  N: Total num of grid points distributed uniformly.
%  ChemTbl_DIR: Target directory where transformed data files are stored.
%  MAX_ITER: Maximun number of global iteration.

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
    mdot_L = mdot_f ; %Fuel stream, Kg/(m^2 * s)
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

    CFL = 0.9;
    max_dT = 1.0;
    max_dY = zeros(K, N);
    
    cl = 0.0; %Upwind coef for i-1
    cm = 0.0; %Upwind coef for i
    cr = 0.0; %Upwind coef for i+1

    global_iter_cnt = 0;
    global_iter_ok = false;
    
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
        global_iter_cnt = max_idx;

        report(0, 'Loading existing data ...');
        fin = fopen(sprintf('../data/iter%d.txt', global_iter_cnt), 'r');
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
        
        %Enforce the velocity inlet B.C. according to mass flux.
        u(PREV, 1) = mdot_L / rho(PREV, 1);
        u(PREV, N) = mdot_R / rho(PREV, N);
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

        %V set to 0 at boundary as no vertical slip physically
        V(PREV, :) = -df(rho(PREV, :) .* u(PREV, :), dz, N) ./ (2 * rho(PREV, :));
        V(PREV,1)=0.0;
        V(PREV,N)=0.0;
        
        %Select initial guess of the eigenvalue
        dVdz = df(V(PREV, :), dz, N);
        lhs1 = dot(rho(PREV, :) .* u(PREV, :), dVdz);
        lhs2 = dot(rho(PREV, :) .* V(PREV, :), V(PREV, :));
        rhs2 = sum(df(mu .* dVdz, dz, N));
        Nbla(PREV) = (rhs2 - lhs1 - lhs2) / N;
    end

    rho(CUR, :) = rho(PREV, :);
    u(CUR, :) = u(PREV, :);
    V(CUR, :) = V(PREV, :);
    Nbla(CUR) = Nbla(PREV);
    T(CUR, :) = T(PREV, :);
    Y(CUR, :, :) = Y(PREV, :, :);

    %==================================Loop=================================
    report(0, 'Main program running ...');
    while(~global_iter_ok && max(T(PREV, :)) > 350 && global_iter_cnt < MAX_ITER)
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

        subplot(3, 6, 5)
        plot(z, squeeze(Y(PREV, no_idx, :)))
        title('Y_{NO}');
        xlabel('z / m');

        subplot(3, 6, 6)
        plot(z, squeeze(Y(PREV, co_idx, :)))
        title('Y_{CO}');
        xlabel('z / m');

        subplot(3, 6, 7)
        plot(z, squeeze(Y(PREV, ch4_idx, :)))
        title('Y_{CH4}');
        xlabel('z / m');

        subplot(3, 6, 8)
        plot(z, squeeze(Y(PREV, o2_idx, :)))
        title('Y_{O2}');
        xlabel('z / m');

        subplot(3, 6, 9)
        plot(z, squeeze(Y(PREV, co2_idx, :)))
        title('Y_{CO2}')
        xlabel('z / m')

        subplot(3, 6, 10)
        plot(z, squeeze(Y(PREV, h2o_idx, :)))
        title('Y_{H2O}')
        xlabel('z / m')

        subplot(3, 6, 11)
        plot(z, squeeze(Y(PREV, h2_idx, :)))
        title('Y_{H2}')
        xlabel('z / m')

        subplot(3, 6, 12)
        plot(z, squeeze(Y(PREV, n2_idx, :)))
        title('Y_{N2}')
        xlabel('z / m')

        subplot(3, 6, 13)
        plot(z, squeeze(RR(ch4_idx, :)))
        title('$$\dot{\omega}_{CH_4}$$','Interpreter','latex');
        xlabel('z / m')
        ylabel('Kg\cdotm^{-3}\cdots^{-1}')

        subplot(3, 6, 14)
        plot(z, squeeze(RR(o2_idx, :)))
        title('$$\dot{\omega}_{O_2}$$','Interpreter','latex');
        xlabel('z / m')
        ylabel('Kg\cdotm^{-3}\cdots^{-1}')

        subplot(3, 6, 15)
        plot(z, squeeze(RR(co2_idx, :)))
        title('$$\dot{\omega}_{CO_2}$$','Interpreter','latex');
        xlabel('z / m')
        ylabel('Kg\cdotm^{-3}\cdots^{-1}')

        subplot(3, 6, 16)
        plot(z, squeeze(RR(h2o_idx, :)))
        title('$$\dot{\omega}_{H_2O}$$','Interpreter','latex');
        xlabel('z / m')
        ylabel('Kg\cdotm^{-3}\cdots^{-1}')

        subplot(3, 6, 17)
        plot(z, squeeze(RR(h2_idx, :)))
        title('$$\dot{\omega}_{H_2}$$','Interpreter','latex');
        xlabel('z / m')
        ylabel('Kg\cdotm^{-3}\cdots^{-1}')

        %saveas(h, sprintf('../pic/iter%d.png', global_iter_cnt));

        %% Update global iteration counter
        global_iter_cnt = global_iter_cnt + 1;
        report(1, sprintf('Iteration %d:', global_iter_cnt));

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
        abs_change_of_Nbla = abs(Nbla(CUR) - Nbla(PREV));
        rel_change_of_Nbla = abs(abs_change_of_Nbla / Nbla(PREV));
        report(2, sprintf('Nbla = %f, abs_err = %f, rel_err = %e', Nbla(CUR), abs_change_of_Nbla, rel_change_of_Nbla));
        global_iter_ok = abs_change_of_Nbla < 1e-2 || rel_change_of_Nbla < 5e-8;

        %% CFL condition
        dt_cfl = CFL * dz / max(abs(u(CUR, :)));

        %% Solve T and Y together
        report(2, 'Solving T and Y equations ...');
        TY_iter_cnt = 0;
        TY_iter_ok = false;
        while(~TY_iter_ok)
            TY_iter_cnt = TY_iter_cnt + 1;
            
            %Update energy source term, reaction rate
            %and physical properties
            for i = 1 : N
                local_T = T(PREV, i);
                set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(PREV, :, i)));
                lambda(i) = thermalConductivity(gas);
                cp(i) = cp_mass(gas);
                D(:, i) = mixDiffCoeffs(gas);
                w = netProdRates(gas); % kmol / (m^3 * s)
                h = enthalpies_RT(gas) * local_T * gasconstant; % J/Kmol
                RS(i) = -dot(w, h); % J / (m^3 * s)
                RR(:, i) = w .* MW; % Kg / (m^3 * s)
            end
            
            %% Choose proper time-step
            for k = 1:K
                for i = 2:N-1
                    max_dY(k, i) = 1e-3 * Y(PREV, k, i);
                     if max_dY(k, i) < 1e-10
                        max_dY(k, i) = 1e-10;
                    end
                end
            end
            dt = dt_cfl;
            for i = 2:N-1
                % Time-step locally determined by Temperature change limit
                dt_Tchem = rho(PREV, i)*cp(i)*max_dT/(abs(RS(i))+1e-20);
                dt = min(dt, dt_Tchem);
                % Time-step locally determined by each Specise change limit
                for k = 1:K
                    dt_Ychem = rho(PREV, i)*max_dY(k, i)/(abs(RR(k, i))+1e-80);
                    dt = min(dt, dt_Ychem);
                end
            end 
            
            %% Solve T
            %Calculate the specise diffusion term
            dYdz = zeros(K, N);
            for k = 1:K
                dYdz(k, :) = df(Y(PREV, k, :), dz, N);
            end
            SpecDiffTerm = zeros(N, 1);
            for i = 2:N-1
                SpecDiffTerm(i) = rho(PREV, i) * cp(i) * dot(D(:, i), dYdz(:, i));
            end
            %Construct Temperature Equation coef matrix
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
            %Construct corresponding rhs
            rhs = zeros(N, 1);
            for i = 2 : N-1
                rhs(i) = rho(PREV, i)*cp(i)*T(PREV, i)+dt*RS(i);
            end
            b = rhs(2:N-1);
            b(1) = b(1) - coef(2, 1) * T(PREV, 1);
            b(N-2) = b(N-2) - coef(N-1, N) * T(PREV, N);
            %Solve the tri-diagonal sparse matrix
            xT = solTriDiag(A, b);
            %Check constraints
            for i = 2:N-1
                idx = i - 1;
                %No less than 300 K
                if xT(idx) < 300.0
                    xT(idx) = 300.0;
                end
                %No greater than 3000 K
                if xT(idx) > 3000.0
                    xT(idx) = 3000.0;
                end
            end
            %Calculate absolute and relative error
            abs_change_of_T = abs(xT - squeeze(T(PREV, 2:N-1))');
            max_abs_change_of_T = max(abs_change_of_T);
            rel_change_of_T = abs_change_of_T ./ T(PREV, 2:N-1)';
            max_rel_change_of_T = max(rel_change_of_T);
            
            %% Solve  each Y
            xY = zeros(K, N-2);
            max_abs_change_of_Y = zeros(K, 1);
            max_rel_change_of_Y = zeros(K, 1);
            for k=1:K
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
                xY(k, :) = solTriDiag(A, b);
                %Check constraints
                for i = 2:N-1
                    idx = i-1;
                    %No less than 0
                    if xY(k, idx) < 0
                        xY(k, idx) = 0.0;
                    end
                    %No greater than 1.0
                    if xY(k, idx) > 1.0
                        xY(k, idx) = 1.0;
                    end
                end
                %Calculate absolute and relative error 
                max_abs_change_of_Y(k) = max(abs(squeeze(xY(k, :)') - squeeze(Y(PREV, k, 2:N-1))));
                max_rel_change_of_Y(k) = max_abs_change_of_Y(k) / (max(Y(PREV, k, 2:N-1)) + 1e-80);
                global_max_rel_change_of_Y = max(max_rel_change_of_Y);
            end
            
            %% Check convergence
            TY_iter_ok = max_rel_change_of_T < 1e-4 && global_max_rel_change_of_Y < 1e-3;
            report(3, sprintf('dt=%es, Tmax=%fK, dT_max=%eK, dT/T_max=%e, dY/Y_max=%e', dt, max(xT), max_abs_change_of_T, max_rel_change_of_T, global_max_rel_change_of_Y));
            
            %% Update T and Y
            T(PREV, 2:N-1) = xT(:);
            for k = 1:K
                Y(PREV, k, 2:N-1) = xY(k, :);
            end
            
            %% Normalization
            for i = 2 : N-1
                Y(PREV, :, i) = Y(PREV, :, i) / sum(Y(PREV, :, i));
            end
        end
              
        %Update
        report(3, sprintf('Converges after %d iterations!', TY_iter_cnt));
        T(CUR, :) = T(PREV, :);
        Y(CUR, k, :) = relaxation(Y(PREV, k, :), Y(CUR, k, :), 0.5);

        %% Update density
        rho(CUR, 1) = rho(PREV, 1);
        for i = 2:N-1
            rho(CUR, i) = P / (gasconstant * T(CUR, i) * sum(squeeze(Y(CUR, :, i)) ./ MW'));
        end
        rho(CUR, N) = rho(PREV, N);

        %% Save current iteration
        report(2, 'Writing data ...');
        fout = fopen(sprintf('../data/iter%d.txt', global_iter_cnt), 'w');
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
    report(0, sprintf('Main program converges after %d iterations!', global_iter_cnt));

    %% ==================Transform to Z space and output============================
    report(0, 'Transforming to Z space ...');
	for i = 1:N
		local_T = T(CUR, i);
		set(gas, 'T', local_T, 'P', P, 'Y', squeeze(Y(CUR, :, i)));
		D(:, i) = mixDiffCoeffs(gas);
	end
	
    MixFrac = zeros(N, 1);
    for i = 1:N
        local_yf = Y(CUR, ch4_idx, i) + Y(CUR, h2_idx, i);
        local_yo = Y(CUR, o2_idx, i);
        MixFrac(i) = calcZ(local_yf, local_yo);
    end
    dMixFrac = df(MixFrac, dz, N);

    DiffCoef = zeros(N, 1);
    for i = 1:N
        DiffCoef(i) = dot(squeeze(Y(CUR, :, i)), squeeze(D(:, i)));
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
    mail_title = sprintf('mf=%f_mo=%f Done!', mdot_f, mdot_o);
    mail_content = sprintf('%d iterations, Tmax=%fK', global_iter_cnt, max(T(CUR, :)));
    mail_to_me(mail_title, mail_content);
end

%% =================================Helpers================================
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
    s = 40/9; %Mass ratio of oxidizer and fuel at stoichiometry
    Yo_0 = 0.232; %Mixture fraction of oxidizer in air stream
    Yf_0 = 1.0; %Mixture fraction of fuel in gas stream
    ret = (s*Yf-Yo+Yo_0)/(s*Yf_0+Yo_0); %The mixture fraction
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

    initial_grid = domain_length*linspace(0,1, 501);  % Units: m
    tol_ss    = [1.0e-5 1.0e-9];        % [rtol atol] for steady-state problem
    tol_ts    = [1.0e-3 1.0e-9];        % [rtol atol] for time stepping
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

function mail_to_me(subject,content) 
    dizhi = 'cysqyh@163.com'; 
    mima = '&LgbNXr&C*m$5WPh';
    setpref('Internet','E_mail',dizhi); 
    setpref('Internet','SMTP_Server','smtp.163.com'); 
    setpref('Internet','SMTP_Username',dizhi); 
    setpref('Internet','SMTP_Password',mima); 
    props = java.lang.System.getProperties; 
    props.setProperty('mail.smtp.auth','true'); 
    sendmail('yu.cang@sjtu.edu.cn',subject,content);
end