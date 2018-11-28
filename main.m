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

ich4 = speciesIndex(gas, 'CH4');
io2 = speciesIndex(gas, 'O2');
in2 = speciesIndex(gas, 'N2');

NumOfSpecies = nSpecies(gas);
MW = molecularWeights(gas);

mu = 1e-5;

NumOfPnt = 1001;
zL = 0.0;
zR = 1.0;
L = zR - zL;
z = linspace(zL, zR, NumOfPnt);
dz = z(2)-z(1);

T_L = 300.0; %K
T_R = 300.0; %K
Tmax = 2000.0; %K
Tmax_pos = lin_dist(zL, zR, 0.5);
Tcoef = polyfit([zL, zR, Tmax_pos], [T_L, T_R, Tmax],2);

PREV = 1;
CUR = 2;

rho = zeros(2, NumOfPnt);
u = zeros(2, NumOfPnt);
V = zeros(2, NumOfPnt);
T = zeros(2, NumOfPnt);
Y = zeros(2, NumOfSpecies, NumOfPnt);
Nbla = zeros(1, 2);

%% Init
% rho(PREV, :) = linspace(rhoL , rhoR, NumOfPnt);
rho(PREV, :) = ones(NumOfPnt, 1) * 1.0;
rho(CUR, :) = rho(PREV, :);
u(PREV, :) = linspace(uL, uR, NumOfPnt);
Nbla(PREV) = -0.1;
Y(PREV, ich4, :) = linspace(1.0, 0.0, NumOfPnt);
Y(PREV, io2, :) = linspace(0.0, massFraction(air, 'O2'), NumOfPnt);
Y(PREV, in2, :) = linspace(0.0, massFraction(air, 'N2'), NumOfPnt);
V(PREV, :) = df(rho(PREV, :) .* u(PREV, :), dz, NumOfPnt);
for i = 1:NumOfPnt
    T(PREV, i) = polyval(Tcoef, z(i));
    V(PREV, i) = -0.5 * V(PREV, i) / rho(PREV, i);
end
T(CUR, :) = T(PREV, :);
% plot(z, T(PREV, :))
% plot(z, squeeze(Y(PREV, in2, :)))

%% Loop
err = 1.0;
iter_cnt = 0;
while(err > 1e-6)
    iter_cnt = iter_cnt + 1;
    %Solve V
    coef = zeros(NumOfPnt, NumOfPnt);
    coef(1, 1) = - rho(PREV, 1) * u(PREV, 1) / dz - mu / dz^2 + rho(PREV, 1) * V(PREV, 1);
    coef(1, 2) =  rho(PREV, 1) * u(PREV, 1) / dz + 2 * mu / dz^2;
    coef(1, 3) = -mu / dz^2;
    for i = 2 : NumOfPnt-1
        coef(i, i-1) = -0.5 * rho(PREV, i) * u(PREV, i) / dz - mu / dz^2;
        coef(i, i) = rho(PREV, i) * V(PREV, i) + 2 * mu / dz^2;
        coef(i, i+1) = 0.5 * rho(PREV, i) * u(PREV, i) / dz - mu / dz^2;
    end
    coef(NumOfPnt, NumOfPnt-2) = -mu / dz^2;
    coef(NumOfPnt, NumOfPnt-1) = -rho(PREV, NumOfPnt) * u(PREV, NumOfPnt) / dz + 2 * mu / dz^2;
    coef(NumOfPnt, NumOfPnt) = rho(PREV, NumOfPnt) * u(PREV, NumOfPnt) / dz - mu / dz^2 + rho(PREV, NumOfPnt) * V(PREV, NumOfPnt);
    V(CUR, :) = linsolve(coef, -Nbla(PREV)*ones(NumOfPnt, 1));
    
    % Solve u
    u(CUR, 1) = u(PREV, 1);
    for i = 2 : NumOfPnt - 1
        u(CUR, i) = (-2 * rho(PREV, i) * V(CUR, i) * dz + rho(PREV, i-1) * u(PREV, i-1)) / rho(PREV, i);
    end
    u(CUR, NumOfPnt) = u(PREV, NumOfPnt);
    
    % Correct V
    flux = 0.0;
    flux = flux - 2 * rho(PREV, 1) * V(CUR, 1) * (dz/2);
    for i = 2 : NumOfPnt-1
        flux = flux - 2 * rho(PREV, i) * V(CUR, i) * dz;
    end
    flux = flux - 2 * rho(PREV, NumOfPnt) * V(CUR, NumOfPnt) * (dz/2);
    gain_factor = flux / (mdot_R - mdot_L) ;
    V(CUR, :) = V(CUR, :) / gain_factor;
    
    % Correct Nbla
    lhs1 = dot(rho(PREV, :) .* u(CUR, :), df(V(CUR, :), dz, NumOfPnt));
    lhs2 = dot(rho(PREV, :) .* V(CUR, :), V(CUR, :));
    rhs2 = mu * sum(ddf(V(CUR, :), dz, NumOfPnt));
    Nbla(CUR) = (rhs2 - lhs1 - lhs2) / NumOfPnt;
    err = abs(Nbla(CUR) - Nbla(PREV));
    Nbla(CUR) = lin_dist(Nbla(PREV), Nbla(CUR), 0.5);
    
    % Solve T
    coef = zeros(NumOfPnt, NumOfPnt);
    rhs = zeros(NumOfPnt, 1);
    for i = 1:NumOfPnt
        local_T = T(PREV, i);
        setTemperature(gas, local_T);
        setMassFractions(gas, Y(PREV, :, i));
        w = netProdRates(gas);
        h = enthalpies_RT(gas) * local_T * gasconstant;                                                      
        rhs(i) = -dot(w, h) ;
    end
    
    
    % Sovle Y_k
    % TODO
    
    % Update density
    for i = 1:NumOfPnt
        tmp = sum(Y(CUR, :, i) ./ MW);
        rho(CUR, i) = P / (gasconstant * T(CUR, i) * tmp);
    end
    
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
