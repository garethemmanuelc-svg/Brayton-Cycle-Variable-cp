clc; clear; close all;

%% ============================================================
%  REAL BRAYTON CYCLE ANALYSIS
%  Variable Cp(T) (NASA) vs Constant Gamma
% ============================================================

%% ------------------ 1. CONSTANTS ----------------------------
R  = 287;              % J/kg-K
gamma = 1.4;
cp_const = gamma*R/(gamma-1);

T1 = 300;              % K
p1 = 1e5;              % Pa

eta_c = 0.88;          % Compressor efficiency
eta_t = 0.90;          % Turbine efficiency

%% ------------------ 2. NASA POLYNOMIAL (AIR, 300–2000 K) ----
% cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
a = [ ...
     3.653
     1.337e-3
    -4.940e-7
     6.458e-11
    -2.917e-15 ];

%% ------------------ 3. RANGES -------------------------------
rp_range  = 5:0.2:50;
TIT_range = [1200 1400 1600 1800];

options = optimoptions('fsolve','Display','off');

%% ------------------- 4. MAIN LOOP (FIXED HEAT INPUT)------------
for i = 1:length(TIT_range)

    for j = 1:length(rp_range)
        rp = rp_range(j);

        %% ========= CONSTANT Cp MODEL (REFERENCE) ========
        % Compressor
        T2g = T1 * (1 + (rp^((gamma-1)/gamma) - 1)/eta_c);

        % Define FIXED heat input using reference TIT
        Qin = cp_const * (TIT_range(i) - T2g);

        % Combustor
        T3g = T2g + Qin / cp_const;

        % Turbine
        T4g = T3g * (1 - eta_t*(1 - rp^(-(gamma-1)/gamma)));

        % Work
        Wc_g = cp_const*(T2g - T1);
        Wt_g = cp_const*(T3g - T4g);
        Wnet_g = Wt_g - Wc_g;
        eta_g = Wnet_g / Qin;

        %% =========== VARIABLE Cp(T) MODEL ============
        % Compressor (isentropic)
        T2s = isentropic_solver(T1, rp, a, R, options);
        h1  = h_T(T1,a,R);
        h2s = h_T(T2s,a,R);

        % Real compressor
        h2 = h1 + (h2s - h1)/eta_c;

        % ---- Solve for T3 from FIXED Qin ----
        funT3 = @(T3) h_T(T3,a,R) - h2 - Qin;
        T3 = fsolve(funT3, TIT_range(i), options);

        % Turbine (isentropic)
        T4s = turbine_solver(T3, rp, a, R, options);
        h3  = h_T(T3,a,R);
        h4s = h_T(T4s,a,R);

        % Real turbine
        h4 = h3 - eta_t*(h3 - h4s);

        % Work
        Wc_var = h2 - h1;
        Wt_var = h3 - h4;
        Wnet_var = Wt_var - Wc_var;
        eta_var = Wnet_var / Qin;

        %% ---- Store results ----
        res(i,j).eta_var = eta_var;
        res(i,j).eta_g   = eta_g;
        res(i,j).Wnet_v  = Wnet_var;
        res(i,j).Wnet_g  = Wnet_g;

        res(i,j).Wc_var  = Wc_var;
        res(i,j).Wt_var  = Wt_var;
        res(i,j).BWR_var = Wc_var / Wt_var;

        res(i,j).Wc_g  = Wc_g;
        res(i,j).Wt_g  = Wt_g;
        res(i,j).BWR_g = Wc_g / Wt_g;

    end
end

% Preallocate
rp_opt_var = zeros(size(TIT_range));
rp_opt_g   = zeros(size(TIT_range));

Wnet_opt_var = zeros(size(TIT_range));
Wnet_opt_g   = zeros(size(TIT_range));

for i = 1:length(TIT_range)

    % ----- Variable Cp -----
    Wnet_v = [res(i,:).Wnet_v];      % extract net work array
    [Wnet_opt_var(i), idx_v] = max(Wnet_v);
    rp_opt_var(i) = rp_range(idx_v);

    % ----- Constant Gamma -----
    Wnet_g = [res(i,:).Wnet_g];
    [Wnet_opt_g(i), idx_g] = max(Wnet_g);
    rp_opt_g(i) = rp_range(idx_g);

end


%% ------------------ 5. PLOTS -------------------------------
figure; hold on; grid on;
for i = 1:length(TIT_range)
    plot(rp_range,[res(i,:).eta_var],'LineWidth',1.5)
end
xlabel('Pressure Ratio'); ylabel('Thermal Efficiency');
legend("Variable Cp | TIT="+TIT_range+" K",'Location','best');
title('Variable Cp Brayton Cycle');
ylim([0.17 0.45])
exportgraphics(gcf,'Variable Cp Brayton Cycle.png','Resolution',600)
exportgraphics(gcf,'Variable Cp Brayton Cycle.pdf','ContentType','vector')


figure; hold on; grid on;
for i = 1:length(TIT_range)
    plot(rp_range,[res(i,:).eta_g],'--','LineWidth',1)
end
xlabel('Pressure Ratio'); ylabel('Thermal Efficiency');
legend("Constant \gamma | TIT="+TIT_range+" K",'Location','best');
title('Constant \gamma Brayton Cycle');
ylim([0.2 0.5])
 exportgraphics(gcf,'Constant gamma Brayton Cycle.png','Resolution',600)
 exportgraphics(gcf,'Constant gamma Brayton Cycle.pdf','ContentType','vector')
drawnow

figure; hold on; grid on;
for i = 1:length(TIT_range)
    plot(rp_range, [res(i,:).Wnet_v]/1000, 'LineWidth', 1.5)
end
xlabel('Pressure Ratio')
ylabel('Net Work Output (kJ/kg)')
legend("Variable Cp | TIT="+TIT_range+" K",'Location','best')
title('Net Work Output vs Pressure Ratio (Variable Cp)')
ylim([50 450])
exportgraphics(gcf,'Net Work Output vs Pressure Ratio Variable Cp .png','Resolution',300)
exportgraphics(gcf,'Net Work Output vs Pressure Ratio Variable Cp .pdf','ContentType','vector')

figure; hold on; grid on;
for i = 1:length(TIT_range)
    plot(rp_range, [res(i,:).Wnet_g]/1000, '--', 'LineWidth', 1.5)
end
xlabel('Pressure Ratio')
ylabel('Net Work Output (kJ/kg)')
legend("Constant \gamma | TIT="+TIT_range+" K",'Location','best')
title('Net Work Output vs Pressure Ratio (Constant \gamma)')
ylim([50 550])
exportgraphics(gcf,'Net Work Output vs Pressure Ratio Constant Gamma .png','Resolution',300)
exportgraphics(gcf,'Net Work Output vs Pressure Ratio Constant Gamma .pdf','ContentType','vector')

figure; hold on; grid on;
for i = 1:length(TIT_range)
    plot(rp_range, [res(i,:).BWR_var], 'LineWidth', 1.5)
end
xlabel('Pressure Ratio')
ylabel('Back Work Ratio')
legend("Variable Cp | TIT="+TIT_range+" K",'Location','best')
title('Back Work Ratio vs Pressure Ratio (Variable Cp)')
exportgraphics(gcf,'Back Work Ratio vs Pressure Ratio Variable Cp.png','Resolution',600)
exportgraphics(gcf,'Back Work Ratio vs Pressure Ratio Variable Cp.pdf','ContentType','vector')

figure; hold on; grid on;
for i = 1:length(TIT_range)
    plot(rp_range, [res(i,:).BWR_g], '--', 'LineWidth', 1.5)
end
xlabel('Pressure Ratio')
ylabel('Back Work Ratio')
legend("Constant \gamma | TIT="+TIT_range+" K",'Location','best')
title('Back Work Ratio vs Pressure Ratio (Constant \gamma)')
exportgraphics(gcf,'Back Work Ratio vs Pressure Ratio Constant Gamma.png','Resolution',600)
exportgraphics(gcf,'Back Work Ratio vs Pressure Ratio Constant Gamma.pdf','ContentType','vector')

figure; hold on; grid on;

for i = 1:length(TIT_range)
    Wnet_v = [res(i,:).Wnet_v];
    Wnet_g = [res(i,:).Wnet_g];

    % Mask very small net work
    idx = Wnet_v > 0.05*max(Wnet_v);

    err = (Wnet_g(idx) - Wnet_v(idx)) ./ Wnet_v(idx) * 100;
    plot(rp_range(idx), err, 'LineWidth', 1.5)
end

xlabel('Pressure Ratio')
ylabel('Net Work Error (%)')
legend("TIT="+TIT_range+" K",'Location','best')
title('Error in Net Work Due to Constant \gamma (Valid Operating Range)')

    exportgraphics(gcf,'Error.png','Resolution',600)
    exportgraphics(gcf,'Error.pdf','ContentType','vector')
xlabel('Pressure Ratio')
ylabel('Net Work Error (%)')
legend("TIT="+TIT_range+" K",'Location','best')
title('Error in Net Work Due to Constant \gamma')
fprintf('\nOPTIMUM PRESSURE RATIO RESULTS\n');
fprintf('---------------------------------------------\n');
fprintf('TIT (K)   rp_opt (Var Cp)   rp_opt (Const γ)\n');
fprintf('---------------------------------------------\n');

for i = 1:length(TIT_range)
    fprintf('%6.0f     %8.2f          %8.2f\n', ...
        TIT_range(i), rp_opt_var(i), rp_opt_g(i));
end
fprintf('---------------------------------------------\n');



%% ============================================================
%                    FUNCTIONS
% ============================================================

function cp = cp_T(T,a,R)
    cp = R*(a(1) + a(2)*T + a(3)*T.^2 + a(4)*T.^3 + a(5)*T.^4);
end

function h = h_T(T,a,R)
    Tref = 300;   % reference temperature (K)
    h = integral(@(x) cp_T(x,a,R), Tref, T);
end

function s = s_T(T,a,R)
    s = integral(@(x) cp_T(x,a,R)./x, 300, T);
end

function T2 = isentropic_solver(T1, rp, a, R, opt)
    fun = @(T2) ...
        (s_T(T2,a,R) - s_T(T1,a,R)) - R*log(rp);
    T2 = fsolve(fun, T1*rp^0.3, opt);
end

function T4 = turbine_solver(T3, rp, a, R, opt)
    fun = @(T4) ...
        (s_T(T4,a,R) - s_T(T3,a,R)) + R*log(rp);
    T4 = fsolve(fun, T3/rp^0.3, opt);
end

function T = T_from_h(h, Tguess, a, R)
    fun = @(T) h_T(T,a,R) - h;
    T = fsolve(fun, Tguess);
end
