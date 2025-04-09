% Constants
global F R T;
F = 96480;  % Faradays constant (C/mol)
R = 8.314;  % Gas constant (J/K/mol)
T = 293;    % Temperature (K)

% Ion concentrations and permeabilities
PK = 4e-9;   % m/s
PNa = 0.12e-9;
PCl = 0.4e-9;
Kin = 400;   % mM
Kout = 10;
Nain = 50;
Naout = 460;
Clin = 40;
Clout = 5;

% 3. Goldman-Hodgkin-Katz equation
function Vrest = GHK_V(PK, PNa, PCl, Kin, Kout, Nain, Naout, Clin, Clout)
    global F R T;
    numerator = PK*Kout + PNa*Naout + PCl*Clin;
    denominator = PK*Kin + PNa*Nain + PCl*Clout;
    Vrest = (R*T/F) * log(numerator / denominator);
end

% Calculate resting potential
Vrest = GHK_V(PK, PNa, PCl, Kin, Kout, Nain, Naout, Clin, Clout);
fprintf('Resting potential: %.2f mV\n', Vrest * 1000);

% Interchange Cin with Cout
Vrest_reversed = GHK_V(PK, PNa, PCl, Kout, Kin, Naout, Nain, Clout, Clin);
fprintf('Resting potential (reversed): %.2f mV\n', Vrest_reversed * 1000);

% Only K permeability non-zero
Vrest_K = GHK_V(PK, 0, 0, Kin, Kout, Nain, Naout, Clin, Clout);
fprintf('Resting potential (K only): %.2f mV\n', Vrest_K * 1000);

% Only Na permeability non-zero
Vrest_Na = GHK_V(0, PNa, 0, Kin, Kout, Nain, Naout, Clin, Clout);
fprintf('Resting potential (Na only): %.2f mV\n', Vrest_Na * 1000);

% GHK current equation
function I = GHK_I(P, z, Cin, Cout, V)
    global F R T;
    if V == 0
        I = P * z * F * (Cin - Cout);
    else
        xi = z * F * V / (R * T);
        I = P * z^2 * F^2 * V / (R * T) * (Cin - Cout * exp(-xi)) / (1 - exp(-xi));
    end
end

% Calculate Im at -70mV and 0mV
V_test = [-0.07, 0];
for i = 1:length(V_test)
    IK = GHK_I(PK, 1, Kin, Kout, V_test(i));
    INa = GHK_I(PNa, 1, Nain, Naout, V_test(i));
    ICl = GHK_I(PCl, -1, Clin, Clout, V_test(i));
    Im = IK + INa + ICl;
    fprintf('Im at %.0f mV: %.2e A/m^2\n', V_test(i)*1000, Im);
end

% I-V relationship
V = -0.08:0.005:0.08;
I = zeros(size(V));
for i = 1:length(V)
    IK = GHK_I(PK, 1, Kin, Kout, V(i));
    INa = GHK_I(PNa, 1, Nain, Naout, V(i));
    ICl = GHK_I(PCl, -1, Clin, Clout, V(i));
    I(i) = IK + INa + ICl;
end

figure;
plot(V*1000, I);
xlabel('Membrane potential (mV)');
ylabel('Current density (A/m^2)');
title('I-V Relationship');

% 4. Spherical membrane compartment
diam = 100e-6;  % 100 μm diameter (soma size)
cm = 0.01;  % F/m^2
area = pi * diam^2;
Cm = cm * area;

dt = 0.1e-3;  % 0.1 ms time step
t = 0:dt:50e-3;  % 50 ms simulation

V = zeros(size(t));
V(1) = -0.05; % Initial Vm = -50 mV

% Calculate Im at Vm = -50 mV
IK_initial = -GHK_I(PK, 1, Kin, Kout, V(1));
INa_initial = -GHK_I(PNa, 1, Nain, Naout, V(1));
ICl_initial = -GHK_I(PCl, -1, Clin, Clout, V(1));
Im_initial = IK_initial + INa_initial + ICl_initial;
fprintf('Im at Vm = -50 mV: %.2e A\n', Im_initial * area);

% Calculate total membrane conductance at -50 mV
g_total = abs(Im_initial / V(1));
fprintf('Total membrane conductance at -50 mV: %.2e S/m^2\n', g_total);

for i = 2:length(t)
    IK = -GHK_I(PK, 1, Kin, Kout, V(i-1));
    INa = -GHK_I(PNa, 1, Nain, Naout, V(i-1));
    ICl = -GHK_I(PCl, -1, Clin, Clout, V(i-1));
    Im = IK + INa + ICl;
    
    I = 0; % Injected current, set to 0 for now
    V(i) = V(i-1) + dt * (Im + I) * area / Cm; % Note: using cm instead of Cm for correct units
end

figure;
plot(t*1000, V*1000);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
title('Membrane potential over time');

% Estimate τm
V_normalized = (V - V(end)) / (V(1) - V(end));
t_63 = interp1(V_normalized, t, 0.37);  % Find time when voltage reaches 63% of final value
tau_m = t_63;
fprintf('Estimated τm: %.2f ms\n', tau_m*1000);

% Plot tau with e^(t/tau)
figure;
plot(t*1000, V_normalized, 'b', t*1000, exp(-t/tau_m), 'r--');
xlabel('Time (ms)');
ylabel('Normalized Voltage');
title('Membrane Potential Decay and Exponential Fit');
legend('Simulated Data', 'Exponential Fit');

% Change permeabilities over time
V_dynamic = zeros(size(t));
V_dynamic(1) = -0.05;

for i = 2:length(t)
    if t(i) < 10e-3
        PK_current = 4e-9;
        PNa_current = 0.12e-9;
    elseif t(i) < 15e-3
        PK_current = 4e-9;
        PNa_current = 6e-9;
    elseif t(i) < 25e-3
        PK_current = 4e-9;
        PNa_current = 0.12e-9;
    elseif t(i) < 30e-3
        PK_current = 40e-9;
        PNa_current = 0.12e-9;
    else
        PK_current = 4e-9;
        PNa_current = 0.12e-9;
    end
    
    IK = -GHK_I(PK_current, 1, Kin, Kout, V_dynamic(i-1));
    INa = -GHK_I(PNa_current, 1, Nain, Naout, V_dynamic(i-1));
    ICl = -GHK_I(PCl, -1, Clin, Clout, V_dynamic(i-1));
    Im = IK + INa + ICl;
        
    I = 0; % Injected current, set to 0 for now
    V_dynamic(i) = V_dynamic(i-1) + dt * (Im + I) *area / Cm; % Note: using cm instead of Cm for correct units
end

figure;
plot(t*1000, V_dynamic*1000);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
title('Membrane potential with changing permeabilities');

% 5. Stochastic voltage-dependent ion channels
function [alpha_m, beta_m] = Na_m_rates(V)
    alpha_m = -10e5 * (V + 0.035) ./ (exp(-(V + 0.035)/0.010) - 1);
    alpha_m(abs(V + 0.035) == 0) = 1e3;  % Handle singularity
    beta_m = 4000 * exp(-(V + 0.060)/0.018);
end

function [alpha_h, beta_h] = Na_h_rates(V)
    alpha_h = 12 * exp(-V/0.020);
    beta_h = 180 ./ (exp(-(V + 0.030)/0.010) + 1);
end

V = -0.08:0.01:0.08;
[alpha_m, beta_m] = Na_m_rates(V);
m_inf = alpha_m ./ (alpha_m + beta_m);
tau_m = 1 ./ (alpha_m + beta_m);

figure;
subplot(2,1,1);
plot(V*1000, m_inf);
xlabel('Membrane potential (mV)');
ylabel('m_∞');
title('m_∞(V)');

subplot(2,1,2);
plot(V*1000, tau_m*1000);
xlabel('Membrane potential (mV)');
ylabel('τ_m (ms)');
title('τ_m(V)');

% Simulate Na-m particle dynamics
function next_state = update_channel_state(state, alpha, beta, dt)
    nch = size(state, 2);
    p01 = rand(1, nch);
    alphadt = repmat(alpha, 1, nch) * dt;
    betadt = repmat(beta, 1, nch) * dt;
    next_state1 = (p01 < alphadt) .* (state == 0);
    next_state0 = (p01 < betadt) .* (state == 1);
    next_state = state + next_state1 - next_state0;
end

t = 0:dt:50e-3;
V_test = -0.08:0.01:0.08;
m_states = zeros(length(V_test), length(t));
m_average = zeros(1, length(V_test));

for v = 1:length(V_test)
    [alpha_m, beta_m] = Na_m_rates(V_test(v));
    m_state = 0;
    for i = 1:length(t)
        m_states(v,i) = m_state;
        m_state = update_channel_state(m_state, alpha_m, beta_m, dt);
    end
    m_average(v) = mean(m_states(v,:));
end

figure;
plot(V_test*1000, m_average);
xlabel('Membrane potential (mV)');
ylabel('<n>');
title('Average degree of openness for Na-m particle');

% Na-h particle
[alpha_h, beta_h] = Na_h_rates(V);
h_inf = alpha_h ./ (alpha_h + beta_h);
tau_h = 1 ./ (alpha_h + beta_h);

figure;
subplot(2,1,1);
plot(V*1000, h_inf);
xlabel('Membrane potential (mV)');
ylabel('h_∞');
title('h_∞(V)');

subplot(2,1,2);
plot(V*1000, tau_h*1000);
xlabel('Membrane potential (mV)');
ylabel('τ_h (ms)');
title('τ_h(V)');

% Simulate m3h Na channel
m3h_states = zeros(length(V_test), length(t));
m3h_average = zeros(1, length(V_test));

for v = 1:length(V_test)
    [alpha_m, beta_m] = Na_m_rates(V_test(v));
    [alpha_h, beta_h] = Na_h_rates(V_test(v));
    m_state = [0, 0, 0];
    h_state = 0;
    for i = 1:length(t)
        m_state = update_channel_state(m_state, alpha_m, beta_m, dt);
        h_state = update_channel_state(h_state, alpha_h, beta_h, dt);
        m3h_states(v,i) = all(m_state) && h_state;
    end
    m3h_average(v) = mean(m3h_states(v,:));
end

figure;
plot(V_test*1000, m3h_average);
xlabel('Membrane potential (mV)');
ylabel('<n>');
title('Average degree of openness for m3h Na channel');