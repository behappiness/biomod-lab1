% Task 3
% Constants
R = 8.314; % J/(K*mol)
F = 96480; % C/mol
T = 293; % K

% Ion properties
P_K = 4e-9; % permeability for K
P_Na = 0.12e-9; % permeability for Na
P_Cl = 0.40e-9; % permeability for Cl

% Concentrations (in mol/L)
C_in_K = 400e-3; % intracellular K concentration
C_out_K = 10e-3; % extracellular K concentration
C_in_Na = 50e-3; % intracellular Na concentration
C_out_Na = 460e-3; % extracellular Na concentration
C_in_Cl = 40e-3; % intracellular Cl concentration
C_out_Cl = 5e-3; % extracellular Cl concentration

% GHK Voltage Equation
function Vm = ghk_voltage(P_K, P_Na, P_Cl, C_in_K, C_out_K, C_in_Na, C_out_Na, C_in_Cl, C_out_Cl, R, F, T)
    numerator = (P_K * C_out_K + P_Na * C_out_Na + P_Cl * C_in_Cl);
    denominator = (P_K * C_in_K + P_Na * C_in_Na + P_Cl * C_out_Cl);
    Vm = (R * T / F) * log(numerator / denominator);
end

% Calculate resting potential
Vm_initial = ghk_voltage(P_K, P_Na, P_Cl, C_in_K, C_out_K, C_in_Na, C_out_Na, C_in_Cl, C_out_Cl, R, F, T);
fprintf('Resting potential (initial): %.2f mV\n', Vm_initial * 1e3);

% Interchange Cin with Cout for all ions
Vm_swapped = ghk_voltage(P_K, P_Na, P_Cl, C_out_K, C_in_K, C_out_Na, C_in_Na, C_out_Cl, C_in_Cl, R, F, T);
fprintf('Resting potential (swapped): %.2f mV\n', Vm_swapped * 1e3);

% Resting potential if only K permeability is non-zero
Vm_K_only = ghk_voltage(P_K, 0, 0, C_in_K, C_out_K, 0, 0, 0, 0, R, F, T);
fprintf('Resting potential (K only): %.2f mV\n', Vm_K_only * 1e3);

% Resting potential if only Na permeability is non-zero
Vm_Na_only = ghk_voltage(0, P_Na, 0, 0, 0, C_in_Na, C_out_Na, 0, 0, R, F, T);
fprintf('Resting potential (Na only): %.2f mV\n', Vm_Na_only * 1e3);

% GHK Current Equation
function I = ghk_current(P, z, C_in, C_out, Vm, R, F, T)
    if(Vm ~= 0)
        exponent = (z * Vm * F) / (R * T);
        I = P * z * F * exponent * (C_in - C_out * exp(-exponent)) / (1 - exp(-exponent));
    else
        I = P * z * F * (C_in - C_out);
    end
    %I = P * F * ((C_out * exp(-Vm * F / (R * T))) - (C_in * exp(Vm * F / (R * T)))) / (exp(-Vm * F / (R * T)) - exp(Vm * F / (R * T)));
end

% Calculate currents at Vm = -70 mV and Vm = 0 mV
Vm_neg70 = -70e-3; % mV to V
Vm_0 = 0; % V

I_K_neg70 = ghk_current(P_K, 1, C_in_K, C_out_K, Vm_neg70, R, F, T);
I_Na_neg70 = ghk_current(P_Na, 1, C_in_Na, C_out_Na, Vm_neg70, R, F, T);
I_Cl_neg70 = ghk_current(P_Cl, -1, C_in_Cl, C_out_Cl, Vm_neg70, R, F, T);
I_total_neg70 = I_K_neg70 + I_Na_neg70 + I_Cl_neg70;

I_K_0 = ghk_current(P_K, 1, C_in_K, C_out_K, Vm_0, R, F, T);
I_Na_0 = ghk_current(P_Na, 1, C_in_Na, C_out_Na, Vm_0, R, F, T);
I_Cl_0 = ghk_current(P_Cl, -1, C_in_Cl, C_out_Cl, Vm_0, R, F, T);
I_total_0 = I_K_0 + I_Na_0 + I_Cl_0;

fprintf('Total current at Vm = -70 mV: %.2f µA\n', I_total_neg70 * 1e6);
fprintf('Total current at Vm = 0 mV: %.2f µA\n', I_total_0 * 1e6);

% I-V relationship
V_range = -80:5:80; % mV
I_values = zeros(size(V_range));

for i = 1:length(V_range)
    I_K = ghk_current(P_K, 1, C_in_K, C_out_K, V_range(i) * 1e-3, R, F, T);
    I_Na = ghk_current(P_Na, 1, C_in_Na, C_out_Na, V_range(i) * 1e-3, R, F, T);
    I_Cl = ghk_current(P_Cl, -1, C_in_Cl, C_out_Cl, V_range(i) * 1e-3, R, F, T);
    I_values(i) = I_K + I_Na + I_Cl;
end

% Plot I-V curve
figure;
plot(V_range, I_values * 1e6); % Convert to µA
xlabel('Membrane Potential (mV)');
ylabel('Total Current (µA)');
title('I-V Relationship');
grid on;

% Task 4

% Initial parameters for soma
diameter_soma = 100e-6; % Diameter of soma in meters
diameter_spine = 1e-6; % Diameter of spine head in meters
membrane_thickness = 100e-10; % Membrane thickness in meters
specific_capacitance = 0.01; % F/m^2

% Calculate surface area and capacitance
surface_area = pi * diameter_soma^2; % Surface area of the sphere
Cm = specific_capacitance * surface_area; % Total capacitance in Farads

% Initial membrane potential
Vm = -50e-3; % Initial membrane potential in Volts

% Time parameters
dt = 0.1e-3; % Time step in seconds
total_time = 50e-3; % Total simulation time in seconds
num_steps = total_time / dt; % Number of time steps

% Preallocate arrays for Vm over time
Vm_array = zeros(num_steps, 1);
Vm_array(1) = Vm;

% Simulation loop
for t = 1:num_steps
    % Calculate ionic currents using GHK current equation
    I_K = ghk_current(P_K, 1, C_in_K, C_out_K, Vm, R, F, T); % K+ current
    I_Na = ghk_current(P_Na, 1, C_in_Na, C_out_Na, Vm, R, F, T); % Na+ current
    I_total = -(I_K + I_Na); % Total current (positive inwards)

    % Update membrane potential using Euler method
    Vm = Vm + (I_total / Cm) * dt; % Vm(t+dt) = Vm(t) + Im/Cm * dt
    Vm_array(t + 1) = Vm; % Store the new Vm value

    % Change permeabilities according to the specified table
    if t * dt < 0.010
        P_K = 4e-9;
        P_Na = 0.12e-9;
    elseif t * dt < 0.015
        P_K = 4e-9;
        P_Na = 6e-9;
    elseif t * dt < 0.025
        P_K = 4e-9;
        P_Na = 0.12e-9;
    elseif t * dt < 0.030
        P_K = 40e-9;
        P_Na = 0.12e-9;
    else
        P_K = 4e-9;
        P_Na = 0.12e-9;
    end
end

% Plot the time course of Vm
time_vector = (0:num_steps) * dt; % Time vector for plotting
figure;
plot(time_vector, Vm_array * 1e3); % Convert to mV
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
title('Time Course of Membrane Potential (Soma)');
grid on;

% Change diameter to spine head size
radius_spine = diameter_spine / 2;
surface_area_spine = 4 * pi * radius_spine^2; % Surface area of spine
Cm_spine = specific_capacitance * surface_area_spine; % Total capacitance for spine

% Reset Vm for spine simulation
Vm_spine = -50e-3; % Initial membrane potential for spine
Vm_array_spine = zeros(num_steps, 1);
Vm_array_spine(1) = Vm_spine;

% Simulation loop for spine
for t = 1:num_steps
    % Calculate ionic currents using GHK current equation for spine
    I_K_spine = ghk_current(P_K, 1, C_in_K, C_out_K, Vm_spine, R, F, T); % K+ current
    I_Na_spine = ghk_current(P_Na, 1, C_in_Na, C_out_Na, Vm_spine, R, F, T); % Na+ current
    I_total_spine = -(I_K_spine + I_Na_spine); % Total current (positive inwards)

    % Update membrane potential using Euler method for spine
    Vm_spine = Vm_spine + (I_total_spine / Cm_spine) * dt; % Vm(t+dt) = Vm(t) + Im/Cm * dt
    Vm_array_spine(t + 1) = Vm_spine; % Store the new Vm value

    % Change permeabilities according to the specified table for spine
    if t * dt < 0.010
        P_K = 4e-9;
        P_Na = 0.12e-9;
    elseif t * dt < 0.015
        P_K = 4e-9;
        P_Na = 6e-9;
    elseif t * dt < 0.025
        P_K = 4e-9;
        P_Na = 0.12e-9;
    elseif t * dt < 0.030
        P_K = 40e-9;
        P_Na = 0.12e-9;
    else
        P_K = 4e-9;
        P_Na = 0.12e-9;
    end
end

% Plot the time course of Vm for spine
figure;
plot(time_vector, Vm_array_spine * 1e3); % Convert to mV
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
title('Time Course of Membrane Potential (Spine Head)');
grid on;

% Function to calculate GHK current
function I = ghk_current2(P, z, C_in, C_out, Vm, R, F, T)
    exponent = (-z * F * Vm) / (R * T);
    I = P * z * F * (C_out * exp(-exponent) - C_in) / (exp(-exponent) - 1);
end



% Constants
dt = 0.1e-3; % Time step in seconds
total_time = 0.050; % Total simulation time in seconds
num_steps = total_time / dt; % Number of time steps

% Voltage range for plotting
Vm_range = -80:10:80; % Voltage from -80 mV to 80 mV
num_voltages = length(Vm_range);

% Preallocate arrays for m_inf and tau_m
m_inf = zeros(num_voltages, 1);
tau_m = zeros(num_voltages, 1);

% Calculate m_inf and tau_m for each voltage
for i = 1:num_voltages
    Vm = Vm_range(i) * 1e-3; % Convert mV to V
    alpha_m = -10^5 * (Vm + 0.035) / (exp(-(Vm + 0.035) / 0.010) - 1);
    beta_m = 4000 * exp(-(Vm + 0.060) / 0.018);
    
    m_inf(i) = alpha_m / (alpha_m + beta_m); % Steady-state activation
    tau_m(i) = 1 / (alpha_m + beta_m); % Time constant for activation
end

% Plot m_inf and tau_m
figure;
subplot(2, 1, 1);
plot(Vm_range, m_inf);
xlabel('Membrane Potential (mV)');
ylabel('m_{\infty}');
title('Steady-State Activation (m_{\infty})');
grid on;

subplot(2, 1, 2);
plot(Vm_range, tau_m * 1e3); % Convert tau_m to ms
xlabel('Membrane Potential (mV)');
ylabel('\tau_m (ms)');
title('Time Constant for Activation (\tau_m)');
grid on;

% Simulate Na-m particle dynamics
num_channels = 1000; % Number of sodium channels
state_m = zeros(1, num_channels); % Initial state (0 = closed, 1 = open)
state_m_array = zeros(num_steps + 1, num_channels); % Store states over time
state_m_array(1, :) = state_m; % Store initial state

% Simulation loop for Na-m particle
for t = 1:num_steps
    for j = 1:num_voltages
        Vm = Vm_range(j) * 1e-3; % Convert mV to V
        alpha_m = -10^5 * (Vm + 0.035) / (exp(-(Vm + 0.035) / 0.010) - 1);
        beta_m = 4000 * exp(-(Vm + 0.060) / 0.018);
        
        % Update state of Na-m particle
        state_m = next_state(state_m, alpha_m, beta_m, dt); % Update current state
        state_m_array(t + 1, :) = state_m; % Store new state
    end
end

% Plot the state of Na-m particle over time
figure;
imagesc((0:num_steps) * dt, Vm_range, state_m_array);
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
title('State of Na-m Particle Over Time');
colorbar;
clim([0 1]); % Color scale from 0 (closed) to 1 (open)
grid on;

% Inactivation dynamics for Na-h particle
state_h = zeros(1, num_channels); % Initial state for inactivation (0 = closed, 1 = open)
state_h_array = zeros(num_steps + 1, num_channels); % Store states over time
state_h_array(1, :) = state_h; % Store initial state

% Simulation loop for Na-h particle
for t = 1:num_steps
    for j = 1:num_voltages
        Vm = Vm_range(j) * 1e-3; % Convert mV to V
        alpha_h = 12 * exp(-Vm / 0.020);
        beta_h = 180 / (exp(-(Vm + 0.030) / 0.010) + 1);
        
        % Update state of Na-h particle
        next_state_h = next_state(state_h, alpha_h, beta_h, dt);
        state_h = next_state_h; % Update current state
        state_h_array(t + 1, :) = state_h; % Store new state
    end
end

% Compose Na channel (m^3h)
state_channel = zeros(num_steps + 1, num_channels); % Store states for Na channel
for t = 1:num_steps
    state_channel(t + 1, :) = state_m_array(t + 1, :).^3 .* state_h_array(t + 1, :); % m^3h channel
end

% Plot the state of the Na channel over time
figure;
imagesc((0:num_steps) * dt, Vm_range, state_channel);
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
title('State of m^3h Sodium Channel Over Time');
colorbar;
clim([0 1]); % Color scale from 0 (closed) to 1 (open)
grid on;

% Function to determine the next state of the channel
function next_state = next_state(state, alpha, beta, dt)
    nch = size(state, 2); % Number of channels
    p01 = rand(1, nch); % Random numbers for state transition
    alphadt = repmat(alpha,1,nch) * dt; % Transition probability for open to closed
    betadt = repmat(beta,1,nch) * dt; % Transition probability for closed to open
    
    % Determine next state based on probabilities
    next_state1 = (p01 < alphadt) .* (state == 0); % Transition from closed to open
    next_state0 = (p01 < betadt) .* (state == 1); % Transition from open to closed
    next_state = state + next_state1 - next_state0; % Update state
end


% Constants
dt = 0.1e-3; % Time step in seconds
total_time = 0.050; % Total simulation time in seconds
num_steps = total_time / dt; % Number of time steps

% Voltage range for plotting
Vm_range = -80:10:80; % Voltage from -80 mV to 80 mV
num_voltages = length(Vm_range);

% Preallocate arrays for average openness <n>
average_openess = zeros(num_voltages, 1);

% Simulate for each voltage
for j = 1:num_voltages
    Vm = Vm_range(j) * 1e-3; % Convert mV to V
    state_m = zeros(1, 1000); % Initial state for Na-m (0 = closed, 1 = open)
    state_m_array = zeros(num_steps + 1, 1000); % Store states over time
    state_m_array(1, :) = state_m; % Store initial state

    % Simulation loop for Na-m particle
    for t = 1:num_steps
        alpha_m = -10^5 * (Vm + 0.035) / (exp(-(Vm + 0.035) / 0.010) - 1);
        beta_m = 4000 * exp(-(Vm + 0.060) / 0.018);
        
        % Update state of Na-m particle
        state_m = next_state(state_m, alpha_m, beta_m, dt); % Update current state
        state_m_array(t + 1, :) = state_m; % Store new state
    end

    % Calculate average degree of openness <n>
    average_openess(j) = mean(state_m_array(end, :)); % Average openness at final time step
end

% Plot <n> versus Vm
figure;
plot(Vm_range, average_openess);
xlabel('Membrane Potential (mV)');
ylabel('<n> (Average Degree of Openness)');
title('Average Degree of Openness vs Membrane Potential');
grid on;

