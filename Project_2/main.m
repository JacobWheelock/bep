%% Initial Conditons
y0(1,1) =  -61.7987;            % V0
y0(2,1) = 0.0529322;            % m0
y0(3,1) = 0.596147;             % h0
%% Input Current
TFinal = 15;

IMag = 23;


%% ODE45
dt=[0,TFinal]; % time of integration in ms
options=odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8 1e-8],'MaxStep',0.01);
[t,y]=ode45(@(t, y) hh_diff_eq(t,y,Id,Idt, y0(1,1)), dt, y0, options);
V = y(:,1); m = y(:,2); h = y(:,3);

%% Find Conductances as a Function of Time
gNABar = 1445;
gNA = gNABar.*m.^3.*h;

%% Find Currents as a Function of Time
% Max Extracellular Concentration (mmol/L)
ECNA = 490; ECK = 20;
% Max Intracellular Concentration (mmol/L)
ICNA = 50; ICK = 400;
% Na and K Potentials
ENA = 25*log(ECNA/ICNA); EK = 25*log(ECK/ICK);
% Current
INA = gNA.*(V - ENA);
IK = gK.*(V - EK);


%% Plot Results
figure()
subplot(2,2,1)
plot(t,V);
title("Membrane Voltage vs. Time")
xlabel("Time (ms)")
ylabel("Membrane Voltage (mV)")

subplot(2,2,2)
plot(t,n);
hold on
plot(t,m);
plot(t,h);
title("Gate Activation Probabilities vs. Time")
xlabel("Time (ms)")
ylabel("Probability")
legend(["n" "m" "h"])

subplot(2,2,3)
plot(t,INA * 10e-3);
hold on
plot(t,IK * 10e-3);
title("Current vs. Time")
xlabel("Time (ms)")
ylabel("Current (mA/cm^2)")
legend(["I_{Na}" "I_K"])

subplot(2,2,4)
plot(Idt,Id);
title("I_d vs. Time")
xlabel("Time (ms)")
ylabel("Input Current (mA/cm^2)")

