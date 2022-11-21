%% Physical Parameters
L = 9;                             % (cm)
axonRad = 7e-4;                    % micron (converted to cm)
NoRWidth = 1e-4;                   % micron (converted to cm)
INGap = 1e-1;                      % mm (10e-3 m)
x = NoRWidth+INGap:(NoRWidth + INGap):(L/2);
x = [flip(x) 0 x];

%% Initial Conditons
V0 = -61.7987; m0 = 0.0529322; h0 = 0.596147;
VN0 = 105.6;
y0 = [V0*ones(size(x)) m0*ones(size(x)) h0*ones(size(x)) VN0*ones(size(x))];

%% Input Current
TFinal = 15;
Idt = 0:0.001:TFinal;
IMag = -23;
Id = zeros(length(Idt),1);
t0 = 1;
t1 = t0 + 0.008;
t2 = t1 + 0.002;
Id(Idt >= t0) = IMag;
Id(Idt >= t1) = (-IMag * 4);
Id(Idt >= t2) = 0;
z = 1e-1;

%% ODE45
dt=[0,TFinal]; % time of integration in ms
options=odeset('RelTol',1e-4,'AbsTol',[1e-8*ones(size(y0))],'MaxStep',0.01);
[t,y]=ode45(@(t, y) hh_diff_eq(t,y,Id,Idt, y0(1,1), z, x), dt, y0, options);
%V = y(:,1); m = y(:,2); h = y(:,3); VN = y(:,4);

%% Find Conductances as a Function of Time
%gNABar = 1445;
%gNA = gNABar.*m.^3.*h;

%% Find Currents as a Function of Time
% Max Extracellular Concentration (mmol/L)
%ECNA = 490; ECK = 20;
% Max Intracellular Concentration (mmol/L)
%ICNA = 50; ICK = 400;
% Na and K Potentials
%ENA = 25*log(ECNA/ICNA); EK = 25*log(ECK/ICK);
% Current
%INA = gNA.*(V - ENA);
%IK = gK.*(V - EK);


%% Plot Results
%figure()
%plot(t,V);
%title("Membrane Voltage vs. Time")
%xlabel("Time (ms)")
%ylabel("Membrane Voltage (mV)")
x(1:45) = -x(1:45);
figure()
h = surf(x, t, y(:,3*length(x)+1:4*length(x)));

xlim([-2 2])
ylim([0 3])
zlim([100 110])
set(h,'LineStyle','none')

