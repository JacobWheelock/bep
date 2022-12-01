%% Physical Parameters
L = 9;                             % (cm)
axonRad = 7e-4;                    % micron (converted to cm)
axonDi = 2*axonRad;                % micron (converted to cm)
NoRWidth = 1e-4;                   % micron (converted to cm)
INGap = 1e-1;                      % mm (converted to cm)
delX = NoRWidth + INGap; 
x = NoRWidth+INGap:(NoRWidth + INGap):(L/2);
x = [flip(x) 0 x];
% Max Conductances (mS/cm2)
gNABar = 1445; gLeak = 128;     %mS/cm^2 
% Na Potential
ENA = 115;                      % mV
% Leak Nernst Potential 
EL = -0.01;                     % mV
% Resistance 
Re = 0.350;                       % kohm-cm
Ri = 0.150;                       % kohm-cm
% Membrane Capacitance (uF/cm^2)
Cm = 2.5;

%% Initial Conditons
V0 = 0; m0 = 0.003; h0 = 0.750;
y0 = [V0*ones(size(x))', m0*ones(size(x))', h0*ones(size(x))'];

%% Input Current
TFinal = 5;
dt = 1e-3;
Idt = 0:dt:TFinal;
IMag = 130;
Id = zeros(height(Idt));
t0 = 2;
t1 = t0 + 0.8;
t2 = t1 + 0.2;
Id(Idt >= t0) = IMag;
Id(Idt >= t1) = (-IMag * 4);
Id(Idt >= t2) = 0;
z = 1e-1;
y = zeros(height(Id), length(x), 3);
y(1,:,:) = y0;

%% Euler's Method
for i = 1:length(Idt)-1
    ven = ((Re * Id(i)))./(4*pi*sqrt(x.^2 + z^2));
    v = y(i,:,1); m = y(i,:,2); h = y(i,:,3);
    vm = v;

    % CRRSS Current
    imi = (gNABar*m.^2.*h.*(vm-ENA)) + (gLeak*(vm-EL));

    % m Gates
    alpham = (97 + (0.363 * vm))./(1 + exp((31 - vm)/5.3));
    betam = alpham./exp((vm - 23.8)/4.17);

    % h Gates
    betah = 15.6./(1 + exp((24 - vm)/10));
    alphah = betah./exp((vm - 5.5)/5);

    % dm dh
    dm = (-1*(alpham + betam).*m + alpham);
    dm(1) = 0; dm(end) = 0;
    dh = (-1*(alphah + betah).*h + alphah);
    dh(1) = 0; dh(end) = 0;
    
    if (sum(isnan(v)) > 0)
        testV = y(i-1,:,1);
        testm = y(i-1,:,2);
        testh = y(i-1,:,3);
        break
    end
    
    % Convolve - try for loop?
    vConv = [0 conv(v, [1 -2 1], 'valid') 0];
    venConv = [0 conv(ven, [1 -2 1], 'valid') 0];

    dv = (1/Cm)*(-imi + ((axonDi * delX)/(4*Ri*NoRWidth))*(((vConv)/(delX^2)) + ((venConv)/(delX^2))));

    y(i+1,:,:) = [(v + dv*dt)', (m + dm*dt)', (h + dh*dt)'];
    

    
end

%% Plot Results
x(1:45) = -x(1:45);
figure()
for i = 1:5:height(y)
    text = sprintf("T = %f",Idt(i));   
    h = plot(x, y(i,:,1));
    hold on
    title(text);
    xlim([-5 5])
    ylim([-100 100])
    hold off
    drawnow
end

close all

