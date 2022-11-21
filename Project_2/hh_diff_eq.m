function dy = hh_diff_eq(t,y, U, Ut, vrest, z)
    %% Set Up
    % Interpolate Input
    Id = interp1(Ut, U, t);

    % Assign Input Parameters
    V=y(1,:); m=y(2,:); h=y(3,:);

    % Max Conductances (mS/cm2)
    gNABar = 1445; gLeak = 128;     %mS/cm^2
   
    % Na Potentials
    ENA = 115;                      % mV

    % Leak Nernst Potential (mV)
    EL = -0.01;                     % mV

    % Resistance 
    Re = 350;                       % ohm-cm
    Ri = 150;                       % ohm-cm

    % Membrane Capacitance (uF/cm^2)
    Cm = 2.5;

    % Physical Lengths
    axonRad = 7;                    % micron (10e-6 m)
    NoRWidth = 1;                   % micron (10e-6 m)
    INGap = 1;                      % mm (10e-3 m)

    ven = (Re * Id)/(4*pi*sqrt())

    % vm = Vm - Vrest
    vN = vin - ven - vrest;

    %% Alphas and Betas 
    % m Gates
    alpham = (97 + 0.363 * vm)/(1 + exp((31 - vm)/5.3));
    betam = alpham/exp((vm - 23.8)/4.17);

    % h Gates
    betah = 15.6/(1 + exp((24 - vm)/10));
    alphah = betah/exp((vm - 5.5)/5);

    %% Differential Equations
    dv = (1/Cm)*(Id - gNABar*m^2*h*(V-ENA) - gLeak*(V-EL));
    dm = (-1*(alpham + betam)*m + alpham);
    dh = (-1*(alphah + betah)*h + alphah);

    dy = [dv; dm; dh];


end