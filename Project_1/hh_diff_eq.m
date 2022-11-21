function dy = hh_diff_eq(t,y, U, Ut, Vrest)
    %% Set Up
    Id = interp1(Ut, U, t);
    % Assign Input Parameters
    V=y(1,1); n=y(2,1); m=y(3,1); h=y(4,1);
    % Max Conductances (mS/cm2)
    gNABar = 120; gKBar = 36; gLeak = 0.3;   
    % Max Extracellular Concentration (mmol/L)
    ECNA = 490; ECK = 20;
    % Max Intracellular Concentration (mmol/L)
    ICNA = 50; ICK = 400;
    % Na and K Potentials
    ENA = 25*log(ECNA/ICNA); EK = 25*log(ECK/ICK);
    % Leak Nernst Potential (mV)
    EL = -50;
    % Membrane Capacitance (uF/cm2)
    Cm = 1;
    % vm = Vm - Vrest
    vm = V - Vrest;
    %% Alphas and Betas 
    alphan = (0.01*(10-vm))/(exp((10-vm)/10) - 1);
    betan = 0.125*exp(-vm/80);
    alpham = (0.1*(25-vm))/(exp((25-vm)/10) - 1);
    betam = 4*exp(-vm/18);
    alphah = 0.07*exp(-vm/20);
    betah = 1/(exp((30-vm)/10) + 1);
    %% Differential Equations
    dv = (1/Cm)*(Id - gNABar*m^3*h*(V-ENA) - gKBar*n^4*(V-EK) - gLeak*(V-EL));
    dm = (-1*(alpham + betam)*m + alpham);
    dn = (-1*(alphan + betan)*n + alphan);
    dh = (-1*(alphah + betah)*h + alphah);

    dy = [dv; dn; dm; dh];


end