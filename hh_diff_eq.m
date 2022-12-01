function dy = hh_diff_eq(t,y, U, Ut, vrest, z, x)
    %% Set Up
    % Interpolate Input
    Id = interp1(Ut, U, t);

    % Assign Input Parameters
    vN = y(1:length(x)); 
    m= y(length(x)+1:2*length(x));
    h= y(2*length(x) + 1:3*length(x));
    
    % Max Conductances (mS/cm2)
    gNABar = 1445; gLeak = 128;     %mS/cm^2
   
    % Na Potential
    ENA = 115;                      % mV

    % Leak Nernst Potential 
    EL = -0.01;                     % mV

    % Resistance 
    Re = 350;                       % ohm-cm
    Ri = 150;                       % ohm-cm

    % Membrane Capacitance (uF/cm^2)
    Cm = 2.5;

    % Physical Lengths
    axonRad = 7e-4;                    % micron (converted to cm)
    axonDi = 2*axonRad;
    NoRWidth = 1e-4;                   % micron (converted to cm)
    INGap = 1e-1;                      % mm (converted to cm)
    delX = NoRWidth + INGap;
    
    ven = ((Re * Id)*ones(size(x)))./(4*pi*sqrt(x.^2 + z^2));

    %vN = vin - ven' - vrest;
    %vN = vin;
    vm = vN - vrest;

    %% Alphas and Betas 
    % m Gates
    alpham = (97 + 0.363 * vm)./(1 + exp((31 - vm)/5.3));
    betam = alpham./exp((vm - 23.8)/4.17);

    % h Gates
    betah = 15.6./(1 + exp((24 - vm)/10));
    alphah = betah./exp((vm - 5.5)/5);

    %% Differential Equations
    dm = (-1*(alpham + betam).*m + alpham);
    dh = (-1*(alphah + betah).*h + alphah);
    imi = gNABar*m.^2.*h.*(vm-ENA) + gLeak*(vm-EL);

    vNConv = conv(vN, [1 -2 1], 'same');
    venConv = conv(ven, [1 -2 1], 'same')';
    disp(vNConv')

    dvN = (1/Cm)*(-imi + ((axonDi * delX)/(4*Ri*NoRWidth))*(((vNConv)/(delX^2)) + ((venConv)/(delX^2))));    

    dy = [dvN; dm; dh];


end