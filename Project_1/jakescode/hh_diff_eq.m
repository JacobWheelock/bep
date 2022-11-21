%% Hodgkin-Huxley Action Potential Model

%  MAIN PROGRAM 
function dy = hh_diff_eq(t,y)
    % Inputs: %need v,m,n,h
    % Outputs   y = vector of derivatives
    %Experimental values
    V_Na = 115;
    V_K = -12;
    V_L = 10.6;  % or should be -50 mv
    g_Na = 120;
    g_K = 36;
    g_L = 0.3;
    C_m = 1e-6;

    % Unload passed variables 
  v=y(1,1); n=y(2,1); m=y(3,1); h=y(4,1); 
    

    %conductances???

    g_Na = g_Na*m^3

    %Currents
    I_Na = g_Na*(m^3)*h*(v-V_Na);
    I_K = g_K*(n^4)*(v-V_K);
    I_L = g_L*(v-V_L);
    
    %alpha and betas
    a_m = 0.1*(25-v)/(exp((25-v)/10)-1);
    b_m =  4*exp(-v/18);
    a_n = 0.01*(10-v)./(exp((10-v)/10)-1);
    b_n = 0.125*exp(-v/80);
    a_h =  0.07*exp(-v/20);
    b_h =  1 ./ (exp((30-v)/10) + 1);
    %Temperature, for now just 1
    T = 6.3;
    k = 3^(0.1*T - 0.63);
    %Initial parameters 
    dy(1,1) = (I_d-I_Na - I_K - I_L)/C_m;     
    dy(2,1) = (-(a_m + b_m)*m + a_m)*k;  
    dy(3,1) = (-(a_n + b_n)*n + a_n)*k; 
    dy(4,1) = (-(a_h + b_h)*h + a_h)*k; 

end

 