%need to find these values
y0(1,1) = -60;    
y0(2,1) = 0.3; 
y0(3,1) = 0.1;
y0(4,1) = 0.7; 
dt=[0,20];  % time of integration in ms 
  
options=odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8 1e-8],'MaxStep',0.01); 
 I_d = 0;
[t,y]=ode45('hh_diff_eq',dt,y0,options); 
 
figure;
plot(t,v)