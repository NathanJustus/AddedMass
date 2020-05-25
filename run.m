clear all; 

% E Kanso, April 22, 2004

tic

% ellipses
a = 10; b = 1;  

% offset distance between joint and tip of ellipse
e = 2;

l = a+e;

% number of panels per ellipse
N = 100;

% discretization of each ellipse
[zcg,tvec,nvec,del] = ellipse(a,b,N);

% Actual mass and inertia of each ellipse
m = a*b/l^2;
II = [(m*(a^2+b^2)/4)/l^2, 0, 0; 0, m, 0; 0, 0, m];

m = 0;
II = [0, 0, 0; 0, 0, 0; 0, 0, 0];

% solve rk4 - fixed time steps
% time span
h = 0.5;  t = 0:h:4*pi;  Nt = length(t);
hhalf = h/2; thalf = hhalf:h:t(end)-hhalf; 
hrk = h/6; 
% --------------------------------------------


% ---------- initialize ode solver -----------
g = zeros(Nt,3);
H = [0; 0; 0];

% shape variables
[th1,th2,th1dot, th2dot] = shape_var(t);
[th1half,th2half,th1dothalf, th2dothalf] = shape_var(thalf);

% connection
A0 = connection(th1(1),th2(1),th1dot(1),th2dot(1),l,zcg,tvec,nvec,del,N,II);
% ----------------------------------------------


% ------------ rk4 fixed steps ------------
tic
for i = 1:Nt-1
    % Afun	
    Ahalf = connection(th1half(i),th2half(i),th1dothalf(i),th2dothalf(i),l,zcg,tvec,nvec,del,N,II);
    A1 = connection(th1(i+1),th2(i+1),th1dot(i+1),th2dot(i+1),l,zcg,tvec,nvec,del,N,II);
   
    % intermediate steps
    K1 = vel_fun(A0,g(i,:));
	K2 = vel_fun(Ahalf,g(i,:) + hhalf*K1);
	K3 = vel_fun(Ahalf,g(i,:) + hhalf*K2);
	K4 = vel_fun(A1,g(i,:) + h*K3);
    
    % new step
    g(i+1,:) = g(i,:) + hrk*(K1+2*(K2+K3)+K4);
   
    % re-initialize
    A0 = A1;
end
toc



