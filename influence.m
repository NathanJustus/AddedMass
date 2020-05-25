function [An,Bt,Phi] = influence(zc,t,n,del,N)
% -----------------
% E Kanso, 14 april 2004



% -----------------INPUT
% 
% zc     position of collocation pts 
% t      components of vectors tangent to panels
% n      components of outward normal vectors 
% del    panel length
%
% zc, t and n  are w.r.t inertial frame
%
% ----------------


% -----------------INTERNAL VARIABLES
% 
% Xrel(i,j) & Yrel(i,j)    coordinates of control pt i (Ci) relative 
%                          to control pt j (Cj) w.r.t inertial frame
%
% Cn(i,j) & Ct(i,j)        normal and tangential coordinates of Ci relative
%                          to Cj w.r.t. a frame attached to the panel j
%                           
% Vn(i,j) & Vt(i,j)         normal and tangential velocities induced 
%                           at Ci due to a constant source distribution 
%                           at panel j,  w.r.t. a frame attached to panel j
%                           
% Vx(i,j) & Vy(i,j)         velocities Vn(i,j) and Vt(i,j)
%                           expressed w.r.t. inertial frame
%
% ----------------


% -----------------OUTPUT 
% 
% An(i,j) & Bt(i,j)       normal and tangential velocities induced at Ci
%                         due to a constant source distribution at panel j 
%                         Expressed w.r.t a frame attached to the panel i
%
% Phi                     potential function
%
% ----------------


% initialize
Xc = zeros(N,N); Yc = zeros(N,N);
tx = zeros(N,N); ty = zeros(N,N);
nx = zeros(N,N); ny = zeros(N,N);
Ds = zeros(N,N);


% assign
Xc = zc(:,1)*ones(1,N); Yc = zc(:,2)*ones(1,N);
tx = t(:,1)*ones(1,N);  ty = t(:,2)*ones(1,N);
nx = n(:,1)*ones(1,N);  ny = n(:,2)*ones(1,N);

Ds = 0.5.*ones(N,1)*del';


% compute Xrel(i,j) and Yrel(i,j) 
XREL = Xc-Xc';  
YREL = Yc-Yc';  

% compute Cn(i,j) and Ct(i,j)
Cn = XREL.*nx' + YREL.*ny';
Ct = XREL.*tx' + YREL.*ty';


% compute Vn(i,j) and Vt(i,j) 
temp1 = Cn.^2;
temp2 = (Ct + Ds).^2;
temp3 = (Ct - Ds).^2;
Vt_Num = temp2 + temp1;
Vt_Den = temp3 + temp1;
Vt = log(Vt_Num./Vt_Den);
Vn_Num = 2.*Cn.*Ds;
Vn_Den = Ct.^2 + temp1 - Ds.^2;
angle = 2.*atan(Vn_Num./Vn_Den); 
Vn = angle + 2*pi*eye(N,N);

% compute Vx(i,j) and Vy(i,j) 
Vx =   Vn.*nx' + Vt.*tx';
Vy =   Vn.*ny' + Vt.*ty';
        
% compute An(i,j) and Bt(i,j) 
An =  Vx.*nx + Vy.*ny;
Bt =  Vx.*tx + Vy.*ty;
        
% compute Phi
Phi = - Ct.*Vt - Cn.*Vn - Ds.*log((temp2 + temp1).*(temp3 + temp1));

% quiver(Xc(:,1),Yc(:,1),Vx(:,1),Vy(:,1),'r'); hold on; 
% N = N/3;
% plot(zc(1:N,1),zc(1:N,2),'g*',zc(1:N,1),zc(1:N,2),'g',...
%     zc(N+1:2*N,1),zc(N+1:2*N,2),'g*',zc(N+1:2*N,1),zc(N+1:2*N,2),'g',...
%     zc(2*N+1:3*N,1),zc(2*N+1:3*N,2),'g*',zc(2*N+1:3*N,1),zc(2*N+1:3*N,2),'g') 

% plot(zc(1:N,1),zc(1:N,2),'g',...
%     zc(N+1:2*N,1),zc(N+1:2*N,2),'g',...
%     zc(2*N+1:3*N,1),zc(2*N+1:3*N,2),'g') 
% 
% plot(zc(2*N+1:2*N+10,1),zc(2*N+1:2*N+10,2),'r*',...
%     zc(2*N+1:2*N+10,1),zc(2*N+1:2*N+10,2),'r') 
% plot(zc(N+1:N+10,1),zc(N+1:N+10,2),'r*',...
%     zc(N+1:N+10,1),zc(N+1:N+10,2),'r') 
% plot(zc(1:10,1),zc(1:10,2),'r*',...
%     zc(1:10,1),zc(1:10,2),'r') 

% axis([-50, 50, -10, 10])