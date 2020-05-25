function [m11,m22,m33,m12,m13,m23] = admass(An,Bt,Phi,zcg,n,Del,l,npts)
% E Kanso, 22 April 2004
% modified 27 April 2004


% -----------------INPUT 
% 
% An(i,j)                 normal velocity induced at Ci due
%                         to a constant source distribution at panel j 
%                         Expressed w.r.t a frame attached to the panel i
%
% phi                     potential function
%
% n                       normal vectors to panels
%
% zcg                     position of the control points relative to c.o.m
%
% del                     length of panels
%
% ----------------


% -----------------INTERNAL VARIABLES
% 
% vfn                   boundary conditions, i.e., normal velocity 
%                       of the fluid at the control points 
% 
% sigma                 source distribution due to a given vfn
%
%------------------


% -----------------OUTPUT 
% 
% m11, m22, m33,        added inertias
% m12, m13, m23       
%
% -----------------

N1 = npts;   N1p1 = N1 + 1;
N2 = 2*npts; N2p1 = N2 + 1;
N3 = 3*npts;

density = 1/pi;
density = 1;

% define bc  

%  velocity due to angular rotation at unit speed
avel = zcg(:,1).*n(:,2) - zcg(:,2).*n(:,1);

xvel = ones(N3,1).*n(:,1);
yvel = ones(N3,1).*n(:,2);


% body 1
vfn_1x = zeros(N3,1); vfn_1y = zeros(N3,1); vfn_1a = zeros(N3,1);

vfn_1x(1:N1,1) = xvel(1:N1,1);  
vfn_1y(1:N1,1) = yvel(1:N1,1);  
vfn_1a(1:N1,1) = avel(1:N1,1);


% body 2
vfn_2x = zeros(N3,1); vfn_2y = zeros(N3,1); vfn_2a = zeros(N3,1);

vfn_2x(N1p1:N2,1) = xvel(N1p1:N2,1);  
vfn_2y(N1p1:N2,1) = yvel(N1p1:N2,1);  
vfn_2a(N1p1:N2,1) = avel(N1p1:N2,1);

% body 3
vfn_3x = zeros(N3,1); vfn_3y = zeros(N3,1); vfn_3a = zeros(N3,1);

vfn_3x(N2p1:N3,1) = xvel(N2p1:N3,1);  
vfn_3y(N2p1:N3,1) = yvel(N2p1:N3,1);  
vfn_3a(N2p1:N3,1) = avel(N2p1:N3,1);


% Source density distribution

inv_An = inv(An);

% body 1
sigma_1x = inv_An*vfn_1x;
sigma_1y = inv_An*vfn_1y;
sigma_1a = inv_An*vfn_1a;


% body 2
sigma_2x = inv_An*vfn_2x;
sigma_2y = inv_An*vfn_2y;
sigma_2a = inv_An*vfn_2a;


% body 3
sigma_3x = inv_An*vfn_3x;
sigma_3y = inv_An*vfn_3y;
sigma_3a = inv_An*vfn_3a;

% % check the circulation
% 
% % body 1
% vft_1x = Bt*sigma_1x; 
% vft_1y = Bt*sigma_1y;
% vft_1a = Bt*sigma_1a;
% 
% G1x = sum(vft_1x.*Del)
% G1y = sum(vft_1y.*Del)
% G1a = sum(vft_1a.*Del)
% 
% % body 2
% vft_2x = Bt*sigma_2x;
% vft_2y = Bt*sigma_2y;
% vft_2a = Bt*sigma_2a;
% 
% G2x = sum(vft_2x.*Del)
% G2y = sum(vft_2y.*Del)
% G2a = sum(vft_2a.*Del)
% 
% % body 3
% vft_3x = Bt*sigma_3x;
% vft_3y = Bt*sigma_3y;
% vft_3a = Bt*sigma_3a;
% 
% G3x = sum(vft_3x.*Del)
% G3y = sum(vft_3y.*Del)
% G3a = sum(vft_3a.*Del)



% compute potential functions 

% body 1
phi_1x = Phi*sigma_1x;
phi_1y = Phi*sigma_1y;
phi_1a = Phi*sigma_1a;

% body 2
phi_2x = Phi*sigma_2x;
phi_2y = Phi*sigma_2y;
phi_2a = Phi*sigma_2a;

% body 3
phi_3x = Phi*sigma_3x;
phi_3y = Phi*sigma_3y;
phi_3a = Phi*sigma_3a;


% compute the added masses

% body 1
m11_xx = density*sum(phi_1x.*vfn_1x.*Del);
m11_yy = density*sum(phi_1y.*vfn_1y.*Del);
m11_aa = density*sum(phi_1a.*vfn_1a.*Del);

m11_xy = density*sum(phi_1x.*vfn_1y.*Del);
m11_xa =  density*sum(phi_1x.*vfn_1a.*Del);
m11_ya =  density*sum(phi_1y.*vfn_1a.*Del);

% body 2
m22_xx = density*sum(phi_2x.*vfn_2x.*Del);
m22_yy = density*sum(phi_2y.*vfn_2y.*Del);
m22_aa = density*sum(phi_2a.*vfn_2a.*Del);

m22_xy = density*sum(phi_2x.*vfn_2y.*Del);
m22_xa =  density*sum(phi_2x.*vfn_2a.*Del);
m22_ya =  density*sum(phi_2y.*vfn_2a.*Del);

% body 3
m33_xx = density*sum(phi_3x.*vfn_3x.*Del);
m33_yy = density*sum(phi_3y.*vfn_3y.*Del);
m33_aa = density*sum(phi_3a.*vfn_3a.*Del);

m33_xy = density*sum(phi_3x.*vfn_3y.*Del);
m33_xa =  density*sum(phi_3x.*vfn_3a.*Del);
m33_ya =  density*sum(phi_3y.*vfn_3a.*Del);

% influence of body 1 on body 2
m12_xx = density*sum(phi_1x.*vfn_2x.*Del);
m12_yy = density*sum(phi_1y.*vfn_2y.*Del);
m12_aa = density*sum(phi_1a.*vfn_2a.*Del);

m12_xy = 0.5*density*(sum(phi_1x.*vfn_2y.*Del) + sum(phi_2x.*vfn_1y.*Del));
m12_xa =  0.5*density*(sum(phi_1x.*vfn_2a.*Del) + sum(phi_2x.*vfn_1a.*Del));
m12_ya =  0.5*density*(sum(phi_1y.*vfn_2a.*Del) + sum(phi_2y.*vfn_1a.*Del));

% influence of body 2 on body 3
m23_xx = density*sum(phi_2x.*vfn_3x.*Del);
m23_yy = density*sum(phi_2y.*vfn_3y.*Del); 
m23_aa = density*sum(phi_2a.*vfn_3a.*Del); 

m23_xy = 0.5*density*(sum(phi_2x.*vfn_3y.*Del) + sum(phi_3x.*vfn_2y.*Del));
m23_xa = 0.5*density*(sum(phi_2x.*vfn_3a.*Del) + sum(phi_3x.*vfn_2a.*Del));
m23_ya = 0.5*density*(sum(phi_2y.*vfn_3a.*Del) + sum(phi_3y.*vfn_2a.*Del));

% influence of body 1 on body 3
m13_xx = density*sum(phi_1x.*vfn_3x.*Del);
m13_yy = density*sum(phi_1y.*vfn_3y.*Del);
m13_aa = density*sum(phi_1a.*vfn_3a.*Del);

m13_xy = 0.5*density*(sum(phi_1x.*vfn_3y.*Del)  + sum(phi_3x.*vfn_1y.*Del));
m13_xa = 0.5*density*(sum(phi_1x.*vfn_3a.*Del) + sum(phi_3x.*vfn_1a.*Del));
m13_ya = 0.5*density*(sum(phi_1y.*vfn_3a.*Del) + sum(phi_3y.*vfn_1a.*Del));


% % non-dimensionlize
m = l^2; j = m*l^2; d = m*l;  % m = density*pi*l^2;

% body 1
m11_xx = m11_xx/m;
m11_yy = m11_yy/m;
m11_aa = m11_aa/j;

m11_xy = m11_xy/m;
m11_xa = m11_xa/d;
m11_ya = m11_ya/d;

% body 2
m22_xx = m22_xx/m;
m22_yy = m22_yy/m;
m22_aa = m22_aa/j;

m22_xy = m22_xy/m;
m22_xa = m22_xa/d;
m22_ya = m22_ya/d;

% body 3
m33_xx = m33_xx/m;
m33_yy = m33_yy/m;
m33_aa = m33_aa/j;

m33_xy = m33_xy/m;
m33_xa = m33_xa/d;
m33_ya = m33_ya/d;

% influence of body 1 on body 2
m12_xx = m12_xx/m;
m12_yy = m12_yy/m;
m12_aa = m12_aa/j;

m12_xy = m12_xy/m;
m12_xa = m12_xa/d;
m12_ya = m12_ya/d;

% influence of body 2 on body 3
m23_xx = m23_xx/m;
m23_yy = m23_yy/m;
m23_aa = m23_aa/j;

m23_xy = m23_xy/m;
m23_xa = m23_xa/d;
m23_ya = m23_ya/d;

% influence of body 1 on body 3
m13_xx = m13_xx/m;
m13_yy = m13_yy/m;
m13_aa = m13_aa/j;

m13_xy = m13_xy/m;
m13_xa = m13_xa/d;
m13_ya = m13_ya/d;



% assign
m11 = [m11_aa, m11_xa, m11_ya;...
        m11_xa, m11_xx, m11_xy;...
        m11_ya, m11_xy, m11_yy];

m22 = [m22_aa, m22_xa, m22_ya;...
        m22_xa, m22_xx, m22_xy;...
        m22_ya, m22_xy, m22_yy];

m33 = [m33_aa, m33_xa, m33_ya;...
        m33_xa, m33_xx, m33_xy;...
        m33_ya, m33_xy, m33_yy];

m12 = [m12_aa, m12_xa, m12_ya;...
        m12_xa, m12_xx, m12_xy;...
        m12_ya, m12_xy, m12_yy];

m23 = [m23_aa, m23_xa, m23_ya;...
        m23_xa, m23_xx, m23_xy;...
        m23_ya, m23_xy, m23_yy];

m13 = [m13_aa, m13_xa, m13_ya;...
        m13_xa, m13_xx, m13_xy;...
        m13_ya, m13_xy, m13_yy];


% call m-file check.m to perform checks
% check

