function [zc,zcg,t,n,del] = threebody(th1,th2,l,zcg1,t1,n1,del1)
% E Kanso, April 22, 2004
% modified April 27, 2004


% -----------------INPUT
%
% zg     position of c.o.m of ellipses
%
% beta   angle of each ellipse relative to inertial frame
% 
% a & b   major and minor axes of the ellipse
%
% npts    total number of panels per ellipse
%
% ----------------

% -----------------OUTPUT
% 
% zc     position of collocation pts w.r.t inertial 
%        frame 
%
% t      components of vectors tangent to panels
%
% n      components of outward normal vectors
%
% del    panel length
%
% ----------------


% ----------------Example to check code
%
%  g = [0,0,0];
% 
%  th1 = 0.2*(1+cos(0)); th2 = 0;
% 
%  a = 10; b = 1; e = 2;
% 
%  npts = 50;
% 
% call ellipse
%
% [zcg1,t1,n1,del1] = ellipse(a,b,npts);
%
% -----------------End example

npts = length(zcg1);
N1 = npts;   N1p1 = N1 + 1;
N2 = 2*npts; N2p1 = N2 + 1;
N  = 3*npts;

% orientation of front and rear ellipses
cth1 = cos(th1); sth1 = sin(th1);
cth2 = cos(th2); sth2 = sin(th2);

% position of c.o.m of ellipses
% l = a + e;

zg(1,1) = 0;              
zg(1,2) = 0;

zg(2,1) = l*(1+cth1);
zg(2,2) = l*sth1;

zg(3,1) = - l*(1+cth2);
zg(3,2) = - l*sth2;

% initialize
xg  = zeros(N,1);  yg = zeros(N,1);
zcg = zeros(N,2);  zc = zeros(N,2);  
t   = zeros(N,2);  n  = zeros(N,2);
Del = zeros(N,1);

% assign

% c.o.m
xg(1:npts,1)    = ones(npts,1)*zg(1,1); 
xg(N1p1:N2,1)   = ones(npts,1)*zg(2,1); 
xg(N2p1:N,1)    = ones(npts,1)*zg(3,1);

yg(1:npts,1)    = ones(npts,1)*zg(1,2); 
yg(N1p1:N2,1)   = ones(npts,1)*zg(2,2); 
yg(N2p1:N,1)    = ones(npts,1)*zg(3,2);


% orientation of front of rear ellipses
orient2(:,1)    = ones(npts,1).*cth1;
orient2(:,2)    = ones(npts,1).*sth1;

orient3(:,1)    = ones(npts,1).*cth2;
orient3(:,2)    = ones(npts,1).*sth2;


% ellipses
zcg(1:npts,1)        = zcg1(:,1);
zcg(1:npts,2)        = zcg1(:,2);

zcg(N1p1:N2,1) = zcg1(:,1).*orient2(:,1) - zcg1(:,2).*orient2(:,2);
zcg(N1p1:N2,2) = zcg1(:,1).*orient2(:,2) + zcg1(:,2).*orient2(:,1);

zcg(N2p1:N,1)   = zcg1(:,1).*orient3(:,1) - zcg1(:,2).*orient3(:,2);
zcg(N2p1:N,2)   = zcg1(:,1).*orient3(:,2) + zcg1(:,2).*orient3(:,1);

zc = [xg,yg] + zcg;


% tangent vectors
t(1:npts,1)    = t1(:,1);
t(1:npts,2)    = t1(:,2);

t(N1p1:N2,1)   = t1(:,1).*orient2(:,1) - t1(:,2).*orient2(:,2);
t(N1p1:N2,2)   = t1(:,1).*orient2(:,2) + t1(:,2).*orient2(:,1);

t(N2p1:N,1)    = t1(:,1).*orient3(:,1) - t1(:,2).*orient3(:,2);
t(N2p1:N,2)    = t1(:,1).*orient3(:,2) + t1(:,2).*orient3(:,1);


% normal vectors
n(:,1)        = -t(:,2);
n(:,2)        =  t(:,1);


% panel length
del(1:npts,:)  = del1;
del(N1p1:N2,:) = del1;
del(N2p1:N,:)  = del1;


