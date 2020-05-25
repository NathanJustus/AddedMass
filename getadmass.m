function [m11,m22,m33,m12,m13,m23] = getadmass(th1,th2,l,zcg1,t1,n1,del1,npts)
% E Kanso, April 22, 2004
% modified April 27, 2004


% -----------------INPUT
% 
% zg             (3,2) matrix -  coordinates of c.o.m of ellipses 
%                                starting with the ellipse in the middle
%
% beta           (3,1) matrix  - orientation of ellipses relative to
%                                inertial frame
%
% b1, b2                       - major and minor length of each ellipse
%
% e                            - offset between joint and ellipse tip
%
% npts                         - total nb of panels per ellipse
%
% -------------------

% ------------------OUTPUT
%
% Added inertia matrices
%
% -------------------

% define geometry
[zc,zcg,t,n,del] = threebody(th1,th2,l,zcg1,t1,n1,del1);

% calculate the influence matrix
[An,Bt,phi] = influence(zc,t,n,del,3*npts);

% calculate the added masses
[m11,m22,m33,m12,m13,m23] = admass(An,Bt,phi,zc,n,del,l,npts);
