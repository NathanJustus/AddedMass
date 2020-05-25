function A = connection(th1,th2,th1dot,th2dot,l,zcg1,t1,n1,del1,npts,II);

% non-dimensional added inertias
[I33a,I11a,I22a,I31,I32,I12] = getadmass(th1,th2,l,zcg1,t1,n1,del1,npts);

    
% motions of 1 and 2 relative to 3
  x1 = [th1,(1+cos(th1)),sin(th1)];
  x2 = [th2,-(1+cos(th2)),-sin(th2)];

% adjoint action
  Adx1inv = adjointinv(x1);
  Adx2inv = adjointinv(x2);

% actual mass
I33 = I33a; %+ II;
I22 = I22a; %+ Adx2inv'*II*Adx2inv;
I11 = I11a; %+ Adx1inv'*II*Adx1inv;

  
% locked moment of inertia
  IIloc = I11 + I22 + I33 + 2*I12 + 2*I31 + 2*I32;
         
% velocities of ellipses 1 and 2 relative to 3 w.r.t B3-fixed frame
  zeta1_temp = [th1dot; 0;  -th1dot]; 
  zeta2_temp = [th2dot; 0; th2dot];

  h1p2  = (I11+I12+I31)*zeta1_temp + (I22+I32+I12)*zeta2_temp;

% connection
  A = inv(IIloc)*h1p2;
