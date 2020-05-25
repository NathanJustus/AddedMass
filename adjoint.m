function Adg = adjoint(g)
% E Kanso, April 6, 2004


% g = (theta,x,y)

theta = g(1);
x = g(2);
y = g(3);


ct = cos(theta);
st = sin(theta);

Adg = [  1, 0, 0;...
         y, ct, -st;...
        -x, st, ct];


