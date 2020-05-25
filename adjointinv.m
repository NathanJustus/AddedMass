function Adginv = adjointinv(g)
% E Kanso, April 22, 2004


% g = (theta,x,y)

theta = g(1);
x = g(2);
y = g(3);

ct = cos(theta);
st = sin(theta);

Adginv = [           1,    0,   0;...
            -y*ct+x*st,   ct,  st;...
             y*st+x*ct,  -st,  ct];
 

