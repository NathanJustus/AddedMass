function gdot = vel_fun(A,g)

% g = (beta,x,y)
% A = (A_angular, Ax, Ay)

% differential equations
gdot = - [ A(1),...
        A(2)*cos(g(1)) - A(3)*sin(g(1)),...
        A(2)*sin(g(1)) + A(3)*cos(g(1))];
