function [th1,th2,th1dot,th2dot] = shape_var(t)
% shape changes
% amp = 1;

% Shape variables
th1 = cos(t);
th2 = sin(t);

th1dot = -sin(t);
th2dot =  cos(t);

