% E Kanso, Nov 5, 2006
% 
% Run the file called run.m to integrate equation (6.7) on page 272 of 
% my JNLS paper
% 
% To animate the result, run the file animate_fish.m 
% 
% To try different shape changes, i.e., to prescribe the controlled shape 
% variables theta1(t) and theta2(t) modify the file shape_var.m
% 
% The files called ellipse.m and threebody.m define the geometric of the 
% three-link fish and the geometric quantities needed to solve for the 
% velocity potential
% 
% The file called influence.m contains the core of the panel method 
% that solves for the velocity potential (you should check your pseudo-code with this)
% 
% The file called admass.m computes the added masses as defined in equations 
% (3.12)-(3.14) of the JNLS paper