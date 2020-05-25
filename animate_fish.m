
% E Kanso, March 23, 2004


N_t = length(t);

% Dimensions of each ellipse
ea = 12; eb = 1; ec = 2;
l = ea + ec;



% Discretization of the main ellipse
% 
%  phi = [phi(1)   phi(2)   ...     phi(25)]
%
%  PHI = [-------------phi-----------------]
%        [-------------phi-----------------]
%                       .
%                       .
%                       .
%        [-------------phi-----------------]
%
% PHI is (N_t,N_s) matrix
%
% M_x, M_y 
%
%---------------------



phi = linspace(0,2*pi,25);
N_s = length(phi);

PHI = ones(N_t,1)*phi;

M_x = ea*cos(PHI);
M_y = eb*sin(PHI);

EYE_x = 0.5*cos(PHI);
EYE_y = 0.5*sin(PHI);

% motion of the main ellipse
xc = l*g(:,2);
yc = l*g(:,3);
beta = g(:,1);

% xc = xc';
% yc = yc';
% beta = beta';

% Generate N_s copies of the motion of the main ellipse
Xc = xc*ones(1,N_s);
Yc = yc*ones(1,N_s);
Beta = beta*ones(1,N_s);


% compute motion of all points on the main ellipse
xM = Xc + M_x.*cos(Beta) - M_y.*sin(Beta);
yM = Yc + M_x.*sin(Beta) + M_y.*cos(Beta);


% Compute motion of c.o.m. of the other two ellipses
ac = l*cos(beta);
as = l*sin(beta);

xA = xc + ac; yA = yc + as;
xB = xc - ac; yB = yc - as;

theta1 = th1';
theta2 = th2';

xc1 = xA + l*cos(beta + theta1); yc1 = yA + l*sin(beta + theta1);
xc2 = xB - l*cos(beta + theta2); yc2 = yB - l*sin(beta + theta2);

% xeye = l + ea - 3;
% yeye = eb*sqrt(1-((ea-3)/ea)^2);
% 
% xcEYE1 = xA + xeye*cos(beta + theta1) - yeye*sin(beta + theta1); 
% ycEYE1 = yA + xeye*sin(beta + theta1) + yeye*cos(beta+theta1);
% 
% xcEYE2 = xA + xeye*cos(beta + theta1) + yeye*sin(beta + theta1); 
% ycEYE2 = yA + xeye*sin(beta + theta1) - yeye*cos(beta+theta1);

% AME 441 experiment
xeye = - (l + ea) + 3;
yeye = eb*sqrt(1-((ea-3)/ea)^2);

xcEYE1 = xB + xeye*cos(beta + theta2) - yeye*sin(beta + theta2); 
ycEYE1 = yB + xeye*sin(beta + theta2) + yeye*cos(beta+theta2);

xcEYE2 = xB + xeye*cos(beta + theta2) + yeye*sin(beta + theta2); 
ycEYE2 = yB + xeye*sin(beta + theta2) - yeye*cos(beta+theta2);



% Generate N_s copies of the motion of ellipse 1
Xc1 = xc1*ones(1,N_s);
Yc1 = yc1*ones(1,N_s);
Beta1 = (beta+theta1)*ones(1,N_s);

% Generate N_s copies of the motion of EYES
XcEYE1 = xcEYE1*ones(1,N_s);
YcEYE1 = ycEYE1*ones(1,N_s);

XcEYE2 = xcEYE2*ones(1,N_s);
YcEYE2 = ycEYE2*ones(1,N_s);


% Generate N_s copies of the motion of ellipse 2
Xc2 = xc2*ones(1,N_s);
Yc2 = yc2*ones(1,N_s);
Beta2 = (beta+theta2)*ones(1,N_s);


% Compute motion of all points on the ellipse 1 
xM1 = Xc1 + M_x.*cos(Beta1) - M_y.*sin(Beta1);
yM1 = Yc1 + M_x.*sin(Beta1) + M_y.*cos(Beta1);

% Compute motion of all points on EYE 
xEYE1 = XcEYE1 + EYE_x.*cos(Beta1) - EYE_y.*sin(Beta1);
yEYE1 = YcEYE1 + EYE_x.*sin(Beta1) + EYE_y.*cos(Beta1);

xEYE2 = XcEYE2 + EYE_x.*cos(Beta1) - EYE_y.*sin(Beta1);
yEYE2 = YcEYE2 + EYE_x.*sin(Beta1) + EYE_y.*cos(Beta1);


% Compute motion of all points on the ellipse 2
xM2 = Xc2 + M_x.*cos(Beta2) - M_y.*sin(Beta2);
yM2 = Yc2 + M_x.*sin(Beta2) + M_y.*cos(Beta2);


npts = N_t;



% make movie
       
mov = avifile('./test.avi','FPS', 25);

figure

set(gca,'AmbientLightColor','green');

% i = 1;

for i=1:1:npts
   h = fill(xM(i,:),yM(i,:),'g',xM1(i,:),yM1(i,:),'g',xEYE1(i,:),yEYE1(i,:),'b',xEYE2(i,:),yEYE2(i,:),'b',xM2(i,:),yM2(i,:),'g');
   axis image
   grid on
   set(gca,'XTick',[],'YTick',[],...
    'XLim',[-150 50],'YLim',[-60 30],...
    'FontSize',12) 
   
  
   F = getframe;
 
  mov = addframe(mov,F);
  
end

mov = close(mov);


% -------------- Forward Movie
       
% figure
% subplot(4,1,1)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(1,:),yM(1,:),'g',...
%        xM1(1,:),yM1(1,:),'g',...
%        xEYE1(1,:),yEYE1(1,:),'b',...
%        xEYE2(1,:),yEYE2(1,:),'b',...
%        xM2(1,:),yM2(1,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-50 75],'YLim',[-30 10],...
%     'FontSize',12) 
%    
%   
% subplot(4,1,2)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(42,:),yM(42,:),'g',...
%        xM1(42,:),yM1(42,:),'g',...
%        xEYE1(42,:),yEYE1(42,:),'b',...
%        xEYE2(42,:),yEYE2(42,:),'b',...
%        xM2(42,:),yM2(42,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-50 75],'YLim',[-30 10],...
%     'FontSize',12) 
% 
% subplot(4,1,3)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(84,:),yM(84,:),'g',...
%        xM1(84,:),yM1(84,:),'g',...
%        xEYE1(84,:),yEYE1(84,:),'b',...
%        xEYE2(84,:),yEYE2(84,:),'b',...
%        xM2(84,:),yM2(84,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-50 75],'YLim',[-30 10],...
%     'FontSize',12) 
% 
% subplot(4,1,4)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(126,:),yM(126,:),'g',...
%        xM1(126,:),yM1(126,:),'g',...
%        xEYE1(126,:),yEYE1(126,:),'b',...
%        xEYE2(126,:),yEYE2(126,:),'b',...
%        xM2(126,:),yM2(126,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-50 75],'YLim',[-40 0],...
%     'FontSize',12) 
% 
% 
% -------------- Turning Movie
       
% figure
% subplot(4,1,1)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(1,:),yM(1,:),'g',...
%        xM1(1,:),yM1(1,:),'g',...
%        xEYE1(1,:),yEYE1(1,:),'b',...
%        xEYE2(1,:),yEYE2(1,:),'b',...
%        xM2(1,:),yM2(1,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-40 40],'YLim',[-30 40],...
%     'FontSize',12) 
%    
%   
% subplot(4,1,2)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(42,:),yM(42,:),'g',...
%        xM1(42,:),yM1(42,:),'g',...
%        xEYE1(42,:),yEYE1(42,:),'b',...
%        xEYE2(42,:),yEYE2(42,:),'b',...
%        xM2(42,:),yM2(42,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-40 40],'YLim',[-30 40],...
%     'FontSize',12) 
% 
% subplot(4,1,3)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(84,:),yM(84,:),'g',...
%        xM1(84,:),yM1(84,:),'g',...
%        xEYE1(84,:),yEYE1(84,:),'b',...
%        xEYE2(84,:),yEYE2(84,:),'b',...
%        xM2(84,:),yM2(84,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-40 40],'YLim',[-30 40],...
%     'FontSize',12) 
% 
% subplot(4,1,4)
% set(gca,'AmbientLightColor','green');
%    h = fill(xM(126,:),yM(126,:),'g',...
%        xM1(126,:),yM1(126,:),'g',...
%        xEYE1(126,:),yEYE1(126,:),'b',...
%        xEYE2(126,:),yEYE2(126,:),'b',...
%        xM2(126,:),yM2(126,:),'g');
%    axis image
%    grid on
%    set(gca,'XTick',[],'YTick',[],...
%     'XLim',[-40 40],'YLim',[-20 50],...
%     'FontSize',12) 
