%%% MIE301 A frame
%%
% Author: Naib Gulam
% Created on: [MM/DD/YYYY]
% 
% Last edited by: Max Kopstein
% Edited on: 11/29/2018 [MM/DD/YY]

%IMPORTANT: MAKE SURE CoefficientMatrix.csv IS IN THE SAME DIRECTORY AS
%           THIS CODE!!!!!

close all; % closes all figures
clear; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
% using axis()
xmin= -1;  % leftmost window edge
xmax= 5;   % rightmost window edge
ymin= -1; % bottom window edge
ymax= 6;  % top window edge

%% Variable Initialization
% Although not every variable used is declared right at the start, most of
% the more important ones are declared here.

couplerPointLength = 0;
couplerPointAngle = pi/4;

scalingFactor = 1.25; % 1 means distances are in metres, scale accordingly
% Link lengths, see diagram for reference
r1 = scalingFactor*0.4;
r2 = scalingFactor*0.3;
r3 = scalingFactor*0.2;
r4 = scalingFactor*0.3;
r5 = scalingFactor*0.1;
r6 = scalingFactor*0.2;
r7 = scalingFactor*0.2;
r8 = scalingFactor*0.1;

linden = 0.27; % Linear density of links [kg/m]

% Link masses (assuming uniform linear densities) [kg]
m1 = linden*r1;
m2 = linden*r2 + linden*r5;     % Combining links 2 and 5 to simplify force analysis
m3 = linden*r3;
m4 = linden*r4 + linden*r8;     % Combining links 4 and 8 to simplify force analysis
m6 = linden*r6;
m7 = linden*r7;

m9 = 0.766; % Combined mass of the desk and geneva mechanisms [kg]
%angle of the base beam (leave at 0)
theta1 = 0;

% Find Equations for Theta2 Theta3 Theta6 Theta7
[theta2_eqn, theta3_eqn] = solveEquation(r1,r2,r3,r4,theta1);
[theta6_eqn, theta7_eqn] = solveSecondEquation(r3,r5,r6,r7,r8);

% independent variable, 0.9999 is there to prevent errors
theta4 = linspace(2.497,pi-acos(2/6),60); %CHANGE THE 2/6!!!!!

omega4 = 0.5;                           % Desired input angular velocity [rad/s]
time = (abs(theta4-theta4(1)))./omega4; % total elapsed time for at each frame
dt = mean(gradient(time));              % Average time step between frames (all time steps between frames should not vary in theory, taking average just eliminates any numerical error)

g = 9.81;   % Acceleration due to gravity (before scaling) [m/s^2]

Coef = csvread('CoefficientMatrix.csv',1,0);        % Matrix of component force and moment coefficients (file only contains rows pertaining to equations for x and y components, rows for angular components are are filled with zeros initially)
x = zeros(length(Coef),length(time)); % Matrix of resultant time varying forces and moments (format is x(i,j) = Force/Moment i at time frame j, e.g. x(1,5) = F14x(5), x(6,10) = F47y(10) )

%% calculate mechanism motion
% step through motion of the mechanism  
%for i=1:length(time)

% function solves the linear vector equation
%[theta2temp, theta3temp] = solveEquation(r1,r2,r3,r4,theta1,theta4(i));

% stores the answer
theta2 = double(subs(theta2_eqn));
theta3 = double(subs(theta3_eqn));

% calculates the the points A,B,C, and D
Ax = zeros(1,length(time));                              
Ay = zeros(1,length(time));

Bx = Ax + r2*cos(theta2);             
By = Ay + r2*sin(theta2);

Cx = Bx + r3*cos(theta3);             
Cy = By + r3*sin(theta3);

Dx = Cx - r4*cos(theta4);             
Dy = Cy - r4*sin(theta4);  


% angles 5 and 8 are equal to 2 and 4 repectively
theta5=theta2;
theta8=theta4;

% stores the answer
theta6 = double(subs(theta6_eqn));
theta7 = double(subs(theta7_eqn));

% calculates the the points E,F,G
Ex = Bx + r5*cos(theta5);             
Ey = By + r5*sin(theta5);

Fx = Ex + r6*cos(theta6);             
Fy = Ey + r6*sin(theta6);

Gx = Fx - r7*cos(theta7);             
Gy = Fy - r7*sin(theta7);

couplerPointX = Fx + couplerPointLength*cos(couplerPointAngle);
couplerPointY = Fy + couplerPointLength*sin(couplerPointAngle);

t=5;

%% Animate motion of system thruogh time
figure(1);
for i=1:length(time);
    
    % clears the previous figure
    clf;
    
    hold on;
    plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','b','LineWidth',3 );
    plot( [Bx(i) Cx(i)], [By(i), Cy(i)],'Color','g','LineWidth',3 );
    plot( [Cx(i) Dx(i)], [Cy(i), Dy(i)],'Color','c','LineWidth',3 );
    plot( [Ax(i) Dx(i)], [Ay(i), Dy(i)],'Color','black','LineWidth',3 );
    
    plot( [Bx(i) Ex(i)], [By(i), Ey(i)],'Color','b','LineWidth',3 );
    plot( [Ex(i) Fx(i)], [Ey(i), Fy(i)],'Color',[255/255,165/255,0],'LineWidth',3 );
    plot( [Fx(i) Gx(i)], [Fy(i), Gy(i)],'Color','r','LineWidth',3 );
    plot( [Gx(i) Cx(i)], [Gy(i), Cy(i)],'Color','c','LineWidth',3 );
    
    % add a couple point to the top point
    plot( [Fx(i) couplerPointX(i)], [Fy(i), couplerPointY(i)],'Color','m','LineWidth',3 );
    line([couplerPointX(1) couplerPointX(i)],[couplerPointY(1) couplerPointY(i)],'Color','red','LineStyle','--')
    
    text(Fx(i)+0.025,Fy(i)-0.0125,'P','FontSize',12,'FontWeight','bold');            % label point F and title it P
    plot(Fx(i),Fy(i),'.k','MarkerSize',15 );
    
    hold on; 
    grid on;                            % add a grid to the figure
    axis equal;                         % make sure the figure is not stretched
    title('A Frame');                   % add a title to the figure
    %axis( [xmin xmax ymin ymax] );      % set figure axis limits
    line([Fx(1) Fx(i)],[Fy(1) Fy(i)],'Color','red','LineStyle','--');
    xlabel('X Position (m)'); 
    ylabel('Y Position (m)');
    pause(0.001);
    
end
text((Dx(i)-Ax(i))/2+1, 2, '1','color','black');            % label link 1
text(((Ex(i)-Ax(i))/2)-3, (Ey(i)-Ay(i))/2,'2','color','b');            % label link 2
text(((Dx(i)-Gx(i))/2)+Gx(i)+2, (Ey(i)-Ay(i))/2,'3','color','b');            % label link 2
text(((Cx(i)-Bx(i))/2)+Bx(i)+1, Cy(i)+2,'4','color','g');            % label link 2

text(((Fx(i)-Ex(i))/2)-3+Ex(i), (Fy(i)-Ey(i))/2+Ey(i),'5','color','r');            % label link 2
text(((Gx(i)-Fx(i))/2)+Fx(i)+2, (Fy(i)-Gy(i))/2+Gy(i),'6','color','r');

velocityOfPointFy = gradient(Fy)./gradient(time);
deployed = find(Fy>=(Fy(end)-13),1);
theta_report = fliplr(theta4.*(180/pi))-110+33.6;

%% Inertia Calculation
% IMPORTANT NOTE:
%   For this section link 8 has been incorporated into link 4, meaning
%   link 4 extends from point D to point G. Similarly, link 5 has been
%   incorporated into link 2 so the at link 2 now extends from point A to
%   point E. This is done to simplify the force analysis. This is possible
%   since links 4 and 8 are both colinear and rigidly connected, as are
%   links 2 and 5.

% Second Note:
%   The desk/geneva assembly itself is being modeled as a mass attached 
%   to point F and is assigned link number 9. It is assumed that any moment
%   resulting from the CG of the desk/geneva assembly being located at some
%   distance from point F will be dealt with by some other mechanism, and 
%   so only the effects of its inertia and gravitational force are
%   considered.

% Absolute position vectors for the CG of each link and the desk/geneva assembly (format is Rg = [Rgx;Rgy])
Rg2 = ([Ex;Ey]-[Ax;Ay])./2;
Rg3 = ([Bx;By]+[Cx;Cy])./2;
Rg4 = ([Dx;Dy]+[Gx;Gy])./2;
Rg6 = ([Ex;Ey]+[Fx;Fy])./2;
Rg7 = ([Gx;Gy]+[Fx;Fy])./2;
Rg9 = [Fx;Fy];

% Differentiating CG position vectors w.r.t time to find CG velocity vectors (i.e. Vg = dRg/dt)
Vg2 = gradient(Rg2,dt);
Vg3 = gradient(Rg3,dt);
Vg4 = gradient(Rg4,dt);
Vg6 = gradient(Rg6,dt);
Vg7 = gradient(Rg7,dt);
Vg9 = gradient(Rg9,dt);

% Differentiating the angle of each link to find their angular velocities (i.e. omega = (d/dt)*(theta))
omega2 = gradient(theta2,dt);
omega3 = gradient(theta3,dt);
%omega4 is constant
omega6 = gradient(theta6,dt);
omega7 = gradient(theta7,dt);

% Differentiating CG velocity vectors w.r.t time to find CG acceleration vectors (i.e. ag = dVg/dt)
ag2 = gradient(Vg2,dt);
ag3 = gradient(Vg3,dt);
ag4 = gradient(Vg4,dt);
ag6 = gradient(Vg6,dt);
ag7 = gradient(Vg7,dt);
ag9 = gradient(Vg9,dt);

% Differentiating the angular velocity of each link to find their angular accelerations (i.e. alpha = (d/dt)*omega)
alpha2 = gradient(omega2,dt);
alpha3 = gradient(omega3,dt);
alpha4 = zeros(2,length(time)); % Since omega4 is constant
alpha6 = gradient(omega6,dt);
alpha7 = gradient(omega7,dt);

% Matrix of component inertias (format is Inrt(i,j) = Component Inertia i at time frame j, e.g. Inrt(4,8) = m7*ag7(1,8), Inrt(12,3) = (m3/12)*(r3^2).*alpha3(3) )
Inrt = [m4*ag4(1,:);
        m4*ag4(2,:)+m4*g;
        ((Gx-Dx)./2)*m4*g;
        
        m2*ag2(1,:);
        m2*ag2(2,:)+m2*g;
        (m2/3)*(r2^2).*alpha2 + ((Ex-Ax)./2)*m2*g;
        
        m3*ag3(1,:);
        m3*ag3(2,:)+m3*g;
        (m3/3)*(r3^2).*alpha3 + ((Cx-Bx)./2)*m3*g;
        
        m7*ag7(1,:);
        m7*ag7(2,:)+m7*g;
        (m7/3)*(r7^2).*alpha7 + ((Gx-Fx)./2)*m7*g;
        
        m6*ag6(1,:);
        m6*ag6(2,:)+m6*g + m9*ag9(2,:)+m9*g;
        (m6/3)*(r6^2).*alpha6 + ((Ex-Fx)./2)*m6*g];
        


%% Force Analysis
% IMPORTANT NOTE:
%   For this section link 8 has been incorporated into link 4, meaning
%   link 4 extends from point D to point G. Similarly, link 5 has been
%   incorporated into link 2 so the at link 2 now extends from point A to
%   point E. This is done to simplify the force analysis. This is possible
%   since links 4 and 8 are both colinear and rigidly connected, as are
%   links 2 and 5.
%
% Second Note:
%   The desk/geneva assembly itself is being modeled as a mass attached 
%   to point F and is assigned link number 9. It is assumed that any moment
%   resulting from the CG of the desk/geneva assembly being located at some
%   distance from point F will be dealt with by some other mechanism, and 
%   so only the effects of its inertia and gravitational force are
%   considered.


% Loop for solving for all component forces and moments in each time frame
for i=1:length(time)
    
    % Start of moment coefficient calculation
    % Coefficients of moment equation for link 4
    Coef(3,3) = Cy(i)-Dy(i);
    Coef(3,4) = -(Cx(i)-Dx(i));
    Coef(3,5) = Gy(i)-Dy(i);
    Coef(3,6) = -(Gx(i)-Dx(i));
    
    % Coefficients of moment equation for link 2
    Coef(6,9) = Ey(i)-Ay(i);
    Coef(6,10) = -(Ex(i)-(i));
    Coef(6,11) = By(i)-Ay(i);
    Coef(6,12) = -(Bx(i)-Ax(i));
    
    % Coefficients of moment equation for link 3
    Coef(9,3) = Cy(i)-By(i);
    Coef(9,4) = -(Cx(i)-Bx(i));
    
    % Coefficients of moment equation for link 7
    Coef(12,5) = Gy(i)-Fy(i);
    Coef(12,6) = -(Gx(i)-Fx(i));
    
    % Coefficients of moment equation for link 6
    Coef(15,9) = Ey(i)-Fy(i);
    Coef(15,10) = -(Ex(i)-Fx(i));
    % End of moment coefficient calculation
    
    % Calculate i'th column of x matrix given Coef matrix and Inrt vector at time frame i
    %disp(i);
    if(i==30)
        disp('');
        disp('');
        %disp(Coef);
        disp('');
        disp('');
        Test = [Coef,Inrt(:,i)];
    end
    x(:,i) = Coef\Inrt(:,i);
    
end

%% Force Analysis Plot Generation
figure(2);
hold on;
yyaxis left;
gca.YColor = [0,0.447,0.741];
plot(theta_report,x(15,:),'Color',[0,0.447,0.741],'linewidth',1);
actuate = find(Fy<=(Fy(end)-0.13),1,'last');
ylabel('Required Torque into Link 2[Nm]');
%text(theta_report(actuate)-1,x(15,actuate)-15, sprintf('Geneva Actuation Complete'));
%plot(theta_report(actuate),x(15,actuate),'.k','MarkerSize',15);


yyaxis right;
gca.YColor = [0.0,0.325,0.098];
plot(theta_report,(x(15,:)./(2.5*0.6)),'linewidth',1);
ylabel('Required User Input Force [N]');
[Peak, PeakIdx] = findpeaks((x(15,:)./(2.5*0.6)));
text(theta_report(PeakIdx)+4, Peak, sprintf('Max Force = %6.3f N', Peak),'Color','k');
plot(theta_report(PeakIdx),Peak,'.k','MarkerSize',15);
%plot(theta_report(actuate),(x(15,actuate)./(2.5*0.6)),'.k','MarkerSize',15);


xlabel('\theta_{2} [degrees]');
grid on;
title('Force Analysis for Chosen Range of Motion');



%% Solve for Theta2 and Theta3 in first loop
function [theta2temp,theta3temp] = solveEquation(r1,r2,r3,r4,theta1)
    syms x1 y1 theta4
    
    eqn1 = r2*cos(x1)+r3*cos(y1)-r4*cos(theta4) == r1*cos(theta1);
    eqn2 = r2*sin(x1)+r3*sin(y1)-r4*sin(theta4) == r1*sin(theta1);
    
    sol = solve([eqn1, eqn2], [x1, y1]);
    theta2temp = real(sol.x1(1));
    theta3temp = real(sol.y1(1));
end

%% Solve for Theta6 and Theta7 in second loop
function [theta6temp,theta7temp] = solveSecondEquation(r3,r5,r6,r7,r8)
    syms x2 y2 theta3 theta8 theta5
    
    eqn1 = r6*cos(x2)-r7*cos(y2) == -r5*cos(theta5)+r8*cos(theta8)+r3*cos(theta3);
    eqn2 = r6*sin(x2)-r7*sin(y2) == -r5*sin(theta5)+r8*sin(theta8)+r3*sin(theta3);
    
    sol = solve([eqn1, eqn2], [x2, y2]);
    theta6temp = real(sol.x2(2));
    theta7temp = real(sol.y2(2));
end