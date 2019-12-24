function doublependulum(l1,l2,m1,m2,theta1,theta2,timespan)
%CTINS_DOUBLEPEND will simulate the motion of a double pendulum with given
%top length (l1), bottom length (l2), top mass (m1), bottom mass (m2), and
%initial angles (theta1 and theta2) over specified time (timespan).

close all

if nargin==6
    timespan=10;
elseif nargin>7
    error('Too Many Input Arguments');
elseif nargin<=5
    error('Not Enough Input Arguments');
end

%Initializing variables
theta1_prime=0;
theta2_prime=0; %initial velocities

y0=[theta1 theta1_prime theta2 theta2_prime]; %Initial conditions

%y(:,1) is theta1, y(:,2) is velocity of m1, y(:,3) is theta2, and
%y(:,4) is velocity of m2.
options=odeset('RelTol',1e-3); %Low tolerance for optimal solving speed.
[t,y]=ode45(@pend,[0,timespan],y0,options,l1,l2,m1,m2);

%Coordinates for m1 and m2 points
x1=l1*sin(y(:,1));
y1=-l1*cos(y(:,1)); %m1
x2=l1*sin(y(:,1))+l2*sin(y(:,3));
y2=-l1*cos(y(:,1))-l2*cos(y(:,3)); %m2

%Visualizing pathway of each mass
figure(1)
%Plot path of m1
plot(x2,y2,'r'); hold on;
%Plot pivot point
plot(0,0,'go')
%Plot path of m2
plot(x1,y1,'b'); hold off;
h=gca;
set(h,'Fontsize',12); %Changes fontsize on axes
legend('Path of M1','Pivot Point','Path of M2');
title('Traced Movement of Double Pendulum','FontSize',20);
xlabel('X Position'); ylabel('Y Position');
%Fit Axis to pendulum radius
axis([-(l1+l2) (l1+l2) -(l1+l2) (l1+l2)]); axis square;

figure(2)
subplot(2,1,1)
plot(t,y(:,1),'b'); hold on; %Plots Theta1
plot(t,y(:,3),'r'); hold off; %Plots theta2
title('Graph of Theta vs. Time','Fontsize',14);
xlabel('Time (t)'); ylabel('Theta');
legend('\theta_1 (Angle of M1)','\theta_2 (Angle of M2)');
g=gca;
set(g,'Fontsize',9);

subplot(2,1,2)
plot(t,y(:,2),'b'); hold on; %Plots Theta1
plot(t,y(:,4),'r'); hold off; %Plots theta2
title('Graph of Velocity vs. Time','Fontsize',14);
xlabel('Time (t)'); ylabel('Velocity');
legend('\theta_1 Prime (Velocity of M1)','\theta_2 Prime (Velocity of M2)');
c=gca;
set(c,'Fontsize',9);

fram=0;
%Animation of the Double Pendulum
figure(3)
for i=1:length(y)
    %Plot pivot point of pendulum
    plot(0,0,'g.','markersize',25); hold on;
    %Plot m1 point (size related to mass)
    plot(x1(i),y1(i),'b.','markersize',30.*(m1./2));
    %Plot m2 point (size related to mass)
    plot(x2(i),y2(i),'r.','markersize',30.*(m2./2)); hold off
    %Set axis to fit the pendulum
    axis([-(l1+l2) (l1+l2) -(l1+l2) (l1+l2)]); axis square;
    %Line from pivot point to m1
    line([0 x1(i)], [0 y1(i)], 'linewidth',2);
    %Line from m1 to m2
    line([x1(i) x2(i)],[y1(i),y2(i)],'linewidth',2);
    %Graph Appearance details
    xlabel('X Position'); ylabel('Y Position');
    title('Animation of Double Pendulum Motion','FontSize',18);
    g=gca;
    set(g,'FontSize',12);
    %Reference number of frames for movie
    A=rem(i,3);
    if A ~= 0
        fram=fram+1;
        %Obtains the current frame for movie
        F=getframe;
    end
end
%Animates the motion of the double pendulum
movie(F,fram,60);

function yprime = pend(t,y,l1,l2,m1,m2)
g=9.81; %gravitational constant

%Segments of equations needed for velocities of m1 and m2
a = (m1+m2)*l1;
b = m2*l2*cos(y(1)-y(3));
c = m2*l1*cos(y(1)-y(3));
d = m2*l2;
e = -m2*l2*y(4)* y(4)*sin(y(1)-y(3))-g*(m1+m2)*sin(y(1));
f = m2*l1*y(2)*y(2)*sin(y(1)-y(3))-m2*g*sin(y(3));

%Updating array
yprime(1)= y(2); %New theta1
yprime(3)= y(4); %New theta2
yprime(2)= (e*d-b*f)/(a*d-c*b); %New theta1_prime (velocity m1)
yprime(4)= (a*f-c*e)/(a*d-c*b); %New theta2_prime (velocity m2)

%Orients yprime for ODE solver
yprime=yprime';
