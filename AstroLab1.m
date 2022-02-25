%Bryan Rathke
%MAE 4410 Astrodynamics
%Lab 1

clc;clear;close all

GM=398600.5;        %Earth Standard Gravitonal Parameter km^3/s^2

a_(1)=10000;        %Orbital parameters to be simulated
e_(1)=.5;
t_(1)=10000;
a_(2)=15000;
e_(2)=.1;
t_(2)=19000;
a_(3)=30000;
e_(3)=.9;
t_(3)=58000;

for i=1:3               %Loop to test different orbits
    a=a_(i); e=e_(i);

    Rp=a*(1-e);         %Orbit radius at perigee
    Ra=a*(1+e);         %Orbit radius at apogee

    Vp=(GM*( (2/Rp)-(1/a) ) )^.5;  %Velocity at perigee

    x(1)=Rp;            %Initial conditions
    xdot(1)=0;
    y(1)=0;
    ydot(1)=Vp;
    t(1)=0;


    for n=1:t_(i)                       %Euler's method
        dt=1;
        t(n+1)=t(n)+dt;
        R=( x(n)^2 + y(n)^2 )^.5;       %Orbit radius
        Rh(n)=R;                        %Radius history
        xddot(n)=-GM*x(n)/R^3;          %x double dot (acceleration)
        xdot(n+1)=xdot(n) + xddot(n)*dt;%x dot (velocity)
        x(n+1)=x(n) + xdot(n)*dt;       %x (position)
        yddot(n)=-GM*y(n)/R^3;          %y double dot (acceleration)
        ydot(n+1)=ydot(n) + yddot(n)*dt;%y dot (velocity)
        y(n+1)=y(n) + ydot(n)*dt;       %y (position)
    end
    
    z0=[Rp; 0; 0; Vp];              %ODE45 Initial conditions
    time=[1,t_(i)];                 %ODE45 time vector
    [T,Z]=ode45(@ohsofun,time,z0);  %ODE45 magic
    
    figure(i)                       %plots on plots on plots
    txt=["a=10,000  e=0.5","a=15,000  e=0.1","a=30,000  e=0.9"];
    subplot(3,2,1)
    plot(x,y)
    title("Euler")
    axis equal 
    
    subplot(3,2,2)
    plot(Z(:,1),Z(:,2))
    title("ODE45")
    axis equal
    
    subplot(3,2,[3,4,5,6])
    plot(x,y)
    hold
    plot(Z(:,1),Z(:,2))
    title(txt(i))
    legend("Euler","ODE45")
    axis equal
    
    
    txt(i)                  %Print Ra values for comparison.
    Ra_Euler = min(x)
    Ra_ODE45 = min(Z(:,1))
    RadE=Ra_Euler+Ra
    RadO=Ra_ODE45+Ra
    
end
    
    function dz = ohsofun(time,z)   %ODE45 function thingy, very exciting.
    R=( z(1)^2 + z(2)^2 )^.5;       %Radius of the orbit
    GM=398600.5;                    %aka mu
    dz=[z(3); z(4); -GM*z(1)/R^3; -GM*z(2)/R^3];    %The first derivative.
    end
    
