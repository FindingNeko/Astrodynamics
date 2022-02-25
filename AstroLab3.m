% 
%                 Lab 3 - Orbital Determination from Tracking Station Data
%                                 MAE 4410 Astrodynamics            
%                                      Bryan Rathke
%                    

clc;clear all;close all
%Inputs
% L - Latitude [deg]
% LST [deg]
% H - Height above sea level
% rho - range [km]
% rho dot - [km/s]
% Az - Azimuth [deg]
% Az dot [rad/s]
% El - Elevation [deg]
% El dot [rad/s]

%constants
a_earth=6378.145; %Earth semi-major axis [km]
e_earth=0.08182;  %Earth eccentricity
Re=6378.145;      %Earth radius
mu=398600.5;
omega_earth=[0;0;0.00007292115856];

L_(1)=30;
LST_(1)=45;
H_(1)=0.145;
rho_(1)=637.8145;
rho_dot_(1)=0;
Az_(1)=30;
Az_dot_(1)=0;
Az_dot_(1)=pi/180*Az_dot_(1);
El_(1)=90;
El_dot_(1)=10/0.001239446309;

L_(2)=30;
LST_(2)=262.5;
H_(2)=2;
rho_(2)=1000;
rho_dot_(2)=0;
Az_(2)=45;
Az_dot_(2)=0;
Az_dot_(2)=pi/180*Az_dot_(2);
El_(2)=90;
El_dot_(2)=0.015;

L_(3)=38.83;
LST_(3)=294.122;
H_(3)=1.839;
rho_(3)=430;
rho_dot_(3)=.5;
Az_(3)=201.1;
Az_dot_(3)=2;
Az_dot_(3)=pi/180*Az_dot_(3);
El_(3)=87;
El_dot_(3)=1;
El_dot_(3)=pi/180*El_dot_(3);

ii=0;
cc=1;
for ii = 1:3
%     L=input('Input Latitude [deg]:');
%     LST=input('Input Local Sidereal Time [deg]:');
%     H=input('Input Height Above Sea Level [km]:');
%     rho=input('Input Range [km]:');
%     rho_dot=input('Input Rate of Change of Range [km/s]:');
%     Az=input('Azimuth [km]:');
%     Az_dot=input('Rate of Change of Azimuth [rad/s]:');
%     El=input('Input Elevation [deg]:');
%     El_dot=input('Input Rate of Change of Elevation [rad/s]:');
    
    fprintf('\n\nSatellite tracking software: Case %0.0f\n\n',ii)

    L=L_(ii);
    LST=LST_(ii);
    H=H_(ii);
    rho=rho_(ii);
    rho_dot=rho_dot_(ii);
    Az=Az_(ii);
    Az_dot=Az_dot_(ii);
    El=El_(ii);
    El_dot=El_dot_(ii);
    
    
    rho_sez=[-rho*cosd(El)*cosd(Az);...     %find rho sez
            rho*cosd(El)*sind(Az);...
            rho*sind(El)];                  
    rho_dot_sez=[-rho_dot*cosd(El)*cosd(Az)+rho*sind(El)*(El_dot)*cosd(Az)+...
                rho*cosd(El)*sind(Az)*(Az_dot);...
                rho_dot*cosd(El)*sind(Az)-rho*sind(El)*El_dot*sind(Az)+...
                rho*cosd(El)*cosd(Az)*Az_dot;...    %Find rho_dot sez
                rho_dot*sind(El)+rho*cosd(El)*El_dot];
                
    Rsquig=[sind(L)*cosd(LST),-sind(LST),cosd(L)*cosd(LST);...
            sind(L)*sind(LST),cosd(LST),cosd(L)*sind(LST);... 
            -cosd(L),0,sind(L)];   %SEZ to IJK transformation matrix
    rho_ijk=Rsquig*rho_sez;        %find rho ijk
    rho_dot_ijk=Rsquig*rho_dot_sez;%find rho_dot ijk
        
    x=abs((a_earth/sqrt(1-e_earth^2*sind(L)^2))+H)*cosd(L); %ground station
    z=abs((a_earth*(1-e_earth^2))/sqrt(1-e_earth^2*sind(L)^2)+H)*sind(L);
    Rsite_ijk=[x*cosd(LST);x*sind(LST);z];      %calculate station position
    
    R_ijk=rho_ijk+Rsite_ijk;       % calculate position vector R_ijk
    V_ijk=rho_dot_ijk + cross(omega_earth,R_ijk); %veloctiy V_ijk
    

% R=input('Enter position vector in [xi yj zk] format:');
% V=input('Enter velocity vector in [xi yj zk] format:');

R=R_ijk;
V=V_ijk;
fprintf('Position vector [R]: (%0.4f i + %0.4f j + %0.4f k) km\n',R_ijk(1),R_ijk(2),R_ijk(3))
fprintf('Velocity vector [V]: (%0.4f i + %0.4f j + %0.4f k) km/s\n',V_ijk(1),V_ijk(2),V_ijk(3))

H=cross(R,V);                              %Calculate specific momentum
E=(1/mu)*((norm(V)^2-(mu/norm(R)))*R-(dot(R,V)*V)); %Eccentricity vector
N=cross([0;0;1],H);                        %Node vector

eps=(norm(V)^2 /2)-mu/norm(R);             %Epsilon, specific mech e

a=-mu/(2*eps);                              %Semi-major axis
e=norm(E);                                  %eccentricity
I=acosd(H(3)/norm(H));                      %inclination

fprintf('Semi-Major axis [a]: %0.7f km\n',a);
fprintf('Eccentricity [e]: %0.7f\n',e)
fprintf('Inclination (i): %0.7f degrees\n',I)


is_equatorial=0;                        %initializing some variables
is_circular=0;
is_elliptical=0;
is_parabolic=0;
is_hyperbolic=0;

if I>0.001                  %if so inclined
    O=acosd(N(1)/norm(N));  %find OMEGA (right ascension of ascending node)
    if N(2)<0               %quad check Nj<0
        O=360-O;            %adjust OMEGA based on quad check
    end
    fprintf('Right ascenscion of the ascending node [OMEGA]: %0.7f degrees\n',O)
    if e<0.001              %if circular
        is_circular=1;      %yep
        u=acosd(dot(N,R)/(norm(N)*norm(R)));  %Argument of latitude

        if R(3)<0           %quad check
            u=360-u;        %adjust u based on quad check
        end
        fprintf('Argument of latitude [u]: %0.7f degrees\n',u)
    else                    %if not circular
        uu=acosd(dot(N,E)/(norm(N)*norm(E)));  %omega argument of perigee
        if E(3)<=0              %quad check
            uu=360-uu;          %adjust omega based on quad check
        end
        fprintf('Argument of perigee [omega]: %0.7f degrees\n',uu)
    end
end
if e>=0.001           %if not circular
    nu=acosd(dot(E,R)/(norm(E)*norm(R)));  %True anomaly
    if dot(R,V)<0       %quad check
        nu=360-nu;      %adjust true anomaly based on quad check
    end

    fprintf('True anomaly [nu]: %0.7f degrees\n',nu)
    
    if e>0.999 && e<1.001    %if parabolic
        is_parabolic=1;
    end
    if e>=1.001             %if hyperbolic
        is_hyperbolic=1;
    end
    if e<=0.999             %if elliptical
        is_elliptical=1;
    end
    PI=O+uu;                %BIG PI logitude of perigee
end
if I<0.001 || I>179.999    %if not so inclined
    is_equatorial=1;
    if e>0.001              %if not circular
        l=PI+nu;            %true longitude
    end
    if e<0.001              %if circular
        l=acosd(R(1)/R);    %true longitude
    end
    fprintf('True longitude [l]: %0.7f degrees\n',l)

end
if is_equatorial==1
    disp('Orbit is equatorial')
end
if is_circular==1
    disp('Orbit is circular')
end
if is_elliptical==1
    disp('Orbit is elliptical')
end
if is_parabolic==1
    disp('Trajectory is parabolic')
end
if is_hyperbolic==1
    disp('Tajectory is hyperbolic')
end
        

    
    cc=cc+1;
%    ii=input('Enter 1 to check another case.') 
     
end


