% 
%                   Lab 4 - Satellite Tracking and Acquisition Software
%                                 MAE 4410 Astrodynamics            
%                                      Bryan Rathke
%                    

%
% Explanation of problem 2:
% For this problem, tracking station POGO at Thule air base, Greenland was
% selected to track the satellite. This site was selected because the
% position of the satellite places it above the local horizon at the LST
% at which the satellite will be sighted.
% 

%Bryan Rathke
%Lab 3 - COEs from radar tracking station data
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

L_(1)=76.53;
Long_(1)=-64.7;
GST_(1)=53/3;
Al_(1)=1.207;
rho_(1)=3000;
rho_dot_(1)=6;
Az_(1)=7.5;
Az_dot_(1)=1;
Az_dot_(1)=pi/180*Az_dot_(1);
El_(1)=85;
El_dot_(1)=0.01;
El_dot_(1)=pi/180*El_dot_(1);
TOF_(1)=14*24;

L_(2)=38.8;
Long_(2)=-104.54;
GST_(2)=22;
Al_(2)=1.915;
rho_(2)=2121.4180;
rho_dot_(2)=-3.32040;
Az_(2)=350;
Az_dot_(2)=-0.07653;
Az_dot_(2)=pi/180*Az_dot_(2);
El_(2)=35.3507;
El_dot_(2)=0.20367;
El_dot_(2)=pi/180*El_dot_(2);
TOF_(2)=96;

ii=0;
cc=1;
for ii = 1:2
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

    GST=GST_(ii)*15;    
    L=L_(ii);
    Long=Long_(ii);
    LST=GST+Long;
    %LST=LST_(ii);
    Al=Al_(ii);
    rho=rho_(ii);
    rho_dot=rho_dot_(ii);
    Az=Az_(ii);
    Az_dot=Az_dot_(ii);
    El=El_(ii);
    El_dot=El_dot_(ii);
    
    fprintf('Case %0.0f initial tracking station data:\n',ii)
    fprintf('GST [Greenwich Sidereal Time]: %0.3f degrees\n',GST)
    fprintf('LST [Local Sidereal Time]: %0.3f degrees\n',LST)
    fprintf('Latitude: %0.7f degrees\n',L)
    fprintf('Longitude: %0.7f degrees\n',Long)
    fprintf('Altitude: %0.3f km\n\n',Al)
    disp('Satellite sighting data from first observation:')
    fprintf('Range: %0.4f km\n',rho)
    fprintf('Azimuth: %0.7f degrees\n',Az)
    fprintf('Elevation: %0.7f degrees\n',El)
    fprintf('Range Rate of Change: %0.7f km/s\n',rho_dot)
    fprintf('Azimuth Rate of Change: %0.7f degrees/s\n',Az_dot)
    fprintf('Elevation Rate of Change: %0.7f degrees/s\n\n',El_dot)
    
    
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
        
    x=abs((a_earth/sqrt(1-e_earth^2*sind(L)^2))+Al)*cosd(L); %ground station
    z=abs((a_earth*(1-e_earth^2))/sqrt(1-e_earth^2*sind(L)^2)+Al)*sind(L);
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
        

    TOF=TOF_(ii);           %Perform TOF calculations
    fprintf('\nAfter %0.3f hours, Satellite''s new position is:\n',TOF)
    
    nu_r=degtorad(nu);      
    
    E=acos((e+cos(nu_r))/(1+e*cos(nu_r)));      %Find E
    if nu>180
        E=(2*pi)-E;                             %HP check
    end
    
    E_n=pi;                                     %Guess new E
    M_o=(E)-e*sin(E);                           %Find M_o
    M=(mu/a^3)^.5 *(TOF*3600)+M_o;              %Find M
    M_n=E_n - e*sin(E_n);                       %Find M_n
    while M>2*pi                                %M is less than 2pi
        M=M-2*pi;
    end
    while abs(M-M_n)>=0.001
                                            %Newtons magic loop
    E_n=E_n+((M-M_n)/(1-e*cos(E_n)));
    M_n=E_n - e*sin(E_n);  
    end
    
    nu_f=acosd((cos(E_n)-e)/(1-e*cos(E_n)));    %Find nu final
    if E_n>pi
        nu_f=360-nu_f;                          %HP check nu final
    end
    
    fprintf('True Anomaly: %0.7f degrees.',nu_f)
    
    r=a*(1-e*cos(E_n));
    
    p=norm(H)^2 /mu;
    
    R_pqw=[r*cosd(nu_f);r*sind(nu_f);0];
    V_pqw=(mu/p)^.5 .* [-sind(nu_f);(e+cosd(nu_f));0];
                      %PQW to IJK transformation matrix        
    PQWtoIJK=[cosd(O)*cosd(uu)-sind(O)*sind(uu)*cosd(I),...
        -cosd(O)*sind(uu)-sind(O)*cosd(uu)*cosd(I),...
        sind(O)*sind(I);...
        sind(O)*cosd(uu)+cosd(O)*sind(uu)*cosd(I),...
        -sind(O)*sind(uu)+cosd(O)*cosd(uu)*cosd(I),...
        -cosd(O)*sind(I);...
        sind(uu)*sind(I),cosd(uu)*sind(I),cosd(I)];
    
    R_f=PQWtoIJK*R_pqw;
    V_f=PQWtoIJK*V_pqw;
    
    fprintf('\nR = %0.7f i + %0.7f j %0.7f k km\n',R_f(1),R_f(2),R_f(3))
    fprintf('V = %0.7f i + %0.7f j + %0.7f k km/s\n',V_f(1),V_f(2),V_f(3))
   
    if ii==2            %Use this site to find the satellite
    GST=mod(((GST_(ii)+TOF)*15),360);    
    L=76.531111;
    Long=-68.703056;
    LST=GST+Long;
    Al=0.112;
    x=abs((a_earth/sqrt(1-e_earth^2*sind(L)^2))+Al)*cosd(L); %ground station
    z=abs((a_earth*(1-e_earth^2))/sqrt(1-e_earth^2*sind(L)^2)+Al)*sind(L);
    Rsite_ijk=[x*cosd(LST);x*sind(LST);z];      %calculate station position
    
    Rsquig=[sind(L)*cosd(LST),-sind(LST),cosd(L)*cosd(LST);...
            sind(L)*sind(LST),cosd(LST),cosd(L)*sind(LST);... 
            -cosd(L),0,sind(L)];   %SEZ to IJK transformation matrix
    end    
    
    
    rho_ijk=R_f-Rsite_ijk;              %Find new range-site information
    rho_sez=transp(Rsquig)*rho_ijk;
    
    rho_f=norm(rho_sez);
    El_f=asind(rho_sez(3)/rho_f);
    Az_f=acosd(-rho_sez(1)/(rho_f*cosd(El_f)));
    
    fprintf('\nCase %0.0f final tracking station data:\n',ii)
    fprintf('GST [Greenwich Sidereal Time]: %0.3f degrees\n',GST)
    fprintf('LST [Local Sidereal Time]: %0.3f degrees\n',LST)
    fprintf('Latitude: %0.7f degrees\n',L)
    fprintf('Longitude: %0.7f degrees\n',Long)
    fprintf('Altitude: %0.3f km\n\n',Al)
    disp('Satellite sighting information from final tracking station:')
    fprintf('Range: %0.4f km\n',rho_f)
    fprintf('Azimuth: %0.7f degrees\n',Az_f)
    fprintf('Elevation: %0.7f degrees\n',El_f)
    
    cc=cc+1;
%    ii=input('Enter 1 to check another case.') 

     
end


