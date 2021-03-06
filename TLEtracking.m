
%		Satellite Tracking using TLE Set Data
%              MAE 4410 Astrodynamics
%                  Bryan Rathke
% 


clc;clear all;close all

%constants
a_earth=6378.145; %Earth semi-major axis [km]
e_earth=0.08182; %Earth eccentricity
Re=6378.145;      %Earth radius
mu=398600.5;
omega_earth=[0;0;0.00007292115856];

%Case one
% ISS (ZARYA)
% 1 25544U 98067A 17105.87060019 .00003156 00000-0 54969-4 0 9997
% 2 25544 51.6434 348.9165 0006992 60.4905 51.3661 15.54103949 52067

meanm_(1)=15.54103949*2*pi/(3600*24);         %Mean motion converted to rad/sec
a_(1)= (mu/meanm_(1)^2)^(1/3);
e_(1)= 0.0006992; 
i_(1)= 51.6434;
OMEGA_(1)= 348.9165; 
omega_(1)= 60.4905; 
meana_(1)=51.3661; 
ndot_(1)=0.00003156;

%Case two
% EO-1 (Earth Observing)
% 1 26619U 00075A   17111.67947533  .00000840  00000-0  16733-3 0  9996
% 2 26619  97.8397 143.5394 0008529 285.0815  74.9399 14.64376714875263

meanm_(2)=14.643767148*2*pi/(3600*24);         %Mean motion converted to rad/sec
a_(2)= (mu/meanm(2)^2)^(1/3);
e_(2)= 0.0008529;
i_(2)= 97.8397;
OMEGA_(2)= 143.5394; 
omega_(2)= 285.0815; 
meana_(2)=74.9399; 
ndot_(2)=0.00000840;

%Case three
% METEOSAT-8 (MSG-1)
% 1 27509U 02040B   17111.93370946  .00000136  00000-0  00000-0 0  9991
% 2 27509   4.7447  58.5438 0000844 164.9477   4.3783  1.00281083 53734

meanm_(3)=1.00281083*2*pi/(3600*24);        %Mean motion converted to rad/sec
a_(3)= (mu/meanm_(3)^2)^(1/3);
e_(3)= 0.0000844; 
i_(3)= 4.7447;
OMEGA_(3)= 58.5438; 
omega_(3)= 164.9477; 
meana_(3)=4.3783; 
ndot_(3)=0.00000136;

%Case four
% COSMOS 2361 (Russian)
% 1 25590U 98076A   17111.70709955  .00000046  00000-0  32887-4 0  9996
% 2 25590  82.9339 157.7326 0032539  43.7054  95.6061 13.72931187917920

meanm_(4)=13.729311879*2*pi/(3600*24);        %Mean motion converted to rad/sec

a_(4)= (mu/meanm_(4)^2)^(1/3);
e_(4)= 0.0032539;
i_(4)= 82.9339;
OMEGA_(4)= 157.7326; 
omega_(4)= 43.7054; 
meana_(4)=95.6061; 
ndot_(4)=0.00000046;

ii=0; 
cc=1; 
for ii = 1:4

   meanm=meanm_(ii); 
   a=a_(ii);
   e=e_(ii);
   i=i_(ii);
   OMEGA=OMEGA_(ii); 
   omega=omega_(ii); 
   meana=meana_(ii); 
   ndot=ndot_(ii);

   I=i;
   O=OMEGA; 
   uu=omega;

   E_n=pi;                                %Guess new E
   M_o=meana*2*pi/360;                    %Find M_o 
   M=M_o;                                 %Find M 
   M_n=E_n - e*sin(E_n);                    %Find M_n
   while abs(M-M_n)>=0.001
                                         %Newtons magic loop
   E_n=E_n+((M-M_n)/(1-e*cos(E_n))); 
   M_n=E_n - e*sin(E_n);

   end

   nu_i=acosd((cos(E_n)-e)/(1-e*cos(E_n))); %Find nu initial

   if E_n>pi
      nu_i=360-nu_i;                      %HP check nu
   end

   fprintf('Case %0.0f:\n\n',ii)
   fprintf('Initial orbital elements:\n') 
   fprintf('Semi-major axis: %0.7f km\n',a) 
   fprintf('Eccentricity: %0.7f\n',e) 
   fprintf('Inclination: %0.7f degrees\n',i) 
   fprintf('RAAN: %0.7f degrees\n',OMEGA) 
   fprintf('Argument of Perigee: %0.7f degrees\n',omega) 
   fprintf('True Anomaly: %0.7f degrees.\n\n',nu_i)

   th=3/2;
   p=a*(1-e^2); 
   Repo= (Re/p)^2;

   ndot=ndot*2*pi/(3600*24)^2; 
   n=sqrt(mu/a^3);

J2=.00108262668; 
nbar=(1+(th*J2*Repo*sqrt(1-e^2)*(1-(th*sind(i)^2))))*n*360/(2*pi); 
Odot_J2=(-th*J2*Repo*cosd(i))*nbar; 
wdot_J2=(th*J2*Repo*(2-((5/2)*sind(i)^2)))*nbar; 
adot=-(a*2*ndot)/(3*n);
edot=-(1-e)*2*ndot/(3*n);

TOF=2*7*24*3600;
a_f=a+adot*TOF;
e_f=e+edot*TOF;
O=OMEGA+Odot_J2*TOF; 
w=omega+wdot_J2*TOF;

E_n=pi;                                 %Guess new E 
M=n*(TOF)+ ndot*(TOF)^2 + M_o;           %Find new M
M_n=E_n - e*sin(E_n);                      %Find new Mn
M=mod(M,2*pi);
while abs(M-M_n)>=0.001                        %Newtons magic loop
E_n=E_n+((M-M_n)/(1-e*cos(E_n)));
M_n=E_n - e*sin(E_n); 
end

nu_f=acosd((cos(E_n)-e)/(1-e*cos(E_n))); %Find nu final
if E_n>pi
    nu_f=360-nu_f;                       %HP check nu final
end

uu=w; 

r=a_f*(1-e_f*cos(E_n));

R_pqw=[r*cosd(nu_f);r*sind(nu_f);0];
V_pqw=(mu/p)^.5 .* [-sind(nu_f);(e+cosd(nu_f));0];
                %PQW to IJK transformation matrix 
PQWtoIJK=[cosd(0)*cosd(uu)-sind(0)*sind(uu)*cosd(I),...
    -cosd(0)*sind(uu)-sind(0)*cosd(uu)*cosd(I),... 
    sind(0)*sind(I);...
    sind(0)*cosd(uu)+cosd(0)*sind(uu)*cosd(I),...
    -sind(0)*sind(uu)+cosd(0)*cosd(uu)*cosd(I),...
    -cosd(0)*sind(I);...
    sind(uu)*sind(I),cosd(uu)*sind(I),cosd(I)];

R_f=PQWtoIJK*R_pqw; 
V_f=PQWtoIJK*V_pqw;

disp('Final orbital elements:')
fprintf('Semi-major axis: %0.7f km\n',a_f)
fprintf('Eccentricity: %0.7f\n',e_f)
fprintf('Inclination: %0.7f degrees\n',i)
fprintf('RAAN: %0.7f degrees\n',O)
fprintf('Argument of Perigee: %0.7f degrees\n',w)
fprintf('True Anomaly: %0.7f degrees.\n',nu_f)
fprintf('\nR = %0.7f i + %0.7f j %0.7f k km\n',R_f(1),R_f(2),R_f(3)) 
fprintf('V = %0.7f i + %0.7f j + %0.7f k km/s\n' ,V_f(1),V_f(2),V_f(3))

		  %Use following observation site to find the satellite
GST=330;
L=38.8;
Long=-104.54; 
LST=GST+Long; 
Al=1.915;
x=abs((a_earth/sqrt(1-e_earth^2*sind(L)^2))+A1)*cosd(L); %ground station 
z=abs((a_earth*(1-e_earth^2))/sqrt(1-e_earth^2*sind(L)^2)+A1)*sind(L); 
Rsite_ijk=[x*cosd(LST);x*sind(LST);z];  %calculate station position

Rsquig=[sind(L)*cosd(LST),-sind(LST),cosd(L)*cosd(LST);...
       sind(L)*sind(LST),cosd(LST),cosd(L)*sind(LST);...
       -cosd(L),0,sind(L)]; %SEZ to IJK transformation matrix

rho_ijk=R_f-Rsite_ijk;          %Find new range-site information 
rho_sez=transp(Rsquig)*rho_ijk;
  
rho_f=norm(rho_sez);
El_f=asind(rho_sez(3)/rho_f);
Az_f=acosd(-rho_sez(1)/(rho_f*cosd(El_f)));

fprintf('\nCase %0.0f final tracking station data:\n',ii) 
fprintf('GST [Greenwich Sidereal Time]: %0.3f degrees\n',GST) 
fprintf('LST [Local Sidereal Time]: %0.3f degrees\n',LST) 
fprintf('Latitude: %0.7f degrees\n',L)
fprintf('Longitude: %0.7f degrees\n',Long)
fprintf('Altitude: %0.3f km\n\n',A1)
disp('Satellite sighting information from final tracking station:') 
fprintf('Range: %0.4f km\n',rho_f)
fprintf('Azimuth: %0.7f degrees\n',Az_f)
fprintf('Elevation: %0.7f degrees\n\n',El_f) 
fprintf('****************End of case %0.0f******************\n\n\n',ii) 
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
cc=cc+1;

%
%   ii=input('Enter 1 to check another case.')
%

%
% The following is the output of the above program when executed in MATLAB
%

% Case 1:

% Initial orbital elements: 
% Semi-major axis: 6782.8959206 km
% Eccentricity: 0.0006992 
% Inclination: 51.6434000 degrees 
% RAAN: 348.9165000 degrees
% Argument of Perigee: 60.4905000 degrees
% True Anomaly: 51.4287219 degrees.

% Final orbital elements:
% Semi-major axis: 6782.7673597 km
% Eccentricity: 0.0006803
% Inclination: 51.6434000 degrees
% RAAN: 279.1145668 degrees
% Argument of Perigee: 112.5385387 degrees 
% True Anomaly: 260.3363094 degrees.

% R = 1973.6967325 i + -6380.9111492 j 1185.2814888 k km 
% V = 4.3063080 i + 2.4258649 j + 5.8585568 k km/s

% Case 1 final tracking station data:
% GST [Greenwich Sidereal Time]: 330.000 degrees 
% LST [Local Sidereal Time]: 225.460 degrees 
% Latitude: 38.8000000 degrees
% Longitude: -104.5400000 degrees
% Altitude: 1.915 km

% Satellite sighting information from final tracking station: 
% Range: 6759.1753 km
% Azimuth: 100.0051817 degrees
% Elevation: -27.9053658 degrees

% ****************End of case 1 ******************

% Case 2:

% Initial orbital elements: 
% Semi-major axis: 7057.2145275 km
% Eccentricity: 0.0008529 
% Inclination: 97.8397000 degrees 
% RAAN: 143.5394000 degrees
% Argument of Perigee: 285.0815000 degrees
% True Anomaly: 75.0766435 degrees.

% Final orbital elements:
% Semi-major axis: 7057.1767444 km
% Eccentricity: 0.0008476
% Inclination: 97.8397000 degrees
% RAAN: 156.8845936 degrees
% Argument of Perigee: 240.7136734 degrees 
% True Anomaly: 80.2523044 degrees.

% R = -5278.9533357 i + 1594.2814587 j -4402.3061539 k km 
% V = -4.0457774 i + 2.5922824 j + 5.7800891 k km/s

% Case 2 final tracking station data:
% GST [Greenwich Sidereal Time]: 330.000 degrees 
% LST [Local Sidereal Time]: 225.460 degrees 
% Latitude: 38.8000000 degrees
% Longitude: -104.5400000 degrees
% Altitude: 1.915 km

% Satellite sighting information from final tracking station: 
% Range: 9992.1337 km
% Azimuth: 135.7941890 degrees
% Elevation: -45.5263765 degrees

% ****************End of case 2 ******************

% Case 3:

% Initial orbital elements:
% Semi-major axis: 42162.1276477 km
% Eccentricity: 0.0000844
% Inclination: 4.7447000 degrees
% RAAN: 58.5438000 degrees
% Argument of Perigee: 164.9477000 degrees 
% True Anomaly: 4.3934917 degrees.

% Final orbital elements:
% Semi-major axis: 42161.5939698 km
% Eccentricity: 0.0000717
% Inclination: 4.7447000 degrees
% RAAN: 58.3566075 degrees
% Argument of Perigee: 165.3201595 degrees 
% True Anomaly: 18.6560092 degrees.

% R = -19584.3216843 i + -37333.0178963 j -241.8079963 k km 
% V = 2.7143977 i + -1.4223812 j + -0.2537384 k km/s

% Case 3 final tracking station data:
% GST [Greenwich Sidereal Time]: 330.000 degrees 
% LST [Local Sidereal Time]: 225.460 degrees 
% Latitude: 38.8000000 degrees
% Longitude: -104.5400000 degrees
% Altitude: 1.915 km

% Satellite sighting information from final tracking station: 
% Range: 37658.0960 km
% Azimuth: 154.3381773 degrees
% Elevation: 41.4330916 degrees

% ****************End of case 3 ******************

% Case 4:

% Initial orbital elements: 
% Semi-major axis: 7367.2037665 km
% Eccentricity: 0.0032539 
% Inclination: 82.9339000 degrees 
% RAAN: 157.7326000 degrees
% Argument of Perigee: 43.7054000 degrees
% True Anomaly: 95.9770363 degrees.

% Final orbital elements:
% Semi-major axis: 7367.2014627 km 
% Eccentricity: 0.0032536
% Inclination: 82.9339000 degrees
% RAAN: 147.3777168 degrees
% Argument of Perigee: 4.8017781 degrees 
% True Anomaly: 171.4262554 degrees.

% R = 6179.2048379 i + -4026.1678163 j 482.5227165 k km 
% V = 0.8944084 i + 0.4960324 j + -7.2602810 k km/s

% Case 4 final tracking station data:
% GST [Greenwich Sidereal Time]: 330.000 degrees 
% LST [Local Sidereal Time]: 225.460 degrees 
% Latitude: 38.8000000 degrees
% Longitude: -104.5400000 degrees
% Altitude: 1.915 km

% Satellite sighting information from final tracking station: 
% Range: 10294.1272 km
% Azimuth: 79.6926994 degrees
% Elevation: -44.4638655 degrees

% ****************End of case 4*******************

% End of output.

%Analysis:
%
% Part 2: Comparing results with OMEGA dot and omega dot figures
% For case 1, RAAN is changing at -4.98 degrees per day, which agrees
% closely with the value on the graph for an altitude of 400km and i=51. 
% Argument of Perigee is changing at 3.72 degrees per day, the value on
% the chart is much higher, possibly due to chart perigee being too close.

% For case 2, RAAN is changing at 0.95 degrees per day, which agrees
% closely with the value on the graph for an altitude of 600km and i=98. 
% Argument of perigee is changing at -3.17 degrees per day, which agrees 
% closely with the value on the graph at i=98 and apogee above 500km.

% For case 3, RAAN is changing at -.013 degrees per day, which agrees
% with the trend towards less precession at higher altitudes. Since case 
% 3 has an altitude several times higher than the max charted value, it 
% should have a very low rate of precession, which is what is seen here. 
% Similarly, argument of perigee is changing very slowly, 0.027 degrees
% per day. The value on the chart is again much higher, likely due to the 
% fact that this orbit is not actually coming within 100km of the planet.

% For case 4, RAAN is changing at -0.739 degrees per day, which agrees
% closely with the value on the graph for an altitude of 1000km and i=83. 
% Argument of perigee is changing at -2.78 degrees per day, which agrees 
% closely with the value on the graph for i=83 and apogee above 1000km.

% The only satellite visible from Schriever AFB at the time specified is 
% METEOSAT-8. It is at 41.4 degrees elevation, above the 20 degrees
% necessary for tracking, and the satellite is at a range of 37658km, 
% which places in just inside the trackable range.

% All other studied satellites are below the horizon, and can not be 
% tracked from Schriever AFB.