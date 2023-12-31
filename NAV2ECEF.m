function [rinexobs,drel]=NAV2ECI(NAV,rinexobs,sAt)
%---------------------------------------------
Findobs = find(rinexobs.data(:,rinexobs.col.PRN)==sAt);
t = rinexobs.data(Findobs,rinexobs.col.Temission);
global GM we
GM=3.986004418e14;
we=7.2921151467e-5;
%----------------------DATUM------------------------
% WGS-84(1984)
global A B E2 GAST xp yp c
c = 2.99792458e8;
A=6378137;	  B=6356752.3142;  E2=(A^2-B^2)/(A^2);
% Greenwich True Sidereal Time at Longitude 0.0°
GAST =0 ;
xp=0;yp=0;
% S=Rotation(2,-xp)*Rotation(1,-yp)*Rotation(3,-GAST);  %ECEF TO ECI
%------------------------------------------
PRN = NAV.data(:,NAV.col.PRN);
Find=find(PRN==sAt);
Mean_Anom=NAV.data(Find,NAV.col.Mean_Anom);
TOE=NAV.data(Find,NAV.col.TOE);
a=NAV.data(Find,NAV.col.a);
deltan=NAV.data(Find,NAV.col.deltan);
Right_Ascen_at_reference_time=NAV.data(Find,NAV.col.Right_Ascen_at_reference_time);
Rate_of_Right_Ascen=NAV.data(Find,NAV.col.Rate_of_Right_Ascen);
Argument_of_Perigee=NAV.data(Find,NAV.col.Argument_of_Perigee);
Cuc=NAV.data(Find,NAV.col.Cuc);
Cus=NAV.data(Find,NAV.col.Cus);
i0=NAV.data(Find,NAV.col.i0);
idot=NAV.data(Find,NAV.col.idot);
Cic=NAV.data(Find,NAV.col.Cic);
Cis=NAV.data(Find,NAV.col.Cis);
e=NAV.data(Find,NAV.col.e);
Crc=NAV.data(Find,NAV.col.Crc);
Crs=NAV.data(Find,NAV.col.Crs);
%% Start
for j=1:length(t)
    MIN=(TOE-t(j)).^2;
    nnn=find(MIN==min(MIN));
    Nnearest(1,j)=nnn(1);
end
j = Nnearest;
tk=t-TOE(j);
for kk = 1:length(tk)
    if tk(kk)>302400
        tk(kk) = tk(kk)-604800;
    elseif tk(kk)< -302400
        tk(kk) = tk(kk)+604800;
    end
end

n=(sqrt(GM./(a(j).^3)) + deltan(j));
Mk=Mean_Anom(j)+n.*tk;
omegak=Right_Ascen_at_reference_time(j)+(Rate_of_Right_Ascen(j)-we).*tk-we.*TOE(j);
omegak = wrapTo2Pi(omegak+GAST);
Mk = wrapTo2Pi(Mk);
E0=Mk;
for iii=1:10
    Ek=Mk+e(j).*sin(E0);
    E0=Ek;
end
drel = -2*((sqrt(GM*a(j)))/c^2).*e(j).*sin(Ek);
True_Anom=2*atan( tan(Ek/2).*sqrt(1+e(j))./sqrt(1-e(j)) );
% True_Anom=atan( (sqrt(1-e(j).^2).*sin(Ek))./( cos(Ek)-e(j) ) );
w=Argument_of_Perigee(j)+Cuc(j).*cos(2*(Argument_of_Perigee(j)+ True_Anom)) +Cus(j).*sin(2*(Argument_of_Perigee(j)+ True_Anom))+True_Anom;
inclination=i0(j)+idot(j).*tk+Cic(j).*cos(2*(Argument_of_Perigee(j)+ True_Anom)) +Cis(j).*sin(2*(Argument_of_Perigee(j)+ True_Anom));
rk=a(j).*(1-e(j).*cos(Ek))+Crc(j).*cos(2*(Argument_of_Perigee(j)+ True_Anom)) +Crs(j).*sin(2*(Argument_of_Perigee(j)+ True_Anom));
r_orb=[[rk]';zeros(1,max(size(Ek)));zeros(1,max(size(Ek)))];
r_ct=[];
for K=1:max(size(Ek))
    j=Nnearest(K);
    RRR=Rotation(3,-omegak(K))*Rotation(1,-inclination(K))*Rotation(3,-w(K)+True_Anom(K));
    RR=Rotation(3,-omegak(K))*Rotation(1,-inclination(K))*Rotation(3,-w(K));
    r_ct(:,K)=RR*r_orb(:,K);
    VxVyVz(:,K)=(n(K).*(a(j)^2).*(1/rk(K)).*( sqrt(1-e(j)^2)*cos(Ek(K))*RRR(:,2)-RRR(:,1)*sin(Ek(K)) ));
end
rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat) = [r_ct'];
rinexobs.VxVyVz = VxVyVz';

%--------------------------------------------------------------------
end
function [R]=Rotation(i,j)
%Rotation matrix
if i==3
    R = [cos(j) sin(j) 0 ; -sin(j) cos(j) 0 ; 0 0 1];
else if i==2
        R = [cos(j) 0 -sin(j) ; 0 1 0 ; sin(j) 0 cos(j)];
    else if i==1
            R = [1 0 0 ; 0 cos(j) sin(j) ; 0 -sin(j) cos(j)];
        end
    end
end
end
