function [u,w,un,ut,veff1,veff2] = xfilm_interp(p1,p2,xfilm);
 % input p1 and p2 as the parameter from ramp up and down
 % xfilm is the data acquisition time
 % output u and w from actual data acquisition

veff1= polyval(p1,xfilm(:,1));
veff2= polyval(p2,xfilm(:,2));

k1=0.3;
k2=0.3;

a=veff1.^2 - k1.^2 * veff2.^2;

b=veff2.^2 - k2.^2 * veff1.^2;

index1 = find(a<0);
index2 = find(b<0);

un = abs((sqrt((veff1.^2 - k1.^2 * veff2.^2) / (1 - k1.^2 * k2.^2)))) ;
ut = abs((sqrt((veff2.^2 - k2.^2 * veff1.^2) / (1 - k1.^2 * k2.^2)))) ;
%ut=1;

un (index1)=0;%-un(index1); 
ut (index2)=0;%-ut(index2);





u = 0.5 * sqrt(2.) * (ut + un);
w = 0.5 * sqrt(2.) * (ut - un);





end