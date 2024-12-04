function [p1,p2,air_speed,xfilm_bin_ave] = xfilm_calibration(xfilm_name,scan_name,start_min,len,order)
rhoa=1.12;
len_scan = len*60*100;
mili_sec=str2num(scan_name(end-6:end-5));
last_num=str2num(scan_name(end-4));
if last_num >= 5
    mili_sec = mili_sec+1;
else
    mili_sec = mili_sec;
end
num_of_mili_sec=100-mili_sec;
sec=str2num(scan_name(end-9:end-7))+1;
num_of_sec = 60 -sec;
min=str2num(scan_name(end-11:end-10))+1;
num_of_min = start_min - min;
offset =  num_of_min * 60 * 100 + num_of_sec * 100 + num_of_mili_sec;

pressure_probe=readtable(scan_name,'numheaderlines',6);


pressure_probe=pressure_probe(:,13:30);
pressure_probe=table2array(pressure_probe);
pressure_probe=pressure_probe(offset+1:offset+len_scan,:);
d1  = pressure_probe(:,17);
s1  = pressure_probe(:,18);
difference = d1-s1;
origin =  mean(difference(1:18000));
ref = 4*1.12/2;
pressure_offset = ref - origin;

for j = 1 : len*60 % nummber of sec
    mean_dynamic(j) = mean(pressure_probe((j-1)*100+1:j*100,17))
    mean_static (j) = mean(pressure_probe((j-1)*100+1:j*100,18))
end

air_speed  = sqrt(2*(mean_dynamic-mean_static)/rhoa);

%air_speed = sqrt(2*(pressure_probe(:,17)-pressure_probe(:,18))/rhoa);
%air_speed = switch_air(air_speed); % turn on the switch
% direction judgement:
% if the air_speed derived from Bernoulli Equation is imaginary, then it is
% upwind
% if the air speed derived from Bernoulli Equation is real
% then it is downwind
%for i = 1 : size(air_speed_com,1)
 %   if(imag(air_speed_com(i))~=0)
  %      air_speed(i)=-abs(air_speed_com(i));
   % else
    %    air_speed(i)=abs(air_speed_com(i));
    %end   


% now load hotfilm

fs=1000;
xfilm_1st=readtable(xfilm_name,'numheaderlines',1);
xfilm_2nd=table2array(xfilm_1st(:,2:3));
mili_sec=str2num(xfilm_name(end-9:end-7));
last_num=str2num(xfilm_name(end-6));
if last_num >= 5
    mili_sec = mili_sec+1;
else
    mili_sec = mili_sec;
end
num_of_mili_sec=fs-mili_sec;
sec=str2num(xfilm_name(end-12:end-11))+1;
num_of_sec = 60 -sec;
min=str2num(xfilm_name(end-14:end-13))+1;
num_of_min = start_min - min;
offset =  num_of_min * 60 * fs + num_of_sec * fs + num_of_mili_sec;
%
xfilm_length = len * 60 * fs;

xfilm_3rd=xfilm_2nd(offset+1:offset+xfilm_length,:);


 
bin=1000;
for j = 1:2
for i = 0:bin:size(xfilm_3rd,1)-bin
     
    xfilm_bin_ave(i/bin+1,j)=mean(xfilm_3rd(i+1:i+bin,j));
end
end

bin=100;

%for j = 1
%for i = 0:bin:size(air_speed)-bin
     
 %   air_speed_bin_ave(i/bin+1,j)=mean(air_speed(i+1:i+bin,j));
%end
%end
air_speed = effective_velocity(air_speed);

air_speed_bin_ave = air_speed;



clear effective_velocity
p1 = polyfit(xfilm_bin_ave(:,1),air_speed_bin_ave,order);
p2 = polyfit(xfilm_bin_ave(:,2),air_speed_bin_ave,order);

end

