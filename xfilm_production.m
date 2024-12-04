function [u,w,xfilm,un,ut,veff1,veff2] = xfilm_production(xfilm_name,start_min,len,p1,p2,f)

[x_offset,e_offset] = phase_correction(f);


fs = 1000;
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

offset= offset-x_offset;

xfilm_length = len * 60 * fs;


xfilm_raw = readtable(xfilm_name,'numheaderlines',1);
xfilm = table2array(xfilm_raw(offset+1:offset+xfilm_length,2:3));


[u,w,un,ut,veff1,veff2] = xfilm_interp(p1,p2,xfilm);



end