clear all
close all
%
addpath(genpath('/Users/peisen/Documents/NRL Ldeo'))

% addpath(genpath('/Users/peisen/Documents/workalholic'))
addpath(genpath('/Users/peisen/Documents/NRL Ldeo'))
xfilm_name = 'u6_stream_20221019_202411.028496.csv'; % calibration hotfilm file
scan_name = 'Scan - 2022-10-19T202405.708.csv'
start_min=31;
len = 13;
order =2;
[p1,p2,air_speed,xfilm_bin_ave,xfilm,d,s] = xfilm_calibration_right(xfilm_name,scan_name,start_min,len,order);

%
air_speed=(air_speed);
plot(real(air_speed));


%%
figure

subplot(3,1,1)

plot(linspace(0,13*60,length(d)),d,'Linewidth',2.0);
ylabel 'Pa'
subplot(3,1,2)
plot(linspace(0,13*60,length(d)),s,'Linewidth',2.0);
ylabel 'Pa'

subplot(3,1,3)
plot(linspace(0,13*60,length(d)),air_speed,'Linewidth',2.0);
ylabel 'm s^{-1}'

for i = 1:3
   subplot(3,1,i)
    set(gca,'fontsize',23,'fontname','Times New Roman')
    grid minor
    xlim([100 700])
    if i <3
        ylim([0 300])
    else
        ylim([4 12])
        xlabel 'Seconds'
    end
    xticks([100 400 700])
    xticklabels({'0','300','600'})
    

    


end




%%
title 'calibration air speed'

xfilm_name = 'u6_stream_20221019_213249.663609.csv'; % hotfilm production file
load('/Users/peisentan/Documents/NRL Ldeo/1.2 Hz 4cm/15 Hz/wave.mat')
start_min =  35; % start min of the data production
len  = 25;    %  running length of the data production


%%
[u,w,xfilm,un,ut,veff1,veff2] = xfilm_production(xfilm_name,start_min,len,p1,p2,1.2);
u_seg = u(end-3*60*1000+1:end);

w_seg = w(end-3*60*1000+1:end);
[u_seg,w_seg] = rotate(u_seg,w_seg);

% here to add some plots
xfilm_seg = xfilm(end-3*60*1000+1:end,:);

t= linspace(0,180,length(u_seg));                 
subplot(2,2,1)
plot(t,xfilm_seg(:,1),'LineWidth',1,'Color','k');
ylabel 'Volt'
subplot(2,2,2)

plot(t,xfilm_seg(:,2),'LineWidth',1,'Color','k');
ylim([1.5 1.9])
ylabel 'Volt'

subplot(2,2,3)

plot(t,u_seg,'LineWidth',1,'Color','k');
ylim([2 12])
ylabel 'm s^{-1}'

subplot(2,2,4)

plot(t,w_seg,'LineWidth',1,'Color','k')
ylabel 'm s^{-1}'
ylim([-4 4])
hold on

for i = 1:4
    subplot(2,2,i)
    xlim([0 8])
    set(gca,'fontsize',20,'fontname','Times New Roman')
    grid minor
    xlabel 'Seconds'      
   
end




%%
u(u>100) =NaN;
w(w>100 | w<-100)=NaN;
w=fillmissing(w,'linear');
u=fillmissing(u,'linear');
plot(u)
yyaxis right
plot(w)
L = 5*60*1000;
w=-w;
%%


for i = 1:5
  eval(['u',num2str(i),'= u((',num2str(i),'-1)*L+1:',num2str(i),'*L);']);
   eval(['w',num2str(i),'= w((',num2str(i),'-1)*L+1:',num2str(i),'*L);']);
  eval(['[u',num2str(i),',w',num2str(i),']=rotate(u',num2str(i),',w',num2str(i),');']);
   %eval(['u((',num2str(i),'-1)*L+1:',num2str(i),'*L)=u',num2str(i),';']);
   %eval(['w((',num2str(i),'-1)*L+1:',num2str(i),'*L)=w',num2str(i),';']);

end 
%[u,w]=rotate(u,w);
%%
mean_u(1)=mean(u1);
mean_u(2)=mean(u2);
mean_u(3)=mean(u3);
mean_u(4)=mean(u4);
mean_u(5)=mean(u5);
color = 'kgbmr'
for i = 1:5

eval(['w_prime = w', num2str(i),' - mean(w', num2str(i),');']);
eval(['u_prime = u', num2str(i),' - mean(u', num2str(i),');']);
tau_EC(i) = (-mean(u_prime(1+60*1000:end).*w_prime(1+60*1000:end)));
ustar(i) = sqrt(tau_EC(i));

uw = spectf(u_prime(1+60*1000:end), w_prime(1+60*1000:end),.001);
co_spectrum = -real(uw(:,4));
tau_co_spectrum(i) = trapz(uw(:,1),co_spectrum);
semilogx(uw(:,1)*height(i)/mean_u(i),co_spectrum.*uw(:,1)/ustar(i)/ustar(i),'linewidth',2.0,'color',color(i));

hold on

%legend('0.45 m','0.35 m','0.25 m','0.15 m','0.10 m')%,'0.05 m')
end
set(gca,'fontsize',23,'fontname','Times New Roman')
%%
xlabel 'f (Hz)'
ylabel '\it{Co-specm(uw) *f }'
grid minor
set(gcf,'position',[   756   419   817   528]);
%%
Suu = spectf(u5,.001);
loglog(Suu(:,1),Suu(:,2))
hold on

Sww = spectf(w5,.001);
loglog(Sww(:,1),Sww(:,2));
legend('S_{uu}','S_{ww}')
xlabel 'Freq.'
ylabel 'Spectrum Density (m^2/s)'
set(gca,'fontsize',23,'fontname','Times New Roman')

%%
f=polyfit([ .25 .15 .10],tau_co_spectrum(3:end),1);
tau_profile = polyval(f,[0:.01:.50]);
figure
plot(tau_profile,[0:.01:.50],'--')
hold on
scatter(tau_EC,[.45 .35 .25 .15 .10],100,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel 'stress (N/m^2)'
ylabel 'height (m)'
set(gca,'fontsize',23,'fontname','Times New Roman')
sqrt(tau_profile(1))
%%
figure
plot(u);
yyaxis right
plot(w)
title 'U and W in data production'
figure
plot(u);
yyaxis right
plot(w)
title 'U and W in data production'
%%
tau = -fliplr(tau)
f=polyfit([ .10 .15 .25 .35 .45],tau,1);
tau_profile = polyval(f,[0:.01:.50]);
figure
plot(tau_profile,[0:.01:.50],'--')
hold on
scatter(tau,[.10 .15 .25 .35 .45]);
%%

mean_u(1)=mean(u1);
mean_u(2)=mean(u2);
mean_u(3)=mean(u3);
mean_u(4)=mean(u4);
mean_u(5)=mean(u5);
%%
clear wave*
f_wave =1.2;
load('/Users/peisentan/Documents/NRL Ldeo/1.2 Hz 4cm/15 Hz/wave.mat')

num_heights = 5;
sample_len = 5*60*1000;
first_one_min = 60*1000;
wave_freq = 1.2; 
rhoa=1.12;
wavewire=wave_signal;
phase_long=[9:18:351];
clear u_along_phase
clear w_along_phase
wavewire=bandpass(wavewire,[1.1 1.3],20);
%% w self check and make sure no drifting

bin =50;

close all
subplot(3,1,1)
plot(linspace(0,5,length(w5)),u5,'linewidth',1,'color','b');
xlim([3.48 3.58])
ylim([2 12])
yyaxis right

plot(linspace(0,5,length(wavewire(end-60*20*5:end))),wavewire(end-60*20*5:end),'linewidth',1.0,'linestyle','--','color','k');

subplot(3,1,2)
plot(linspace(0,5,length(w5)),fliplr(w5),'linewidth',1,'color','r');
xlim([3.48 3.58])
ylim([-3 3])
yyaxis right
plot(linspace(0,5,length(wavewire(end-60*20*5:end))),wavewire(end-60*20*5:end),'linewidth',1.0,'linestyle','--','color','k')
%%


%% 

bin =50;

for i = 1: num_heights
    eval(['u_seg = u',num2str(i),';']);
    eval(['w_seg = w',num2str(i),';']);
    seg_wave1 = wavewire((i-1) * 20*5*60+1:i * 20*5*60);
    % get rid of first 10 and last 10 seconds


    % U_TILDE
    [long,short,phase,short_collection_tot,short_collection] = conditional_phase_ave_NRL_hotfilm(seg_wave1,u_seg,50);

    vel=short-mean(short); % i remove the mean u at each height,

    %vel = flipwind(short);
    % and i should try removing mean u at each zeta
    plot(vel,'linewidth',2,'Color','r','marker','^','markersize',12);
    
    u_tilde = vel;
    u_along_phase(num_heights+1-i,:) = u_tilde;

    % W_TILDE
    [long,short,phase,short_collection_tot,short_collection] = conditional_phase_ave_NRL_hotfilm(seg_wave1,w_seg,50);
   
    
    %short = fliplr(short);
    vel=short-mean(short); % i remove the mean u at each height,
    %vel = flipwind(short);
    % and i should try removing mean u at each zeta
    vel_along_phase=vel;
    plot(vel,'linewidth',2,'Color','b','marker','^','markersize',12);
    w_tilde = vel_along_phase;
    tau_wave_coherent(num_heights+1-i) = -mean((u_tilde.*w_tilde))
    hold on
    w_along_phase(num_heights+1-i,:) = vel_along_phase;

    [tke,tau_turb,u_prime,w_prime] =  turbulence_component(u_seg,w_seg,seg_wave1, 50);

    tke_along_phase(num_heights+1-i,:) = tke;
    tau_turb_along_phase(num_heights+1-i,:) = tau_turb;
    %stress1 (num_heights+1-i) = -mean(u_prime.*)

end

yyaxis right
plot(long)

for i = 1:num_heights
    tau_turb_ave_in_height(i) = mean(tau_turb_along_phase(i,:));
end 
%%
plot(tau_wave_coherent,z.*k_wave,'linewidth',2.0,'color','B');hold on
plot(tau_turb_ave_in_height,z.*k_wave,'linewidth',2.0,'color','r')

ylabel '\it kz'
set(gca,'xscale','linear','fontsize',23,'fontname','Times New Roman')
grid minor

xlabel 'N/m^2'
% this is a warning to flip wind!!
%% calculate the residueal

close all
figure(1)
scatter(-tau_turb_ave_in_height,[.10 .15 .25 .35 .45],100,'markeredgecolor','r','markerfacecolor','r')
set(gca,'fontsize',23,'fontname','Times New Roman')
xlabel 'turbulence stress'
ylabel 'Heights (m)'

figure(2)
scatter(-tau_wave_coherent,[.10 .15 .25 .35 .45],100,'markeredgecolor','r','markerfacecolor','r')
set(gca,'fontsize',23,'fontname','Times New Roman')
xlabel 'wave coherent stress'
ylabel 'Heights (m)'

%%




close all
eta = long - mean(long);
eta = flipwind(eta)
for i = 1:5
     w_along_phase(i,:) = flipwind(w_along_phase(i,:));
     u_along_phase(i,:) = flipwind(u_along_phase(i,:));
     tke_along_phase(i,:) = flipwind(tke_along_phase(i,:));
     tau_turb_along_phase(i,:) = flipwind(tau_turb_along_phase(i,:));
 
end 
%%
z =[ .10, .15,.25,.35,.45];
k = double(dispersion(1.2*2*pi,.70));
% establishing the new coordinate
clear zeta
for i = 1: num_heights
    [zeta(i,:)] = decay_coordinate(k,z(i),eta);
end



% start interpolation

clear p2_total
res = .01;

h2= [0:res:45]/100;

for j = 1: 20 % interpolate on each bin

    
    h1 = z;%eta(:,j);
    disp(h1)
    p1 = u_along_phase(:,j);
   
    p2 = spline(h1,p1,h2);
    p2_total(:,j)=p2;



end
%
%

% next step is figure out how to block out the dots that are not in the range of those bent lines.
% remember, each phase bin has their own thresh hold , so identify those
% bins respectively
[X,Y]=meshgrid(phase_long,h2);

for i = 1:20
   minimum = .10;%zeta(1,i);
   maximum = .45;%zeta(end,i);
   level2_height = Y(:,i);
   index =find(level2_height>maximum | level2_height<minimum);
   p2_total(index,i) =NaN;


end
%


p2_total_demean = nan(size(p2_total));

for i= 1:size(p2_total,1)

    mean_on_zeta(i) = nanmean(p2_total(i,:));
    p2_total_demean(i,:) = p2_total(i,:);% - mean_on_zeta(i);

end


close all
scatter3(X(:),Y(:),p2_total_demean(:),1000,p2_total_demean(:),'filled','marker','s')
clim([-.1 .1])

view(2);
colormap(bluewhitered(128))
set(gca,'fontsize',23,'fontname','Times New Roman');
ylabel 'Height above water surface (m)'
xlabel 'Phase'
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})
set(gcf,'position',[   560   242   698   705])
c=colorbar;
set(c,'location','southout')
hold on
xlim([0 360])
set(gcf,'position',[ 560   242   610   705])
clim auto

plot(phase_long,zeros(size(eta)),'k','linewidth',2.0)


%% tke visualization
close all
z =[ .05, .10, .15,.25,.35,.45];
k = double(dispersion(1.2*2*pi,.70));
% establishing the new coordinate
clear zeta
for i = 1: num_heights
    [zeta(i,:)] = decay_coordinate(k,z(i),eta);
end



% start interpolation

clear p2_total
res = .01;

h2= [0:res:45]/100;

for j = 1: 20 % interpolate on each bin

    
    h1 = zeta(:,j);
    disp(h1)
    p1 = tke_along_phase(:,j);
   
    p2 = spline(h1,p1,h2);
    p2_total(:,j)=p2;



end
%
%

% next step is figure out how to block out the dots that are not in the range of those bent lines.
% remember, each phase bin has their own thresh hold , so identify those
% bins respectively
[X,Y]=meshgrid(phase_long,h2);

for i = 1:20
   minimum = zeta(1,i);
   maximum = zeta(end,i);
   level2_height = Y(:,i);
   index =find(level2_height>maximum | level2_height<minimum);
   p2_total(index,i) =NaN;


end
%


p2_total_demean = nan(size(p2_total));

for i= 1:size(p2_total,1)

    mean_on_zeta(i) = nanmean(p2_total(i,:));
    p2_total_demean(i,:) = p2_total(i,:) ;%- mean_on_zeta(i);

end


close all
scatter3(X(:),Y(:),p2_total_demean(:),1000,p2_total_demean(:),'filled','marker','s')
clim([-.1 .1])

view(2);
colormap(bluewhitered(128))
set(gca,'fontsize',23,'fontname','Times New Roman');
ylabel 'Height above water surface (m)'
xlabel 'Phase'
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})
set(gcf,'position',[   560   242   698   705])
c=colorbar;
set(c,'location','southout')
hold on
xlim([0 360])
set(gcf,'position',[ 560   242   610   705])
clim auto

plot(phase_long,zeros(size(eta)),'k','linewidth',2.0)



















%%
velocity_u_along_phase =u_along_phase;
velocity_w_along_phase = w_along_phase;
z =[.10, .15,.25,.35,.45];
k = double(dispersion(1.2*2*pi,.70));
% establishing the new coordinate


clear zeta
for i = 1: num_heights
    [zeta(i,:)] = decay_coordinate(k,z(i),eta);
end



% start interpolation

clear p2_total
res = .01;

h2= [0:res:45]/100;

for j = 1: 20 % interpolate on each bin

    
    h1 = zeta(:,j);
    disp(h1)
    p1 = u_along_phase(:,j);
   
    p2 = spline(h1,p1,h2);
    p2_total(:,j)=p2;



end
%
%

% next step is figure out how to block out the dots that are not in the range of those bent lines.
% remember, each phase bin has their own thresh hold , so identify those
% bins respectively
[X,Y]=meshgrid(phase_long,h2);

for i = 1:20
   minimum = zeta(1,i);
   maximum = zeta(end,i);
   level2_height = Y(:,i);
   index =find(level2_height>maximum | level2_height<minimum);
   p2_total(index,i) =NaN;


end
%


p2_total_demean = nan(size(p2_total));

for i= 1:size(p2_total,1)

    mean_on_zeta(i) = nanmean(p2_total(i,:));
    p2_total_demean(i,:) = p2_total(i,:);% - mean_on_zeta(i);

end


close all
subplot(1,2,1)
scatter3(X(:),Y(:),p2_total_demean(:),1000,p2_total_demean(:),'filled','marker','s')
clim([-1 1])


view(2);
colormap(bluewhitered(128))
set(gca,'fontsize',23,'fontname','Times New Roman');
ylabel '\zeta (m)'
xlabel 'Phase'
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})
set(gcf,'position',[   560   242   698   705])
c=colorbar;
set(c,'location','southout')
hold on
xlim([0 360])
ylim([0 .50])

set(gcf,'position',[ 560   242   610   705])

plot(phase_long,zeros(size(eta)),'k','linewidth',2.0)
c.Label.String = '$$\tilde{u}  (ms^{-1})$$ ';
c.Label.Interpreter='Latex';
%
subplot(1,2,2)

clear p2_total
res = .01;

h2= [0:res:45]/100;

for j = 1: 20 % interpolate on each bin

    
    h1 = zeta(:,j);
    disp(h1)
    p1 = w_along_phase(:,j);
   
    p2 = spline(h1,p1,h2);
    p2_total(:,j)=p2;



end
%
%

% next step is figure out how to block out the dots that are not in the range of those bent lines.
% remember, each phase bin has their own thresh hold , so identify those
% bins respectively
[X,Y]=meshgrid(phase_long,h2);

for i = 1:20
   minimum = zeta(1,i);
   maximum = zeta(end,i);
   level2_height = Y(:,i);
   index =find(level2_height>maximum | level2_height<minimum);
   p2_total(index,i) =NaN;


end
%


p2_total_demean = nan(size(p2_total));

for i= 1:size(p2_total,1)

    mean_on_zeta(i) = nanmean(p2_total(i,:));
    p2_total_demean(i,:) = p2_total(i,:);% - mean_on_zeta(i);

end



scatter3(X(:),Y(:),p2_total_demean(:),1000,p2_total_demean(:),'filled','marker','s')
clim([-1 1])

ylim([0 .50])
view(2);
colormap(bluewhitered(128))
set(gca,'fontsize',23,'fontname','Times New Roman');
ylabel 'Height above water surface (m)'
xlabel 'Phase'
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})
set(gcf,'position',[   560   242   698   705])
c=colorbar;
set(c,'location','southout')
hold on
xlim([0 360])
set(gcf,'position',[ 560   242   610   705])
grid minor
plot(phase_long,zeros(size(eta)),'k','linewidth',2.0)
c.Label.String = '$$\tilde{w}  (ms^{-1})$$ ';
c.Label.Interpreter='Latex';
set(gcf,'position',[         177         242        1460         705]);
set(gca,'yticklabel',[])
grid minor
ylabel('')


%% tke visualization
close all
eta = long - mean(long);
z =[ .05, .10, .15,.25,.35,.45];
k = double(dispersion(1.2*2*pi,.70));
% establishing the new coordinate
clear zeta
for i = 1: num_heights
    [zeta(i,:)] = decay_coordinate(k,z(i),eta);
end



% start interpolation

clear p2_total
res = .01;

h2= [0:res:45]/100;

for j = 1: 20 % interpolate on each bin

    
    h1 = zeta(:,j);
    disp(h1)
    p1 = u_along_phase(:,j);
   
    p2 = spline(h1,p1,h2);
    p2_total(:,j)=p2;



end
%
%

% next step is figure out how to block out the dots that are not in the range of those bent lines.
% remember, each phase bin has their own thresh hold , so identify those
% bins respectively
[X,Y]=meshgrid(phase_long,h2);

for i = 1:20
   minimum = zeta(1,i);
   maximum = zeta(end,i);
   level2_height = Y(:,i);
   index =find(level2_height>maximum | level2_height<minimum);
   p2_total(index,i) =NaN;


end
%


p2_total_demean = nan(size(p2_total));

for i= 1:size(p2_total,1)

    mean_on_zeta(i) = nanmean(p2_total(i,:));
    p2_total_demean(i,:) = p2_total(i,:) ;%- mean_on_zeta(i);

end


close all
scatter3(X(:),Y(:),p2_total_demean(:),1000,p2_total_demean(:),'filled','marker','s')
clim([-.1 .1])

view(2);
colormap(bluewhitered(128))
set(gca,'fontsize',23,'fontname','Times New Roman');
ylabel 'Height above water surface (m)'
xlabel 'Phase'
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})
set(gcf,'position',[   560   242   698   705])
c=colorbar;
set(c,'location','southout')
hold on
xlim([0 360])
set(gcf,'position',[ 560   242   610   705])
clim auto

plot(phase_long,zeros(size(eta)),'k','linewidth',2.0)
