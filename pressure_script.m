clear all
close all
%%
f=1.2;

[x_offset,e_offset] = phase_correction(f);
% in this case we have 40 min of average
start_min=35 % declare the starting minutes
scan_name='Scan - 2022-10-19T213244.203.CSV';
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
offset =  num_of_min * 60 * 100 + num_of_sec * 100 + num_of_mili_sec-e_offset;
%
pressure_probe=readtable(scan_name,'numheaderlines',6);
pressure_probe=pressure_probe(:,13:30);
pressure_probe=table2array(pressure_probe);
% you need to declare offset value here
length = 25*60*100; %but i will only extract the first 30 minutes

pressure_probe=pressure_probe(offset+1:offset+length,:);
%
color_c='ykgbmr';
subplot(3,1,1)
for i =1:3:16
 Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
   disp(i)
  plot((pressure_probe(:,i)),color_c((i+2)/3));
    hold on
end
legend('1','2','3','4','5','6')
title 'Elliot'
subplot(3,1,2)
for i =2:3:17
    plot((pressure_probe(:,i)),color_c((i+1)/3));
    hold on
end
title 'Dynamic'
legend('1','2','3','4','5','6')
subplot(3,1,3)
for i =3:3:18
    plot((pressure_probe(:,i)),color_c((i)/3))
   
    hold on
end
legend('1','2','3','4','5','6')
title 'Static'
% looks good

%%
bin=5;
for j = 1:18
for i = 0:5:size(pressure_probe,1)-5
     
    pressure_bin_ave(i/5+1,j)=mean(pressure_probe(i+1:i+5,j));
end
end 





%
load('wave.mat')
%%
close all
% this is just a linear plot
fs=20;
color_c='ykgbmr'
wavewire=bandpass(wave_signal,[1.1 1.3],20); % i will only take the first 30 min of the 40 min wave wire data. % the full length of wave
% wire data is reserved for hotfilm
phase_long=[9:18:351];

phase_long = [9:18:351];
for i = 1:3:16
[long,short,phase,short_wave_tot] = conditional_phase_ave_NRL(wavewire,pressure_bin_ave(:,i));
yyaxis left
press=short-mean(short); % demean
pressure_along_phase = press;

%pressure_along_phase(1)=press(20);
%pressure_along_phase(2:20)=press(1:19);
plot(phase_long,pressure_along_phase,'linewidth',2,'Color',color_c((i+2)/3),'marker','^');
hold on
pressure_array((i+2)/3,:)=short-mean(short);
end
yyaxis right
plot(phase_long,long,'linewidth',2,'Color','k','linestyle','--');

%%
addpath(genpath('/Users/peisen/Documents/NRL Ldeo'))
eta = long - mean(long);
eta = flipwind(eta);
%
for i = 1:6
   pressure_array(i,:) = flipwind(pressure_array(i,:));
end
%


%%
k_wave =double(dispersion(1.2*2*pi,.70));
z=[5:5:30]/100;
for i = 1:6
    [zeta(i,:)] = decay_coordinate(k_wave,z(i),eta);
end
clear p2_total
res = .1;

h2= [0:res:35]/100;
for j = 1: 20 % interpolate on each bin

    
    h1 = z;%zeta(:,j);
    disp(h1)
    p1 = pressure_array(:,j);
   
    p2 = spline(h1,p1,h2);
    p2_total(:,j)=p2;



end

% p2 is the pressure array after spline interpolation
%
close all
[X,Y]=meshgrid(phase_long,h2);



% next step is figure out how to block out the dots that are not in the range of those bent lines.
% remember, each phase bin has their own thresh hold , so identify those
% bins respectively

for i = 1:20
   minimum = .05;%zeta(1,i);
   maximum = .30;%;zeta(end,i);
   level2_height = Y(:,i);
   index =find(level2_height>maximum | level2_height<minimum);
   p2_total(index,i) =NaN;


end
%
figure
scatter3(X(:),Y(:),p2_total(:),1000,p2_total(:),'filled','marker','s')
view(2);
clim([-8 8])
colormap(jet)
set(gca,'fontsize',23,'fontname','Times New Roman');
ylabel 'Height above water surface (m)'
xlabel 'Phase'
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})
set(gcf,'position',[      560   298   594   649])
c=colorbar;
set(c,'location','southout')
hold on
plot(phase_long,zeros(size(eta)),'k','linewidth',2.0)
c.Label.String = 'wave-induced pressure (Pascal)';
colormap(bluewhitered(128));
xlim([9 351])
%%
f_wave= 1.2;

z =[.05:.05:.30];
k_wave =double(dispersion(f_wave*2*pi,.70));
alpha=zeros(1,6);
x=z;
for i = 1:20
    y = pressure_array(:,i);
    x =  z - eta(i);
    f=fit(x',y,'exp1');
    p0(i) = f.a;
    alpha(i)=-f.b/k_wave;
end
%%

%%


% form stress <deta/dx * p> , a first examination of the form stress on the z plane with only 6 dots
 k = double(dispersion(1.2*2*pi,.70));
 clear tau_form
 lambda  = 2 * pi./k;
 
dx = lambda/20;

deta = zeros(size(eta));

for i = 2:19
    deta(i) = eta(i+1)  - eta(i-1);
end

deta(1) = eta(2) -  eta(20);
deta(20) = eta(1) -eta(19);


deta_dx =  deta./dx;
% now we calculate the form stress at each zeta 


tau_form_discrete  = zeros(size(zeta));
%%
z_new =[0:.05:.30];
for i = 1:7
   
    relaxation_factor(i,:) = deta_dx .* exp(-k_wave.*z_new(i));
    plot(relaxation_factor(i,:)); hold on

end
legend('0','5 cm','10 cm','15 cm','20 cm','25 cm','30 cm')
set(gca,'fontsize',20,'fontname','times new roman')
title 'relaxation factor'
%%



for i = 1:7
    tau_form_discrete_1 =  relaxation_factor(i,:) .* (pressure_array_new(i,:));
    tau_form_discrete(i,:) = tau_form_discrete_1;
    tau_form(i) =  nanmean(tau_form_discrete(i,:));
end 

plot(tau_form,[0:.05:.30],'-*','markersize',20,'linewidth',2.0,'color','k')
set(gca,'fontsize',23,'fontname','times')
xlabel '\tau_f'
ylabel 'zeta'
%%

tau_form_surface = p0.* deta_dx/2;
%%
% now lets take a look at the discrete p * deta/dx
%scatter3(X(:),Y(:),tau_form_discrete(:),1000,tau_form_discrete(:),'filled','marker','s')


% plot the form stress profile and lets take a look!

%xlabel([ 'p $\frac{d\eta}{dx}$'],'interpreter','latex')
%set(gca,'fontsize',23,'fontname','Times New Roman')
%grid minor
%ylabel '\zeta'

%% quick check on the interpolation/extrapolation of pressure 
% toward surface
 pressure_array_new = [p0;pressure_array];






for i = 1:20

plot(pressure_array_new(:,i),[0, .05, .10,.15,.20,.25,.30],'-o','markersize',10,'linewidth',2.0);
hold on


end



set(gca,'fontsize',24,'fontname','Times New Roman')


%%
figure
plot(p0)
hold on
yyaxis right
plot(deta_dx)
%%
figure
subplot(2,1,1);
plot(eta,'--','Color','k','linewidth',2.0);

hold on
plot(deta_dx);
hold on
yyaxis right
plot(p0)
subplot(2,1,2)
plot(p0.*deta_dx);

%%
scatter3(X(:),Y(:),p2_total(:),1000,p2_total(:),'filled','marker','s')
view(2);
clim([-1.5 1.5])
colormap(jet)
set(gca,'fontsize',23,'fontname','Times New Roman');
ylabel '\zeta (m)'
xlabel 'Phase'
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})
set(gcf,'position',[   560   242   698   705])
c=colorbar;
set(c,'location','southout')
hold on
plot(phase_long,zeros(size(eta)),'k','linewidth',2.0)
c.Label.String = '$$\tilde{p}$$ (Pa)';
c.Label.Interpreter='latex';

colormap(bluewhitered(128));
clim([-1 1])
xlim([4 356])

set(gcf,"Position", [639   242   619   705])

%%

figure




subplot(2,1,2)
yyaxis right
plot(phase_long,eta,'--','Color','k','linewidth',2.0);
hold on
yticks([])
yyaxis left
plot(phase_long,p0.*deta_dx,'Color','b','linewidth',2.0);

plot(phase_long,zeros(size(p0.*deta_dx)),'linestyle',':','Color','b','linewidth',2.0);

set(gcf,"Position", [639   242   619   705])
set(gca,'fontsize',23,'fontname','Times New Roman')
grid minor
ylabel('$\it{P_0 \frac{d\eta}{dx}} (N m^{-2})$','Interpreter','Latex')
set(gcf,'position',[198   439   619   423]);
xticks([9 180 351])
xticklabels({'0','\pi','2\pi'})

xlim([4 356])
%%
l=legend('$\it{P_0 \frac{d\eta}{dx}}$','zero reference','wave shape')
set(l,'Interpreter','latex','Orientation','horizontal','location','Southout')
%%











