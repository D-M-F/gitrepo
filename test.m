%%
x = 0.01:0.001:.3;
y = exp(0.01./x).*(1-0.00001./(x.^3));
plot(x,y)
ylim([0 2])
%% 

%%
t = [5770 288];
h = 6.63*10^(-34);
c = 3*10^8;
k = 1.3807*10^(-23);
lambdai = 1:5000;
lambda = lambdai*10^(-9);

figure(1)
subplot(2,1,1)
i = 1;
T = t(i);
b = (2*h*c^2)./((lambda.^5).*(exp((h*c)./(lambda*k*T))-1))*10^(-9);
plot(lambdai,b)
xlabel('Wavelength (nm)')
ylabel('Intensity (W*m^{-2}*nm^{-1}*sr^{-1})')
title('Sun(T_{S}=5770K)')

subplot(2,1,2)
lambdai = 1:50000;
lambda = lambdai*10^(-9);
i = 2;
T = t(i);
b = (2*h*c^2)./((lambda.^5).*(exp((h*c)./(lambda*k*T))-1))*10^(-9);
plot(lambdai,b)
xlabel('Wavelength (nm)')
ylabel('Intensity (W*m^{-2}*nm^{-1}*sr^{-1})')
title('Earth(T_{E}=288K)')
% % tspan = 1:200;
% % p0 = [0.3,0.6,0.7,0.2,0.5];
% % e = 0.05;
% % traj = CML(tspan,p0,e);
% % 
% % for i = 1:5
% %     subplot(5,1,i);
% %     x = traj(:,i);
% %     plot(tspan',x)
% % end
%%
vlst = zeros(1,11);
for i = 0:10
    B = .1;
    c = i/10;
    traj = CHenonM(1:10000,rand(1,4),c,B);
    tmp = xcorrelation(DTS(re_traj(traj(1001:10000,:),10,1),9,205),'hvg','uni','xco');
    vlst(i+1) = tmp(1,2);
end
plot(0:0.1:1,vlst)
xlabel('c')
ylabel('estimated synchronization')
%%
% x=-10:1:10;
% y=-10:1:10;
% for xx=x
%     for yy=y
%         z(xx+11,yy+11)=100-xx^2-yy^2;
%     end
% end
% figure(1)
% pcolor(x,y,z);
% shading interp;
% colorbar;colormap(jet);
% xlabel('X');

% ylabel('Y');
% 
%%
% tspan = 1:2500;
% pulsespan = 1000:1600;
% c = 0.5;
% B = 0.8; %B>=0.4难看出
% p0 = rand(1,4);
% traj = CHM_pulse(tspan,pulsespan,p0,c,B);
% tmp = xcorrelation(DTS(re_traj(traj,10,1),9,309),'hvg','uni','syn'); %数据从309开始到2500-9-308=2183
% plot(309:2183,tmp)
tspan = 1:2500;
pulsespan = 1000:1600;
c = 0.5;
p0 = rand(1,4);
for bb = 0:11
    subplot(3,4,bb+1);
    B = bb/10;
    traj = CHM_pulse(tspan,pulsespan,p0,c,B);
    tmp = xcorrelation(DTS(re_traj(traj,10,1),9,309),'hvg','uni','syn'); %数据从309开始到2500-9-308=2183
    plot(309:2183,tmp)
    xlabel('c')
    ylabel('estimated synchronization')
    title(sprintf('B=%.1f',B))
end
%%
addpath('D:\学习\本研\statistics\本研材料\AQI');
M = readmatrix('D:\学习\本研\statistics\本研材料\AQI\北京.xlsx','Sheet','空气情况','Range','D32:I1856');
% 1825*6, 1825=5*(90+91+92+92)
m = [0 90 181 273 365];
wlst = zeros(1,20);
for i=1:20
    year = floor((i-1)/4);
    season = mod(i,4);
    if season == 0
        season = 4;
    end
    st = 365*year+m(season)+1;
    ed = 365*year+m(season+1);
    ts = M(st:ed,:);
%     wlst(i) = avg_edge_overlap(VG(ts,'hvg','bi'))
    tmp = Interlayer_mutual_info(VG(ts,'hvg','bi'));
    wlst(i) = tmp(1,2);
end
plot(1:20,wlst)

% windowSize = 15; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% plot(1:1825,filter(b,a,M(:,1)))
% datetime(,,):datetime(,,)
%%
% d,refd为列
tmp = size(refd);
len = tmp(1);
aux = 0;
itp_index = [];

for i = 1:len
    if strcmp(d(i-aux),refd(i)) == 0
        itp_index(end+1) = i;
        aux = aux+1;
    end
end

t_complete = 1:len;
t_origin = setdiff(t_complete,itp_index);
%%
tts = ["12-2","3-5","6-8","9-11"];
for i = 1:4
    jpd_mat0 = zeros(101,101);
    for k = 0:4:16
        j = k+i;
        vg = VG([shanghai0(91*(j-1)+1:91*j,5),shanghai0(91*(j-1)+1:91*j,6)],'hvg','bi');
        jpd_mat = plot_jpd(vg);
        [~,M] = size(jpd_mat);
        jpd_mat0(1:M,1:M) = jpd_mat0(1:M,1:M)+jpd_mat;
    end
    jpd_mat0 = jpd_mat0/5;
    M = 21;
    subplot(2,2,i);
    pcolor(0.5:1:M-0.5,0.5:1:M-0.5,jpd_mat0(1:21,1:21));
%     shading interp;
    colorbar;colormap(jet);
    caxis([0 0.1]);
    xlabel('k_{NO_2}');
    ylabel('k_{O_3}');
    title(tts(i));
end
%%
fh = figure('Visible','off');
set(gca,'position',[0.1,0.1,0.8,0.8])
set(gcf,'unit','normalized','position',[0.1,0.1,0.9,0.75])
plot(1:1917,shanghai0(:,1))
upper = 500;
line([31 31],[0 upper],'linestyle','--','Color','r')
line([396 396],[0 upper],'linestyle','--','Color','r')
line([761 761],[0 upper],'linestyle','--','Color','r')
line([1127 1127],[0 upper],'linestyle','--','Color','r')
line([1492 1492],[0 upper],'linestyle','--','Color','r')
line([1857 1857],[0 upper],'linestyle','--','Color','r')
legend('上海');
xlabel('day')
xlim([-10 1928])
title('上海PM2.5(2013/12/2-2019/3/2)');
set(fh,'Visible','on')
print2eps shanghaiPM25;
%%
% !!!imi
jmin = 0;
jmax = 80;
jstep = 0.01;
len = (jmax-jmin)+1;
nc = zeros(1000*len,1);
lc = zeros(1000*len,1);
steps = ["delay=0,AEO","delay=0,Pearson","delay=1,AEO","delay=1,Pearson","delay=2,AEO","delay=2,Pearson"];
figure(1)
step = 1;

for delay = 0:2
    for j = jmin:1:jmax 
        tsn = step*3000;
        ts = zeros(tsn,2);
        ts(1,1) = 0.3;
        ts(1,2) = 0.5;
        for i = 2:tsn
            ts(i,1) = ts(i-1,1)*(3.78-3.78*ts(i-1,1));
            ts(i,2) = ts(i-1,2)*(3.77-3.77*ts(i-1,2)-j*jstep*ts(i-1,1));
        end
        ts = ts(1:step:end,:);
        ts = ts(911-delay:2000,:);
        for i = 1:1000
            nc(i+(j-jmin)*1000) = avg_edge_overlap(VG([ts(i:i+90,1),ts(i+delay:i+90+delay,2)],'hvg','bi'));
%                 tmp = Interlayer_mutual_info(VG([ts(i:i+90,1),ts(i+delay:i+90+delay,2)],'hvg','bi'));
%                 nc(i+(j-jmin)*1000) = tmp(1,2);
            lc(i+(j-jmin)*1000) = corr(ts(i:i+90,1),ts(i+delay:i+90+delay,2),'Type','Pearson');
        end
    end
%         plot(1:len,nc)
%         hold on
    plot(1:len*1000,nc,1:len*1000,lc)
    hold on
end
%     nc = reshape(nc,[1000 len]);
%     nc = mean(nc,1);

legend(steps)
%%
% scales = [1 2 3 4 5 6 8 10 12 15 16 20 24 25 30 40 48 50 60 80 100 120 160 200];
scales = 1:40;
tsn = 10600;
ts = zeros(tsn,2);
slst = zeros(1,40);
rep = 5;
for j = 1:rep
    ts(1:10,1) = rand(10,1);
    ts(1:10,2) = rand(10,1);
    for i = 11:tsn
        ts(i,1) = ts(i-1,1)*(3.78-3.78*ts(i-1,1));%-0.01*ts(i-1,2));
        ts(i,2) = ts(i-1,2)*(3.77-3.77*ts(i-1,2)-0.1*ts(i-10,1));
    end
    tsx = ts(1001:end,1);
    tsy = ts(1001:end,2);
    for n = 1:40
        scale = scales(n);
        num = floor((tsn-1000)/scale);
        tmpx = tsx(1:num*scale);
        tmpy = tsy(1:num*scale);
        tmpx = reshape(tmpx,[scale num]);
        tmpy = reshape(tmpy,[scale num]);
        tsx1 = tmpx(1,:);
        tsy1 = tmpy(1,:);
%         slst(n) = slst(n)+avg_edge_overlap(VG([tsx1',tsy1'],'hvg','bi'));
        vg = VG([tsx1',tsy1'],'hvg','bi');
        tmp1 = Interlayer_mutual_info(vg);
        tmp2 = joint_entropy(vg);
        imi = tmp1./tmp2;
        slst(n) = slst(n) + imi(1,2);
    end
end
slst = slst/rep;
scatter(scales,slst)
xlabel('time scale')
ylabel('AEO')
%%
% for i = 1:5000
%     ts(1,1) = randn();
%     ts(1:tau,2) = randn(tau,1);
%     nc(i) = avg_edge_overlap(VG([ts(i:i+90,1),ts(i:i+90,2)],'hvg','bi'));
%     lc(i) = corr(ts(i:i+90,1),ts(i:i+90,2),'Type','Pearson');
% end
% plot(1:5000,nc,1:5000,lc)
% pearson
% spearman
% hvg vg
% plot(1:5000,nc/100,1:5000,lc/100)
% legend(["nc","lc"])
%%
%bj2017
%temp
witype = 15;
pn = 6;
ptype = ["PM2.5","PM10","SO2","CO","NO2","O3"];
t = tiledlayout(2,3);
ini_c = beijingwi(:,witype);
cmin = min(ini_c);
cmax = max(ini_c);
for i = 0:4
    nexttile
    z = beijing0(46+365*i:410+365*i,pn);
    c = beijingwi(1+365*i:365+365*i,witype);
%     c = nanjing(1+365*i:365+365*i,45);
    x = beijingt(1+365*i:365+365*i);
    y = beijingr(1+365*i:365+365*i);
    scatter3(x,y,z,[],c)
    ylim([0 15])
    xlim([-10 40])
    zlim([0 500])
    xlabel('temperature (℃)')
    ylabel('radiation (hours/day)')
    zlabel([ptype(pn),' (\mug/m^3)']) % μg/m3(CO为mg/m3)
    colorbar
    caxis([cmin cmax])% 4:.74-.86 15:.68-.76
    title(string(2014+i))
end
t.TileSpacing = 'compact';
t.Padding = 'compact';
%%
t = tiledlayout(2,3);
for i = 0:4
    nexttile
    X = [beijing0(46+365*i:410+365*i,[1 5 6]),beijingt(1+365*i:365+365*i),beijingr(1+365*i:365+365*i)];
    [~,score,~] = pca(X);
    x = score(:,1);
    y = score(:,2);
    z = score(:,3);
    c = beijingwi(1+365*i:365+365*i,30);
    scatter3(x,y,z,[],c)
    colorbar
    xlim([-100 200])
    ylim([-50 150])
    zlim([-50 50])
    caxis([0.15 0.5])
end
t.TileSpacing = 'compact';
t.Padding = 'compact';
%%
ts1 = lorenz63(0.01:0.01:31.45,[0 .99 0],[10,28,8/3]);
ts2 = lorenz63(0.01:0.01:40,ts1(end,:),[10,28,0]);
ts3 = lorenz63(0.01:0.01:31.45,ts2(end,:),[10,28,3]);
ts = [ts1;ts2;ts3];
n = 5000;
ts = randn(n+90,2);
nc = zeros(n,1);
c = zeros(n,1);
for i = 1:n
    nc(i) = avg_edge_overlap(VG([ts(i:i+90,1),ts(i:i+90,2)],'vg','bi'));
    c(i) = corr(ts(i:i+90,1),ts(i:i+90,2),'Type','Pearson');
end
nc = mean(nc);
% pearson
% spearman
% hvg vg
plot(1:n,nc*ones(1,n),1:n,c)
legend(["nc","c"])
%%
ts = lorenz63(0.01:0.01:50,[0 1.01 0],[10,28,8/3]);
nc = zeros(5000,1);
c = zeros(5000,1);
for i = 1:5000
    nc(i) = avg_edge_overlap(VG([ts(1:i,1),ts(1:i,2),ts(1:i,3)],'vg','bi'));
    c(i) = corr(ts(1:i,1),ts(1:i,2),'Type','Pearson');
end
% pearson
% spearman
% hvg vg
plot(1:5000,nc,1:5000,c)
legend(["nc","c"])
%%
tn = 10000;
a = Hmap(tn,rand(1,2),[1.4 0.3]);
plot(1:tn,a(:,1),1:tn,a(:,2))
%%
B = 0.3;
w = 90;
tspan = 1:5000+w;
m = 1;

nc = zeros(5000,1);
lc = zeros(5000,1);
n = 1:5090;
c = 0.3*sin(n*pi/250)+0.3;
for j = 1:m
    p0 = rand(1,4);
    ts = CHenonM(tspan,p0,c,B)+0.5*randn(5090,2);
    %     pulsespan = 2045:3045;
    %     ts = CHM_pulse(tspan,pulsespan,p0,cp,B);
    for i = 1:5000
        nc(i) = nc(i)+avg_edge_overlap(VG([ts(i:i+w,1),ts(i:i+w,2)],'hvg','bi'));
        lc(i) = lc(i)+corr(ts(i:i+w,1),ts(i:i+w,2),'Type','Pearson');
    end
end
nc = nc/m;
lc = lc/m;
figure(1)
plot(1:5000,nc,1:5000,lc)
% legend(["0",".1",".2",".3",".4",".5",".6",".7",".8",".9","1"])
% plot(1:5000,nc,1:5000,c)
legend(["nc","lc"])
%%
w = 90;
nc = zeros(1916-w,1);
c = zeros(1916-w,1);
for j = 1:1
    ts = beijing0(1:1916,[1 6]);
    for i = 1:1916-w
        nc(i) = avg_edge_overlap(VG([ts(i:i+w,1),ts(i:i+w,2)],'hvg','bi'));
%         tmp = Interlayer_mutual_info(VG([ts(i:i+w,1),ts(i:i+w,2)],'hvg','bi'));
%         nc(i) = tmp(1,2);
        c(i) = corr(ts(i:i+w,1),ts(i:i+w,2),'Type','Pearson');
    end
    plot(1:1916-w,c,1:1916-w,nc)
    hold on
end

%%
a = VG([beijing0(1:1800,6)],'hvg','bi');
a = sum(a,2);
a = reshape(a,[30 60]);
a = mean(a,1);
plot(1:60,a)
%%
w2 = 159;
t2),15),16:350,running_mean(tmp(88+365*(i-1):452+365*(i-1),3),15))
  
h = legend(["PM2.5-NO2","PM2.5-O3","NO2-O3"]);



set(h,'Fontsize',5);
    ylim([0 .6])
    ylabel('synchronisation')
    xlim([0 366])
    title(tts(i))
end
t.TileSpacing = 'compact';
t.Padding = 'compact';
% x = (w2+rn+1:len*24-w2-10+1-rn)/24;

%% syn
cd 'D:\学习\本研\statistics\2021刘升葳\HVGtest';
w2 = 729;
beijing2016all = [beijing2015(end-w2+1:end,:);beijing2016;beijing2017(1:w2+9,:)];
traj1 = beijing2016all(:,[6 9]);
tmp1 = xcorrelation(DTS(re_traj(traj1,10,1),9,w2),'hvg','uni','syn');
traj2 = beijing2016all(:,[6 11]);
tmp2 = xcorrelation(DTS(re_traj(traj2,10,1),9,w2),'hvg','uni','syn');
traj3 = beijing2016all(:,[9 11]);
tmp3 = xcorrelation(DTS(re_traj(traj3,10,1),9,w2),'hvg','uni','syn');
tmp = [tmp1,tmp2,tmp3];
%%
m = matfile('re2016.mat','Writable',true);
m.bj2016syn = tmp;
%%
tts = "bj2016";
len = 366;
rn = 7;
% x = (w2+rn+1:len*24-w2-10+1-rn)/24;
x = (rn+1:len*24-rn)/24;
for tt = 1:3
    plot(x,running_mean(tmp(:,tt),rn));
%     plot(x,tmp(:,tt))
%     plot(x,bandpass(tmp(:,tt),[1/30 1/3],24))
    hold on
end
h = legend(["PM2.5-NO2","PM2.5-O3","NO2-O3"]);
ylim([0 .6])
xlim
[0 len+1])
ylabel('synchronization')

title(tts)
%%
cd 'D:\学习\本研\statistics\2021刘升葳\HVGtest'
%%
%%
wlen = 84;%)*2+1
wi = zeros(8760-wlen*2,6);
aqi = highpass(beijing2016,1/7,24);
for i = 1:8760-wlen*2
    vg1 = VG([aqi(i:i+wlen,6),aqi(i:i+wlen,9)],'hvg','bi');
    wi(i,1) = avg_edge_overlap(vg1);
    imi = Interlayer_mutual_info(vg1);
    nfac = joint_entropy(vg1);
    nimi = imi./nfac;
    wi(i,4) = nimi(1,2);
    vg2 = VG([aqi(i:i+wlen,6),aqi(i:i+wlen,11)],'hvg','bi');
    wi(i,2) = avg_edge_overlap(vg2);
    imi = Interlayer_mutual_info(vg2);
    nfac = joint_entropy(vg2);
    nimi = imi./nfac;
    wi(i,5) = nimi(1,2);
    vg3 = VG([aqi(i:i+wlen,9),aqi(i:i+wlen,11)],'hvg','bi');
    wi(i,3) = avg_edge_overlap(vg3);
    imi = Interlayer_mutual_info(vg3);
    nfac = joint_entropy(vg3);
    nimi = imi./nfac;
    wi(i,6) = nimi(1,2);
end
%% scaling
scalem = 24;
wimean = zeros(scalem,6);
wistd = zeros(scalem,6);
for newscale=1:24
    wiscale = winHVGana(beijing2016,newscale,'avg');
    wimean(newscale,:) = mean(wiscale,1);
    wistd(newscale,:) = std(wiscale,1);
end
%%
ind = 4;
errorbar(1:24,wimean(:,ind),wistd(:,ind))
%%
m = matfile('211006tmp.mat','Writable',true);
m.bj2016wi721 = wi;
% m.bj2016wih = wih;
% m.bj2016wil = wil;
%%
%% wi
aqi = [beijing2015(end-360+1:end,:);beijing2016;beijing2017(1:360,:)];
pm25 = aqi(:,1);

o2 = aqi(:,2);
o3 = aqi(:,3);
ts = [pm25,no2,o3];
freq = 1;
wi = [winHVGana(ts(:,[1 2]),1,'no'),winHVGana(ts(:,[1 3]),1,'no'),winHVGana(ts(:,[2 3]),1,'no')];
% tsh = highpass(ts,freq,24);

% wih = [winHVGana(ts(:,[1 2]),1,'no'),winHVGana(ts(:,[1 3]),1,'no'),winHVGana(ts(:,[2 3]),1,'no')];
% tsl = lowpass(ts,freq,24);
% wil = [winHVGana(ts(:,[1 2]),1,'no'),winHVGana(ts(:,[1 3]),1,'no'),winHVGana(ts(:,[2 3]),1,'no')];
% [len,~] = size(wi);
% m.bj2016wi = wi;
%% fc-wlen filter dependence
T = [0
    25 1 3 7];
tnum = length(T);
hours = [0,744,696,744,720,744,720,744,744,720,744,720,744];
% xcm = zeros(12,2,tnum,6);
xcm = zeros(2,tnum,6);

a
i
= [beijing2015(end-360+1:end,:);beijing2016;beijing2017(1:360,:)];
p25 = aqi(:,1);
no2 = aqi(:,2);
o3 = a
i(:,3);
ts = [pm25,no2,o3];
wi

= [wnVana(ts(:,[1 2]),1,'no'),winHVGana(ts(:,[1 3]),1,'no'),winHVGana(ts(:,[2 3]),1,'no')];
for i=:
    eth(T)
    fre = 1/T(i);
    tsh
     
    hghass(ts,freq,24);
    wih=
    
    
    in
    VGana(tsh(:,[1 2]),1,'no'),winHVGana(tsh(:,[1 3]),1,'no'),winHVGana(tsh(:,[2 3]),1,'no')];    
    s
    =lowpass(ts,freq,24);
    wi
     
    =
    [winHVGana(tsl(:,[1 2]),1,'no'),winHVGana(ts(:,[1 3]),1,'no'),winHVGana(ts(:,[2 3]),1,'no')];
    
    o
    r j=1:12
     
    i
    (sum(hours(1:j))+1):sum(hours(1:j+1));
       for k=1:6
            [hxcm,~] = xcorr(wih(ti,k),wi(ti,k),0,'normalized');
            [lxcm,~] = xcorr(wil(ti,k),wi(ti,k),0,'normalized');
    
            xcm(j,:,i,k) = [hxcm,lxcm];
       end    end
end
%%
t

ype = 1;
colors = ['r','g','b','y'];
fo
r
tn = 1:tnum
  
lot(1:12,xcm(:,1,tnum,type),1:12,xcm(:,2,tnum,type),'color',colors(tnum))
    hold on
end

 normal plot
% wi = bj2016wi;
% wil = bj2016wimin;
% wih = bj2016wimax;
tl = tiledlayout(2,3);
wlen = 0;
rn = 0;
t = (wlen+1+rn:8784-wlen-rn)/24;
yl = ["AEO","AEO","AEO","IMI","IMI","IMI"];
tts = ["PM2.5-NO2","PM2.5-O3","NO2-O3","PM2.5-NO2","PM2.5-O3","NO2-O3"];
ord = [1 3 5 2 4 6];
for ind = 1:6
    nexttile
    original = running_mean(bj2016wi(:,ord(ind)),rn);
    upper = running_mean(wimax(:,3*ind-2),rn);
    lower = running_mean(wimin(:,3*ind-2),rn);
    plot(t,original)%,t,upper,t,lower)
%     errorbar(1:20,bj2016wimean(:,ind),bj2016wistd(:,ind))
    xlim([0 367])
%     xlim([0 21])
%     xlabel('sampling interval(hour)')
    ylabel(yl(ind))
    title(tts(ind))
%     legend(["original","upper","lower"])
%     legend(["original","highpass","lowpass"])
end
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
%% surrogate
wis = zeros(len,18,19);
wimin = zeros(len,18);
wimax = zeros(len,18);%pm25s-no2...;pm25-no2s...;pm25s-no2s...
for i=1:19
    [pm25s,~] = generate_iAAFT(pm25,0);
    [no2s,~] = generate_iAAFT(no2,0);
    [o3s,~] = generate_iAAFT(o3,0);

    wis(:,[1 10],i) = winHVGana([pm25s,no2],1,'no');
    wis(:,[2 11],i) = winHVGana([pm25,no2s],1,'no');
    wis(:,[3 12],i) = winHVGana([pm25s,no2s],1,'no');
    
    wis(:,[4 13],i) = winHVGana([pm25s,o3],1,'no');
    wis(:,[5 14],i) = winHVGana([pm25,o3s],1,'no');
    wis(:,[6 15],i) = winHVGana([pm25s,o3s],1,'no');
    
    wis(:,[7 16],i) = winHVGana([no2s,o3],1,'no');
    wis(:,[8 17],i) = winHVGana([no2,o3s],1,'no');
    wis(:,[9 18],i) = winHVGana([no2s,o3s],1,'no');
    
end
for j=1:len
    for k=1:18
        wimin(j,k) = min(wis(j,k,:));
        wimax(j,k) = max(wis(j,k,:));
    end
end

m = matfile('211006tmp.mat','Writable',true);

m.bj2016wimin = wimin;
m.bj2016wimax = wimax;

%%
t = (wlen+1:8784-wlen)/24;
t = t';
t = [0;t;367];
%%
wlen = 84;
ind = 1;
plot(t,wi(:,ind))%,t,wimin(:,ind),t,wimax(:,ind))
%%
w0 = (wimin(1,:)+wimax(1,:))/2;
we = (wimin(end,:)+wimax(end,:))/2;
wimin = [w0;wimin;we];
wimax = [w0;wimax;we];
%%
wimin = wimin(2:end-1,:);
wimax = wimax(2:end-1,:);
%%
patch([t fliplr(t)],[wimin(:,ind) fliplr(wimax(:,ind))],'y','edgealpha',0,'facealpha',0.5);
% h = fill([t fliplr(t)],[sin(t) fliplr(cos(t))],'r');
% h = fill(wimin(:,ind),fliplr(wimax(:,ind)),'r');
xlim([0 367])
%%
x = linspace(0, 10, 100);
y1 = sin(x);
y2 = cos(x);
% plot(x, y1, x, y2, 'LineWidth', 1.5)
area = fill([x fliplr(x)],[y1 fliplr(y2)],[0.93333, 0.83529, 0.82353],'edgealpha', '0', 'facealpha', '.5');
%%
t = 1:8784;
x = beijing2016(:,11






= sin(t.^2/4000000);
plot(t,x,t,generate_iAAFT(x))
%%

ax_pos=zeros(12,2);


gure(1)
f

i=1:12
    t = 720*(i-1)+1:720*i;
 
    x = bj2016syn(t,1);
    y = bj2016syn(t,3);
    [c,lags] = xcorr(x,y,'normalized');
    [m,p] = max(c);
    max_pos(i,:) = [m,lags(p)];
    plot(lags,c)
%     scatter()
    hold on

    end
% scatt
r(max_pos(:,1),max_pos(:,2))
legend(string(1:12))

%

lim([-60 60])







































































































