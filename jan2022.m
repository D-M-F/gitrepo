%%
% shanghai2016(1417:1440,:) = [];
wi0 = FTtest(beijing2017,0);
[tsrow,tscol] = size(wi0);
wimin = zeros(tsrow,tscol);
wimax = zeros(tsrow,tscol);
wis = zeros(tsrow,tscol,19);
%%
rep = 19;
for i = 1:rep
    wii = FTtest(beijing2017,'FT');
    wis(:,:,i) = wii;
end
for j = 1:tsrow
    for k = 1:tscol
        wimin(j,k) = min(wis(j,k,:));
        wimax(j,k) = max(wis(j,k,:));
    end
end
%%
x = 1:8760;
y = beijing2017(:,6);
y1 = FT(y);
plot(x,y,x,y1)
%%
m = matfile('20210113ht.mat','Writable',true);
m.bj2017wiht = wi0;
m.bj2017wiminht = wimin;
m.bj2017wimaxht = wimax;
%%
wlen = 84;
num = 6;
x = [1+wlen:8760-wlen]/24;
p = plot(x,wi0(:,num),'r');%,'b',x,wimin(:,num),'r--',x,wimax(:,num),'r--');
p(1).LineWidth = 2;
hold on
f = fill([x,fliplr(x)],[fliplr(wimin(:,num)'),wimax(:,num)'],'b');
set(f,'edgealpha',0,'facealpha',0.3)
xlim([0 365])

%%
tl = tiledlayout(2,3);
wlen = 84;
rn = 0;
x = [1+wlen+rn:8760-wlen-rn]/24;
yl = ["AEO","AEO","AEO","IMI","IMI","IMI"];
tts = ["PM2.5-NO2","PM2.5-O3","NO2-O3","PM2.5-NO2","PM2.5-O3","NO2-O3"];
% ord = [1 3 5 2 4 6];
for ind = 1:6
    nexttile

    original = running_mean(wi0(:,ind),rn);
    upper = running_mean(wimax(:,ind),rn)';
    lower = running_mean(wimin(:,ind),rn)';

    p = plot(x,original,'r');%,'b',x,wimin(:,num),'r--',x,wimax(:,num),'r--');
    p(1).LineWidth = 2;
    hold on
    f = fill([x,fliplr(x)],[fliplr(lower),upper],'b');
    set(f,'edgealpha',0,'facealpha',0.3)
    xlim([0 366])
    ylabel(yl(ind))
    title(tts(ind))
end
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
%%
%% syn
% cd 'D:\学习\本研\statistics\2021刘升葳\HVGtest';
w2 = 729;
bj2016_365 = beijing2016;
bj2016_365(1417:1440,:) = [];
beijing2016all = [beijing2015(end-w2+1:end,:);bj2016_365;beijing2017(1:w2+9,:)];
traj1 = beijing2016all(:,[6 9]);
tmp1 = xcorrelation(DTS(re_traj(traj1,10,1),9,w2),'hvg','uni','syn');
traj2 = beijing2016all(:,[6 11]);
tmp2 = xcorrelation(DTS(re_traj(traj2,10,1),9,w2),'hvg','uni','syn');
traj3 = beijing2016all(:,[9 11]);
tmp3 = xcorrelation(DTS(re_traj(traj3,10,1),9,w2),'hvg','uni','syn');
tmp = [tmp1,tmp2,tmp3];
%%
[tsrow,tscol] = size(tmp);
tmpmin = zeros(tsrow,tscol);
tmpmax = zeros(tsrow,tscol);
tmps = zeros(tsrow,tscol,19);

rep = 19;
for i = 1:rep
    tmpi = FTtest1(beijing2015,beijing2016,beijing2017,'FT');
    tmps(:,:,i) = tmpi;
end
for j = 1:tsrow
    for k = 1:tscol
        tmpmin(j,k) = min(tmps(j,k,:));
        tmpmax(j,k) = max(tmps(j,k,:));
    end
end
%%
m = matfile('20210113ht.mat.mat','Writable',true);
m.bj2016syn_365 = tmp;
%%
tts = "bj2016";
len = 365;
rn = 15;
x = (rn+1:len*24-rn)/24;
for tt = 1:1
    original = running_mean(tmp(:,tt),rn);
    upper = running_mean(tmpmax(:,tt),rn)';
    lower = running_mean(tmpmin(:,tt),rn)';

    p = plot(x,original);
%     p(1).LineWidth = 2;
    hold on
    f = fill([x,fliplr(x)],[fliplr(lower),upper],'b');
    set(f,'edgealpha',0,'facealpha',0.3)
end
h = legend(["PM2.5-NO2","PM2.5-O3","NO2-O3"]);
ylim([0 .6])
xlim([0 len])
ylabel('synchronization')

title(tts)