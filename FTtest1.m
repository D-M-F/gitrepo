function tmp = FTtest1(v1,v2,v3,mode)

v3(1417:1440,:) = [];
w2 = 729;
beijing2016all = [v1(end-w2+1:end,:);v2;v3(1:w2+9,:)];

if strcmp(mode,'FT')
    beijing2016all = FT(beijing2016all);
end

traj1 = beijing2016all(:,[6 9]);
tmp1 = xcorrelation(DTS(re_traj(traj1,10,1),9,w2),'hvg','uni','syn');
traj2 = beijing2016all(:,[6 11]);
tmp2 = xcorrelation(DTS(re_traj(traj2,10,1),9,w2),'hvg','uni','syn');
traj3 = beijing2016all(:,[9 11]);
tmp3 = xcorrelation(DTS(re_traj(traj3,10,1),9,w2),'hvg','uni','syn');
tmp = [tmp1,tmp2,tmp3];

end
