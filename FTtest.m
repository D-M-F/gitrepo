function wi = FTtest(source,mode)

aqi = source;

if strcmp(mode,'FT')
    aqi = FT(aqi);
end
wlen = 84;
wi = zeros(8760-wlen*2,6);
for i = 1:8760-2*wlen
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

end
