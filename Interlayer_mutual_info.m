function I_mat = Interlayer_mutual_info(vg)
% 返回dim*dim 矩阵
[~,~,dim] = size(vg);
I_mat = zeros(dim,dim);
pd_mat = pd_degree(vg); %dim行，行向量的分量值从1开始
[~,dmax] = size(pd_mat);
jpd_mat = jpd_degree(vg); %a,b,pda,pdb
for a = 1:dim-1
    for b = a+1:dim
        Iab = 0;
        for ka = 1:dmax
            for kb = 1:dmax
                if jpd_mat(a,b,ka,kb) ~= 0 
                    Iab = Iab + jpd_mat(a,b,ka,kb)*log((jpd_mat(a,b,ka,kb))/(pd_mat(a,ka)*pd_mat(b,kb)));
                end
            end
        end
        I_mat(a,b) = Iab;
    end     
end
I_mat = I_mat + I_mat';

end