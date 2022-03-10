function prs = FT(ts)
% Fs = ?;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = ?;             % Length of signal
%t = (0:L-1)*T;        % Time vector
[L,cols] = size(ts); %列向量
if rem(L,2) == 1
    ts(L+1,:) = zeros(1,cols);
    L = L+1;
end
p2 = fft(ts);
p2_arg = angle(p2);
p2_abs1 = abs(p2(2:L/2,:));
% p2_arg1 = p2_arg(2:L/2);
% p2_arg1 = p2_arg1(randperm(L/2-1)); %DC???
p2_arg1 = repelem(2*pi*rand(L/2-1,1),1,cols); %repeation
p2(2:L/2,:) = p2_abs1.*exp(1i*p2_arg1);
p2((L/2+2):end,:) = conj(flipud(p2(2:L/2,:)));

prs = ifft(p2);

end