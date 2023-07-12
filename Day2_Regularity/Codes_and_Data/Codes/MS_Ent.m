function [MSE] = MS_Ent(data,tau,m,r)
%
% inputs     - x, single column time seres
%            - tau, scale factor
%            - m, length of vectors to be compared
%            - R, radius for accepting matches (as a proportion of the
%                 standard deviation)
%
% output     - MSE, Multiscale Entropy
%
% Remarks
% This code finds the MSE of a data series using the
% methods described by - Wu, Shuen-De, et al. 2014. "Analysis of complex 
% time series using refined composite multiscale entropy." Physics Letters 
% A. 378, 1369-1374.
%
% Created 20170828 by Will Denton (21denton@gmail.com) and adjusted by PC.
% Raffalt 20190904
x=data;
R = r*std(x);
N = length(x);

%Coarse-graining for MSE and derivatives
y_tau_kj = zeros(length(1:tau),length(1:N/tau));
for j = 1:N/tau
    for k = 1:tau
        try
            y_tau_kj(k,j) = 1/tau*sum(x((j-1)*tau+k:j*tau+k-1));
        catch
            y_tau_kj(k,j) = NaN;
        end
    end
end

%Multiscale Entropy (MSE)
MSE = Samp_Ent(y_tau_kj(1,~isnan(y_tau_kj(1,:))),m,R);

end

function [SE,sum_nm,sum_nm1] = Samp_Ent(data,m,r)
% [SE,sum_nm,sum_nm1] = Samp_Ent(data,m,r)
% This is a faster version of the previous code - Samp_En.m

% inputs     - data, single column time seres
%            - m, length of vectors to be compared
%            - R, radius for accepting matches (as a proportion of the
%                 standard deviation)

% output     - SE, sample entropy
%            - sum_nm, total number of matches for vector length m
%            - sum_nm1, total number of matches for vector length m+1
%
% Remarks
% This code finds the sample entropy of a data series using the method
% described by - Richman, J.S., Moorman, J.R., 2000. "Physiological
% time-series analysis using approximate entropy and sample entropy."
% Am. J. Physiol. Heart Circ. Physiol. 278, H2039–H2049.
%

% J McCamley May, 2016
% W Denton August, 2017 (Made count total number of matches for each vector length, necessary for CMSE and RCMSE)

N = length(data);
dij=zeros(N-m,m+1);
Bm=zeros(N-m,1);
Am=zeros(N-m,1);
sum_nm = 0;
sum_nm1 = 0;
for i = 1:N-m
    for k = 1:m+1
        dij(:,k) = abs(data(1+k-1:N-m+k-1)-data(i+k-1));
    end
    dj = max(dij(:,1:m),[],2);
    dj1 = max(dij,[],2);
    d = find(dj<=r);
    d1 = find(dj1<=r);
    nm = length(d)-1; % subtract the self match
    sum_nm = sum_nm+nm;
    Bm(i) = nm/(N-m);
    nm1 = length(d1)-1; % subtract the self match
    sum_nm1 = sum_nm1+nm1;
    Am(i) = nm1/(N-m);
end
Bmr = sum(Bm)/(N-m);
Amr = sum(Am)/(N-m);
SE = -log(Amr/Bmr);
end

