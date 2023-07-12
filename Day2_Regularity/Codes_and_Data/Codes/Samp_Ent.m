function [SE,Amr,Bmr] = Samp_Ent(data,m,R)
% SE = Samp_Ent(data,m,R)
% This is a faster version of the previous code - Samp_En.m

% inputs     - data, single column time seres
%            - m, length of vectors to be compared
%            - R, radius for accepting matches (as a proportion of the
%                 standard deviation)

% output     - SE, sample entropy
%
% Remarks
% This code finds the sample entropy of a data series using the method
% described by - Richman, J.S., Moorman, J.R., 2000. "Physiological
% time-series analysis using approximate entropy and sample entropy."
% Am. J. Physiol. Heart Circ. Physiol. 278, H2039–H2049.
%

% J McCamley May, 2016

% Define r as R times the standard deviation
r = R * std(data);
N = length(data);

dij=zeros(N-m,m+1);
dj=zeros(N-m,1);
dj1=zeros(N-m,1);
Bm=zeros(N-m,1);
Am=zeros(N-m,1);

for i = 1:N-m
    for k = 1:m+1
        dij(:,k) = abs(data(1+k-1:N-m+k-1)-data(i+k-1));
    end
    dj = max(dij(:,1:m),[],2);
    dj1 = max(dij,[],2);
    d = find(dj<=r);
    d1 = find(dj1<=r);
    nm = length(d)-1; % subtract the self match
    Bm(i) = nm/(N-m);
    nm1 = length(d1)-1; % subtract the self match
    Am(i) = nm1/(N-m);
end
Bmr = sum(Bm)/(N-m);
Amr = sum(Am)/(N-m);

SE = -log(Amr/Bmr);
end
