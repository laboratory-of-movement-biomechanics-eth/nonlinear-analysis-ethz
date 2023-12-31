%% Section 1. Simulating Data from Ilhen 2012
% Feel free to load your own times series and give it a try!

load fractaldata.mat
x = multifractal;
plot(multifractal)

%% Step 1: Perform MFDFA on time series
% Pre-define the inout parameters. Details can be found in the MFDFA1 function
scale = [16, 32, 64, 128, 256]; % Note the increase by power of 2

% DFA looks at q = 2, with MFDFA we can look at negative q-orders giving more
% information about the structure of the data.

q2 = -2:0.1:2; % Second statistical moment
q3 = -3:0.1:3; % Third statistical moment
q5 = -5:0.1:5; % Fifth statistical moment

% Polynomial order for detrending
m = 1; 

% Flag for output plot (1 = on, 0 = off)
plotOption = 0; 

% The next line of code runs the actual MFDFA function. An explanation of
% the steps in the analysis can be seen in slides 11-19 of "Mulitfractals - Part 3". 

figure(2);
[Hq, tq, hq, Dq, Fq] = MFDFA1(x, scale,q3 , m, plotOption);

%% Step 2: Compute Width
W = max(hq)- min(hq);

%% Step 3:  Use IAAFT to generate N surrogate

% Our interest in constructing 95% confidence intervals. Hence, we will
% need at least 40 surrogates (e.g., 1/40 = 0.025, 39/40 = 0.975). 
n_surrogates = 40;

% construct n_surrogates total surrogates for each input signal
surrogate_x = iaaft(x, n_surrogates);

%% Step 4: Perform MFDFA on each N surrogates

for i = 1:n_surrogates
   [Hqs, tqs, hqs(:,i), Dqs, Fqs] = MFDFA1(surrogate_x(:,i), scale,q3 , m, plotOption);
end

%% Step 5:  Compute WidthSurrogate:
Ws = max(hqs)- min(hqs);
% Ws = Ws';

%% Step 6: Given N observation of Ws, Compute 2.5th and 97.5TH quantiles
% construct 95% confidence intervals as the 2.5% and 97.5% quantiles at
% each width, Ws. Hence we will length(Ws) confidence intervals in all, one
% for each of the original spectral width(w) estimated in Section !
% surrogate_cis = zeros(1,length(Ws));

ci_ws = quantile(Ws, [0.025 0.975]);
%% Step 7: Compare W to CI(Ws), reject the H0= linear system.
msg1 = sprintf('The observed width, W = %d',W);
msg1
msg2 = sprintf('The CI(Ws) is between %d and %d', ci_ws(1), ci_ws(2));
msg2
is_significant = W < ci_ws(1) | W > ci_ws(2);
if is_significant 
    msg3 = sprintf('We reject H0 = linear system.');
else 
    msg3 = sprintf('We fail to reject H0 = linear system.');
    
end
msg3
%Ws goes from 0.5 to 0.8... and W= 1.4923 
%Reject the null hypothesis (=linear system)
