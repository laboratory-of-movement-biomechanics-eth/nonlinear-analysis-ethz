%% Phase space reconstruction

%% time series of pendulum
t = linspace(0, 2*pi * 10, 634);
x = sin(t)';
y = cos(t)';

% plot time series
plot(t,x, t,y)

% plot time series
plot(x,y)


%% autocorrelation
lag = 200;
c = xcorr(x,lag,'coeff');
plot(-lag:lag, c)
grid on


tau = 16;

%% false nearest neighbours

fnn(x,6,tau)


%% embedding vector
m = 2;
N = length(x) - (m-1)*tau;

xe = zeros(N,m);
for j = 1:m
   xe(:,j) = x([1:N]+tau*(j-1));
end

plot(xe(:,1), xe(:,2))


%% create external function 'embed.m'
xe = embed(x,2,tau);


%% Lorenz data
clear
x = load('lorenz.csv');

plot(x)

%% phase space
plot3(x(:,1),x(:,2),x(:,3))


%% autocorrelation
lag = 500;
c = xcorr(x(:,1),lag,'coeff');
plot(-lag:lag, c)
grid on

% correlation time
line([-500 500],[1 1]./exp(1))


% false nearest neighbours
fnn(x(:,1))


% phase space reconstruction
tau = 4;

plot3(x(1:end-2*tau,1),x(1+tau:end-tau,1),x(1+2*tau:end,1))

%% use mutual information instead of ACF
mi(x(:,1),10,50)


%% Recurrence plot

% RP for a 3-dimensional system
whos x
D = squareform(pdist(x));

% show distance matrix
imagesc(D)

% recurrence matrix
imagesc(D < quantile(D(:),.1))

%% using CRP toolbox
crp2(x(1:end,:)) % for multidimensional vectors

crp(x(1:3:end,1),3,5) % for univariate time series

R = crp(x(1:3:end,1),3,5, 'nogui'); % store RP in variable R
imagesc(R)
colormap([1 1 1; 0 0 0])


% some test models
N = 1000;
x = rand(N,1);
x = sin(linspace(0,2*pi * 10,N));
x = sin(linspace(0,2*pi * 10,N)) + sin(linspace(0,2*pi * 32.4,N));

x = rand(1,1);
for i = 2:N
    x(i) = .99 * x(i-1) + randn(1,1);
end

plot(x)
crp(x)


% some data
x = load('tianmen.csv');
x = load('qunf.csv');
x = load('lonar.csv');

t = x(1,1):5:x(end,1);
xi = interp1(x(:,1),x(:,2),t);

plot(t,xi)


%% Recurrence quantification

% for noise
N = 1000;
x = rand(N,1);

X = crp(x,1,1,.1,'rr');
[l set_l] = dl(X);

bins = 1:20;
hl = hist(set_l,bins);

bar(bins,hl)

% determinism
sum(set_l(set_l > 1)) / sum(set_l)


% determinism as inline function
det = @(x) sum(x(x > 1)) ./ sum(x)
det(set_l)



crqa(x)


% for logistic map
N = 10000;
a = linspace(3.4,4,N)';
x = .451;
for i = 2:N
    x(i,1) = a(i) * x(i-1) * (1 - x(i-1));
end

plot(a,x,'k.')


%% RP
crp([a(1:2000), x(1:2000)])

crp([a(3000:end), x(3000:end)])


%% RQA
crqa([a, x],1,1,.1)

% windowed RQA
w = 200; % window size
ws = 100; % window step

d = zeros((length(x)-w)/ws,1); % vector for determinism
lmean = zeros((length(x)-w)/ws,1); % vector for mean line length

cnt = 1;

% loop through all windows
setall_l = [];
for i = 1:ws:length(x)-w
   x_win = x(i:i+w);
   X = crp(x_win, 1, 1, .1, 'euc', 'sil'); % RP
   [l set_l] = dl(X); % diagonal lines
   setall_l = [setall_l; set_l];
   d(cnt) = det(set_l); % determinism
   lmean(cnt) = mean(set_l(set_l>1));
   cnt = cnt + 1;
end

plot(a(1:ws:length(x)-w), d)

plot(a(round(w/2)+(1:ws:length(x)-w)), d)

plot(a(w+(1:ws:length(x)-w)), d)



subplot(211)
plot(a,x,'k.')
xlim([a(1) a(end)])

subplot(212)
plot(a(round(w/2)+(1:ws:length(x)-w)), d)


% significance test
X = crp(x, 1, 1, .1, 'euc', 'sil'); % RP of all
[l setall2_l] = dl(X); % set of diagonal lines

bootstat = bootstrp(10, @mean, setall_l(setall_l>1));
l_CI = quantile(bootstat, .95);

clf
plot(a(1:ws:length(x)-w), lmean)
line([a(1) a(end)], [l_CI l_CI], 'color', [1 .3 0])

bootstat = bootstrp(100, det, setall_l);
det_CI = quantile(bootstat, .95);

clf
plot(a(1:ws:length(x)-w), d)
line([a(1) a(end)], [det_CI det_CI], 'color', [1 .3 0])

