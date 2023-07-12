%% CRP
t = (1:1000)*2*pi/67;
a = sin(t);
b = sin(.01*t.^2);

crp(a,b,3,12,'distance')



%% JRP

% separate RP for a and b, using 10% RR
Ra = crp(a,3,12,.1,'rr');
Rb = crp(b,3,12,.1,'rr');

% the JRP
J = Ra .* Rb;

imagesc(J)


jrp(a,b,3,12,'rr')


%% use modified variable
c = a.^2 - .5;

plot(t,a,t,c)

% linear correlation
corrcoef(a,c)

% JRP
Ra = crp(a,3,12,.1,'rr');
Rc = crp(c,3,12,.1,'rr');
JR = Ra .* Rc;

subplot(131)
imagesc(Ra)
subplot(132)
imagesc(Rc)
subplot(133)
imagesc(JR)


mean(Ra(:))
mean(JR(:))

imagesc(JR)


jrp(a,c,3,12,'fan')



%% Phasesync

N = 2000;
phaseshift = 2*cos(linspace(0,2*pi*10,N))';
ampshift = .8+.4*sin(linspace(0,2*pi*15,N))';
t = linspace(0,2*pi*44,N)';
x = sin(t+phaseshift);
y = ampshift.*cos(t+phaseshift);

plot(t,x, t,y)
xlabel('Time'), ylabel('x_1, x_2')

corrcoef(x,y)

%% tau-RR
maxW = 800;
Rx = taucrp(x,2,12,.1,'rr',maxW);
Ry = taucrp(y,3,12,.1,'rr',maxW);

tauRRx = mean(Rx');
tauRRy = mean(Ry');

plot(0:maxW,tauRRx(maxW+1:end), 0:maxW,tauRRy(maxW+1:end))

tau=20
corr(tauRRx(maxW+1+tau:end)', tauRRy(maxW+1+tau:end)')


%% CENOGRID example
x = flipud(load('cenogrid.txt')); 

t = x(1,1):-0.005:x(end,1); % new time axis
x_ = interp1(x(:,1), x(:,2), t); % interpolate to new time axis
x = [-1000*t' x_'];

plot(x(:,1),x(:,2))



y = load('laskar_67Ma.txt'); % reverse order to start with youngest age



% show RP
toi = [12700:13420]; % time of interest

crp(x(toi,:))

% calculate tauRR
maxW = 120;
Rx = taucrp(x(toi,2),4,3,.1,'rr',maxW);
Ry = taucrp(y(toi,4),3,2,.1,'rr',maxW);

tauRRx = mean(Rx');
tauRRy = mean(Ry');
plot(5*(0:maxW),tauRRx(maxW+1:end), 5*(0:maxW),tauRRy(maxW+1:end))
xlabel('Lag (ka)')


% CPR measure
tau = 6;
corr(tauRRx(tau+maxW+1:end)', tauRRy(tau+maxW+1:end)', 'Type', 'Pearson')


% filter out trends
[a b] = butter(3,.005);
xs = filtfilt(a,b, x(toi,2));

plot(x(toi,2)-xs)

crp(x(toi,2)-xs,4,3,.1,'rr')

Rx = taucrp(x(toi,2)-xs,4,3,.1,'rr',maxW);
tauRRx = mean(Rx');
plot(5*(0:maxW),tauRRx(maxW+1:end), 5*(0:maxW),tauRRy(maxW+1:end))
xlabel('Lag (ka)')


% windowed analysis
maxW = 120;
w = 500;
ws = 100;

tau = 6;

cpr = zeros(length(x),3);
rho = zeros(length(x),3);
for i = 1:ws:length(x)-w
    Rx = taucrp(x(i:i+w,2),4,3,.1,'rr',maxW);
    Ry = taucrp(y(i:i+w,4),3,2,.1,'rr',maxW);
    Ry2 = taucrp(y(i:i+w,2),3,2,.1,'rr',maxW);
    Ry3 = taucrp(y(i:i+w,3),3,2,.1,'rr',maxW);

    tauRRx = mean(Rx');
    tauRRy = mean(Ry');
    tauRRy2 = mean(Ry2');
    tauRRy3 = mean(Ry3');
    
    cpr(i,1) = corr(tauRRx(tau+maxW+1:end)', tauRRy(tau+maxW+1:end)');
    
    cpr(i,2) = corr(tauRRx(tau+maxW+1:end)', tauRRy2(tau+maxW+1:end)');
    cpr(i,3) = corr(tauRRx(tau+maxW+1:end)', tauRRy3(tau+maxW+1:end)');
    
    rho(i,1) = corr(x(i:i+w,2), y(i:i+w,4));
    rho(i,2) = corr(x(i:i+w,2), y(i:i+w,2));
    rho(i,3) = corr(x(i:i+w,2), y(i:i+w,3));
 
end

plot(x(1:ws:end,1),cpr(1:ws:end,:)')

plot(x(1:ws:end,1),rho(1:ws:end,:)')
legend('Obl','Ecc','Prc')
xlabel('Age (kyr)')







