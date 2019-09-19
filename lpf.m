fc=400; %carrier frequency
f0=400;  %LPF cutoff frequency 
wc=2*pi*fc;   %fc in rad/s
w0=2*pi*f0;   %f0 in rad/s
 
dw = 0.2*pi;    %summation increment for w
 
fw=1000;       %maximum frequency to analyze
ww=2*pi*fw;     %fw in rad/s
w=-ww:dw:ww;    %range of w for summation
 
ts=2000;      %time period for sampling
dt=1/ts;        %summation increment for CTFT 
t=0:dt:.2;     %time for analysis (number of time points = max(t)*Ts)
 
xt = t.*exp(-100*t);     %baseband signal
xt=xt(1:length(t));  %limit signal to time of analysis 
 
XW = FT(t,dt,xt,w);    %compute FT of x(t)
xtr = IFT(w,dw,XW,t);  %reconstruct directly from X(jw) 
 
tic;     %start clock for processing time
 
subplot(3,3,1);              %3 columns, 3 rows of plots, 1st plot
plot(t,real(xt),'r');       %plot the baseband signal on time axis
title('x(t)');         
 
subplot(3,3,2);
plot(w/(2*pi), abs(XW));      %plot CTFT for x(t)
title('X(jw)');
 
 
%_________________________________________________________
%fill in code here
yt=xt.*cos(wc*t);   %modulate the signal with the carrier
subplot(3,3,3);
plot(t,real(yt))                    %plot the transmitted signal y(t)
title('y(t)');
 
YW = FT(t,dt,yt,w);     %compute FT of y(t)
 
subplot(3,3,4);
plot(w/(2*pi), abs(YW));                    %plot Y(jw)
title('Y(jw)');
 
ft= yt.*(cos(wc*t));    %calculate demodulated signal
FW = FT(t,dt,ft,w);    %compute FT of f(t) 
 
subplot(3,3,5);
plot(w/(2*pi), abs(FW))               %Plot F(jw)
title('F(jw)');
 
%reconstruction
%LPF    design LPF that with cutoff f0 and magnitude = 2
fb=0;
wb=2*pi*fb;
frt = yrt.*(cos(wb*t));    
XrW = FT(t,dt,frt,w);
 
 
subplot(3,3,6);
plot(w/(2*pi), abs(XrW))           %plot F(jw) after LPF   
title('Xr(jw)');
 
 
subplot(3,3,7);
xrt = IFT(w,dw,XrW,t)          %compute the reconstructed signal 
plot(t,real(xrt));            %plot reconstrcuted signal
title('Xr(t)');
 
 
toc               %how long did the code take to run?
 
 
%Working functions for FT and IFT
function xt = IFT(w,dw,XW,t)
LengthT = length(t);
xt = zeros(1,LengthT);
for len=1:LengthT
    tt=t(len);
    xt(len) = sum(XW.*exp(j*w*tt))*dw/(2*pi);
end;
end
 
function XW = FT(t,dt,xt,w)
LengthW = length(w);
XW = zeros(1,LengthW);
for len=1:LengthW
    wl=w(len);
    XW(len) = sum(xt.*exp(-j*wl*t))*dt;
end;
end
