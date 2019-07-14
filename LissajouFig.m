function [ Mag,Phase ] = LissajouFig( u,y,t,w,T,TDOM )
%[ Mag,Phase ] = LissajouFig( u,y,t,w )
% plot Lissajou figure with 
%   u = A*cos(w*t)     input siganl
%   y                  outputsignal
%   t                  time [sec]
%   w                  input frequency [Hz]
%   T = [tstart tstop] time interval [sec]
%   TDOM               plot time domain response (=1)

% Created: Daniel Rubin, 20-Aug-2015
% Update: 14-Jan-2016

if nargin<6, TDOM=0; end

dt = t(2)-t(1);
Nf=2.*round(1/w/dt/2); % number of points in 1 cycle

tstart=T(1); tstop=T(2);
if (tstop-tstart)/dt<length(u) && tstop/dt<length(u)
    u = u(floor(tstart/dt):floor(tstop/dt));
    y = y(floor(tstart/dt):floor(tstop/dt));
    t = t(floor(tstart/dt):floor(tstop/dt));
else
    disp('SIM time exceeds number of saved points. Using all points.')
end

%ds = min(Nf/1000,1/dt); % downsample
ds=1;
y = y-mean(y);

hF=figure('Name',['Lissajou w=',num2str(w)]);

if TDOM==1
    set(hF,'Position',[500,400,600,600])
    subplot(3,1,3), 
    %hAx=plotyy(t(end-2*Nf:end),u(end-2*Nf:end),t(end-2*Nf:end),y(end-2*Nf:end));
    hAx=plotyy(t(end-Nf:end),u(end-Nf:end),t(end-Nf:end),y(end-Nf:end));
    xlabel('Time [sec]'), ylabel(hAx(1),'Input'), ylabel(hAx(2),'Output')
    hold on
    umin=find(u==min(u(end-Nf:end))); 
    ymin=find(y==min(y(end-Nf:end)));
    plot(t(umin)*[1 1],[-10^6 10^6],'--b');
    plot(t(ymin)*[1 1],[-10^6 10^6],'--g');
    subplot(3,1,[1 2])
end

% main Lissajou plot
scatter(u(1:ds:end),y(1:ds:end),1,t(1:ds:end),'*');  hold on
xlabel('u');
ylabel('y');
c=colorbar;
ylabel(c,'time [sec]')

npt = round(Nf/20);
p1=polyfit(u(end-npt:end),y(end-npt:end),1);
u11=u(end); u10=u(end-Nf-npt);
y11=polyval(p1,u11); y10=polyval(p1,u10);

p2 = polyfit(u(end-Nf/2-npt:end-Nf/2),y(end-Nf/2-npt:end-Nf/2),1);
u21 = u(end-Nf/2); u20=u(end-Nf/2-npt);
y21 = polyval(p2,u21); y20 = polyval(p2,u20);

quiver([u11 u21],[y11 y21],[u11-u10 u21-u20],[y11-y10 y21-y20],'k','linewidth',3,'AutoScale','off');

A = max(u);
B = max(y);
Mag = 20*log10(B/A);
%u0=min(abs(u)); 
%Bsin=max(y(u==u0));
% y0=min(abs(y(end-Nf:end)));  
% Asin=max(u(y==y0));
% Phase = -180+[asin(Asin/A)*180/pi -asin(Asin/A)*180/pi];
Phase = (asin(y(end-Nf)/B)-asin(u(end-Nf)/A))*180/pi;

title(sprintf('Mag=%.3g db, Phase=%.3g deg',Mag,Phase))

end

