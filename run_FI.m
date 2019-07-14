% correlation method 
clc

% exmaple for a 2 input 2 output system
% T11 = tf(1,[1 4 7]);
% T12 = tf(2,[6 3 1]);
% T21 = tf(1,[1 3 10]);
% T22 = tf(2,[3 3 1]);
% P = [T11 T12 ; T21 T22];

Model='FI_tpl'; % simulink model name
Prefix='test';   % Prefix for template file names
LissajouPlot = 0 ; % Lissajou Plot 
% select input channel
%I=1;    
I=2; 

CSclick=1; % resend from CarSim at first run

Oname={'Out 1','Out 2'};
ny=2; % number of outputs

% Frequencies in rad/s:

frq = logspace(-1,2,20);

ncyc = ceil(interp1([0.1 6 250],[3 8 30],frq)); % number of cycles until steady state

%% generate parameter grid
Parnames = {'Amp','dc1','dc2','a1','a2'}.'; % a1,a2 are some uncertain parameters..
% Mix different amplitudes and dc steering
[a1,a2,a3,a4,a5] = pargrid(...
            [1,3],...          
            [0,1],...              
            [0,1],...
            [2,3],...
            [5,6,8]);
PARGRID = [a1; a2; a3; a4; a5];

PARGRID = [1 0 0 2 5]';

Ncase = size(PARGRID,2);   % number of simulated cases

%% initiate some parameters...
dt_sim = 0.001; % simulation time step
t0=2; % time until input starts
N0=1/dt_sim*t0;

% for initial diagram update
w=0.1; 
Amp=0; 
tstart=3; 
corr_select=zeros(1,ny);

% filter
wc=20; % cut-of frequency, in Hz
ws=70; % stopband frequecny in Hz
Rp = 3; % Passband ripple, in db
Rs = 60; % Stopband attenuation, in db

[n,Wn] = buttord(wc*2*pi,ws*2*pi,Rp,Rs,'s');
[AAF_NUM,AAF_DEN] = butter(n,Wn,'s');

%% RUN!!!
switch I
    case 1, corr_select = [1 0]; Iname='input 1';      
    case 2, corr_select = [0 1]; Iname='input 2'; 
end

for ko=1:ny
    Pname{ko}=sprintf('%s_%g%g',Prefix,ko,I);        
end

nrun=0; % counter for number of total runs
for icase=1:Ncase
    
    pack; % free some memory 
    Amp = PARGRID(1,icase);                     % update sin amplitude
    DC = [PARGRID(2,icase) PARGRID(3,icase)];   % update dc imput
    % update plant parameter
    a1 = PARGRID(4,icase);
    a2 = PARGRID(5,icase);
    % in this toy exmaple the plant is reconstructed here for each run
    T11 = tf(a2+5,[1 a1  7]);
    T12 = tf(a2-3,[1 a1 5]);
    T21 = tf(a2-2,[1 a1 10]);
    T22 = tf(a2,[1 a1  1]);
    P=[T11 T12; T21 T22];
    
    % run all frequencies
    for irun=1:length(frq)
        
        w=frq(irun); % [rad/s]
        par=[];
        tpl=[];
        
        if look([Pname{1},'.tpl'])==1
               
            w_tpl=getfrom([Pname{1},'.tpl'],'w_tpl');
            if any(w_tpl(:,1)==w) %&& ~ismember(w,repfreq);
                parname=sprintf('par_%g',w_tpl(find(w_tpl(:,1)==w),2));
                if ismember(PARGRID(:,icase).',getfrom([Pname{1},'.tpl'],parname).','rows')
                    continue
                else % new set of parameters @ this frequency
                    for ko=1:ny
                        [tpl(:,ko),par]=gettpl(Pname{ko},w);
                    end
                end
            end
            
        end
         
        nrun=nrun+1; % count simulink runs
           
        fprintf('\n****** Case %g/%g : ******\n',icase,Ncase);
        fprintf('\n       ****** w = %g [rad/s] (%g/%g) ********\n\n',w,irun,length(frq));
        
        w_Hz=w/2/pi; 
        T1 = 1/w_Hz; % period time [s]
        
        tstart=t0+ncyc(irun)/w_Hz+10; % number of oscillations until settling, according to table ncyc
        tstop=max([T1*ncyc(irun) ceil(2/T1)*T1])+tstart; % minimal time for stop is ncyc oscillations or 2 sec.
        tstop=roundn(tstop,-3);
 
        set_param(getActiveConfigSet(Model),'StopTime',num2str(tstop),'SimulationMode','accelerator')
        Np=max(2/0.016/dt_sim,10000); % number of point to keep in workspace
        
        sim(Model,[]);
        
        Icc = Ic.signals.values; 
        Iss = Is.signals.values; 
        for ko=1:ny
            corrdata=corrmethod2(Icc(:,ko),Iss(:,ko),w_Hz,Amp,[tstart tstop]); 
            Pk=180/pi*corrdata(2)+1j*20*log10(corrdata(1)); % Nichols form   
            if isempty(tpl)
                Tk = Pk;
            else
                Tk = [Pk ; tpl(:,ko)];
            end
            add2tpl(Pname{ko},w,Tk,[],'r',[PARGRID(:,icase) par]);
        end
        
        if LissajouPlot
            time=InputSignal.time;
            U=InputSignal.signals.values;
            Y=OutputSignal.signals.values(:,2); 
            LissajouFig( U,Y,time,w_Hz,[tstart tstop],1 );
        end
        
    end
           
end

for ko=1:ny
    insert([Pname{ko},'.tpl'],sprintf('from %s to %s',Iname,Oname{ko}),'IO','r');
    insert([Pname{ko},'.tpl'],date,'date','r');
    insert([Pname{ko},'.tpl'],Parnames,'parnames','r');
end

showtpl(Pname{1})
showtpl(Pname{2})


