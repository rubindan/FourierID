function [ G,K ] = corrmethod2( Ic,Is,om,Amp,T )
% ver2: Yc and Ys are recorded 
% [ G ] = corrmethod( y,t,om,Amp,T )
% input:    Ic  - cosine channal
%           Is  - sin channal
%           om  - input frequency [Hz]
%           Amp - sin. amplitude
%           T   - time interval        
% output:   G = [ mag phi om], where
%             mag - magnitude at given frequency [abs]
%             phi - phase at given frequency [rad]
%             om  - given frequency [Hz]
%           K = Coherence. 0<=K<=1

if size(Ic,2)>1, error('SISO systems only!'); end  
if length(T)>1, T = T(2)-T(1); end
Ic=Ic(end)/T; Is=Is(end)/T; % devide by integraion time!
mag = 2*sqrt(Ic^2+Is^2)/Amp;
% phi = -atan2(Is,Ic);
phi = atan2(Ic,Is); % rad
%if phi>0, phi=phi-2*pi; end

G = [mag phi om];



end % function

