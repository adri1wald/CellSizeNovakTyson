
function [value,isterminal,direction] = NovTyson95eventse_cycloPulse(t,y)
evarname; 

global tSTART;
global tleaveM;
global tleaveS;
global pulseStartPostG2;
global pulseLength;

% What events do we want to know about:
% 1) when mass > 0.35
% 2) when DNA > 1.0
% 3) when Taphos/Tatotal > 0.5
% 4) when Taphos/Tatotal < 0.5
% 5) when time elapse in G2 > pulseStartPostG2
% 6) when time elapse in G2 > pulseStartPostG2+pulseLength

% Locate the time when values passes through zero in a 
% increasing direction and stop integration.

% value=[(y(nmass) - 0.35); (y(nDNA) - 1.0); (y(nTaphos) - 0.5); (y(nTaphos) - 0.5); ((t-tleaveM) - 10); ((t-tleaveM) - 20); ((t-tleaveM) - (15+(tSTART-tleaveM))); ((t-tleaveS)-pulseStartPostG2)]; 
% isterminal=[1; 1; 1; 1; 1; 1; 1; 1];  % 1 indicates that solver stops when event is reached
% direction=[+1; +1; +1; -1; +1; +1; +1; +1];

value=[(y(nmass) - 0.35); (y(nDNA) - 1.0); (y(nTaphos) - 0.5); (y(nTaphos) - 0.5); ((t-tleaveM) - 10); ((t-tleaveM) - 20); ((t-tleaveM) - (15+(tSTART-tleaveM))); ((t-tleaveS)-pulseStartPostG2); ((t-(tleaveS+pulseStartPostG2))-pulseLength); ((t-tleaveM) - 3); ((t-tleaveM) - 23); ((t-tleaveM) - 43)];
isterminal=[1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];  % 1 indicates that solver stops when event is reached
direction=[+1; +1; +1; -1; +1; +1; +1; +1; +1; +1; +1; +1];