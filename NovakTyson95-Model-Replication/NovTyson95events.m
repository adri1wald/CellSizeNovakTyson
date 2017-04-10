
function [value,isterminal,direction] = NovTyson95events(t,y)
varname; 

global tSTART;
global tleaveM;

% What events do we want to know about:
% 1) when mass > 0.35
% 2) when DNA > 1.0
% 3) when Taphos/Tatotal > 0.5
% 4) when Taphos/Tatotal < 0.5

% Locate the time when values passes through zero in a 
% increasing direction and stop integration.

value=[(y(nmass) - 0.35); (y(nDNA) - 1.0); (y(nTaphos) - 0.5); (y(nTaphos) - 0.5); ((t-tleaveM) - 10); ((t-tleaveM) - 20); ((t-tleaveM) - (15+(tSTART-tleaveM)))]; 
isterminal=[1; 1; 1; 1; 1; 1; 1];  % 1 indicates that solver stops when event is reached
direction=[+1; +1; +1; -1; +1; +1; +1];