varname; 

global mu;
global kt;
global k1AA;
global kdeg;
global V;
global tSTART;
global tleaveM;

kt = 0.0;
k1AA = 0.025;
kdeg = 0.2;
V = 0.0;

% NEED CODE HERE TO GENERATE INITIAL CONDS
y0 = NovTyson95initconds;
t0=0;
cctime = 180;
tf=7*cctime;

mu = 0.0;

[t,y] = ode15s(@NovTyson95eqns,[t0 tf],y0);

y0 = y(length(t),:);
y0(nmass) = 0.6;
%options=odeset('Events',@NovTyson95events,'RelTol',1e-4,'AbsTol',1e-6);
options = odeset('Events',@NovTyson95events,'InitialStep',1.0);

% Setting variables for the numerical integration
tout=t0; % tout contains the time values integrated over
yout=y0; % yout contains the conc. values for each variable in the ODEs
teout = []; % The "e" means events, te is the time  when the event takes place 
yeout = []; % ye is the y(value of all variables) when the event takes place.
ieout = []; % the "e" refers to the event."i" is the index of event. index ie=1 means one event  

tSTART = 0;
lasttSTART = 0;
tenterS = 0;
tleaveS = 0;
tenterM = 0;
tleaveM = 0;
lasttenterM = 0;
lasttleaveM = 0;

cycles = [];
G1s = [];
Ss = [];
G2s = [];
Ms = [];

mu = 0.00385;

disp('ENTERING SOLVER');

while t0<tf    
    
    [t,y,te,ye,ie]=ode15s(@NovTyson95eqns,[t0 tf],y0,options);
    
    nt=length(t);
    
    tout=[tout;t(2:nt)];  
    yout=[yout;y(2:nt, :)];
    
    teout = [teout;te]; % The "e" means events. "te" is the time  when the event takes place 
  	yeout = [yeout;ye];
  	ieout = [ieout;ie];
    
    % If any of the reset rules are applied (detected by the event in "options"
    % above), we refresh the initial conditions for ODEs 
    y0 = y(nt, :);
      	
     % the "e" refers to the event."i" is the index of event. index ie=1
     % means one event
     
    if ~isempty(ie)
        ie = ie(length(ie));
        te = te(length(te));
        ye = ye(size(ye,1),:);
    end;
 	%if isscalar(ie) == 0
    %	ie = 0;
    %end; 
    
    if ie == 5
        disp('Elapsed 10 minutes from mitosis');
        if (ye(nmass) > 0.35)
            disp('Mass exceeds 0.35 after 10 mins');
            % mass exceeds 0.35 for the first time
            tSTART = te;
            G1s = [G1s; (tSTART+15-tleaveM)];
            k1AA = 0.025;
            kt = 0.1;
        else
            disp('Continue in G1');
            k1AA = 0.003;
            kt = 0;
        end;
    end;
        
    if (ie == 1)
        disp('Mass exceeds 0.35');
        if ((te-tleaveM) >= 10)
            tSTART = te;
            G1s = [G1s; (tSTART+15-tleaveM)];
            k1AA = 0.025;
            kt = 0.1;
        else
            k1AA = 0.003;
            kt = 0.0;
        end;    
    end;
    
    if ie == 6
        disp('Completes 20 mins of high Cdc25 degradation after mitosis');
        kdeg = 0.04;
    end;
        
    if ((ie == 7) & (tSTART ~= lasttSTART))
        disp('Enters S phase');
        tenterS = te;
        kt = 0.1;
        V = 0.0555;
    end;
    
    if ((ie == 2) & (tenterS ~= 0))
        disp('Leaves S phase');
        tleaveS = te;
        Ss = [Ss; (tleaveS-tenterS)];
        kt = 0;
        V = 0;       
    end;
    
    if ((ie==3) & (tSTART ~= lasttSTART))
        disp('Enters M phase');
        % Taphos exceeds 0.5 for the first time
        tenterM = te;
        G2s = [G2s; (tenterM-tleaveS)];
        kt = 0;
    end;
    
    if ((ie==4) & (tenterM ~= lasttenterM))
        % Taphos falls below 0.5
        disp('Leaves M phase');
        tleaveM = te;
        disp(tleaveM);
        Ms = [Ms; (tleaveM-tenterM)];
        cycles = [cycles; (tleaveM-lasttleaveM)];
        lasttenterM = tenterM;
        lasttleaveM = tleaveM;
        lasttSTART = tSTART;
        y0(1,nmass) = y0(1,nmass)/2;
        y0(1,nDNA) = 0.0;
        k1AA = 0.025;
        kdeg = 0.2;
        kt = 0;
        V=0;
    end;
    
    % After every event takes place, reinitialize the start time (t0) to
    % time step at event (t(nt)) and begin integration again
	t0 = t(nt);
    
	if t0 >= tf
		break;
	end

end

%% PLotting
% Here we plot the graphs of the concentrations of the variables over
% time, and reproduce Chen(2004) figure 2 (NB: only the plot of variables
% mass,BUD,ORI and SPN)

divs = [];

for i=1:(size(yout,1)-1)
    if (round(yout(i,size(yout,2)-1)) == 1.0) & (round(yout(i+1,size(yout,2)-1)) == 0.0)
        divs = [divs; i];
    end
end;

divs = divs(2:end);

yd = yout((divs(1)+1):divs(end),:);
td = tout((divs(1)+1):divs(end),:);

td = td - td(1);

yd = yd(td <= 560,:);
td = td(td <= 560);

figure(1);
plot(td,yd(:,nmass),'b','LineWidth',2); set(gca,'FontSize',18,'Position',[.05 .05 .9 .9]);
hold on
plot(td,yd(:,nDNA),'r','LineWidth',2)
hold on
plot(td,yd(:,nCdc25active),'g','LineWidth',2)
hold on
plot(td,yd(:,nCdc25total),'LineWidth',2)
hold on
plot(td,yd(:,nMPFactive),'LineWidth',2)
hold on
plot(td,yd(:,nCdc13free)+yd(:,nMPFactive)+yd(:,nMPFphos2)+yd(:,nMPFinactive)+yd(:,nMPFphos0),'LineWidth',2)
hold on
plot(td,yd(:,nWee1active),'LineWidth',2)
hold on
legend('MASS','DNA','Cdc25active','Cdc25total','MPFactive','Cdc13total','Wee1active')
title('WT ODE Model (NT95) Simulations')
xlabel('Time')
ylabel('Concentration')