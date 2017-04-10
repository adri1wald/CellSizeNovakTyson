function [tout, yout, G1s, Ss, G2s, Ms, divmass, cycles,errcode,kout] = novaktyson95e(k)

    evarname;

    %global mu;
    global kt;
    global k1AA;
    global kdeg;
    global ksynscale;
    global V;
    global tSTART;  
    global tleaveM;
    global tleaveS;
    global pulseStartPostG2;
    global negPulse;
    global pulseLength;
    
    errcode = 0;

    pulseStartPostG2 = 0;
    negPulse = 23;
    pulseLength = 20;
    ksynscale = 1.0;
    pulse_yet = 0;
    pulse_over = 0;
    pulse_time = 0;
    pulse_end = 0;

    kt = k(42);
    k1AA = k(1);
    kdeg = k(45);
    V = 0.0;

    % NEED CODE HERE TO GENERATE INITIAL CONDS
    y0 = NovTyson95initcondse;
    t0=0;
    cctime = 180;
    tf=7*cctime;
    
    %mu = 0.0;

    %[t,y] = ode15s(@(t,y) NovTyson95eqnsc(t,y,k),[t0 tf],y0);

    %y0 = y(length(t),:);
    y0(nmass) = 0.5;

    % UNCOMMENT THIS LINE FOR ORIGINAL
    %options = odeset('Events',@NovTyson95eventse,'InitialStep',1.0);

    % UNCOMMENT THIS LINE FOR CYCLO PULSING EXPERIMENTS
    options = odeset('Events',@NovTyson95eventse_cycloPulse,'InitialStep',1.0);
    
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
    divmass = [];

    disp('ENTERING SOLVER');

    while t0<tf

        try 
            [t,y,te,ye,ie]=ode15s(@(t,y) NovTyson95eqnse(t,y,k),[t0 tf],y0,options);
        catch
            disp(k);
            disp('Bad parameter set:');
            errcode = 1;
            kout = k;
            return;
        end;
        
        kout = k;
            
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
        
        if ((ie == 10) & (negPulse == 3) & (t > 2*cctime) & (~pulse_yet))
            disp('cycloheximide pulse HERE');
            ksynscale = 0.0;
            kt = 0.0;
            k1AA = 0.0;
            pulse_yet = 1;
            pulse_time = te;
            disp(pulse_time);
        end;
        
        if ie == 5
            disp('Elapsed 10 minutes from mitosis');
            if (ye(nmass) > 0.35)
                %disp('Mass exceeds 0.35 after 10 mins');
                disp('Passes START after 10 mins');
                % mass exceeds 0.35 for the first time
                tSTART = te;
                G1s = [G1s; (tSTART+15-tleaveM)];
                k1AA = k(1);
                kt = k(43);
            else
                disp('Continue in G1');
                k1AA = k(47);
                kt = k(42);
            end;
        end;
        
        % (ie==9) & (tleaveS ~= 0) & (t > 2*cctime) & (pulse_yet)
        if ((ie == 11) & (negPulse > 0))
            if ((~pulse_over) & (t > 2*cctime) & (pulse_yet))
                disp('reset Cdc25 synthesis rate');
                ksynscale = 1.0;
                kt = k(43);
                k1AA = k(1);
                pulse_over = 1;
                pulse_end = te;
            end;
            
            if ((t > 2*cctime) & (~pulse_yet))
                disp('cycloheximide pulse HERE');
                ksynscale = 0.0;
                kt = 0.0;
                k1AA = 0.0;
                pulse_yet = 1;
                pulse_time = te;
                disp(pulse_time);                
            end;
        end;
    
        if (ie == 1)
            if ((te-tleaveM) >= 10)
                disp('Passes START after 10 mins');
                tSTART = te;
                G1s = [G1s; (tSTART+15-tleaveM)];
                k1AA = k(1);
                kt = k(43);
            else
                k1AA = k(47);
                kt = k(42);
            end;
        end;
    
        if ((ie == 6) | ((ie == 9) & (tleaveS == 0)))
            disp('Completes 20 mins of high Cdc25 degradation after mitosis');
            kdeg = k(46);
        end;
    
        if ((ie == 7) & (tSTART ~= lasttSTART))
            disp('Enters S phase');
            tenterS = te;
            disp(tenterS);
            kt = k(43);
            V = k(44);
        end;
    
        if ((ie == 2) & (tenterS ~= 0))
            disp('Leaves S phase');
            tleaveS = te;
            lasttleaveS = tleaveS;
            Ss = [Ss; (tleaveS-tenterS)];
            kt = k(42);
            V = 0.0;
            if ((pulseStartPostG2 == 0) & (t > 2*cctime) & (~pulse_yet) & (negPulse == 0))
                disp('cycloheximide pulse HERE');
                ksynscale = 0.0;
                kt = 0.0;
                k1AA = 0.0;
                pulse_yet = 1;
                pulse_time = te;
                disp(pulse_time);
            end;
        end;
    
        if ((ie==3) & (tSTART ~= lasttSTART))
            disp('Enters M phase');
            % Taphos exceeds 0.5 for the first time
            tenterM = te;
            G2s = [G2s; (tenterM-tleaveS)];
            kt = k(42);
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
            y0(1,nDNA) = 0.0;
            divmass = [divmass; y0(1,nmass)];
            y0(1,nmass) = y0(1,nmass)/2;
            k1AA = k(1);
            kdeg = k(45);
            kt = k(42);
            V=0.0;
        end;

        if ((ie==8) & (tleaveS ~= 0) & (t > 2*cctime) & (~pulse_yet) & (negPulse == 0))
            % set size-dependent synthesis rate of Cdc25 to ksynscale*(it's
            % current value)
            disp('cycloheximide pulse HERE');
            ksynscale = 0.0;
            kt = 0.0;
            k1AA = 0.0;
            pulse_yet = 1;
            pulse_time = te;
        end;
        
        if ((ie==9) & (tleaveS ~= 0) & (t > 2*cctime) & (pulse_yet) & (~pulse_over) & (negPulse == 0))
            % reset size-dependent synthesis and growth rate to 100%
            disp('reset Cdc25 synthesis rate');
            ksynscale = 1.0;
            kt = k(42);
            k1AA = k(1);
            pulse_over = 1;
            pulse_end = te;
        end;
    
        if ((ie == 12) & (t > 2*cctime) & (pulse_yet) & (~pulse_over) & (negPulse == 23))
            disp('reset Cdc25 synthesis rate');
            ksynscale = 1.0;
            kt = k(43);
            k1AA = k(1);
            pulse_over = 1;
            pulse_end = te;            
        end
        
        % After every event takes place, reinitialize the start time (t0) to
        % time step at event (t(nt)) and begin integration again
        t0 = t(nt);
    
        if t0 >= tf
            break;
        end;
    
    end;
    
    divs = [];
    
    for i=1:(size(yout,1)-1)
        if (round(yout(i,size(yout,2)-1)) == 1.0) & (round(yout(i+1,size(yout,2)-1)) == 0.0)
            divs = [divs; i];
        end
    end;
    
    divs = divs(2:end);
    
    disp('numel(divs)');
    disp(numel(divs));
    
    if (numel(divs) < 2)
        disp('Bad parameter set:');
        errcode = 1;
        kout = k;
        return;
    end;
    
    %yd = yout((divs(1)+1):divs(end),:);
    %td = tout((divs(1)+1):divs(end),:);
    
    yd = yout(1:divs(end),:);
    td = tout(1:divs(end),:);
    
    %pulse_time = pulse_time - td(1);
    %td = td - td(1);
    
    disp(pulse_time);
    disp(pulse_end);
    
    %yd = yd(td <= 560,:);
    %td = td(td <= 560);
        
    figure(2);
    plot(td,yd(:,nmass),'b','LineWidth',2); set(gca,'FontSize',18,'Position',[.05 .05 .9 .9]);
    hold on;
    plot(td,yd(:,nDNA),'r','LineWidth',2);
    hold on;
    plot(td,yd(:,nCdc25active),'g','LineWidth',2);
    hold on;
    plot(td,yd(:,nCdc25total),'LineWidth',2);
    hold on;
    plot(td,yd(:,nMPFactive),'LineWidth',2);
    hold on;
    plot(td,yd(:,nCdc13free)+yd(:,nMPFactive)+yd(:,nMPFinactive)+yd(:,nMPFphos2)+yd(:,nMPFphos0),'LineWidth',2);
    hold on;
    plot(td,yd(:,nWee1active),'LineWidth',2);
    hold on;
    %plot(td,yd(:,nX),'LineWidth',2);
    %hold on;
    %plot(td,yd(:,nW),'LineWidth',2);
    %hold on;
    %plot(td,yd(:,nTaphos),'LineWidth',2);
    %hold on;
    if (pulse_over == 1)
        line([pulse_time pulse_time],[0 1.5],'LineWidth',2,'LineStyle','--','Color','k');
        hold on;
    end;
    h = legend('MASS','DNA','Cdc25active','Cdc25total','MPFactive','Cdc13total','Wee1active');
    %h = legend('MASS','DNA','Cdc25active','W','TAphos');
    h.FontSize = 18;
    title('Wild-Type ODE Model (SC15) Simulations');
    xlabel('Time');
    ylabel('Concentration');
    hold off;
    drawnow;



end