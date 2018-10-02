%% Spark.m
% Simulates a Calcium Spark

% Input: concentration - desired Calcium concentration in uM. plotthis -
% can be 0 or 1 indicating whether or not to plot graphs

function F = Spark(concentration, plotthis)

    plot_num = 1;
    xrange = 0.5; % um 
    yrange = 0.5; % um
    dt = 1; % time step
    r_inf = 0.015; %um
    r_move = sqrt(223*4*dt*10^-6); % um/us
    uptake = 1; % uptake, 12 ions/us = 0.6uM/us, so 6 is 0.3uM, us.
    Iter = 1e6; % 1 second simulation

    % Workflow of converting from ions to concentration
    % Volume (10^-15) * Avogadro's# (10^23) * [Ca2+](10^-6 M) = 10^2*[Ca2+] = #ions
    volume = 2*xrange*2*yrange*0.012; % um3 == 10^-15 L
    AvoConstant = 6.022;
    ConCoeff = volume*AvoConstant*100;
    num_ions = round(concentration*ConCoeff);

    F = zeros(Iter,1); % Spark frequency vector
    times = 1:Iter; % time vector for data analysis
    OpenRyR = zeros(Iter,1); % vector holding number of open RyRs over time

    % Create directory for organization
    Dest_Directory = sprintf('SparkVideo_%i', concentration);
    mkdir(Dest_Directory);

    % intialization of ions in the SR
    Xions = rand(num_ions,1)*2*xrange-xrange;
    Yions = rand(num_ions,1)*2*yrange-yrange;
    LocIons = zeros(num_ions,1); % Location of ions, 0 for SR & 1 for Cytosol

    % Initialize RyRs and their states
    num_RyR = 30;
    [posRyR, grid] = SingleCluster(num_RyR, r_inf);
    state = zeros(num_RyR,1);

    tau_open = 5000; % changed to 5 ms
    tau_inactive = 100000; % changed to 100 ms

    % included probability function from Zana's Thesis. Follows from previous simulation
    %% Begin Probability function to evaluate the probability function to see if a firing occurs or not
    Pmax = 0.35; %maximum probability
    Pmin = 1e-6;%minimum probability
    hillcoeff = 5; %hill coefficient (how steply the curve changes)
    Phalf = 2; %the half way point (where the curve starts changing concavity)

    % Stochastic Probability Function
    fprob = @(n_ions) Pmin+Pmax*n_ions.^hillcoeff./(n_ions.^hillcoeff+Phalf.^hillcoeff);

    % Deterministic Probability Function
    %fprob = @(n_ions) 1*(n_ions>5);

    state(1) = 1; % force open one RyR 

    %% Main Algorithm
    for time = 1:Iter

        % Move ions by random walk
        p = rand(num_ions,1)*(2*pi);
        Xions = Xions + r_move*cos(p);
        Yions = Yions + r_move*sin(p);

        % Boundary Conditions
        Xions(Xions<=-xrange) = -xrange + abs(-xrange-Xions(Xions<=-xrange));
        Xions(Xions>=xrange) = xrange - abs(Xions(Xions>=xrange)-xrange);

        Yions(Yions<=-yrange) = -yrange + abs(-yrange-Yions(Yions<=-yrange));
        Yions(Yions>=yrange) = yrange - abs(Yions(Yions>=yrange)-(yrange));

        % Determine number of ions in the cytosol, along with their exact
        % indices. Do a random uptake of ions back into SR

        SumCyto = sum(LocIons==1);
        FindLocation = find(LocIons==1);

        if SumCyto > uptake
            R = randperm(SumCyto, uptake);
            Remove = FindLocation(R);
            LocIons(Remove) = 0;
        else
            LocIons(FindLocation) = 0;
        end

        % generate random number for each RyR
        rand_num = rand(num_RyR,1);

        % Begin checking for state of each RyR
        for i = 1:num_RyR

            ind_C = LocIons==1; % indicates the location of ions in the cytosol to ind_C
            dist = (posRyR(i,1)-Xions(ind_C)).^2 + (posRyR(i,2)-Yions(ind_C)).^2 < r_inf^2;
            sum_ions = sum(dist(:));

            if state(i) == 0

                %tau_fire = 0.000001; % units should be in Probability per microsecond
                tau_fire = fprob(sum_ions);

                if rand_num(i) < (tau_fire)*dt
                    state(i) = 1;
                    %t_open(i) = time;
                end

            elseif state(i) == 1

                if rand_num(i) < (1/tau_open)*dt
                    %t_inactive(i) = time;
                    state(i) = -1;
                end

            elseif state(i) == -1

                if rand_num(i) < (1/tau_inactive)*dt
                    state(i) = 0;
                end

            end

            % When RyR is open, move ions from SR to Cytosol
            if state(i) == 1
                ind_SR = LocIons==0;
                move = (Xions - posRyR(i,1)).^2 + (Yions - posRyR(i,2)).^2 < r_inf^2 & ind_SR;
                LocIons(move) = 1;
            end
            
        end

        if mod(time,1) == 0
            SRIonsx = Xions(LocIons==0);
            SRIonsy = Yions(LocIons==0);
            CytosolIonsx = Xions(LocIons==1);
            CytosolIonsy = Yions(LocIons==1);

            OpenRyR(time) = sum(state==1);
            sumC = length(CytosolIonsx);
            F(time) = sumC/ConCoeff;

            f = figure('visible', 'off');
            plot(SRIonsx, SRIonsy + 3*yrange, 'r.', 'DisplayName', 'Calcium Ion in SR');
            hold on
            plot(CytosolIonsx, CytosolIonsy, 'b.', 'DisplayName', 'Calcium Ion in Cytosol');
            hold on
            plot([-xrange,xrange,xrange, -xrange,-xrange],[-yrange,-yrange,yrange, yrange, -yrange],'k-','linewidth', 5, 'DisplayName', 'Sarcoplasmic Reticulum');
            hold on
            plot(posRyR(:,1), posRyR(:,2), 'ko', 'MarkerSize', 3, 'DisplayName', 'RyR');
            hold on
            plot([-xrange,xrange,xrange, -xrange,-xrange],[-yrange,-yrange,yrange, yrange, -yrange]+3*yrange, 'Color', [0.3 0.3 0.6],'linewidth', 5, 'DisplayName', 'Cytosol');
            title(['Calcium Spark, Time = ' num2str(time) ' us'] , 'FontSize', 24);
            hold on
            axis([-xrange xrange -yrange 4*yrange]);
            axis equal
            %lgd = legend('show', 'location', 'northeastoutside');
            f_name = sprintf('test%.4d.png',plot_num);
            plot_num = plot_num+1;
            saveas(f, [pwd '/' Dest_Directory '/' f_name]);
            close(f);
        end
        
    end

    if plotthis == 1
        g = figure('visible', 'off');
        gname = 'OpenRyRs.png';
        plot(times, OpenRyR, 'linewidth', 2);
        xlabel('time (us)');
        ylabel('Number of open RyRs');
        saveas(g, [pwd '/' Dest_Directory '/' gname]);
        close(g);

        h = figure('visible', 'off');
        hname = 'CytosolicConcentration.png';
        plot(times, F, 'linewidth', 2);
        xlabel('time (us)');
        ylabel('[Ca^2^+] in Cytosol (uM)');
        saveas(h, [pwd '/' Dest_Directory '/' hname]);
        close(h);
    end

end
