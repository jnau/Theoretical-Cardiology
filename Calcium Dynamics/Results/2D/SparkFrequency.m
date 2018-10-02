%% SparkFrequency.m
% Automated version to run Spark.m, in order to find average number of
% sparks for specific parameter changes. Matlab script with some
% modifications in GNU Octave

% pcase - parameter case: 0, 1, 2, vchange- variable change, vchange2 - secondary variable change if needed
% run - indicates which run, plotthis - indicates to plot

%% Initializations
function F = SparkFrequency(pcase, vchange, run, plotthis,  vchange2)

    % Initialize States
    xrange = 0.25; % um
    yrange = 0.25; % um
    dt = 1; % time step
    r_inf = 0.015; %um
    r_move = sqrt(223*4*dt*10^-6); % um/us
    Iter = 1e6; % us
    concentration = zeros(Iter,1);
    times = 1:Iter;
    Open_RyR = zeros(Iter,1);
    volume = 2*xrange*2*yrange*0.012;
    AvoConstant = 6.022;
    ConCoeff = volume*AvoConstant*100;
    
    switch pcase
        
        case 0 % To test change in [Ca2+] in SR
            num_ions = round(vchange*ConCoeff);
            uptake = 1; % uptake 0.3uM/us
            uptake_time = 1;
            num_RyR = 30;
            % num_RyR = vchange2;
            p_name = 'SRCon';
            text_name = sprintf('%s_%i_%i.txt', p_name,vchange,run);

        case 1
            num_ions = round(250*ConCoeff); % 250
            %num_ions = round(vchange2*7.2);
            num_RyR = vchange;
            uptake = 1; % uptake 0.3uM/us
            uptake_time = 1;
            p_name = 'NumRyR';
            text_name = sprintf('%s_%i_%i.txt', p_name,vchange,run);

        case 2
            num_ions = round(250*ConCoeff); % 250
            num_RyR = 30;
            uptake = vchange;
            uptake_time = vchange2;
            p_name = 'UptakeStrength';
            text_name = sprintf('%s_%i_%i.txt', p_name, vchange,run);
    end

    Xions = rand(num_ions,1)*2*xrange-xrange;
    Yions = rand(num_ions,1)*2*yrange-yrange;
    LocIons = zeros(num_ions,1); % Location of ions, 0 for SR & 1 for Cytosol
    [posRyR, grid] = SingleCluster(num_RyR, r_inf);
    state = zeros(num_RyR,1);

    %% Begin Probability function to evaluate the probability function to see if a firing occurs or not
    Pmax = 0.35; %maximum probability
    Pmin = 1e-6;%minimum probability
    hillcoeff = 5; %hill coefficient (how steply the curve changes)
    Phalf = 2; %the half way point (where the curve starts changing concavity)

    % Stochastic Probability Function
    fprob = @(n_ions) Pmin+Pmax*n_ions.^hillcoeff./(n_ions.^hillcoeff+Phalf.^hillcoeff);

    % Deterministic Probability Function
    %fprob = @(n_ions) 1*(n_ions>5);

    tau_open = 5000; % 5ms
    tau_inactive = 100000; % 100ms

    state(1) = 1;
    t1 = 0;
    t2 = zeros(Iter, 1);
    
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

        % Uptake
        if mod(time,uptake_time)==0

            % find how many ions in cytosol and their exact indices
            SumCyto = sum(LocIons==1);
            FindLocation = find(LocIons==1);

            if SumCyto > uptake
                R = randperm(SumCyto, uptake);
                Remove = FindLocation(R);
                LocIons(Remove) = 0;
            else
                LocIons(FindLocation) = 0;
            end
            
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
                % if sum_ions > 2 || rand_num(i) < tau_fire*dt
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

            if state(i) == 1
                ind_SR = LocIons==0;
                move = (Xions - posRyR(i,1)).^2 + (Yions - posRyR(i,2)).^2 < r_inf^2 & ind_SR;
                LocIons(move) = 1;
            end

        end


        if mod(time,1) == 0

            CytosolIonsx = Xions(LocIons==1);
            Open_RyR(time) = sum(state==1);
            sumC = length(CytosolIonsx);
            concentration(time) = sumC/ConCoeff;
            
            if Open_RyR(time) > 4.5 && t1 ==0
                t1 = time;
            elseif Open_RyR(time) < 4.5 && t1 ~= 0
                t2(time) = 1;
                t1 = 0;
            end
            
        end
        
    end

    SumRyR = sum(Open_RyR);
    SumSparks = sum(t2);
    F = SumSparks/Iter;

    if plotthis == 1

        g = figure("visible", "off");
        g_name = sprintf('OpenRyR_%i_%i.png', vchange, run);
        plot(times, Open_RyR, '-');
        ylabel('Number of Open RyRs');
        xlabel('Time (us)');
        title('OpenRyR over Time');
        axis([1 Iter 0 num_RyR]);
        saveas(g, [pwd '/' p_name '/' g_name]);
        close(g);

        h = figure("visible", "off");
        h_name = sprintf('[Ca2+]_%i_%i.png', vchange, run);
        plot(times, concentration, '-');
        ylabel('Concentration of Calcium Ions in Cytosol (uM)');
        xlabel('Time (us)');
        title('Concentration vs Time');
        saveas(h, [pwd '/' p_name '/' h_name]);
        close(h);

    end

    save ("-text", text_name, "F", "SumSparks", "SumRyR");

endfunction
