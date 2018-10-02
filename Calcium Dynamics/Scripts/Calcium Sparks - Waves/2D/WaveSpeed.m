%% WaveSpeed.m
% Automates CaWave.m and finds the average propagation wave speed for
% specific parameter changes. Matlab script, with some changes in Octave

% Input: pcase - either 0,1,2,3 to indicate first variable change, vchange
% - denotes the variable change, run - rum number, mainly for
% organizational purposes.

function R = WaveSpeed(pcase, vchange, run)

xrange = 2.5; % um
yrange = 0.5; % um
dt = 1;
b = 0; %breaking point
r_move = sqrt(223*4*dt*10^-6); % um/us
uptake = 1; % uptake, 12 ions/us = 0.6uM/us, so 6 is 0.3uM, us.
iter = 1e6; %us
concentration = zeros(iter,1);
time_vector = 1:iter;
Open_RyR = zeros(iter,1);

switch pcase
    
    case 0 % Test Sensitivity of RyR
        Con = 500;
        dcru = 0.5;
        r_influence = 0.015*vchange; % um
        num_RyR = 30;
        p_name = 'Sensitivity';
        text_name = sprintf('%s_%i_%i', p_name,vchange,run);
        
    case 1 % Number of RyRs per CRU
        Con = 500;
        dcru = 0.5;
        r_influence = 0.015; % um
        num_RyR = vchange;
        p_name = 'NumRyR';
        text_name = sprintf('%s_%i_%i', p_name,vchange,run);
        
    case 2 % Distance between CRUs
        Con = 500;
        dcru = vchange;
        r_influence = 0.015; % um
        num_RyR = 30;
        p_name = 'dCRU';
        text_name = sprintf('%s_%i_%i', p_name,vchange,run);
        
    case 3 % [Ca2+] in the SR
        Con = vchange;
        dcru = 0.5;
        r_influence = 0.015; % um
        num_RyR = 30;
        p_name = 'SR[Ca2+]';
        text_name = sprintf('%s_%i_%i', p_name,vchange,run);
end

volume = 2*xrange*2*yrange*0.012;
AvoConstant = 6.022;
ConCoeff = volume*AvoConstant*100;
num_ions = Con*ConCoeff;
Xions = rand(num_ions,1)*2*xrange-xrange;
Yions = rand(num_ions,1)*2*yrange-yrange;
LocIons = zeros(num_ions,1); % Location of ions, 0 for SR & 1 for Cytosol
num_cluster = 9; % must be odd
%vector_clusters = zeros(num_cluster,1);
n = (num_cluster-1)/2;
vector_clusters = -n:1:n;

% generate multiple clusters
Identify_cluster = [];

for j = 1:num_cluster
    posRyR = CableCluster(num_RyR, r_influence, (vector_clusters(j)*dcru));
    Identify_cluster = [Identify_cluster; posRyR];
end

RuntoLength = length(Identify_cluster);
most_negative = (min(vector_clusters))*dcru;
most_positive = (max(vector_clusters))*dcru;

stopping_point = find((Identify_cluster(:,1)==most_positive) & (Identify_cluster(:,2)==0));
force_state = find((Identify_cluster(:,1)==most_negative) & (Identify_cluster(:,2)==0));

state = zeros(RuntoLength,1);
state(force_state) = 1;

tau_open = 5000;
tau_inactive = 100000;

for time = 1:iter
    
    p = rand(num_ions,1)*(2*pi);
    Xions = Xions + r_move*cos(p);
    Yions = Yions + r_move*sin(p);
    
    % Boundary Conditions
    Xions(Xions<=-xrange) = -xrange + abs(-xrange-Xions(Xions<=-xrange));
    Xions(Xions>=xrange) = xrange - abs(Xions(Xions>=xrange)-xrange);
    
    Yions(Yions<=-yrange) = -yrange + abs(-yrange-Yions(Yions<=-yrange));
    Yions(Yions>=yrange) = yrange - abs(Yions(Yions>=yrange)-(yrange));
    
    
    
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
    
    
    % generate random number for each RyR
    rand_num = rand(RuntoLength,1);
    
    
    % Begin checking for state of each RyR
    for i = 1:RuntoLength
        
        ind_C = LocIons==1; % indicates the location of ions in the cytosol to ind_C
        dist = (Identify_cluster(i,1)-Xions(ind_C)).^2 + (Identify_cluster(i,2)-Yions(ind_C)).^2 < r_influence^2;
        sum_ions = sum(dist(:));
        
        if state(i) == 0
            
            % NEW MODIFICATION:
            % we had a deterministic opening. Need to make it stochastic.
            
            tau_fire = 0.000001; % units should be in Probability per microsecond
            if sum_ions > 2 || rand_num(i) < (tau_fire)*dt
                
                state(i) = 1;
                %t_open(i) = time;
                
            end
            
        elseif state(i) == 1;
            
            if rand_num(i) < (1/tau_open)*dt
                state(i) = -1;
            end
            
        elseif state(i) == -1
            
            if rand_num(i) < (1/tau_inactive)*dt
                state(i) = 0;
            end
            
        end
        
        
        if state(i) == 1
            ind_SR = LocIons==0;
            move = (Xions - Identify_cluster(i,1)).^2 + (Yions - Identify_cluster(i,2)).^2 < r_influence^2 & ind_SR;
            LocIons(move) = 1;
            
        end
        
        if state(stopping_point)==1
            b = 1;
            break
        end
        
    end
    
    
    if mod(time,1) == 0
        
        SRIonsx = Xions(LocIons==0);
        SRIonsy = Yions(LocIons==0);
        CytosolIonsx = Xions(LocIons==1);
        CytosolIonsy = Yions(LocIons==1);
        
        Open_RyR(time) = sum(state(state==1));
        sumC = length(CytosolIonsx);
        concentration(time) = sumC/ConCoeff;
        
        
        if b == 1
            
            S = time;
            save ("-text", text_name, "S");
            
            f = figure('visible', 'off');
            plot(SRIonsx, SRIonsy + 3*yrange, 'r.');
            hold on
            plot(CytosolIonsx, CytosolIonsy, 'b.');
            hold on
            plot([-xrange,xrange,xrange, -xrange,-xrange],[-yrange,-yrange,yrange, yrange, -yrange],'k-','linewidth', 3);
            hold on
            plot(Identify_cluster(:,1), Identify_cluster(:,2), 'ko', 'MarkerSize', 3);
            hold on
            plot([-xrange,xrange,xrange, -xrange,-xrange],[-yrange,-yrange,yrange, yrange, -yrange]+3*yrange,'k-','linewidth', 3);
            title(['Wave Speed Time t=' num2str(time)]);
            hold on
            axis([-xrange xrange -yrange 4*yrange]);
            axis equal
            f_name = sprintf('Break_%s',text_name);
            saveas(f,f_name,'png');
            close(f);
            
            break
            
        end
        
    end
    
end

endfunction
