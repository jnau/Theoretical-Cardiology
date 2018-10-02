%% CaWave.m
% Simulates a Calcium Wave in 2D, using Sparks.m as a base, we implemented
% multiple CRUs in order to simulate a wave.

%% Initializations
concentration = 500;
plot_num = 1;
xrange = 2.5; % um
yrange = 0.5; % um
dt = 1;
dcru = 0.5;
r_influence = 0.015; % um
r_move = sqrt(223*4*dt*10^-6); % um/us
uptake = 1; % uptake, 12 ions/us = 0.6uM/us, so 6 is 0.3uM, us.

Dest_Directory = sprintf('WaveVideo_%i', concentration); % Destination directory
mkdir(Dest_Directory);

% Workflow of converting from ions to concentration
% Volume (10^-15) * Avogadro's# (10^23) * [Ca2+](10^-6 M) = 10^2*[Ca2+] = #ions

volume = 2*xrange*2*yrange*0.012; % um3 == 10^-15 L
AvoConstant = 6.022;
ConCoeff = volume*AvoConstant*100;
num_ions = round(concentration*ConCoeff);

iter = 1e5; %us
concentration = zeros(iter,1);
time_vector = 1:iter;
Open_RyR = zeros(iter,1);

Xions = rand(num_ions,1)*2*xrange-xrange;
Yions = rand(num_ions,1)*2*yrange-yrange;
LocIons = zeros(num_ions,1); % Location of ions, 0 for SR & 1 for Cytosol
num_cluster = 9; % must be odd
%vector_clusters = zeros(num_cluster,1);
n = (num_cluster-1)/2;
vector_clusters = -n:1:n;

num_RyR = 30;

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

%% Main Algorithm
for time = 1:iter
    
    p = rand(num_ions,1)*(2*pi);
    Xions = Xions + r_move*cos(p);
    Yions = Yions + r_move*sin(p);
    
    % Boundary Conditions
    Xions(Xions<=-xrange) = -xrange + abs(-xrange-Xions(Xions<=-xrange));
    Xions(Xions>=xrange) = xrange - abs(Xions(Xions>=xrange)-xrange);
    Yions(Yions<=-yrange) = -yrange + abs(-yrange-Yions(Yions<=-yrange));
    Yions(Yions>=yrange) = yrange - abs(Yions(Yions>=yrange)-(yrange));

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
            
            tau_fire = 0.000001; % units should be in Probability per microsecond
            if sum_ions > 2 || rand_num(i) < (tau_fire)*dt
                state(i) = 1;
                %t_open(i) = time;
            end
            
        elseif state(i) == 1
            
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
    end
    
    if mod(time,1) == 0
        
        SRIonsx = Xions(LocIons==0);
        SRIonsy = Yions(LocIons==0);
        CytosolIonsx = Xions(LocIons==1);
        CytosolIonsy = Yions(LocIons==1);
        
        Open_RyR(time) = sum(state==1);
        sumC = length(CytosolIonsx);
        concentration(time) = sumC/ConCoeff;
        
        f = figure('visible', 'off');
        plot(SRIonsx, SRIonsy + 3*yrange, 'r.','DisplayName', 'Calcium Ion in SR', 'MarkerSize', 3);
        hold on
        plot(CytosolIonsx, CytosolIonsy, 'b.', 'DisplayName', 'Calcium Ion in Cytosol', 'MarkerSize', 3);
        hold on
        plot([-xrange,xrange,xrange, -xrange,-xrange],[-yrange,-yrange,yrange, yrange, -yrange],'k-','linewidth',  7, 'DisplayName', 'Sarcoplasmic Reticulum');
        hold on
        plot(Identify_cluster(:,1), Identify_cluster(:,2), 'ko', 'MarkerSize', 3, 'DisplayName', 'RyR');
        hold on
        plot([-xrange,xrange,xrange, -xrange,-xrange],[-yrange,-yrange,yrange, yrange, -yrange]+3*yrange,'k-','Color', [0.3 0.3 0.6],'linewidth', 7, 'DisplayName', 'Cytosol');
        title(['Calcium Wave, Time = ' num2str(time) ' us'] , 'FontSize', 14);
        hold on
        axis([-xrange xrange -yrange 4*yrange]);
        axis equal
        f_name = sprintf('test%.4d.png',plot_num);
        plot_num = plot_num+1;
        saveas(f, [pwd '/' Dest_Directory '/' f_name]);
        close(f);
        
    end
    
end