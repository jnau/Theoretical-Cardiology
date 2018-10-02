%CalciumSpark_Simulation
% This script will simulate a calcium spark using the Shannon-Bers 2011
% Model

% Initialization of constants
xrange = 0.5; % um
yrange = 0.5; % um
plotthis = 1;
plot_num = 1;

dt = 1; % time step        
r_inf = 0.015; %um
Ca_move = sqrt(223*4*dt*10^-6); % um/us based on free Ca2+ diffusion coeff
uptake = 1;
init_calcium = 800; % uM
total_b = 123; % uM
iter = 1e4;
times = 1:iter;

Calcium_Concentration = zeros(iter,1);
Change_trpn = zeros(iter,1);
Change_Calcium = zeros(iter,1);
TrpnBuff_Concentration = zeros(iter,1);
TroponinBuff = 0;
   
Dest_Directory = sprintf('SparkVideo_%i', init_calcium);
mkdir(Dest_Directory);

k2_forward = 0.0001; % uM-1 us-1 converted from 100uM-1s-1
k2_backward = 0.0001; % us-1, converted from 100 s-1

% Convert concentration of B into species.
% For now we are making a big assumption, let us convert B into ions
% because there should be enough buffers to bind to the number of ions.

ConCoeff = 2*xrange*2*yrange*0.012*6.022*100;
num_ions = round(init_calcium*ConCoeff);

% Now randomly place them into their respective domains.
% We had LocIons == 0 for SR ions, LocIons == 1 for Cytosol, now LocIons ==
% 2 for bound ions to these buffers. 
Xions = rand(num_ions,1)*2*xrange-xrange;
Yions = rand(num_ions,1)*2*yrange-yrange;
LocIons = zeros(num_ions,1); 

% Initialize RyRs
num_RyR = 50; 
[posRyR, grid] = SingleCluster(num_RyR, r_inf);
OpenRyR = zeros(iter,1);
state = zeros(num_RyR,1);
Num_States = zeros(num_RyR,4);
state(1) = 1;

% RyR gating constants
k0Ca = 1e-8; % uM-2 us-1
kiCa = 5e-7; % uM-1us-1
SR_max = 15;
SR_min = 1;
EC_50SR = 450; %uM
H = 2.5;
kim = 5e-6; % us-1
kom = 6e-5; %us-1

k_CaSR = @(c) SR_max - ((SR_max-SR_min)/(1 + (EC_50SR/c)^H));
k0_SR = @(c) (k0Ca/k_CaSR(c^2));
ki_SR = @(c) kiCa*k_CaSR(c);

for time = 1:iter
    
    % After initialization, we perform random walks for ions and buffers.
    ion_move = rand(num_ions,1)*(2*pi);
    Xions(Xions~=2) = Xions(Xions~=2) + Ca_move*cos(ion_move);
    Yions(Yions~=2) = Yions(Yions~=2) + Ca_move*sin(ion_move);

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
    
    SumCyto = sum(LocIons==1);
    FindLocation = find(LocIons==1);
    
    Ca_Con = SumCyto/ConCoeff;
    d_TroponinBuff = k2_forward*Ca_Con*(total_b - TroponinBuff) - k2_backward*TroponinBuff;
    TroponinBuff = TroponinBuff+d_TroponinBuff;
    
    TrpnBuff_Concentration(time)= TroponinBuff;
    d_Caions = round(-TroponinBuff*ConCoeff);
    
    Change_Calcium(time) = -d_TroponinBuff;
    Change_trpn(time) = d_TroponinBuff;
    
    % if-else statement to add/remove ions
    if d_Caions < 0
        % remove ions
        if abs(d_Caions) < SumCyto
            R = randperm(SumCyto, abs(d_Caions));
            LocIons(FindLocation(R))= 2;
            %display('Removing...');
        else
            LocIons(FindLocation)=2;
            
        end
        
    elseif d_Caions > 0
        
        BoundIons = sum(LocIons==2);
        R = randperm(BoundIons, abs(d_Caions));
        LocIons(FindLocation(R))=1;
        %display('Unbinding...');
    end
    
    
    % generate random number for each RyR
    rand_num = rand(num_RyR,1);
    
    % Begin checking for state of each RyR
    for i = 1:num_RyR
        
        ind_C = LocIons==1; % indicates the location of ions in the cytosol to ind_C
        SumSR = sum(LocIons==0);
        dist = (posRyR(i,1)-Xions(ind_C)).^2 + (posRyR(i,2)-Yions(ind_C)).^2 < r_inf^2;
        sum_ions = sum(dist(:));
        %sum_ions = sum(LocIons==1);

        SR_calcium = SumSR/ConCoeff;
        ConCyto = sum_ions/ConCoeff;

        % Shannon model
        if state(i) == 0
            
            if (rand_num(i) < abs((k0_SR(SR_calcium)*ConCyto^2)-kom)*dt)
                state(i) = 1;
            elseif rand_num(i) < (kim -(ki_SR(SR_calcium)*ConCyto))*dt
                state(i) = 3;
            end
            
        elseif state(i) == 1
            
            if rand_num(i) < abs((ki_SR(SR_calcium)*ConCyto)-kim)*dt
                state(i) = 2;
            elseif rand_num(i) < ((k0_SR(SR_calcium)*ConCyto^2)-kom)*dt
                state(i) = 0;
            end
                
        elseif state(i) == 2
            
            if rand_num(i) < abs(kom-(k0_SR(SR_calcium)*ConCyto^2))*dt
                state(i) = 3;
            elseif rand_num(i) < ((ki_SR(SR_calcium)*ConCyto)-kim)*dt
                state(i) = 1;
            end
                    
        elseif state(i) == 3
            
            if rand_num(i) < abs(kim -(ki_SR(SR_calcium)*ConCyto))*dt
                state(i) = 0;
            elseif rand_num(i) < (kom-(k0_SR(SR_calcium)*ConCyto^2))*dt
                state(i) = 2;
            end
            
        end
        
        
        if state(i) == 1
            
            ind_SR = LocIons==0;
            move = (Xions - posRyR(i,1)).^2 + (Yions - posRyR(i,2)).^2 < r_inf^2 & ind_SR;
            LocIons(move) = 1;
            Num_States(i,2) = Num_States(i,2)+1;
            
        elseif state(i) == 0
            Num_States(i,1) = Num_States(i,1)+1;
        elseif state(i) == 2
            Num_States(i,3) = Num_States(i,3)+1;  
        elseif state(i) == 3
            Num_States(i,4) = Num_States(i,4)+1;
        end
        
    end
    
    %if mod(time,1) == 0
     if time == 1   
        SRIonsx = Xions(LocIons==0);
        SRIonsy = Yions(LocIons==0);
        CytosolIonsx = Xions(LocIons==1);
        CytosolIonsy = Yions(LocIons==1);
        BoundIonsx = Xions(LocIons==2);
        BoundIonsy = Yions(LocIons==2);
        
        f = figure('visible', 'off');
        plot(SRIonsx, SRIonsy + 3*yrange, 'r.', 'DisplayName', 'Calcium Ion in SR');
        hold on
        plot(CytosolIonsx, CytosolIonsy, 'b.', 'DisplayName', 'Calcium Ion in Cytosol');
        hold on
       % plot(BoundIonsx, BoundIonsy, 'gx');
        %hold on
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
    
    if mod(time,1) == 0
        
        SRIonsx = Xions(LocIons==0);
        SRIonsy = Yions(LocIons==0);
        CytosolIonsx = Xions(LocIons==1);
        CytosolIonsy = Yions(LocIons==1);
        BoundIonsx = Xions(LocIons==2);
        BoundIonsy = Yions(LocIons==2);
        
        OpenRyR(time) = sum(state==1);
        sumC = length(CytosolIonsx);
        Calcium_Concentration(time) = sumC/ConCoeff;

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
		plot(times, Calcium_Concentration, 'linewidth', 2);
		xlabel('time (us)');
		ylabel('[Ca^2^+] in Cytosol (uM)');
		saveas(h, [pwd '/' Dest_Directory '/' hname]);
		close(h);
        
        
        l = figure('visible', 'off');
		lname = 'BufferRelationship.png';
		plot(times, Calcium_Concentration, 'linewidth', 2);
        hold on
        plot(times, TrpnBuff_Concentration, 'linewidth', 2);
        hold off
		xlabel('time (us)');
		ylabel('solutions');
		saveas(l, [pwd '/' Dest_Directory '/' lname]);
		close(l);
        
	
	end
