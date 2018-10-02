%% OneDim_Diffusion.m
% Pre-simulation work. Illustrates diffusion implemented using random walks
% in 1D.

%% Initialization
numions = 100; % number of ions
boundary = 100; % boundary condition
Total_time = 1000; % time; arbitrary unit
dt = 0.1; % time step
ion_position = zeros(numions,1); % vector holding position of each ion
cytosol = NaN; % if ions are moved to the cytosol, we indicate them as such
t = zeros(numions,1); % time vector
%size(x);

%% Main Algorithm
for i = 1:dt:Total_time
    
    walk = randi([-2 2],numions,1); %slow diffusion; can move length of 2 steps
    %walk = randi([-10 10],a,1); %fast diffusion; can move length of 10 steps
    ion_position = ion_position+walk;
    
    for j = 1:numions
        
        if ion_position(j) >= boundary
            diffx = ion_position(j)-boundary;
            ion_position(j) = boundary - diffx;
        elseif ion_position(j) <= -boundary
            diffx = abs(ion_position(j)-(-boundary));
            ion_position(j) = -boundary + diffx;
        end
        
        
     %{ implementing a cytosol portion to simulation
        % Initialize number of openings and evenly distribute them
        nopen = 2;
        Opening = boundary/nopen;
       
        % If ion lands in opening, move them to cytosol.
        if ion_position(j) == Opening
            cytosol = [cytosol; ion_position(j)];
            ion_position(j) = NaN;
        elseif ion_position(j) == -Opening
            cytosol = [cytosol; ion_position(j)];
            ion_position(j) = NaN;
        end
    
        %}
    end
    
   %{
    num_cytosol = length(cytosol); % number of ions in cytosol
    
    % If there are ions in the cytosol, implement random walks on them.
    if num_cytosol > 1
        cwalk = randi([-2 2],num_cytosol,1 ); % walk for cytosolic ions
        cytosol = cytosol+cwalk;
    end
    %}
    
    % plot every time step. can be modified for efficiency
    plot(ion_position,0, 'ro');
    hold on
    %plot(cytosol,2, 'bo');
    histogram(ion_position, 10);
    axis([-100 100 -10 50]);
    axis square
    pause(0.1);
    hold off
    
end