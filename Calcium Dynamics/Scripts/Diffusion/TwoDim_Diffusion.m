%% TwoDim_Diffusion.m
% Pre-simulation work. Simulates diffusion using random walks in 2D.

%% Initialization
num_ions = 1;
x = zeros(num_ions,1); % vector holding positions of ions along x-axis
y = zeros(num_ions,1); % vector holding positions of ions along y-axis
radius = 1; % length of walk; using polar coordinates, therefore specify radius

% Specify domain for simulation
Positive_Boundary = 5;
Negative_Boundary = -5;
time = 5; % total simulation time; arbitrary units

% Color legend to illustrate 2D walk step-by-step
twod_diagram = 1; % indicates plot for 2D walk
red = [1,0,0];
pink = [255, 250, 250]/255;
colors_p = [linspace(red(1), pink(1), time)', linspace(red(2), pink(2), time)', linspace(red(3), pink(3), time)'];
C = time;

%% Main Algorithm
for i = 1:time
    
    p = rand(num_ions,1)*(2*pi); % random theta
    x1 = x; % dummy vector
    y1 = y; % dummy vector
    
    % polar coordinates
    x = x + radius*cos(p);
    y = y + radius*sin(p);
    
    for j = 1:num_ions
        
        if x(j) >= Positive_Boundary
            diffx = abs(x(j)-Positive_Boundary);
            x(j) = Positive_Boundary - diffx;
            
        elseif x(j) <= Negative_Boundary
            diffx = abs(x(j)-Negative_Boundary);
            x(j) = Negative_Boundary + diffx;
            
        elseif y(j) >= Positive_Boundary
            diffy = abs(y(j)-Positive_Boundary);
            y(j) = Positive_Boundary - diffy;
            
        elseif y(j) <= Negative_Boundary
            diffy = abs(y(j)-Negative_Boundary);
            y(j) = Negative_Boundary + diffy;
        end
        
        x = round(x);
        y = round(y);
        
        x1= round(x1);
        y1 = round(y1);

       if twod_diagram == 1
            plot(0,0,'o', 'MarkerSize', 20, 'MarkerEdgeColor', 'k' , 'MarkerFaceColor', 'k');
            hold on
            plot(x,y, 'o', 'MarkerFaceColor',colors_p(C,:),'markersize',20, 'MarkerEdgeColor', 'k');
            hold on
            line([x x1], [y y1], 'linewidth', 3, 'LineStyle', '--');
            C = C-1;
            %axis([Bminus Bplus Bminus Bplus]);
            box off
            axis off
       end
       
    end
end