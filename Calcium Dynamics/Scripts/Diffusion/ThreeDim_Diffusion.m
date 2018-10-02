%% ThreeDim_Diffusion.m
% 3D simulation of random walk diffusion

%% Initialization
NumIons = 1; % number of atoms
Time = 5; % iterations
Boundary = 10; % boundary
%step = zeros(N,3,T);
x = zeros(NumIons,1); % x position of ions
y = zeros(NumIons,1); % y position of ions
z = zeros(NumIons,1); % z position of ions
r = 6; % radius length; spherical coordinates

% color legend for 3D illustration of random walk; for poster
red = [1,0,0];
pink = [255, 250, 250]/255;
colors_p = [linspace(red(1), pink(1), Time)', linspace(red(2), pink(2), Time)', linspace(red(3), pink(3), Time)'];
C = Time;

%% Main Algorithm
for i = 1:Time
    
    t = rand(NumIons,1)*(2*pi); % random theta
    p = rand(NumIons,1)*(pi); % random rho
    
    x1 = x; % dummy 
    y1 = y; % dummy
    z1 = z; % dummy
    
    % Random walk implemented using spherical coordinates
    x = x + r.*cos(t).*sin(p);
    y = y + r.*sin(t).*sin(p);
    z = z + r.*cos(p);
   
    % Enforcing boundary conditions; for ions to stay within specified
    % domain, reflective boundaries
    for j = 1:NumIons
        
        % diffx, diffy, diffz denote the distance from boundary to current
        % position of ion and will reflect it according to this distance.
        
        if x(j) >= Boundary
            diffx = abs(x(j)-Boundary);
            x(j) = Boundary - diffx;
        elseif x(j) <= -Boundary
            diffx = abs(x(j)-(-Boundary));
            x(j) = -Boundary + diffx;
        elseif y(j) >= Boundary
            diffy = abs(y(j)-Boundary);
            y(j) = Boundary- diffy;
        elseif y <= -Boundary
            diffy = abs(y-(-Boundary));
            y = -Boundary + diffy;
        elseif z(j) >= Boundary
            diffz = abs(z(j)-Boundary);
            z(j) = Boundary-diffz;
        elseif z(j) <= -Boundary
            diffz = abs(z(j)-(-Boundary));
            z(j) = -Boundary + diffz;
        end
        
    end
    
    % plot diagram
    x = round(x);
    y = round(y);
    z = round(z);
    x1 = round(x1);
    y1 = round(y1);
    z1 = round(z1);
    %plot3(x,y,z,'.', 'Markersize', '10', 'MarkerEdgeColor', 'r');
    plot3(0,0,0, 'o', 'MarkerSize', 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    hold on
    plot3(x,y,z, 'o', 'Markersize', 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor',colors_p(C,:));
    hold on
    line([x x1], [y y1], [z z1], 'linewidth', 3, 'LineStyle', '--');
    axis([-Boundary Boundary -Boundary Boundary -Boundary Boundary]);
    axis off
    C = C-1;
    ax = gca;
    get(ax, 'fontSize');
    set(ax, 'fontSize', 12);
    fh = figure(1);
    set(fh, 'color', 'white');
    F = getframe;
    
end

movie(F)
