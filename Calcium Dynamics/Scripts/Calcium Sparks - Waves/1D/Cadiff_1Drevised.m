%% Cadiff_1Drevised.m
% Simulates CaSparks and Waves in 1D

format long e
plot_num = 0;

%% initialization of variables
a = 10000; % number of atoms, will become +1 due to uniform spacing.
B = 100; % Boundary
N = 10000; % number of Iteration
inc = (2*B)/a; % increment between atoms
x = -B:inc:B; %a atoms equally spanning boundary
a = length(x);
c = zeros(a,1);
c(1:end) = NaN;
x = round(x);

%% Begin Probability function to evaluate the probability function to see if a firing occurs or not
Pmax = 0.3; %maximum probability
Pmin = 0.0001;%minimum probability % start zero
hillcoeff = 10; %hill coefficient (how steply the curve changes)
Phalf = 100; %the half way point (where the curve starts changing concavity)

% Stochastic Probability Function
fprob = @(n_ions) Pmin+Pmax*n_ions.^hillcoeff./(n_ions.^hillcoeff+Phalf.^hillcoeff);

% Deterministic Probability Function
%fprob = @(n_ions) 1*(n_ions>5); 

%% Determining number of openings from SR to Cytosol, and positions of CRUs
nopen = 20; % not including 0; 
div = (2*B)/nopen; % divisor

%Allocating position of CRUs
cru = zeros((nopen+1),1); % to be number of openings + location at 0.
len = length(cru);
cru(1:end) = -B:div:B;

%% initialization of states
% 0 = closed; 1 = open; -1 = inactive. REFER to Markov Model
s = zeros((nopen+1),1);
s(1) = 1;
t_inactive = zeros((nopen+1),1);
t_open = zeros((nopen+1),1);
t_open(1) = 1;

tau_open = 50; % time constant for how long open state will be
tau_inactive = 200; % time constant for how long inactive state will be

%% main algorithm
for i = 1:N
    
    % check states % consider increments 0 closed, 1-5 open,.. consider changing to get rid of ifelse statements 
    nc = 0;
    for j = 1:len
        
        prand = rand(1);
        nc = sum(c>(cru(j)-(div)) & c<(cru(j)+(div)));
  
        if s(j) == 0

            actual_prob = fprob(nc);

            if prand < actual_prob
                s(j) = 1;
                t_open(j) = i;
            end
            
        elseif s(j) == 1;
            
            if i >= (t_open(j) + tau_open)
                t_inactive(j) = i;
                s(j) = -1;
            end
            
        elseif s(j) == -1
            
            if i >= t_inactive(j) + tau_inactive
                s(j) = 0; 
            end

        end

    end

    % random walk simulation
    w = randi([-1 1], a,1)';
    x = x+w;

    cl = length(c);
    cw = randi([-1 1], cl,1);
    c = c+cw;

    
    % boundary conditions
    x(x<-B) = -B + abs(-B-x(x<-B));
    x(x>B) = B - abs(x(x>B)-B);

    c(c<-B) = -B + abs(-B-c(c<-B));
    c(c>B) = B - abs(c(c>B)-B);
    
    % find ions at CRUs
    num_ions = zeros(len,1);
    
    for j = 1:len
        
        num_ions(j) = sum(x==cru(j)); 
        % change position of x to NaN
        % if position is NaN, make c to have an ion at that index.

        if s(j) == 1 && num_ions(j)>0
            c(x==cru(j))= x(x==cru(j));
            x(x==cru(j)) = NaN;
        end
    end
    
    % uptake portion.
    % Please change this section of the script in order to alter the uptake strength.
    
    % uptake < 0.01 for 1% uptake, should see calcium waves
    % uptake < 0.1 for 10% uptake, should mimic an Over Active SERCA pump
    % uptake = 0 for 0% uptake, should see a deterministic wave because there is no uptake.


    for k = 1:a
        
        uptake = rand(1);
        
        if uptake< 0.01

            if isnan(c(k))==false
                
                x(k) = c(k);
                c(k) = NaN;
            
            end
        end
    end
    
    % keeping track of the states. When they are close, open or inactive. 
    sclose = cru.*(s==0);
    sopen = cru.*(s>0);
    sinactive = cru.*(s<0);
    
    % plot simulation
    if mod(i,200)==0
        %figure(1);
        f = figure('visible','on'); % either close the figure or reuse.
		plot(x,0, 'ro');
		hold on
		plot(c,a/8, 'bo');
		hold on
		plot(cru, 0, 'ks');
        hold on
        plot(sclose, a/10, 'rs');
        hold on
        plot(sinactive, a/10, 'bs');
        hold on
        plot(sopen, a/10, 'gs');
        hold on
       
		title(['time t=' num2str(i)]);
		%legend('Ca Ions in SR', 'Ca Ions in Cytosol', 'CRUs');
		axis([-110 110 -0.5 a/4]);
		axis square
		histogram(c,'BinWidth', div);
        %histogram(x,'BinWidth', 25);
		%pause(0.01);
		hold off
        f_name = sprintf('fsave%.4d',plot_num);
        plot_num = plot_num+1; 
        saveas(f,f_name,'png');
        close(f);
    end
    
end

