close all;
clear all;
% Different Domain
Dest_Directory = sprintf('Images');
plot_num = 1;
mkdir(Dest_Directory);

IonInfo = load('IonInfo.txt');
numIons = IonInfo(1,1);
TotalTime = IonInfo(1,2);
numRyRs = IonInfo(1,3);
sr_xrange = IonInfo(2,1);
sr_yrange = IonInfo(2,2);
c_xrange = IonInfo(2,3);
c_yrange = IonInfo(3,1);
sizex = IonInfo(3,2);
sizey = IonInfo(3,3);

dx = 1/sizex;
dy = 1/sizey;

% fopen
SR_file = fopen('LocalSR.dat');
%Cyto_file = fopen('Cyto_concentration.dat');
RyR_location = fopen('RyRLocation.dat');
local_cyto = fopen('LocalCytosol.dat');
bufferfile = fopen('Buffer.dat');
dyadFile = fopen('dyad.dat');
StateFile = fopen('State.dat');

%SR_con = fread(SR_file, [sizex*sizey TotalTime], 'double');
SR_con = fread(SR_file, [1 TotalTime], 'double');
%Cyto_con = fread(Cyto_file, [sizex*sizey TotalTime], 'double');
ryr = fread(RyR_location,[numRyRs 2],'double');
cytochange = fread(local_cyto, [1 TotalTime], 'double');
buffer = fread(bufferfile, [1 TotalTime], 'double');
dyad = fread(dyadFile, [1 TotalTime], 'double');
openryr = dlmread('OpenRyR.txt');
state = load('StateFreq.txt');
States = fread(StateFile, [numRyRs TotalTime], 'int');


%Grid_SR = zeros(sizex,sizey);
%Grid_Cyto = zeros(sizex,sizey);
%Max_SR = max(max(SR_con));
%Max_Cyto = max(max(Cyto_con));
%zax = ones(length(ryr), 1)*(Max_Cyto)+1;
%t1 = 0;
%t2 = zeros(TotalTime,1);

 timevector = 1:TotalTime;
% %SumSparks = sum(t2);
% %Freq = SumSparks/TotalTime;
% 
FigureN = figure('visible', 'off');
f_name = sprintf('Concentration_Buffer.png');
plot(timevector/1000, cytochange, 'Linewidth', 2, 'DisplayName', '[Ca^{2+}]_{Cytosol}');
hold on
plot(timevector/1000, buffer, 'Linewidth', 2, 'DisplayName', '[CaT]');
hold on
plot(timevector/1000, dyad, 'Linewidth', 2, 'DisplayName', '[Ca^{2+}]_{cleft}');
hold on
plot(timevector/1000, SR_con, 'Linewidth', 2, 'DisplayName', '[Ca^{2+}]_{SR}')
title('Change in concentration over time');
xlabel('Time (ms)');
ylabel('Concentration (\mu M)');
legend show
saveas(FigureN, [pwd '/' Dest_Directory '/' f_name]);
close(FigureN);

FigureN = figure('visible', 'off');
f_name = sprintf('OpenRyR.png');
plot(timevector/1000, openryr, 'Linewidth', 2);
xlabel('Time (ms)');
ylabel('Number of Open RyRs');
title('Number of Open RyRs over time')
saveas(FigureN, [pwd '/' Dest_Directory '/' f_name]);
close(FigureN);

FigureN = figure('visible', 'off');
f_name = sprintf('States.png');
plot(timevector/1000, state(:,1), 'Linewidth', 2, 'DisplayName', 'State 0 [Closed]');
hold on
plot(timevector/1000, state(:,2), 'Linewidth', 2, 'DisplayName', 'State 1 [Open]');
hold on
plot(timevector/1000, state(:,3), 'Linewidth', 2, 'DisplayName', 'State 2 [Inactive]');
hold on
plot(timevector/1000, state(:,4), 'Linewidth', 2, 'DisplayName', 'State 3 [Inactive]');
title('Change of State Over Time');
xlabel('Time (ms)');
ylabel('Number of RyRs at State');
legend show
saveas(FigureN, [pwd '/' Dest_Directory '/' f_name]);
close(FigureN);
% %system(['ffmpeg -i ' Dest_Directory '\test%04d.png -c:v libx265 -preset placebo -y -r 10 ' Dest_Directory '\movie.mp4']);
