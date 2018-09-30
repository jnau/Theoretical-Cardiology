% plotv.m
% Given .dat, .txt file obtained from cell.cc, generate a movie using ffmpeg.
% Scripted by D. Sato, Ph.D, Edited by Jessica Au

%clear all;
tic
XN=400;%x-dim
YN=400;%y-dim
vmax=10;%v range
vmin=-80;%v range

if (YN/XN<720/1280)
intpres=XN/1280;
else
intpres=YN/720;
end

%%make colormap
c=colormap(hsv(48000));
c(32001:48000,:)=[];
c=flipud(c);
clen=(length(c)-1);

%%find data length
fid = fopen('vca50.dat','r');
fseek(fid,0,'eof');
t1=ftell(fid)/1/XN/YN;
fclose(fid);
vmm=load('vmmca50.txt');
vmmmax=vmm(:,1);
vmmmin=vmm(:,2);
t2=length(vmm);
t=min([t1 t2]);

%open file and prepare for movie
fidv = fopen('vca50.dat','r');
mov = VideoWriter('vtmp.avi','Uncompressed AVI');
open(mov);
[x,y]=meshgrid(1:XN,1:YN);
[x2,y2]=meshgrid(1:intpres:XN,1:intpres:YN);
g = ones(3,3)*1/9;

et = clock;
st = clock;%start time

for frm=1:t % frame
    
    [vz,count] = fread(fidv,[YN,XN],'uchar');
    
    if (etime(clock, et)>5)
        %[num2str(frm/t*100) '% done  ETA ' num2str(etime(clock, st)/frm*(t-frm)) ' sec']
        fprintf('%0.1f %% done  ETA %0.1f sec\n',(frm/t*100) , (etime(clock, st)/frm*(t-frm)));
        et = clock;
    end
    
    vz=vz*(vmmmax(frm)-vmmmin(frm))/255+vmmmin(frm);
    vz = interp2(x,y,vz,x2,y2,'cubic');
    vz=filter2(g, vz,'valid');
    vz=filter2(g, vz,'valid');
    
    %normalize
    vz(vz>vmax)=vmax;
    vz(vz<vmin)=vmin;
    vz=(vz-vmin)./(vmax-vmin)*clen+1;
    
    rgb = ind2rgb(uint16(flipud(vz)),c);
    
    F = im2frame(rgb);
    writeVideo(mov,F);

end

fclose(fidv);
close(mov);
comm=['start /wait /low ffmpeg -i vtmp.avi -vcodec libx264 -preset placebo -y vol.mp4'];
system(comm);
delete('vtmp.avi');

toc
close all;
