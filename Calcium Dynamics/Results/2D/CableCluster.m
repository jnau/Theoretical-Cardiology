% CableCluster.m
% Similar to SingleCluster.m, this script generates a 1D cable of CRUs


function pos = CableCluster(n, r_inf, loc)

  gridsize=100;
  grid=zeros(gridsize);
  grid(gridsize/2,gridsize/2)=1;

  prob_g=2.5*10^-2;
  prob_del=0.001;%prob to remove isolated ryr
  prob_r=0.999;

  while (sum(sum(grid))~=n)
  grid1=circshift(grid,[0 1]);
  grid2=circshift(grid,[0 -1]);
  grid3=circshift(grid,[1 0]);
  grid4=circshift(grid,[-1 0]);
  grid5=circshift(grid,[-1 -1]);
  grid6=circshift(grid,[-1 1]);

  %growth
  prob_g2=(grid1+grid2+grid3+grid4+grid5+grid6)*prob_g;
  r=rand(gridsize);
  grid(r<prob_g2)=1;

  grid1=circshift(grid,[0 1]);
  grid2=circshift(grid,[0 -1]);
  grid3=circshift(grid,[1 0]);
  grid4=circshift(grid,[-1 0]);
  grid5=circshift(grid,[-1 -1]);
  grid6=circshift(grid,[-1 1]);
  prob_del2=(6-grid1+grid2+grid3+grid4+grid5+grid6)*prob_del;

  r=rand(gridsize);
  grid(r<prob_del2)=0;

  r=rand(gridsize);
  grid(r>prob_r)=0;

  grid(gridsize/2,gridsize/2)=1;

  while (sum(sum(grid))>n)
    r=rand(gridsize);
    grid(r>prob_r)=0;
    grid(gridsize/2,gridsize/2)=1;
  end
  sum(sum(grid));
  end
  
  cnt=1;
  for lpx=1:gridsize
    for lpy=1:gridsize
      if grid(lpx,lpy)==1
        if mod(lpy,2)==1
          pos(cnt,1)=2*r_inf*(lpx-gridsize/2)+(r_inf+loc);
        else
          pos(cnt,1)=2*r_inf*(lpx-gridsize/2)+loc;
        end
        pos(cnt,2)=sqrt(3)*r_inf*(lpy-gridsize/2);
        cnt=cnt+1;
      end
    end
  end
    

end