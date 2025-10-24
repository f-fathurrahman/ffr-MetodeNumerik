% update magnetic fields

current_time  = current_time + dt/2;
    
Hx = Chxh.*Hx+Chxey.*(Ey(1:nxp1,1:ny,2:nzp1)-Ey(1:nxp1,1:ny,1:nz)) ...
    + Chxez.*(Ez(1:nxp1,2:nyp1,1:nz)-Ez(1:nxp1,1:ny,1:nz)); 
                  
Hy = Chyh.*Hy+Chyez.*(Ez(2:nxp1,1:nyp1,1:nz)-Ez(1:nx,1:nyp1,1:nz)) ...
    + Chyex.*(Ex(1:nx,1:nyp1,2:nzp1)-Ex(1:nx,1:nyp1,1:nz)); 
               
Hz = Chzh.*Hz+Chzex.*(Ex(1:nx,2:nyp1,1:nzp1)-Ex(1:nx,1:ny,1:nzp1))  ...
    + Chzey.*(Ey(2:nxp1,1:ny,1:nzp1)-Ey(1:nx,1:ny,1:nzp1)); 
