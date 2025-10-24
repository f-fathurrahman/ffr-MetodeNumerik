% update magnetic fields
current_time  = current_time + dt/2;
    
if is_TEz
Hz(1:nx,1:ny) = Chzh(1:nx,1:ny).* Hz(1:nx,1:ny) ...
    + Chzex(1:nx,1:ny) .* (Ex(1:nx,2:nyp1)-Ex(1:nx,1:ny))  ...
    + Chzey(1:nx,1:ny) .*(Ey(2:nxp1,1:ny)-Ey(1:nx,1:ny)); 
end

if is_TMz
    Hx(1:nxp1,1:ny) = Chxh(1:nxp1,1:ny) .* Hx(1:nxp1,1:ny) ...
        + Chxez(1:nxp1,1:ny) .* (Ez(1:nxp1,2:nyp1)-Ez(1:nxp1,1:ny)); 
 
    Hy(1:nx,1:nyp1) = Chyh(1:nx,1:nyp1) .* Hy(1:nx,1:nyp1) ...
        + Chyez(1:nx,1:nyp1) .* (Ez(2:nxp1,1:nyp1)-Ez(1:nx,1:nyp1)); 
end               
