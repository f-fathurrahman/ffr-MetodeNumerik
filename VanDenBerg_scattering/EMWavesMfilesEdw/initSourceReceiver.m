function input = initSourceReceiver(input)
global nDIM;  

if nDIM == 1;      %------------------------------------------------------- 
            
   input.xS = -170;                                   % source position
   
   input.NR = 1;                                      % receiver positions
   input.xR(1,1) = -150;                                
   input.xR(1,2) =  150;             
                         
elseif nDIM == 2;  %-------------------------------------------------------  

   input.xS(1) = -170;                                % source position
   input.xS(2) =    0; 
   
   input.NR = 180;                                    % receiver positions
   input.rcvr_phi(1:input.NR) = (1:input.NR) * 2*pi/input.NR;
   input.xR(1,1:input.NR)     = 150 * cos(input.rcvr_phi);
   input.xR(2,1:input.NR)     = 150 * sin(input.rcvr_phi);
   
elseif nDIM == 3;  %-------------------------------------------------------

   input.xS(1) = -170;                                % source position
   input.xS(2) =    0; 
   input.xS(3) =    0;  
   
   input.NR = 180;                                    % receiver positions
   input.rcvr_phi(1:input.NR) = (1:input.NR) * 2*pi/input.NR;
   input.xR(1,1:input.NR)     = 150 * cos(input.rcvr_phi);
   input.xR(2,1:input.NR)     = 150 * sin(input.rcvr_phi);
   input.xR(3,1:input.NR)     = 0;
  
end % if
