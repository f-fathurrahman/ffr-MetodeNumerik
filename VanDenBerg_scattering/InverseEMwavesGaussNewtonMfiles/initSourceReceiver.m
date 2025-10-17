function input = initSourceReceiver(input)
                         
                 input.NS  = 40;                      % source positions
 input.src_phi(1:input.NS) = (1:input.NS) * 2*pi/input.NS;
    input.xS(1,1:input.NS) = 170 * cos(input.src_phi);
    input.xS(2,1:input.NS) = 170 * sin(input.src_phi);
   
                 input.NR  = 50;                      % receiver positions
input.rcvr_phi(1:input.NR) = (1:input.NR) * 2*pi/input.NR;
    input.xR(1,1:input.NR) = 150 * cos(input.rcvr_phi);
    input.xR(2,1:input.NR) = 150 * sin(input.rcvr_phi);  
