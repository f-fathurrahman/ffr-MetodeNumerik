function [] = saveRfl(Rfl)
global it       

if it == 0;   Rfl_0{1}   = real(Rfl{1});   
              Rfl_0{2}   = real(Rfl{2});   save Contrast0   Rfl_0;    end
          
if it == 32;  Rfl_32{1}  = real(Rfl{1});  
              Rfl_32{2}  = real(Rfl{2});   save Contrast32  Rfl_32;   end
          
if it == 64;  Rfl_64{1}  = real(Rfl{1});  
              Rfl_64{2}  = real(Rfl{2});   save Contrast64  Rfl_64;   end
          
if it == 128; Rfl_128{1} = real(Rfl{1});        
              Rfl_128{2} = real(Rfl{2});   save Contrast128 Rfl_128;  end
          
if it == 256; Rfl_256{1} = real(Rfl{1});  
              Rfl_256{2} = real(Rfl{2});   save Contrast256 Rfl_256;  end