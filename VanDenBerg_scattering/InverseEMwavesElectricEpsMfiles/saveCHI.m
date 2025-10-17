function [] = saveCHI(CHI)
global it       

if it ==   0;  CHI_0   = real(CHI);  save Contrast0   CHI_0;    end  
if it ==   1;  CHI_1   = real(CHI);  save Contrast1   CHI_1;    end   
if it ==   4;  CHI_4   = real(CHI);  save Contrast4   CHI_4;    end
if it ==   8;  CHI_8   = real(CHI);  save Contrast8   CHI_8;    end   
if it ==  16;  CHI_16  = real(CHI);  save Contrast16  CHI_16;   end       
if it ==  32;  CHI_32  = real(CHI);  save Contrast32  CHI_32;   end
if it ==  64;  CHI_64  = real(CHI);  save Contrast64  CHI_64;   end   
if it == 128;  CHI_128 = real(CHI);  save Contrast128 CHI_128;  end  
if it == 256;  CHI_256 = real(CHI);  save Contrast256 CHI_256;  end 
   