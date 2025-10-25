%WUBSTEST generate label data
%   IFISS scriptfile: DJS; 29 December 2024
% Copyright (c) 2024 D.J. Silvester
global DELTA
DELTA=0;
diary wubstest32.txt
fprintf('Generating Wubs problem test data set ... ')
[RID,message]=fopen('wubs_results.txt','w');
fprintf(RID,'%33s\n','%------------ Wubs stretched grid32 label results');
test=0;
prandtlset=[1000];  %<--- Prandtl numbers
npr=length(prandtlset);
rayleighset=[2.87e9,2.89e9];  %<--- Rayleigh numbers
nra=length(rayleighset);
for dra=1:nra   %<--- Rayleigh numbers
   for dpr = 1:npr
      test=test+1;
      fprintf(['\n [%g]'],test)
      ra=rayleighset(dra);
      fprintf(['\n Rayleigh number is %g'],ra)
      pr=prandtlset(dpr);
      fprintf(['\n Prandtl number is %g\n'],pr)
      testproblem=['B-NS42_test',num2str(test)];
      batchfile=[testproblem,'_batch.m'];
      gohome, cd batchfiles
%---------- write the batch file
      [FID,message]=fopen(batchfile,'w');
      fprintf(FID,'2%%\n1%%\n0.051%%\n2%%\n32%%\n8%%\n2%%\n1.5%%\n1%%\n1%%\n%g%%\n%g%%\n',ra,pr);
      fprintf(FID,'0%%\n1e13%%\n10%%\n2e-5%%\n2%%\n30%%\n0%%\n0%%\n13%%\n');
      fprintf(FID,'%53s\n','%---------- data file for Wubs test problem');
      fclose(FID);
%---------- execute the batch file and compute the label
      batchmode(testproblem)
      load stabtrBouss_end
      kk=KE.*KE*0.051;
      figure(100+test),subplot(121)
      plot(time(120:2000),kk(120:2000),'.-r'), axis square
      subplot(122)
      plot(time(1200:2000),kk(1200:2000),'.-r'), axis square
%---------- check for oscillating solution
      meanE=mean(kk(1200:2000));
      maxE=max(kk(1200:2000));
      minE=min(kk(1200:2000));
      testE=(maxE-minE)/meanE;
      if testE<5e-4, label=0;
      else, label=1; end
%---------- check for incomplete time integration
      if 2001==length(time), flag=0;
      else, flag=1; end
%---------- write the result to the file
      fprintf(RID,'%g,%g,%g,%g,%g,%g\n',ra,pr,label,meanE,testE,flag);
   end
end
fclose(RID);
fprintf('\nAll done')
diary off
