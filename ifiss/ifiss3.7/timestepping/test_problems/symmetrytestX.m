%SYMMETRYTEST generate label data
%   IFISS scriptfile: DJS; 24 December 2024
% Copyright (c) 2024 D.J. Silvester
global DELTA
diary flowtest5.txt
fprintf('Generating symmetric step flow test data set ... ')
[RID,message]=fopen('flow_results.txt','w');
fprintf(RID,'%33s\n','%------------ grid5 label results');
test=0;
deltaset=[1e-16];  %<--- perturbations to inflow
% [1e-16,1e-14,1e-12,1e-10,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3];
ndelta=length(deltaset);
for nu=[0.0047,0.0045]; %[0.00443:0.00001:0.00468]  %<--- Reynolds numbers
   for dk = 1:ndelta
      test=test+1;
      fprintf(['\n [%g]'],test)
      fprintf(['\n Reynolds number %g'],1/nu)
      delta=deltaset(dk); DELTA = delta;
      fprintf(['\n perturbation magnitude is %g\n'],DELTA)
      testproblem=['T-NS_test',num2str(test)];
      batchfile=[testproblem,'_batch.m'];
      gohome, cd batchfiles
%---------- write the batch file
      [FID,message]=fopen(batchfile,'w');
      fprintf(FID,'6\n16\n5\n1\n4\n%g\n%g\n1e15\n6\n3e-6\n2\n2500\n0\n13\n',nu,delta);
      fprintf(FID,'%53s\n','%---------- data file for symmetric step test problem');
      fclose(FID);
%---------- execute the batch file and compute the label
      batchmode(testproblem)
      load unsteadyrun
      load step_unsteadyflow
      [w]=step_vorticityplot(qmethod,U(:,end),By,Bx,G,xy,x,y,bound,1,12);
      [label,it,ib] = symstep_bdryvorticity(domain,qmethod,xy,bound,w,100+test,'or-');
%---------- check for incomplete time integration
      if nmax*200==length(time), flag=0;
      else, flag=1; end
%---------- write the result to the file
      fprintf(RID,'%g,%g,%g,%g,%g,%g\n',nu,delta,label,it,ib,flag);
   end
end
fclose(RID);
fprintf('\nAll done')
diary off
