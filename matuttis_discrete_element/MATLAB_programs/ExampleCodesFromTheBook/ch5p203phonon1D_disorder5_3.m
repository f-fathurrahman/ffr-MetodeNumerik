% PURPOSE: Computation of the numerical
%  dispersion relation for a linear chain 
%  in one dimension with disorder (p. 203)
% GEOMETRY:
%          +---+---+---+---+---+
% mass     1   2   3   4   5   6
% spring 1   2   3   4   5   6   7
% DIFFERENCE to ch5p201phonon1D_simplified_dos:
%  The wave bynber must be computed via fast
%  Fourier transform, and the points must be
%  reordered accordingly
% CAVEAT: If the disorder (constant A1,A2) is 
%  increased to unphysically large values 
%  (definitely larger than 0.3), some
%  observables which should be real become
%  constant and the results make no sense any
%  more
% REVISION HISTORY: 21-May-2014 H.-G. Matuttis
clear
format compact

n=200
clf
randn('seed',4)
A1=.3 % disorder for the mass 
A2=.15 % disorder for the force constant 
mass=ones(n+1,1)+A1*(randn(n+1,1)-.5);
k=ones(n+1,1)+A2*(rand(n+1,1)-.5);
k(n+1)=k(1);

% set up the dynamic matrix
for i=1:n
  ip=i+1;
  if (ip>n) % periodicity
    ip=ip-n;
  end

  D(i,i)=(k(i)+k(ip))/mass(i);
  D(ip,i)=-(k(i))/sqrt(mass(i)*mass(ip));
  D(i,ip)=-(k(ip))/sqrt(mass(i)*mass(ip));
end

% diagonalization of the dynamic matrix
[U,Deigval]=eig(D);

% computation of the wave number
% via fast Fourier transform
for i=1:n
  absfft=abs(fft(U(:,i)));
  [f,j]=max(absfft(1:n/2));
  kvec(i)=2*pi*j/n;
end


fullk1=[kvec'];
fullomega=[sqrt(diag(Deigval))'];% sqrt(diag(Deigval))'];


plot(fullk1,fullomega,'+')
xlabel('k')
ylabel('\omega')
axis tight
axis([0 pi 0 2.2])


return