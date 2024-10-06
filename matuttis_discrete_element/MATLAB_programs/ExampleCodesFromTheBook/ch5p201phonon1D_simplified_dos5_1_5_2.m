clear % 
format compact
% PURPOSE: Computation of the numerical
%  dispersion relation for a linear chain in one 
%  dimension (p. 199), and computation of the 
%  corresponding density of states (p. 200)
% GEOMETRY:
%          +---+---+---+---+---+
% mass     1   2   3   4   5   6
% spring 1   2   3   4   5   6   7
% REVISION HISTORY: 21-May-2014 H.-G. Matuttis

n=100 % number of masses
mass=ones(n+1,1);
k=ones(n+1,1); % force constant
k(n+1)=k(1); % periodicity

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

% normalization of the k-vectors/ wave-numbers
% with distribution to the correct Brillouin zones
kvec(1:2:n)=pi*[1:2:n]/n; 
kvec(2:2:n)=-pi*[2:2:n]/n; 

fullk1=[kvec'];
fullomega=[sqrt(diag(Deigval))'];% sqrt(diag(Deigval))'];


clf

% plot the dispersion relation
axes('position',[ 0.07  .45 .67 .5 ])
plot(fullk1,fullomega,'+')
xlabel('k')
ylabel('\omega')
axis tight
axis([-pi pi 0 2.2])
axes('position',[ 0.78  .45 .15 .5 ])

% computation of the density of states as
% histogram of the omega-values
[n,omega]=hist(fullomega,[0:.10:2.05]);

% addition of the density of states to the
% dispersion relation
barh(omega,n)
axis([0 1.15*max(n) 0 2.2])
xlabel('A(\omega)')
a=ylabel('\omega')
set(gca,'YAxisLocation','right')




return
