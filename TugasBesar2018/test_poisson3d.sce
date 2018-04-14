exec("poisson3d.sce",-1)

function u = gauss(x,y,z,sigma,R0)
  dx2 = (x - R0(1))^2
  dy2 = (y - R0(2))^2
  dz2 = (z - R0(3))^2
  dr2 = dx2 + dy2 + dz2
  C = sqrt(2*%pi*sigma^2)^3
  u = exp(-dr2/(2*sigma^2))/C
endfunction

Nx = 35
Ny = 35
Nz = 35
x = linspace(-8,8,Nx)
y = linspace(-8,8,Ny)
z = linspace(-8,8,Nz)

u0 = zeros(Nx,Ny,Nz)
// set BC, not necessary actually
u0(1,:,:) = 0.0
u0(:,1,:) = 0.0
u0(:,:,1) = 0.0
u0(Nx,:,:) = 0.0
u0(:,Ny,:) = 0.0
u0(:,:,Nz) = 0.0

// calculate array for RHS
sigma1 = 0.75
sigma2 = 0.50
R0 = [0.0,0.0,0.0]
rho = zeros(Nx,Ny,Nz)

t1 = timer()
for k = 1:Nz
for j = 1:Ny
for i = 1:Nx
  g1 = gauss(x(i),y(j),z(k),sigma1,R0)
  g2 = gauss(x(i),y(j),z(k),sigma2,R0)
  rho(i,j,k) = g2 - g1
end
end
end
t2 = timer()

dVol = 16.0^3/((Nx-1)*(Ny-1)*(Nz-1))
printf("integRho = %18.10f\n", sum(rho)*dVol)

VHartree = poisson3d( u0, x,y,z, Nx,Ny,Nz, 1e-5, -4.0*%pi*rho )

EHartree = 0.5*sum(VHartree.*rho)*dVol
EHartree_analytic = ((1/sigma1+1/sigma2)/2-sqrt(2)/sqrt(sigma1^2+sigma2^2))/sqrt(%pi)
