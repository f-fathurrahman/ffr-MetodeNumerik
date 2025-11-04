% updating electric field components 
% associated with the inductors

for ind = 1:number_of_inductors
    fi = inductors(ind).field_indices;
    switch (inductors(ind).direction(1))
    case 'x' 
        Ex(fi) = Ex(fi) + inductors(ind).Cexj ...
            .* inductors(ind).Jix;
        inductors(ind).Jix = inductors(ind).Jix ...
            + inductors(ind).Cjex .* Ex(fi);
    case 'y'
        Ey(fi) = Ey(fi) + inductors(ind).Ceyj ...
            .* inductors(ind).Jiy;
        inductors(ind).Jiy = inductors(ind).Jiy ...
            + inductors(ind).Cjey .* Ey(fi);
    case 'z'
        Ez(fi) = Ez(fi) + inductors(ind).Cezj ...
            .* inductors(ind).Jiz;
        inductors(ind).Jiz = inductors(ind).Jiz ...
            + inductors(ind).Cjez .* Ez(fi);
    end
end
