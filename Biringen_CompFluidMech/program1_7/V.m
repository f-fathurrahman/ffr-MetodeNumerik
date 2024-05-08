function [V] = V(XI,YI);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF (4.2.11) *****

global I M G X Y;

V = 0.0;
for J=1:M
    if (J ~= I)
        XX = XI - X(J);
        V = V - G(J)*XX/(XX^2+(YI-Y(J))^2);
    end
end
    