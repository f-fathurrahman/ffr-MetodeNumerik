function [RES] = UV(T,VEC);

global I M G X Y;

XI = VEC(1);
UU = VEC(2);
YI = VEC(3);
VV = VEC(4);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF (4.2.10) *****

U = 0.0;
for J=1:M
    if (J ~= I)
        YY = YI - Y(J);
        U = U + G(J)*YY/((XI-X(J))^2+YY^2);
    end
end


% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF (4.2.11) *****

V = 0.0;
for J=1:M
    if (J ~= I)
        XX = XI - X(J);
        V = V - G(J)*XX/(XX^2+(YI-Y(J))^2);
    end
end

RES = [U; U; V; V];
    