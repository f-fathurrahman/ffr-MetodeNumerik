function [U] = U(XI,YI);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF (4.2.10) *****

global I M G X Y;

U = 0.0;
for J=1:M
    if (J ~= I)
        YY = YI - Y(J);
        U = U + G(J)*YY/((XI-X(J))^2+YY^2);
    end
end
    