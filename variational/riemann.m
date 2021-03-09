classdef riemann

methods (Static, Access = private)
  
  function v = trapz2(x1, x2, F)
    sy = zeros(size(x1));
    for i = 1:length(sy)
      sy(i) = trapz(x2, F(i,:));
    end;
    v = trapz(x1, sy);
  end % function

end % methods


methods (Static)

  function v = trpzd(data)
    p = data.points;
    F = data.values;
    n = p.dim;
    switch n
    case 1
      v = trapz(p.x, F);
    case 2
      v = riemann.trapz2(p.x, p.y, F);
    otherwise
      v = [];
      disp('Error in arguments for trpzd')
    end % case
  end

end % methods

end % classdef

