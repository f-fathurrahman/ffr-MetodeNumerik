classdef riemann

methods (Static, Access = private)
  
  function v = trapz2(x1, x2, F)
    sy = zeros(size(x1));
    for i = 1:length(sy)
      sy(i) = trapz(x2, F(i,:));
    end;
    v = trapz(x1, sy);
  end % function trapz2

  function v = trapz3(x1, x2, x3, F)
    sz = zeros(length(x1), length(x2));
    for i = 1:length(x1)
      for j = 1:length(x2)
        sz(i,j) = trapz(x2, F(i,j,:));
      end
    end
    sy = zeros(size(x1));
    for i = 1:length(x1)
      sy(i) = trapz(x2, sz(i,:));
    end
    v = trapz(x1, sy);
  end % function trapz3

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
    case 3
      v = riemann.trapz3(p.x, p.y, p.z, F);
    otherwise
      v = [];
      disp('Error in arguments for trpzd')
    end % case
  end % function trpzd

  function v = mean_value(data)
    p = data.points;
    F = data.values;
    n = p.dim;
    switch n
    case 1
      v = mean(F)*p.measure;
    case 2
      v = mean(mean(F))*p.measure;
    case 3
      v = mean(mean(mean(F)))*p.measure;
    otherwise
      v = [];
      disp('Error in arguments for mean_value')
    end % case
  end % function mean_value

end % methods

end % classdef riemann

