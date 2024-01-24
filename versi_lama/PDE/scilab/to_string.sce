function str = to_string(i)
  if i < 10
    str = "000" + string(i)
  elseif i < 100
    str = "00" + string(i)
  elseif i < 1000
    str = "0" + string(i)
  else
    str = string(i)
  end
endfunction