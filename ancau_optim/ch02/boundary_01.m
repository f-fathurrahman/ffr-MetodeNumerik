clear all; close all;

xa = 0.0;
tau = 0.381966;
xb = tau;

fa = my_fun_01(xa);
fb = my_fun_01(xb);

fprintf('xa=%f  xb=%f  fa=%f  fb=%f\n', xa, xb, fa, fb);

if fa < fb
    left_bound = xa;
    right_bound = xb;
end

while fa > fb
    x1 = xb;
    f1 = fb;
    xb = (1 + tau)*x1 + tau*xa;
    fb = my_fun_01(xb);
    fprintf('xa=%f  xb=%f  fa=%f  fb=%f\n', xa, xb, fa, fb);
    if f1 >= fb
        xa = x1;
        fa = f1;
    elseif f1 < fb
        left_bound = xa;
        right_bound = xb;
    end
end

