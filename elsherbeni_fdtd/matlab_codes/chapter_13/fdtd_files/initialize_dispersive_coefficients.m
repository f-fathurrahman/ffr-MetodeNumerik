% initialize Lorentz updating coefficients

for ind = 1:number_of_lorentz
    
    material = material_types(lorentz(ind).material);
   
    n_elements = size(lorentz(ind).Ex_indices,2);
    lorentz(number_of_lorentz).Qxp = zeros(1, n_elements); 
    lorentz(number_of_lorentz).Qx  = zeros(1, n_elements); 

    n_elements = size(lorentz(ind).Ey_indices,2);
    lorentz(number_of_lorentz).Qyp = zeros(1, n_elements); 
    lorentz(number_of_lorentz).Qy  = zeros(1, n_elements); 

    n_elements = size(lorentz(ind).Ez_indices,2);
    lorentz(number_of_lorentz).Qzp = zeros(1, n_elements); 
    lorentz(number_of_lorentz).Qz  = zeros(1, n_elements); 

    psi = eps_0 * material.lorentz_A ...
        * (material.lorentz_eps_s - material.eps_r) ...
        * material.lorentz_omega^2;

    den = material.lorentz_delta * dt + 1;
    lorentz(number_of_lorentz).Cqq =  (2 - dt^2 ...
        * material.lorentz_omega^2)/den;
    lorentz(number_of_lorentz).Cqqm =  ...
        (material.lorentz_delta * dt - 1)/den;
    lorentz(number_of_lorentz).Cqe =  ...
        dt^2 * psi / den;                    
end

