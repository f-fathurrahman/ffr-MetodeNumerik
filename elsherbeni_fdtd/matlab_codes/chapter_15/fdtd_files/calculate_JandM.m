% Calculate J and M on the imaginary farfiled surface
if number_of_farfield_frequencies == 0
    return;
end
j = sqrt(-1);

tmyxp(1,1,:,:) =  0.5 * (Ez (ui,lj:uj-1,lk:uk-1) + Ez (ui,lj+1:uj,lk:uk-1));
tmzxp(1,1,:,:) = -0.5 * (Ey (ui,lj:uj-1,lk:uk-1) + Ey (ui,lj:uj-1,lk+1:uk));
tmxyp(1,:,1,:) = -0.5 * (Ez (li:ui-1,uj,lk:uk-1) + Ez (li+1:ui,uj,lk:uk-1));
tmzyp(1,:,1,:) =  0.5 * (Ex (li:ui-1,uj,lk:uk-1) + Ex (li:ui-1,uj,lk+1:uk));
tmxzp(1,:,:,1) =  0.5 * (Ey (li:ui-1,lj:uj-1,uk) + Ey (li+1:ui,lj:uj-1,uk));
tmyzp(1,:,:,1) = -0.5 * (Ex (li:ui-1,lj:uj-1,uk) + Ex (li:ui-1,lj+1:uj,uk));

tjyxp(1,1,:,:) = -(txp_h1 * (Hz(ui,lj:uj-1,lk:uk-1)   + Hz(ui,lj:uj-1,lk+1:uk)) ... 
                 + txp_h2 * (Hz(ui-1,lj:uj-1,lk:uk-1) + Hz(ui-1,lj:uj-1,lk+1:uk)));

tjzxp(1,1,:,:) = txp_h1 * (Hy(ui,lj:uj-1,lk:uk-1)   + Hy(ui,lj+1:uj,lk:uk-1)) ... 
               + txp_h2 * (Hy(ui-1,lj:uj-1,lk:uk-1) + Hy(ui-1,lj+1:uj,lk:uk-1));

tjzyp(1,:,1,:) = -(typ_h1 * (Hx (li:ui-1,uj,lk:uk-1)  + Hx(li+1:ui,uj,lk:uk-1)) ... 
                 + typ_h2 * (Hx(li:ui-1,uj-1,lk:uk-1) + Hx(li+1:ui,uj-1,lk:uk-1)));

tjxyp(1,:,1,:) = typ_h1 * (Hz(li:ui-1,uj,lk:uk-1)   + Hz(li:ui-1,uj,lk+1:uk)) ... 
               + typ_h2 * (Hz(li:ui-1,uj-1,lk:uk-1) + Hz(li:ui-1,uj-1,lk+1:uk));

tjyzp(1,:,:,1) = tzp_h1 * (Hx(li:ui-1,lj:uj-1,uk)   + Hx(li+1:ui,lj:uj-1,uk)) ... 
               + tzp_h2 * (Hx(li:ui-1,lj:uj-1,uk-1) + Hx (li+1:ui,lj:uj-1,uk-1));

tjxzp(1,:,:,1) = -(tzp_h1 * (Hy(li:ui-1,lj:uj-1,uk)   + Hy(li:ui-1,lj+1:uj,uk)) ... 
                 + tzp_h2 * (Hy(li:ui-1,lj:uj-1,uk-1) + Hy(li:ui-1,lj+1:uj,uk-1)));

tmyxn(1,1,:,:) = -0.5 * (Ez (li,lj:uj-1,lk:uk-1) + Ez (li,lj+1:uj,lk:uk-1));
tmzxn(1,1,:,:) =  0.5 * (Ey (li,lj:uj-1,lk:uk-1) + Ey (li,lj:uj-1,lk+1:uk));

tmxyn(1,:,1,:) =  0.5 * (Ez (li:ui-1,lj,lk:uk-1) + Ez (li+1:ui,lj,lk:uk-1));
tmzyn(1,:,1,:) = -0.5 * (Ex (li:ui-1,lj,lk:uk-1) + Ex (li:ui-1,lj,lk+1:uk));

tmxzn(1,:,:,1) = -0.5 * (Ey (li:ui-1,lj:uj-1,lk) + Ey (li+1:ui,lj:uj-1,lk));
tmyzn(1,:,:,1) =  0.5 * (Ex (li:ui-1,lj:uj-1,lk) + Ex (li:ui-1,lj+1:uj,lk));

tjyxn(1,1,:,:) = txn_h1 * (Hz(li,lj:uj-1,lk:uk-1)   + Hz(li,lj:uj-1,lk+1:uk)) ... 
               + txn_h2 * (Hz(li-1,lj:uj-1,lk:uk-1) + Hz(li-1,lj:uj-1,lk+1:uk));

tjzxn(1,1,:,:) = -(txn_h1 * (Hy(li,lj:uj-1,lk:uk-1)   + Hy(li,lj+1:uj,lk:uk-1)) ... 
                 + txn_h2 * (Hy(li-1,lj:uj-1,lk:uk-1) + Hy(li-1,lj+1:uj,lk:uk-1)));

tjzyn(1,:,1,:) = tyn_h1 * (Hx(li:ui-1,lj,lk:uk-1)   + Hx(li+1:ui,lj,lk:uk-1)) ... 
               + tyn_h2 * (Hx(li:ui-1,lj-1,lk:uk-1) + Hx(li+1:ui,lj-1,lk:uk-1));

tjxyn(1,:,1,:) = -(tyn_h1 * (Hz(li:ui-1,lj,lk:uk-1)   + Hz(li:ui-1,lj,lk+1:uk)) ... 
                 + tyn_h2 * (Hz(li:ui-1,lj-1,lk:uk-1) + Hz(li:ui-1,lj-1,lk+1:uk)));

tjyzn(1,:,:,1) = -(tzn_h1 * (Hx(li:ui-1,lj:uj-1,lk)   + Hx(li+1:ui,lj:uj-1,lk)) ... 
                 + tzn_h2 * (Hx(li:ui-1,lj:uj-1,lk-1) + Hx (li+1:ui,lj:uj-1,lk-1)));

tjxzn(1,:,:,1) = tzn_h1 * (Hy(li:ui-1,lj:uj-1,lk)   + Hy(li:ui-1,lj+1:uj,lk)) ... 
               + tzn_h2 * (Hy(li:ui-1,lj:uj-1,lk-1) + Hy(li:ui-1,lj+1:uj,lk-1));

% fourier transform
for mi=1:number_of_farfield_frequencies
    exp_h = dt * exp(-j*farfield_w(mi)*(time_step-0.5)*dt);
    cjxyp(mi,:,:,:) = cjxyp(mi,:,:,:) + exp_h * tjxyp(1,:,:,:) ; 
    cjxzp(mi,:,:,:) = cjxzp(mi,:,:,:) + exp_h * tjxzp(1,:,:,:); 
    cjyxp(mi,:,:,:) = cjyxp(mi,:,:,:) + exp_h * tjyxp(1,:,:,:); 
    cjyzp(mi,:,:,:) = cjyzp(mi,:,:,:) + exp_h * tjyzp(1,:,:,:); 
    cjzxp(mi,:,:,:) = cjzxp(mi,:,:,:) + exp_h * tjzxp(1,:,:,:); 
    cjzyp(mi,:,:,:) = cjzyp(mi,:,:,:) + exp_h * tjzyp(1,:,:,:); 

    cjxyn(mi,:,:,:) = cjxyn(mi,:,:,:) + exp_h * tjxyn(1,:,:,:); 
    cjxzn(mi,:,:,:) = cjxzn(mi,:,:,:) + exp_h * tjxzn(1,:,:,:); 
    cjyxn(mi,:,:,:) = cjyxn(mi,:,:,:) + exp_h * tjyxn(1,:,:,:); 
    cjyzn(mi,:,:,:) = cjyzn(mi,:,:,:) + exp_h * tjyzn(1,:,:,:); 
    cjzxn(mi,:,:,:) = cjzxn(mi,:,:,:) + exp_h * tjzxn(1,:,:,:); 
    cjzyn(mi,:,:,:) = cjzyn(mi,:,:,:) + exp_h * tjzyn(1,:,:,:); 

    exp_e = dt * exp(-j*farfield_w(mi)*time_step*dt);
    
    cmxyp(mi,:,:,:) = cmxyp(mi,:,:,:) + exp_e * tmxyp(1,:,:,:); 
    cmxzp(mi,:,:,:) = cmxzp(mi,:,:,:) + exp_e * tmxzp(1,:,:,:); 
    cmyxp(mi,:,:,:) = cmyxp(mi,:,:,:) + exp_e * tmyxp(1,:,:,:); 
    cmyzp(mi,:,:,:) = cmyzp(mi,:,:,:) + exp_e * tmyzp(1,:,:,:); 
    cmzxp(mi,:,:,:) = cmzxp(mi,:,:,:) + exp_e * tmzxp(1,:,:,:); 
    cmzyp(mi,:,:,:) = cmzyp(mi,:,:,:) + exp_e * tmzyp(1,:,:,:); 

    cmxyn(mi,:,:,:) = cmxyn(mi,:,:,:) + exp_e * tmxyn(1,:,:,:); 
    cmxzn(mi,:,:,:) = cmxzn(mi,:,:,:) + exp_e * tmxzn(1,:,:,:); 
    cmyxn(mi,:,:,:) = cmyxn(mi,:,:,:) + exp_e * tmyxn(1,:,:,:); 
    cmyzn(mi,:,:,:) = cmyzn(mi,:,:,:) + exp_e * tmyzn(1,:,:,:); 
    cmzxn(mi,:,:,:) = cmzxn(mi,:,:,:) + exp_e * tmzxn(1,:,:,:); 
    cmzyn(mi,:,:,:) = cmzyn(mi,:,:,:) + exp_e * tmzyn(1,:,:,:); 
end 
