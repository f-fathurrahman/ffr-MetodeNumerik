% initialize farfield arrays
number_of_farfield_frequencies = size(farfield.frequencies,2);

if number_of_farfield_frequencies == 0
    return;
end
nc_farbuffer = farfield.number_of_cells_from_outer_boundary; 
li = nc_farbuffer + 1;
lj = nc_farbuffer + 1;
lk = nc_farbuffer + 1;
ui = nx - nc_farbuffer+1;
uj = ny - nc_farbuffer+1;
uk = nz - nc_farbuffer+1;

farfield_w = 2*pi*farfield.frequencies;

tjxyp = zeros(1,ui-li,1,uk-lk);
tjxzp = zeros(1,ui-li,uj-lj,1);
tjyxp = zeros(1,1,uj-lj,uk-lk);
tjyzp = zeros(1,ui-li,uj-lj,1);
tjzxp = zeros(1,1,uj-lj,uk-lk);
tjzyp = zeros(1,ui-li,1,uk-lk);
tjxyn = zeros(1,ui-li,1,uk-lk);
tjxzn = zeros(1,ui-li,uj-lj,1);
tjyxn = zeros(1,1,uj-lj,uk-lk);
tjyzn = zeros(1,ui-li,uj-lj,1);
tjzxn = zeros(1,1,uj-lj,uk-lk);
tjzyn = zeros(1,ui-li,1,uk-lk);
tmxyp = zeros(1,ui-li,1,uk-lk);
tmxzp = zeros(1,ui-li,uj-lj,1);
tmyxp = zeros(1,1,uj-lj,uk-lk);
tmyzp = zeros(1,ui-li,uj-lj,1);
tmzxp = zeros(1,1,uj-lj,uk-lk);
tmzyp = zeros(1,ui-li,1,uk-lk);
tmxyn = zeros(1,ui-li,1,uk-lk);
tmxzn = zeros(1,ui-li,uj-lj,1);
tmyxn = zeros(1,1,uj-lj,uk-lk);
tmyzn = zeros(1,ui-li,uj-lj,1);
tmzxn = zeros(1,1,uj-lj,uk-lk);
tmzyn = zeros(1,ui-li,1,uk-lk);

cjxyp = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
cjxzp = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cjyxp = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cjyzp = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cjzxp = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cjzyp = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
cjxyn = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
cjxzn = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cjyxn = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cjyzn = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cjzxn = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cjzyn = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
cmxyp = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
cmxzp = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cmyxp = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cmyzp = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cmzxp = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cmzyp = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
cmxyn = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
cmxzn = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cmyxn = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cmyzn = zeros(number_of_farfield_frequencies,ui-li,uj-lj,1);
cmzxn = zeros(number_of_farfield_frequencies,1,uj-lj,uk-lk);
cmzyn = zeros(number_of_farfield_frequencies,ui-li,1,uk-lk);
