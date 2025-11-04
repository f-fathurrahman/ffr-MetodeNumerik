function [node_coordinates_e, node_coordinates_h] ...
    = insert_subregion_into_domain ...
    (node_coordinates_e, node_coordinates_h, subregion, base_cell_size)

% adjust subregion cell size
subregion_length = ...
    subregion.region_end - subregion.region_start;
n_subregion_cells = ...
    round(subregion_length/subregion.cell_size);
subregion.cell_size = subregion_length/n_subregion_cells;

bcs  = base_cell_size;
srcs = subregion.cell_size;

n_nodes_e = size(node_coordinates_e,2); 

% transition before the subregion (indicated with bsr)
% adjust transition length
tr_start_position = (subregion.region_start - ...
    subregion.transition_length);
[mval, I] = min(abs(node_coordinates_e - tr_start_position));
I = find(node_coordinates_e == node_coordinates_e(I(1)));
tr_start_node = I;
tr_length = subregion.region_start ...
    - node_coordinates_e(tr_start_node);

% find number of cells for transition
R = (tr_length + bcs)/(tr_length + srcs);
number_of_transition_cells = ...
    floor(log10(bcs/srcs)/log10(R))-1;

% find transition cell sizes
transition_cell_sizes = ...
    srcs*R.^(1:number_of_transition_cells);

% adjust cell sizes so that total length is equal to
% transition length
transition_cell_sizes = transition_cell_sizes ...
    * (tr_length / sum(transition_cell_sizes));

bsr_transition_cell_sizes = fliplr(transition_cell_sizes);
bsr_tr_start_node = tr_start_node;

% transition after the subregion (indicated with asr)
% adjust transition length
tr_end_position = (subregion.region_end + ...
    subregion.transition_length);
[mval,I] = min(abs(node_coordinates_e - tr_end_position));
I = find(node_coordinates_e == node_coordinates_e(I(1)));
tr_end_node = I(1);
tr_length = node_coordinates_e(tr_end_node) ...
    - subregion.region_end;

% find number of cells for transition
R = (tr_length + bcs)/(tr_length+srcs);
number_of_transition_cells = ...
    floor(log10(bcs/srcs)/log10(R))-1;

% find transition cell sizes
transition_cell_sizes = ...
    srcs*R.^(1:number_of_transition_cells);

% adjust cell sizes so that total length is equal to
% transition length
transition_cell_sizes = transition_cell_sizes ...
    * (tr_length / sum(transition_cell_sizes));

asr_transition_cell_sizes = transition_cell_sizes;
asr_tr_end_node = tr_end_node;

% create cell size array for subregion
subregion_cell_sizes_e = ...
    srcs * ones(1, n_subregion_cells);

% adjust node coordinates
sr_all_cs_e = [bsr_transition_cell_sizes ...
    subregion_cell_sizes_e asr_transition_cell_sizes];

nodes_before_subregion_e = node_coordinates_e(1:bsr_tr_start_node);
nodes_after_subregion_e = node_coordinates_e(asr_tr_end_node:end);
a_node = node_coordinates_e(bsr_tr_start_node);

n_subregion_cells = size(sr_all_cs_e,2);
subregion_node_coordinates_e = zeros(1, n_subregion_cells-1);

for si = 1:n_subregion_cells-1
    a_node = a_node + sr_all_cs_e(si);
    subregion_node_coordinates_e(si) = a_node;
end

node_coordinates_e = [nodes_before_subregion_e ...
    subregion_node_coordinates_e nodes_after_subregion_e];

n_nodes = size(node_coordinates_e, 2);

% magnetic field components are placed at the centers of the cells
node_coordinates_h = ...
    0.5*(node_coordinates_e(1:n_nodes-1)+node_coordinates_e(2:n_nodes));

node_coordinates_h = [node_coordinates_e(1)-bcs/2 node_coordinates_h ...
    node_coordinates_e(n_nodes)+bcs/2];
                   
