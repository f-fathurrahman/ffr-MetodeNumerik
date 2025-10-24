% update electric fields at the CPML regions
if is_any_side_cpml == false
    return;
end
if is_TEz
    update_electric_fields_for_CPML_2d_TEz;
end
if is_TMz
    update_electric_fields_for_CPML_2d_TMz;
end
