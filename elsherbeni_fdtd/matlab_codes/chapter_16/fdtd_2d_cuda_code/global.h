#pragma once
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#define TILE_SIZE 16

#include <windows.h>
#include <GL/glew.h>
#include <cuda_gl_interop.h>

extern char input_file_name[200];
extern char output_file_name[200];

int copyFdtdArraysToGpuMemory();
bool check_device();
bool fdtd_time_marching_loop_on_gpu();
bool fdtdIterationOnGpu();
bool fdtdIterationOnCpu();
bool fetchResultsFromGpuMemory();
bool deallocateArrays();
bool deallocateCudaArrays();
void createColormapOnGpu();

void runIterationAndDisplay();
bool check_device();
bool read_project_parameters_from_file(int, char **);
bool initialize_sources();
bool saveSampledFieldsToFile();

extern int time_step;
extern int maximum_threads_per_block;

extern int computation_platform;
extern int number_of_time_steps;
extern int plotting_step;
extern bool is_TEz;
extern bool is_TMz;
extern bool show_Hz;
extern bool show_Ez;

extern float dt;
extern float dx;
extern float dy;
extern int nx;
extern int ny;
extern int nxx;
extern int nyy;
extern int number_of_cells;
extern int number_of_cells_with_pads;
extern float domain_min_x;
extern float domain_min_y;
extern float domain_max_x;
extern float domain_max_y;

extern int number_of_circles;
extern float* circles_center_x;
extern float* circles_center_y;
extern float* circles_radius;

extern int number_of_rectangles;
extern float* rectangles_min_x;
extern float* rectangles_min_y;
extern float* rectangles_max_x;
extern float* rectangles_max_y;

extern float *Ex;
extern float *Ey;
extern float *Ez;
extern float *Hx;
extern float *Hy;
extern float *Hz;

extern float *Cexe;
extern float *Cexhz;
extern float *Ceye;
extern float *Ceyhz;
extern float *Ceze;
extern float *Cezhy;
extern float *Cezhx;

extern float *Chxh;
extern float *Chxez;
extern float *Chyh;
extern float *Chyez;
extern float *Chzh;
extern float *Chzex;
extern float *Chzey;

extern int cpml_1d_array_size_xnyn;
extern int cpml_1d_array_size_xpyp;
extern int cpml_2d_array_size_xn;
extern int cpml_2d_array_size_yn;
extern int cpml_2d_array_size_xp;
extern int cpml_2d_array_size_yp; 

extern bool is_anyy_side_cpml;
extern bool is_cpml_xn;
extern bool is_cpml_xp;
extern bool is_cpml_yn;
extern bool is_cpml_yp;

extern int n_cpml_xn;
extern int n_cpml_xp;
extern int n_cpml_yn;
extern int n_cpml_yp;

extern float* cpml_a_ex_xn;
extern float* cpml_b_ex_xn;
extern float* cpml_a_ex_xp;
extern float* cpml_b_ex_xp;
extern float* cpml_a_ey_yn;
extern float* cpml_b_ey_yn;
extern float* cpml_a_ey_yp;
extern float* cpml_b_ey_yp;

extern float* cpml_a_mx_xn;
extern float* cpml_b_mx_xn;
extern float* cpml_a_mx_xp;
extern float* cpml_b_mx_xp;
extern float* cpml_a_my_yn;
extern float* cpml_b_my_yn;
extern float* cpml_a_my_yp;
extern float* cpml_b_my_yp;

extern float* Psi_eyx_xn;
extern float* CPsi_eyx_xn;
extern float* Psi_eyx_xp;
extern float* CPsi_eyx_xp;
extern float* Psi_exy_yn;
extern float* CPsi_exy_yn;
extern float* Psi_exy_yp;
extern float* CPsi_exy_yp;

extern float* Psi_ezx_xn;
extern float* CPsi_ezx_xn;
extern float* Psi_ezx_xp;
extern float* CPsi_ezx_xp;
extern float* Psi_ezy_yn;
extern float* CPsi_ezy_yn;
extern float* Psi_ezy_yp;
extern float* CPsi_ezy_yp;

extern float* Psi_hzx_xn;
extern float* CPsi_hzx_xn;
extern float* Psi_hzx_xp;
extern float* CPsi_hzx_xp;
extern float* Psi_hzy_yn;
extern float* CPsi_hzy_yn;
extern float* Psi_hzy_yp;
extern float* CPsi_hzy_yp;

extern float* Psi_hyx_xn;
extern float* CPsi_hyx_xn;
extern float* Psi_hyx_xp;
extern float* CPsi_hyx_xp;
extern float* Psi_hxy_yn;
extern float* CPsi_hxy_yn;
extern float* Psi_hxy_yp;
extern float* CPsi_hxy_yp;

extern int number_of_sampled_electric_fields;
extern int number_of_sampled_magnetic_fields;
extern int number_of_impressed_M;
extern int number_of_impressed_J;

extern char*  sampled_electric_fields_component;
extern int*   sampled_electric_fields_is;
extern int*   sampled_electric_fields_js;
extern float* sampled_electric_fields_sampled_value;

extern char*  sampled_magnetic_fields_component;
extern int*   sampled_magnetic_fields_is;
extern int*   sampled_magnetic_fields_js;
extern float* sampled_magnetic_fields_sampled_value;

extern int*   number_of_cej_components;
extern int    total_number_of_cej_components;
extern char*  impressed_J_direction;
extern int*   impressed_J_is;
extern int*   impressed_J_ie;
extern int*   impressed_J_js;
extern int*   impressed_J_je;
extern float*   impressed_J_min_x;
extern float*   impressed_J_min_y;
extern float*   impressed_J_max_x;
extern float*   impressed_J_max_y;
extern float* impressed_J_waveform;
extern float* impressed_J_Cej;
extern int*   J_indices;

extern int*   number_of_chm_components;
extern int    total_number_of_chm_components;
extern char*  impressed_M_direction;
extern int*   impressed_M_is;
extern int*   impressed_M_ie;
extern int*   impressed_M_js;
extern int*   impressed_M_je;
extern float*   impressed_M_min_x;
extern float*   impressed_M_min_y;
extern float*   impressed_M_max_x;
extern float*   impressed_M_max_y;
extern float* impressed_M_waveform;
extern float* impressed_M_Chm;
extern int*   M_indices;

// parameters on device

extern float *dvEx;
extern float *dvEy;
extern float *dvEz;
extern float *dvHx;
extern float *dvHy;
extern float *dvHz;

extern float *dvCexe;
extern float *dvCexhz;
extern float *dvCeye;
extern float *dvCeyhz;
extern float *dvCeze;
extern float *dvCezhy;
extern float *dvCezhx;

extern float *dvChxh;
extern float *dvChxez;
extern float *dvChyh;
extern float *dvChyez;
extern float *dvChzh;
extern float *dvChzex;
extern float *dvChzey;

extern float* dvimpressed_J_Cej;
extern int*   dvJ_indices;

extern float* dvimpressed_M_Chm;
extern int*   dvM_indices;

extern char*  dvsampled_electric_fields_component;
extern int*   dvsampled_electric_fields_is;
extern int*   dvsampled_electric_fields_js;
extern float* dvsampled_electric_fields_sampled_value;

extern char*  dvsampled_magnetic_fields_component;
extern int*   dvsampled_magnetic_fields_is;
extern int*   dvsampled_magnetic_fields_js;
extern float* dvsampled_magnetic_fields_sampled_value;

extern float* dvPsi_eyx_xn;
extern float* dvCPsi_eyx_xn;
extern float* dvPsi_eyx_xp;
extern float* dvCPsi_eyx_xp;
extern float* dvPsi_exy_yn;
extern float* dvCPsi_exy_yn;
extern float* dvPsi_exy_yp;
extern float* dvCPsi_exy_yp;

extern float* dvPsi_ezx_xn;
extern float* dvCPsi_ezx_xn;
extern float* dvPsi_ezx_xp;
extern float* dvCPsi_ezx_xp;
extern float* dvPsi_ezy_yn;
extern float* dvCPsi_ezy_yn;
extern float* dvPsi_ezy_yp;
extern float* dvCPsi_ezy_yp;

extern float* dvPsi_hzx_xn;
extern float* dvCPsi_hzx_xn;
extern float* dvPsi_hzx_xp;
extern float* dvCPsi_hzx_xp;
extern float* dvPsi_hzy_yn;
extern float* dvCPsi_hzy_yn;
extern float* dvPsi_hzy_yp;
extern float* dvCPsi_hzy_yp;

extern float* dvPsi_hyx_xn;
extern float* dvCPsi_hyx_xn;
extern float* dvPsi_hyx_xp;
extern float* dvCPsi_hyx_xp;
extern float* dvPsi_hxy_yn;
extern float* dvCPsi_hxy_yn;
extern float* dvPsi_hxy_yp;
extern float* dvCPsi_hxy_yp;

extern float* dvminimum_field_value;
extern float* dvmaximum_field_value;

extern int cpml_shift_xp;
extern int cpml_shift_yp;

extern float global_min_field;
extern float global_max_field;

extern GLuint objects_display_list;

#endif
