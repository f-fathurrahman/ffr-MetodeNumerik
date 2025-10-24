#include <windows.h>
#include <GL/glew.h>

char input_file_name[200];
char output_file_name[200];

int maximum_threads_per_block;

int time_step = 0;

int computation_platform;
int number_of_time_steps;
int plotting_step;

bool is_TEz;
bool is_TMz;
bool show_Hz;
bool show_Ez;

float dt;
float dx;
float dy;
int nx;
int ny;
int nxx;
int nyy;
int number_of_cells;
int number_of_cells_with_pads;
float domain_min_x;
float domain_min_y;
float domain_max_x;
float domain_max_y;

int number_of_circles;
float* circles_center_x;
float* circles_center_y;
float* circles_radius;

int number_of_rectangles;
float* rectangles_min_x;
float* rectangles_min_y;
float* rectangles_max_x;
float* rectangles_max_y;

float *Ex;
float *Ey;
float *Ez;
float *Hx;
float *Hy;
float *Hz;

float *Cexe;
float *Cexhz;
float *Ceye;
float *Ceyhz;
float *Ceze;
float *Cezhy;
float *Cezhx;

float *Chxh;
float *Chxez;
float *Chyh;
float *Chyez;
float *Chzh;
float *Chzex;
float *Chzey;

int cpml_1d_array_size_xnyn;
int cpml_1d_array_size_xpyp;
int cpml_2d_array_size_xn;
int cpml_2d_array_size_yn;
int cpml_2d_array_size_xp;
int cpml_2d_array_size_yp; 

bool is_anyy_side_cpml;
bool is_cpml_xn;
bool is_cpml_xp;
bool is_cpml_yn;
bool is_cpml_yp;

int n_cpml_xn;
int n_cpml_xp;
int n_cpml_yn;
int n_cpml_yp;

float* cpml_a_ex_xn;
float* cpml_b_ex_xn;
float* cpml_a_ex_xp;
float* cpml_b_ex_xp;
float* cpml_a_ey_yn;
float* cpml_b_ey_yn;
float* cpml_a_ey_yp;
float* cpml_b_ey_yp;

float* cpml_a_mx_xn;
float* cpml_b_mx_xn;
float* cpml_a_mx_xp;
float* cpml_b_mx_xp;
float* cpml_a_my_yn;
float* cpml_b_my_yn;
float* cpml_a_my_yp;
float* cpml_b_my_yp;

float* Psi_eyx_xn;
float* CPsi_eyx_xn;
float* Psi_eyx_xp;
float* CPsi_eyx_xp;
float* Psi_exy_yn;
float* CPsi_exy_yn;
float* Psi_exy_yp;
float* CPsi_exy_yp;

float* Psi_ezx_xn;
float* CPsi_ezx_xn;
float* Psi_ezx_xp;
float* CPsi_ezx_xp;
float* Psi_ezy_yn;
float* CPsi_ezy_yn;
float* Psi_ezy_yp;
float* CPsi_ezy_yp;

float* Psi_hzx_xn;
float* CPsi_hzx_xn;
float* Psi_hzx_xp;
float* CPsi_hzx_xp;
float* Psi_hzy_yn;
float* CPsi_hzy_yn;
float* Psi_hzy_yp;
float* CPsi_hzy_yp;

float* Psi_hyx_xn;
float* CPsi_hyx_xn;
float* Psi_hyx_xp;
float* CPsi_hyx_xp;
float* Psi_hxy_yn;
float* CPsi_hxy_yn;
float* Psi_hxy_yp;
float* CPsi_hxy_yp;

int number_of_sampled_electric_fields;
int number_of_sampled_magnetic_fields;
int number_of_impressed_M;
int number_of_impressed_J;

char*  sampled_electric_fields_component;
int*   sampled_electric_fields_is;
int*   sampled_electric_fields_js;
float* sampled_electric_fields_sampled_value;

char*  sampled_magnetic_fields_component;
int*   sampled_magnetic_fields_is;
int*   sampled_magnetic_fields_js;
float* sampled_magnetic_fields_sampled_value;

int*   number_of_cej_components;
int    total_number_of_cej_components;
char*  impressed_J_direction;
int*   impressed_J_is;
int*   impressed_J_ie;
int*   impressed_J_js;
int*   impressed_J_je;
float*   impressed_J_min_x;
float*   impressed_J_min_y;
float*   impressed_J_max_x;
float*   impressed_J_max_y;

float* impressed_J_waveform;
float* impressed_J_Cej;
int*   J_indices;

int*   number_of_chm_components;
int    total_number_of_chm_components;
char*  impressed_M_direction;
int*   impressed_M_is;
int*   impressed_M_ie;
int*   impressed_M_js;
int*   impressed_M_je;
float*   impressed_M_min_x;
float*   impressed_M_min_y;
float*   impressed_M_max_x;
float*   impressed_M_max_y;
float* impressed_M_waveform;
float* impressed_M_Chm;
int*   M_indices;

// parameters on device

float *dvEx;
float *dvEy;
float *dvEz;
float *dvHx;
float *dvHy;
float *dvHz;

float *dvCexe;
float *dvCexhz;
float *dvCeye;
float *dvCeyhz;
float *dvCeze;
float *dvCezhy;
float *dvCezhx;

float *dvChxh;
float *dvChxez;
float *dvChyh;
float *dvChyez;
float *dvChzh;
float *dvChzex;
float *dvChzey;

float* dvimpressed_J_Cej;
int*   dvJ_indices;

float* dvimpressed_M_Chm;
int*   dvM_indices;

char*  dvsampled_electric_fields_component;
int*   dvsampled_electric_fields_is;
int*   dvsampled_electric_fields_js;
float* dvsampled_electric_fields_sampled_value;

char*  dvsampled_magnetic_fields_component;
int*   dvsampled_magnetic_fields_is;
int*   dvsampled_magnetic_fields_js;
float* dvsampled_magnetic_fields_sampled_value;

float* dvPsi_eyx_xn;
float* dvCPsi_eyx_xn;
float* dvPsi_eyx_xp;
float* dvCPsi_eyx_xp;
float* dvPsi_exy_yn;
float* dvCPsi_exy_yn;
float* dvPsi_exy_yp;
float* dvCPsi_exy_yp;

float* dvPsi_ezx_xn;
float* dvCPsi_ezx_xn;
float* dvPsi_ezx_xp;
float* dvCPsi_ezx_xp;
float* dvPsi_ezy_yn;
float* dvCPsi_ezy_yn;
float* dvPsi_ezy_yp;
float* dvCPsi_ezy_yp;

float* dvPsi_hzx_xn;
float* dvCPsi_hzx_xn;
float* dvPsi_hzx_xp;
float* dvCPsi_hzx_xp;
float* dvPsi_hzy_yn;
float* dvCPsi_hzy_yn;
float* dvPsi_hzy_yp;
float* dvCPsi_hzy_yp;

float* dvPsi_hyx_xn;
float* dvCPsi_hyx_xn;
float* dvPsi_hyx_xp;
float* dvCPsi_hyx_xp;
float* dvPsi_hxy_yn;
float* dvCPsi_hxy_yn;
float* dvPsi_hxy_yp;
float* dvCPsi_hxy_yp;
float* dvminimum_field_value;
float* dvmaximum_field_value;

int cpml_shift_xp;
int cpml_shift_yp;

float global_min_field = 1e9;
float global_max_field = -1e9;

GLuint objects_display_list;

