#include <stdio.h>
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>

#include "global.h"

using namespace std;


bool read_project_parameters_from_file(int argc, char **argv)
{

	int size_int   = sizeof(int);
	int size_char  = sizeof(char);
	int size_float = sizeof(float);

	ifstream input_file;
	if (argc != 2)
	{
		cout<<"Please rerun using an input file name: "<< argv[0] <<" <filename>\n";
		return false;
	}
	else 
	{
		strcpy(input_file_name, argv[1]);

		input_file.open(input_file_name, ios::in | ios::binary);
		if (!input_file.is_open())
		{
			cout<<"File <" << input_file_name << "> can not be opened! \n";
			return false;
		}
		else 
		{
			input_file.read((char*)&computation_platform, size_int);
			input_file.read((char*)&number_of_time_steps, size_int);
			input_file.read((char*)&plotting_step, size_int);
			input_file.read((char*)&dt, size_float);
			input_file.read((char*)&dx, size_float);
			input_file.read((char*)&dy, size_float);
			input_file.read((char*)&nx, size_int);
			input_file.read((char*)&ny, size_int);
			input_file.read((char*)&nxx, size_int);
			input_file.read((char*)&nyy, size_int);
			input_file.read((char*)&domain_min_x, size_float);
			input_file.read((char*)&domain_min_y, size_float);
			input_file.read((char*)&domain_max_x, size_float);
			input_file.read((char*)&domain_max_y, size_float);

			input_file.read((char*)&number_of_circles, size_int);
			circles_center_x = (float*) malloc(number_of_circles*size_float);
			circles_center_y = (float*) malloc(number_of_circles*size_float);
			circles_radius   = (float*) malloc(number_of_circles*size_float);
			for (int i=0;i<number_of_circles;i++)
			{
				input_file.read((char*)&circles_center_x[i], size_float);
				input_file.read((char*)&circles_center_y[i], size_int);
				input_file.read((char*)&circles_radius[i], size_int);
			}

			input_file.read((char*)&number_of_rectangles, size_int);
			rectangles_min_x = (float*) malloc(number_of_rectangles*size_float);
			rectangles_min_y = (float*) malloc(number_of_rectangles*size_float);
			rectangles_max_x = (float*) malloc(number_of_rectangles*size_float);
			rectangles_max_y = (float*) malloc(number_of_rectangles*size_float);
			for (int i=0;i<number_of_rectangles;i++)
			{
				input_file.read((char*)&rectangles_min_x[i], size_float);
				input_file.read((char*)&rectangles_min_y[i], size_float);
				input_file.read((char*)&rectangles_max_x[i], size_float);
				input_file.read((char*)&rectangles_max_y[i], size_float);
			}

			input_file.read((char*)&is_TEz, size_char);
			input_file.read((char*)&is_TMz, size_char);

			number_of_cells = nxx*nyy;
			int array_size  = number_of_cells * size_float; 
			number_of_cells_with_pads  = number_of_cells + 2*nxx; 
			int array_size_with_pads   = number_of_cells_with_pads * size_float; 

			Ex = (float*) calloc(number_of_cells_with_pads, size_float);
			Ey = (float*) calloc(number_of_cells_with_pads, size_float);
			Ez = (float*) calloc(number_of_cells_with_pads, size_float);
			Hx = (float*) calloc(number_of_cells_with_pads, size_float);
			Hy = (float*) calloc(number_of_cells_with_pads, size_float);
			Hz = (float*) calloc(number_of_cells_with_pads, size_float);

			Cexe  = (float*) calloc(number_of_cells_with_pads, size_float);
			Cexhz = (float*) calloc(number_of_cells_with_pads, size_float);
			Ceye  = (float*) calloc(number_of_cells_with_pads, size_float);
			Ceyhz = (float*) calloc(number_of_cells_with_pads, size_float);
			Ceze  = (float*) calloc(number_of_cells_with_pads, size_float);
			Cezhy = (float*) calloc(number_of_cells_with_pads, size_float);
			Cezhx = (float*) calloc(number_of_cells_with_pads, size_float);

			Chxh  = (float*) calloc(number_of_cells_with_pads, size_float);
			Chxez = (float*) calloc(number_of_cells_with_pads, size_float);
			Chyh  = (float*) calloc(number_of_cells_with_pads, size_float);
			Chyez = (float*) calloc(number_of_cells_with_pads, size_float);
			Chzh  = (float*) calloc(number_of_cells_with_pads, size_float);
			Chzex = (float*) calloc(number_of_cells_with_pads, size_float);
			Chzey = (float*) calloc(number_of_cells_with_pads, size_float);

			input_file.read((char*)Ex, array_size_with_pads);
			input_file.read((char*)Ey, array_size_with_pads);
			input_file.read((char*)Ez, array_size_with_pads);
			input_file.read((char*)Hx, array_size_with_pads);
			input_file.read((char*)Hy, array_size_with_pads);
			input_file.read((char*)Hz, array_size_with_pads);

			input_file.read((char*)Cexe,  array_size_with_pads);
			input_file.read((char*)Cexhz, array_size_with_pads);
			input_file.read((char*)Ceye,  array_size_with_pads);
			input_file.read((char*)Ceyhz, array_size_with_pads);
			input_file.read((char*)Ceze,  array_size_with_pads);
			input_file.read((char*)Cezhy, array_size_with_pads);
			input_file.read((char*)Cezhx, array_size_with_pads);

			input_file.read((char*)Chxh,  array_size_with_pads);
			input_file.read((char*)Chxez, array_size_with_pads);
			input_file.read((char*)Chyh,  array_size_with_pads);
			input_file.read((char*)Chyez, array_size_with_pads);
			input_file.read((char*)Chzh,  array_size_with_pads);
			input_file.read((char*)Chzex, array_size_with_pads);
			input_file.read((char*)Chzey, array_size_with_pads);

			input_file.read((char*)&is_anyy_side_cpml, size_char);
			input_file.read((char*)&is_cpml_xn, size_char);
			input_file.read((char*)&is_cpml_xp, size_char);
			input_file.read((char*)&is_cpml_yn, size_char);
			input_file.read((char*)&is_cpml_yp, size_char);

			if (is_anyy_side_cpml)
			{
				input_file.read((char*)&n_cpml_xn, size_int);
				input_file.read((char*)&n_cpml_xp, size_int);
				input_file.read((char*)&n_cpml_yn, size_int);
				input_file.read((char*)&n_cpml_yp, size_int);
			}

			cpml_1d_array_size_xnyn   = size_float*TILE_SIZE; 
			cpml_1d_array_size_xpyp   = size_float*TILE_SIZE*2; 
			cpml_2d_array_size_xn = size_float*(nyy+2)*TILE_SIZE; 
			cpml_2d_array_size_yn = size_float*nxx*TILE_SIZE; 
			cpml_2d_array_size_xp = size_float*(nyy+2)*TILE_SIZE*2; 
			cpml_2d_array_size_yp = size_float*nxx*TILE_SIZE*2; 

			if (is_cpml_xn)
			{
				cpml_a_ex_xn = (float*) malloc(cpml_1d_array_size_xnyn);
				cpml_b_ex_xn = (float*) malloc(cpml_1d_array_size_xnyn);
				cpml_a_mx_xn = (float*) malloc(cpml_1d_array_size_xnyn);
				cpml_b_mx_xn = (float*) malloc(cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_a_ex_xn, cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_b_ex_xn, cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_a_mx_xn, cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_b_mx_xn, cpml_1d_array_size_xnyn);

				if (is_TEz)
				{
					Psi_hzx_xn  = (float*) malloc(cpml_2d_array_size_xn);
					CPsi_hzx_xn = (float*) malloc(cpml_2d_array_size_xn);
					Psi_eyx_xn  = (float*) malloc(cpml_2d_array_size_xn);
					CPsi_eyx_xn = (float*) malloc(cpml_2d_array_size_xn);
					input_file.read((char*)Psi_hzx_xn,  cpml_2d_array_size_xn);
					input_file.read((char*)CPsi_hzx_xn, cpml_2d_array_size_xn);
					input_file.read((char*)Psi_eyx_xn,  cpml_2d_array_size_xn);
					input_file.read((char*)CPsi_eyx_xn, cpml_2d_array_size_xn);
				}
				if (is_TMz)
				{
					Psi_ezx_xn  = (float*) malloc(cpml_2d_array_size_xn);
					CPsi_ezx_xn = (float*) malloc(cpml_2d_array_size_xn);
					Psi_hyx_xn  = (float*) malloc(cpml_2d_array_size_xn);
					CPsi_hyx_xn = (float*) malloc(cpml_2d_array_size_xn);
					input_file.read((char*)Psi_ezx_xn,  cpml_2d_array_size_xn);
					input_file.read((char*)CPsi_ezx_xn, cpml_2d_array_size_xn);
					input_file.read((char*)Psi_hyx_xn,  cpml_2d_array_size_xn);
					input_file.read((char*)CPsi_hyx_xn, cpml_2d_array_size_xn);
				}
			}

			if (is_cpml_xp)
			{
				cpml_a_ex_xp = (float*) malloc(cpml_1d_array_size_xpyp);
				cpml_b_ex_xp = (float*) malloc(cpml_1d_array_size_xpyp);
				cpml_a_mx_xp = (float*) malloc(cpml_1d_array_size_xpyp);
				cpml_b_mx_xp = (float*) malloc(cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_a_ex_xp, cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_b_ex_xp, cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_a_mx_xp, cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_b_mx_xp, cpml_1d_array_size_xpyp);

				if (is_TEz)
				{
					Psi_hzx_xp  = (float*) malloc(cpml_2d_array_size_xp);
					CPsi_hzx_xp = (float*) malloc(cpml_2d_array_size_xp);
					Psi_eyx_xp  = (float*) malloc(cpml_2d_array_size_xp);
					CPsi_eyx_xp = (float*) malloc(cpml_2d_array_size_xp);
					input_file.read((char*)Psi_hzx_xp,  cpml_2d_array_size_xp);
					input_file.read((char*)CPsi_hzx_xp, cpml_2d_array_size_xp);
					input_file.read((char*)Psi_eyx_xp,  cpml_2d_array_size_xp);
					input_file.read((char*)CPsi_eyx_xp, cpml_2d_array_size_xp);
				}
				if (is_TMz)
				{
					Psi_ezx_xp  = (float*) malloc(cpml_2d_array_size_xp);
					CPsi_ezx_xp = (float*) malloc(cpml_2d_array_size_xp);
					Psi_hyx_xp  = (float*) malloc(cpml_2d_array_size_xp);
					CPsi_hyx_xp = (float*) malloc(cpml_2d_array_size_xp);
					input_file.read((char*)Psi_ezx_xp,  cpml_2d_array_size_xp);
					input_file.read((char*)CPsi_ezx_xp, cpml_2d_array_size_xp);
					input_file.read((char*)Psi_hyx_xp,  cpml_2d_array_size_xp);
					input_file.read((char*)CPsi_hyx_xp, cpml_2d_array_size_xp);
				}
				cpml_shift_xp = nxx - TILE_SIZE*2;
			}

			if (is_cpml_yn)
			{
				cpml_a_ey_yn = (float*) malloc(cpml_1d_array_size_xnyn);
				cpml_b_ey_yn = (float*) malloc(cpml_1d_array_size_xnyn);
				cpml_a_my_yn = (float*) malloc(cpml_1d_array_size_xnyn);
				cpml_b_my_yn = (float*) malloc(cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_a_ey_yn, cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_b_ey_yn, cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_a_my_yn, cpml_1d_array_size_xnyn);
				input_file.read((char*)cpml_b_my_yn, cpml_1d_array_size_xnyn);
				if (is_TEz)
				{
					Psi_hzy_yn  = (float*) malloc(cpml_2d_array_size_yn);
					CPsi_hzy_yn = (float*) malloc(cpml_2d_array_size_yn);
					Psi_exy_yn  = (float*) malloc(cpml_2d_array_size_yn);
					CPsi_exy_yn = (float*) malloc(cpml_2d_array_size_yn);
					input_file.read((char*)Psi_hzy_yn,  cpml_2d_array_size_yn);
					input_file.read((char*)CPsi_hzy_yn, cpml_2d_array_size_yn);
					input_file.read((char*)Psi_exy_yn,  cpml_2d_array_size_yn);
					input_file.read((char*)CPsi_exy_yn, cpml_2d_array_size_yn);
				}
				if (is_TMz)
				{
					Psi_hxy_yn  = (float*) malloc(cpml_2d_array_size_yn);
					CPsi_hxy_yn = (float*) malloc(cpml_2d_array_size_yn);
					Psi_ezy_yn  = (float*) malloc(cpml_2d_array_size_yn);
					CPsi_ezy_yn = (float*) malloc(cpml_2d_array_size_yn);
					input_file.read((char*)Psi_hxy_yn, cpml_2d_array_size_yn);
					input_file.read((char*)CPsi_hxy_yn, cpml_2d_array_size_yn);
					input_file.read((char*)Psi_ezy_yn, cpml_2d_array_size_yn);
					input_file.read((char*)CPsi_ezy_yn, cpml_2d_array_size_yn);
				}
			}

			if (is_cpml_yp)
			{
				cpml_a_ey_yp = (float*) malloc(cpml_1d_array_size_xpyp);
				cpml_b_ey_yp = (float*) malloc(cpml_1d_array_size_xpyp);
				cpml_a_my_yp = (float*) malloc(cpml_1d_array_size_xpyp);
				cpml_b_my_yp = (float*) malloc(cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_a_ey_yp, cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_b_ey_yp, cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_a_my_yp, cpml_1d_array_size_xpyp);
				input_file.read((char*)cpml_b_my_yp, cpml_1d_array_size_xpyp);
				if (is_TEz)
				{
					Psi_hzy_yp  = (float*) malloc(cpml_2d_array_size_yp);
					CPsi_hzy_yp = (float*) malloc(cpml_2d_array_size_yp);
					Psi_exy_yp  = (float*) malloc(cpml_2d_array_size_yp);
					CPsi_exy_yp = (float*) malloc(cpml_2d_array_size_yp);
					input_file.read((char*)Psi_hzy_yp,  cpml_2d_array_size_yp);
					input_file.read((char*)CPsi_hzy_yp, cpml_2d_array_size_yp);
					input_file.read((char*)Psi_exy_yp,  cpml_2d_array_size_yp);
					input_file.read((char*)CPsi_exy_yp, cpml_2d_array_size_yp);
				}
				if (is_TMz)
				{
					Psi_hxy_yp  = (float*) malloc(cpml_2d_array_size_yp);
					CPsi_hxy_yp = (float*) malloc(cpml_2d_array_size_yp);
					Psi_ezy_yp  = (float*) malloc(cpml_2d_array_size_yp);
					CPsi_ezy_yp = (float*) malloc(cpml_2d_array_size_yp);
					input_file.read((char*)Psi_hxy_yp,  cpml_2d_array_size_yp);
					input_file.read((char*)CPsi_hxy_yp, cpml_2d_array_size_yp);
					input_file.read((char*)Psi_ezy_yp,  cpml_2d_array_size_yp);
					input_file.read((char*)CPsi_ezy_yp, cpml_2d_array_size_yp);
				}
				cpml_shift_yp = nyy - TILE_SIZE*2;
			}

			input_file.read((char*)&number_of_sampled_electric_fields, size_int);
			input_file.read((char*)&number_of_sampled_magnetic_fields, size_int);
			input_file.read((char*)&number_of_impressed_J, size_int);
			input_file.read((char*)&number_of_impressed_M, size_int);

			sampled_electric_fields_component = (char*) malloc(number_of_sampled_electric_fields*sizeof(char));
			sampled_electric_fields_is = (int*) malloc(number_of_sampled_electric_fields*size_int);
			sampled_electric_fields_js = (int*) malloc(number_of_sampled_electric_fields*size_int);
			sampled_electric_fields_sampled_value = 
				(float*) malloc(number_of_sampled_electric_fields*number_of_time_steps*size_float);

			for (int i=0;i<number_of_sampled_electric_fields;i++)
			{
				input_file.read((char*)&sampled_electric_fields_component[i], sizeof(char));
				input_file.read((char*)&sampled_electric_fields_is[i], size_int);
				input_file.read((char*)&sampled_electric_fields_js[i], size_int);
				input_file.read((char*)&sampled_electric_fields_sampled_value[i*number_of_time_steps], 
					size_float*number_of_time_steps);
			}

			sampled_magnetic_fields_component = (char*) malloc(number_of_sampled_magnetic_fields*size_int);
			sampled_magnetic_fields_is = (int*) malloc(number_of_sampled_magnetic_fields*size_int);
			sampled_magnetic_fields_js = (int*) malloc(number_of_sampled_magnetic_fields*size_int);
			sampled_magnetic_fields_sampled_value = 
				(float*) malloc(number_of_sampled_magnetic_fields*number_of_time_steps*size_float);

			for (int i=0;i<number_of_sampled_magnetic_fields;i++)
			{
				input_file.read((char*)&sampled_magnetic_fields_component[i], sizeof(char));
				input_file.read((char*)&sampled_magnetic_fields_is[i], size_int);
				input_file.read((char*)&sampled_magnetic_fields_js[i], size_int);
				input_file.read((char*)&sampled_magnetic_fields_sampled_value[i*number_of_time_steps], 
					size_float*number_of_time_steps);
			}

			impressed_J_direction = (char*) malloc(number_of_impressed_J*sizeof(char));
			impressed_J_is = (int*) malloc(number_of_impressed_J*size_int);
			impressed_J_ie = (int*) malloc(number_of_impressed_J*size_int);
			impressed_J_js = (int*) malloc(number_of_impressed_J*size_int);
			impressed_J_je = (int*) malloc(number_of_impressed_J*size_int);
			impressed_J_min_x = (float*) malloc(number_of_impressed_J*size_float);
			impressed_J_min_y = (float*) malloc(number_of_impressed_J*size_float);
			impressed_J_max_x = (float*) malloc(number_of_impressed_J*size_float);
			impressed_J_max_y = (float*) malloc(number_of_impressed_J*size_float);
			impressed_J_waveform = (float*) malloc(number_of_impressed_J*number_of_time_steps*size_float);

			number_of_cej_components = (int*) malloc(number_of_impressed_J*size_int);
			input_file.read((char*)number_of_cej_components, number_of_impressed_J*size_int);

			total_number_of_cej_components = 0;
			for (int i=0;i<number_of_impressed_J;i++) total_number_of_cej_components 
				= total_number_of_cej_components + number_of_cej_components[i];

			impressed_J_Cej = (float*) malloc(total_number_of_cej_components*size_float);

			int cej_ind = 0;
			for (int i=0;i<number_of_impressed_J;i++)
			{
				input_file.read((char*)&impressed_J_direction[i], sizeof(char));
				input_file.read((char*)&impressed_J_is[i], size_int);
				input_file.read((char*)&impressed_J_ie[i], size_int);
				input_file.read((char*)&impressed_J_js[i], size_int);
				input_file.read((char*)&impressed_J_je[i], size_int);

				input_file.read((char*)&impressed_J_min_x[i], size_float);
				input_file.read((char*)&impressed_J_min_y[i], size_float);
				input_file.read((char*)&impressed_J_max_x[i], size_float);
				input_file.read((char*)&impressed_J_max_y[i], size_float);

				input_file.read((char*)&impressed_J_waveform[i*number_of_time_steps], 
					size_float*number_of_time_steps);

				input_file.read((char*)&impressed_J_Cej[cej_ind], 
					size_float*number_of_cej_components[i]);
				cej_ind = cej_ind + number_of_cej_components[i];
			}

			impressed_M_direction = (char*) malloc(number_of_impressed_M*sizeof(char));
			impressed_M_is = (int*) malloc(number_of_impressed_M*size_int);
			impressed_M_ie = (int*) malloc(number_of_impressed_M*size_int);
			impressed_M_js = (int*) malloc(number_of_impressed_M*size_int);
			impressed_M_je = (int*) malloc(number_of_impressed_M*size_int);
			impressed_M_min_x = (float*) malloc(number_of_impressed_J*size_float);
			impressed_M_min_y = (float*) malloc(number_of_impressed_J*size_float);
			impressed_M_max_x = (float*) malloc(number_of_impressed_J*size_float);
			impressed_M_max_y = (float*) malloc(number_of_impressed_J*size_float);
			impressed_M_waveform = (float*) malloc(number_of_impressed_M*number_of_time_steps*size_float);

			number_of_chm_components = (int*) malloc(number_of_impressed_M*size_int);
			input_file.read((char*)number_of_chm_components, number_of_impressed_M*size_int);

			total_number_of_chm_components = 0;
			for (int i=0;i<number_of_impressed_M;i++) total_number_of_chm_components 
				= total_number_of_chm_components + number_of_chm_components[i];

			impressed_M_Chm = (float*) malloc(total_number_of_chm_components*size_float);

			int chm_ind = 0;
			for (int i=0;i<number_of_impressed_M;i++)
			{
				input_file.read((char*)&impressed_M_direction[i], sizeof(char));
				input_file.read((char*)&impressed_M_is[i], size_int);
				input_file.read((char*)&impressed_M_ie[i], size_int);
				input_file.read((char*)&impressed_M_js[i], size_int);
				input_file.read((char*)&impressed_M_je[i], size_int);

				input_file.read((char*)&impressed_M_min_x[i], size_float);
				input_file.read((char*)&impressed_M_min_y[i], size_float);
				input_file.read((char*)&impressed_M_max_x[i], size_float);
				input_file.read((char*)&impressed_M_max_y[i], size_float);

				input_file.read((char*)&impressed_M_waveform[i*number_of_time_steps], 
					size_float*number_of_time_steps);

				input_file.read((char*)&impressed_M_Chm[chm_ind], 
					size_float*number_of_chm_components[i]);
				chm_ind = chm_ind + number_of_chm_components[i];
			}

			input_file.read((char*)&show_Hz, size_char);
			input_file.read((char*)&show_Ez, size_char);
		}
		input_file.close();
	}
	return true;
}

bool initialize_sources()
{
	int is, ie, js, je, ind;

	J_indices = (int*) malloc(total_number_of_cej_components*sizeof(int));
	ind = 0;
	for (int i=0;i<number_of_impressed_J;i++)
	{
		is = impressed_J_is[i]-1;
		ie = impressed_J_ie[i]-1;
		js = impressed_J_js[i]-1;
		je = impressed_J_je[i]-1;
		switch (impressed_J_direction[i])
		{
		  case 'x' : 
			  for (int mi=is;mi<ie;mi++)
				  for (int mj=js;mj<=je;mj++)
				  {
						J_indices[ind] = (mj+1)*nxx+mi;
						ind++;
					}
			break;
		  case 'y' : 
			  for (int mi=is;mi<=ie;mi++)
				  for (int mj=js;mj<je;mj++)
					{
						J_indices[ind] = (mj+1)*nxx+mi;
						ind++;
					}
			break;
		  default : 
			  for (int mi=is;mi<=ie;mi++)
				  for (int mj=js;mj<=je;mj++)
					{
						J_indices[ind] = (mj+1)*nxx+mi;
						ind++;
					}

		}
	}

	M_indices = (int*) malloc(total_number_of_chm_components*sizeof(int));
	ind = 0;
	for (int i=0;i<number_of_impressed_M;i++)
	{
		is = impressed_M_is[i]-1;
		ie = impressed_M_ie[i]-1;
		js = impressed_M_js[i]-1;
		je = impressed_M_je[i]-1;
		switch (impressed_M_direction[i])
		{
		  case 'x' : 
			  for (int mi=is;mi<=ie;mi++)
				  for (int mj=js;mj<je;mj++)
					{
						M_indices[ind] = (mj+1)*nxx+mi;
						ind++;
					}
			break;
		  case 'y' : 
			  for (int mi=is;mi<ie;mi++)
				  for (int mj=js;mj<=je;mj++)
					{
						M_indices[ind] = (mj+1)*nxx+mi;
  						ind++;
					}

			break;
		  default : 
			  for (int mi=is;mi<ie;mi++)
				  for (int mj=js;mj<je;mj++)
					{
						M_indices[ind] = (mj+1)*nxx+mi;
						ind++;
					}
		}
	}

	return true;
}

bool deallocateArrays()
{
	free(Ex);
	free(Ey);
	free(Ez);
	free(Hx);
	free(Hy);
	free(Hz);

	free(Cexe);
	free(Cexhz);
	free(Ceye);
	free(Ceyhz);
	free(Ceze);
	free(Cezhy);
	free(Cezhx);

	free(Chxh);
	free(Chxez);
	free(Chyh);
	free(Chyez);
	free(Chzh);
	free(Chzex);
	free(Chzey);

	free(cpml_a_ex_xn);
	free(cpml_b_ex_xn);
	free(cpml_a_ex_xp);
	free(cpml_b_ex_xp);
	free(cpml_a_ey_yn);
	free(cpml_b_ey_yn);
	free(cpml_a_ey_yp);
	free(cpml_b_ey_yp);

	free(cpml_a_mx_xn);
	free(cpml_b_mx_xn);
	free(cpml_a_mx_xp);
	free(cpml_b_mx_xp);
	free(cpml_a_my_yn);
	free(cpml_b_my_yn);
	free(cpml_a_my_yp);
	free(cpml_b_my_yp);

	free(Psi_eyx_xn);
	free(CPsi_eyx_xn);
	free(Psi_eyx_xp);
	free(CPsi_eyx_xp);
	free(Psi_exy_yn);
	free(CPsi_exy_yn);
	free(Psi_exy_yp);
	free(CPsi_exy_yp);

	free(Psi_ezx_xn);
	free(CPsi_ezx_xn);
	free(Psi_ezx_xp);
	free(CPsi_ezx_xp);
	free(Psi_ezy_yn);
	free(CPsi_ezy_yn);
	free(Psi_ezy_yp);
	free(CPsi_ezy_yp);

	free(Psi_hzx_xn);
	free(CPsi_hzx_xn);
	free(Psi_hzx_xp);
	free(CPsi_hzx_xp);
	free(Psi_hzy_yn);
	free(CPsi_hzy_yn);
	free(Psi_hzy_yp);
	free(CPsi_hzy_yp);

	free(Psi_hyx_xn);
	free(CPsi_hyx_xn);
	free(Psi_hyx_xp);
	free(CPsi_hyx_xp);
	free(Psi_hxy_yn);
	free(CPsi_hxy_yn);
	free(Psi_hxy_yp);
	free(CPsi_hxy_yp);

	free(sampled_electric_fields_component);
	free(sampled_electric_fields_is);
	free(sampled_electric_fields_js);

	free(sampled_magnetic_fields_component);
	free(sampled_magnetic_fields_is);
	free(sampled_magnetic_fields_js);

	free(number_of_cej_components);
	free(impressed_J_direction);
	free(impressed_J_is);
	free(impressed_J_ie);
	free(impressed_J_js);
	free(impressed_J_je);
	free(impressed_J_min_x);
	free(impressed_J_min_y);
	free(impressed_J_max_x);
	free(impressed_J_max_y);
	free(impressed_J_waveform);
	free(impressed_J_Cej);
	free(J_indices);

	free(number_of_chm_components);
	free(impressed_M_direction);
	free(impressed_M_is);
	free(impressed_M_ie);
	free(impressed_M_js);
	free(impressed_M_je);
	free(impressed_M_min_x);
	free(impressed_M_min_y);
	free(impressed_M_max_x);
	free(impressed_M_max_y);
	free(impressed_M_waveform);
	free(impressed_M_Chm);
	free(M_indices);

	return true;
}
bool saveSampledFieldsToFile()
{

	strcpy(output_file_name, "result_");
	strncat(output_file_name, input_file_name, strlen(input_file_name)); 

	ofstream output_file;
	output_file.open (output_file_name, ios::out | ios::binary);
	
	if (!output_file.is_open())
	{
		cout<<"File <" << output_file_name << "> can not be opened! \n";
		return false;
	}
	else 
	{
		output_file.write((char*)sampled_electric_fields_sampled_value,number_of_sampled_electric_fields*sizeof(float)*number_of_time_steps);
		output_file.write((char*)sampled_magnetic_fields_sampled_value,number_of_sampled_magnetic_fields*sizeof(float)*number_of_time_steps);
		free(sampled_electric_fields_sampled_value);
		free(sampled_magnetic_fields_sampled_value);
	}

	output_file.close();
	return true;
}

