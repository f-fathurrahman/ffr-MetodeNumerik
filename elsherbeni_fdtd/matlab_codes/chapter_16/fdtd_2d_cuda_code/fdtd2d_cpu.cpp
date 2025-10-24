
#include <stdlib.h>
#include <stdio.h>

#include "global.h"


void update_magnetic_fields_on_cpu_TEz()
{
	for (int i=0;i<number_of_cells-nxx;i++)
	Hz[i] = Chzh[i] * Hz[i] + Chzex[i] * (Ex[i+nxx]-Ex[i])  
					+ Chzey[i] * (Ey[i+1]-Ey[i]); 
}

void update_magnetic_fields_on_cpu_TMz()
{
	for (int i=0;i<number_of_cells-nxx;i++)
		Hx[i] = Chxh[i] * Hx[i] + Chxez[i] * (Ez[i+nxx]-Ez[i]); 

	for (int i=0;i<number_of_cells-nxx;i++)
		Hy[i] = Chyh[i] * Hy[i] + Chyez[i] * (Ez[i+1] - Ez[i]);
}

void update_impressed_magnetic_currents_on_cpu(int time_step)
{
	float *H;
	float M_value;
	int chm_ind=0;
	int field_index;

	for (int i=0;i<number_of_impressed_M;i++)
	{
		if (impressed_M_direction[i]=='x') H = Hx;
		else if (impressed_M_direction[i]=='y')	H = Hy; else H = Hz;
		
		M_value = impressed_M_waveform[number_of_time_steps*i+time_step];
			
		for (int j=0;j<number_of_chm_components[i];j++)
		{
			field_index = M_indices[chm_ind];
			H[field_index] = H[field_index] + M_value * impressed_M_Chm[chm_ind];		
			chm_ind++;
		}
	}
}

void capture_sampled_magnetic_fields_on_cpu(int time_step)
{
	for (int i=0;i<number_of_sampled_magnetic_fields;i++)
	{
		int ci = sampled_magnetic_fields_js[i]*nxx+sampled_magnetic_fields_is[i]-1;
		float value;

		switch (sampled_magnetic_fields_component[i])
		{
		  case 'x' : 
				value = 0.5*(Hx[ci] + Hx[ci-nxx]);
			break;
		  case 'y' : 
				value = 0.5*(Hy[ci] + Hy[ci-1]);
			break;
		  default : 
				value = 0.25*(Hz[ci] + Hz[ci-1] + Hz[ci-nxx] + Hz[ci-nxx-1]);
		}
		sampled_magnetic_fields_sampled_value[number_of_time_steps*i+time_step] = value;
	}
}

void update_electric_fields_on_cpu_TEz()
{
	for (int i=nxx;i<number_of_cells;i++)
		Ex[i] = Cexe[i] * Ex[i] + Cexhz[i] * (Hz[i]-Hz[i-nxx]);

	for (int i=nxx;i<number_of_cells;i++)
		Ey[i] = Ceye[i] * Ey[i] + Ceyhz[i] * (Hz[i]-Hz[i-1]);
}

void update_electric_fields_on_cpu_TMz()
{
	for (int i=nxx;i<number_of_cells;i++)
		Ez[i] = Ceze[i] * Ez[i] + Cezhy[i] * (Hy[i]-Hy[i-1]) + Cezhx[i] * (Hx[i]-Hx[i-nxx]);
}

void update_impressed_electric_currents_on_cpu(int time_step)
{
	float* E;
	float J_value;
	int cej_ind = 0;
	int field_index;

	for (int i=0;i<number_of_impressed_J;i++)
	{
		if (impressed_J_direction[i]=='x') E = Ex;
			else if (impressed_J_direction[i]=='y')	E = Ey; else E = Ez;

		J_value = impressed_J_waveform[number_of_time_steps*i+time_step];

		for (int j=0;j<number_of_cej_components[i];j++)
		{
			field_index = J_indices[cej_ind];
			E[field_index] = E[field_index] + J_value * impressed_J_Cej[cej_ind];
			cej_ind++; 
		}
	}
}

void capture_sampled_electric_fields_on_cpu(int time_step)
{
	for (int i=0;i<number_of_sampled_electric_fields;i++)
	{
		int ci = sampled_electric_fields_js[i]*nxx+sampled_electric_fields_is[i]-1;
		float value;
	
		switch (sampled_electric_fields_component[i])
		{
		  case 'x' : 
				value = 0.5*(Ex[ci] + Ex[ci-1]);
			break;
		  case 'y' : 
				value = 0.5*(Ey[ci] + Ey[ci-nxx]);
			break;
		  default : 
				value = Ez[ci];
		}	
		sampled_electric_fields_sampled_value[number_of_time_steps*i+time_step] = value;
	}
}

void update_magnetic_fields_on_cpu_CPML_TEz_xn()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx;
		ti = (j+1)*TILE_SIZE;
		for (i=0;i<n_cpml_xn;i++)
		{
			Psi_hzx_xn[ti] = cpml_b_mx_xn[i] * Psi_hzx_xn[ti] + cpml_a_mx_xn[i] * (Ey[ci+1]-Ey[ci]);
			Hz[ci] = Hz[ci] + CPsi_hzx_xn[ti] * Psi_hzx_xn[ti];	
			ci++; ti++;
		}
	}
}

void update_magnetic_fields_on_cpu_CPML_TMz_xn()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx;
		ti = (j+1)*TILE_SIZE;
		for (i=0;i<n_cpml_xn;i++)
		{
			Psi_hyx_xn[ti] = cpml_b_mx_xn[i] * Psi_hyx_xn[ti] + cpml_a_mx_xn[i]*(Ez[ci+1] - Ez[ci]); 
			Hy[ci] = Hy[ci] + CPsi_hyx_xn[ti] * Psi_hyx_xn[ti];
			ci++; ti++;
		}
	}
}

void update_magnetic_fields_on_cpu_CPML_TEz_xp()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx+cpml_shift_xp;
		ti = (j+1)*TILE_SIZE*2;
		for (i=0;i<TILE_SIZE*2;i++)
		{
			Psi_hzx_xp[ti] = cpml_b_mx_xp[i] * Psi_hzx_xp[ti] + cpml_a_mx_xp[i] * (Ey[ci+1]-Ey[ci]);
			Hz[ci] = Hz[ci] + CPsi_hzx_xp[ti] * Psi_hzx_xp[ti];	
			ci++; ti++;
		}
	}
}

void update_magnetic_fields_on_cpu_CPML_TMz_xp()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx+cpml_shift_xp;
		ti = (j+1)*TILE_SIZE*2;
		for (i=0;i<TILE_SIZE*2;i++)
		{
			Psi_hyx_xp[ti] = cpml_b_mx_xp[i] * Psi_hyx_xp[ti] + cpml_a_mx_xp[i]*(Ez[ci+1] - Ez[ci]); 
			Hy[ci] = Hy[ci] + CPsi_hyx_xp[ti] * Psi_hyx_xp[ti];
			ci++; ti++;
		}
	}
}

void update_magnetic_fields_on_cpu_CPML_TEz_yn()
{
	int i, j, ci, ti;
	for (j=0;j<n_cpml_yn;j++)
	{
		ci = (j+1)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_hzy_yn[ti] = cpml_b_my_yn[j] * Psi_hzy_yn[ti] + cpml_a_my_yn[j] * (Ex[ci+nxx] - Ex[ci]); 	
			Hz[ci] = Hz[ci] + CPsi_hzy_yn[ti] * Psi_hzy_yn[ti];
			ci++; ti++;
		}
	}
}

void update_magnetic_fields_on_cpu_CPML_TMz_yn()
{
	int i, j, ci, ti;
	for (j=0;j<n_cpml_yn;j++)
	{
		ci = (j+1)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_hxy_yn[ti] = cpml_b_my_yn[j] * Psi_hxy_yn[ti] + cpml_a_my_yn[j] * (Ez[ci+nxx] - Ez[ci]); 	
			Hx[ci] = Hx[ci] + CPsi_hxy_yn[ti] * Psi_hxy_yn[ti];
			ci++; ti++;
		}
	}
}

void update_magnetic_fields_on_cpu_CPML_TEz_yp()
{
	int i, j, ci, ti;
	for (j=0;j<TILE_SIZE*2;j++)
	{
		ci = (j+1+cpml_shift_yp)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_hzy_yp[ti] = cpml_b_my_yp[j] * Psi_hzy_yp[ti] + cpml_a_my_yp[j] * (Ex[ci+nxx] - Ex[ci]); 	
			Hz[ci] = Hz[ci] + CPsi_hzy_yp[ti] * Psi_hzy_yp[ti];
			ci++; ti++;
		}
	}
}
void update_magnetic_fields_on_cpu_CPML_TMz_yp()
{
	int i, j, ci, ti;
	for (j=0;j<TILE_SIZE*2;j++)
	{
		ci = (j+1+cpml_shift_yp)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_hxy_yp[ti] = cpml_b_my_yp[j] * Psi_hxy_yp[ti] + cpml_a_my_yp[j] * (Ez[ci+nxx] - Ez[ci]); 	
			Hx[ci] = Hx[ci] + CPsi_hxy_yp[ti] * Psi_hxy_yp[ti];
			ci++; ti++;
		}
	}
}

void update_magnetic_fields_for_CPML_on_cpu()
{
	if (is_cpml_xn)
	{
		if (is_TEz)
			update_magnetic_fields_on_cpu_CPML_TEz_xn();

		if (is_TMz)
			update_magnetic_fields_on_cpu_CPML_TMz_xn();
	}

	if (is_cpml_xp)
	{
		if (is_TEz)
			update_magnetic_fields_on_cpu_CPML_TEz_xp();

		if (is_TMz)
			update_magnetic_fields_on_cpu_CPML_TMz_xp();
	}

	if (is_cpml_yn)
	{
		if (is_TEz)
			update_magnetic_fields_on_cpu_CPML_TEz_yn();

		if (is_TMz)
			update_magnetic_fields_on_cpu_CPML_TMz_yn();
	}

	if (is_cpml_yp)
	{
		if (is_TEz)
			update_magnetic_fields_on_cpu_CPML_TEz_yp();

		if (is_TMz)
			update_magnetic_fields_on_cpu_CPML_TMz_yp();
	}
}

void update_electric_fields_on_cpu_CPML_TEz_xn()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx;
		ti = (j+1)*TILE_SIZE;
		for (i=0;i<n_cpml_xn;i++)
		{
			Psi_eyx_xn[ti] = cpml_b_ex_xn[i] * Psi_eyx_xn[ti] + cpml_a_ex_xn[i]*(Hz[ci]-Hz[ci-1]); 
			Ey[ci] = Ey[ci] + CPsi_eyx_xn[ti] * Psi_eyx_xn[ti];
			ci++; ti++;
		}
	}
}

void update_electric_fields_on_gpu_CPML_TMz_xn()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx;
		ti = (j+1)*TILE_SIZE;
		for (i=0;i<n_cpml_xn;i++)
		{
			Psi_ezx_xn[ti] = cpml_b_ex_xn[i] * Psi_ezx_xn[ti] + cpml_a_ex_xn[i]*(Hy[ci]-Hy[ci-1]); 
			Ez[ci] = Ez[ci] + CPsi_ezx_xn[ti] * Psi_ezx_xn[ti];
			ci++; ti++;
		}
	}
}

void update_electric_fields_on_cpu_CPML_TEz_xp()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx+cpml_shift_xp;
		ti = (j+1)*TILE_SIZE*2;
		for (i=0;i<TILE_SIZE*2;i++)
		{
			Psi_eyx_xp[ti] = cpml_b_ex_xp[i] * Psi_eyx_xp[ti] + cpml_a_ex_xp[i]*(Hz[ci]-Hz[ci-1]); 
			Ey[ci] = Ey[ci] + CPsi_eyx_xp[ti] * Psi_eyx_xp[ti];
			ci++; ti++;
		}
	}
}

void update_electric_fields_on_cpu_CPML_TMz_xp()
{
	int i, j, ci, ti;
	for (j=0;j<nyy;j++)
	{
		ci = (j+1)*nxx+cpml_shift_xp;
		ti = (j+1)*TILE_SIZE*2;
		for (i=0;i<TILE_SIZE*2;i++)
		{
			Psi_ezx_xp[ti] = cpml_b_ex_xp[i] * Psi_ezx_xp[ti] + cpml_a_ex_xp[i]*(Hy[ci]-Hy[ci-1]); 
			Ez[ci] = Ez[ci] + CPsi_ezx_xp[ti] * Psi_ezx_xp[ti];
			ci++; ti++;
		}
	}
}

void update_electric_fields_on_cpu_CPML_TEz_yn()
{
	int i, j, ci, ti;
	for (j=0;j<n_cpml_yn;j++)
	{
		ci = (j+1)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_exy_yn[ti] = cpml_b_ey_yn[j] * Psi_exy_yn[ti] + cpml_a_ey_yn[j]*(Hz[ci]-Hz[ci-nxx]); 
			Ex[ci] = Ex[ci] + CPsi_exy_yn[ti] * Psi_exy_yn[ti];
			ci++; ti++;
		}
	}
}
void update_electric_fields_on_cpu_CPML_TMz_yn()
{
	int i, j, ci, ti;
	for (j=0;j<n_cpml_yn;j++)
	{
		ci = (j+1)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_ezy_yn[ti] = cpml_b_ey_yn[j] * Psi_ezy_yn[ti] + cpml_a_ey_yn[j]*(Hx[ci]-Hx[ci-nxx]); 
			Ez[ci] = Ez[ci] + CPsi_ezy_yn[ti] * Psi_ezy_yn[ti];
			ci++; ti++;
		}
	}
}

void update_electric_fields_on_cpu_CPML_TEz_yp()
{
	int i, j, ci, ti;
	for (j=0;j<TILE_SIZE*2;j++)
	{
		ci = (j+1+cpml_shift_yp)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_exy_yp[ti] = cpml_b_ey_yp[j] * Psi_exy_yp[ti] + cpml_a_ey_yp[j]*(Hz[ci]-Hz[ci-nxx]); 
			Ex[ci] = Ex[ci] + CPsi_exy_yp[ti] * Psi_exy_yp[ti];
			ci++; ti++;
		}
	}
}

void update_electric_fields_on_cpu_CPML_TMz_yp()
{
	int i, j, ci, ti;
	for (j=0;j<TILE_SIZE*2;j++)
	{
		ci = (j+1+cpml_shift_yp)*nxx;
		ti = j*nxx;
		for (i=0;i<nxx;i++)
		{
			Psi_ezy_yp[ti] = cpml_b_ey_yp[j] * Psi_ezy_yp[ti] + cpml_a_ey_yp[j]*(Hx[ci]-Hx[ci-nxx]); 
			Ez[ci] = Ez[ci] + CPsi_ezy_yp[ti] * Psi_ezy_yp[ti];
			ci++; ti++;
		}
	}
}
void update_electric_fields_for_CPML_on_cpu()
{
	if (is_cpml_xn)
	{
		if (is_TEz)
			update_electric_fields_on_cpu_CPML_TEz_xn();

		if (is_TMz)
			update_electric_fields_on_gpu_CPML_TMz_xn();
	}

	if (is_cpml_xp)
	{
		if (is_TEz)
			update_electric_fields_on_cpu_CPML_TEz_xp();

		if (is_TMz)
			update_electric_fields_on_cpu_CPML_TMz_xp();
	}

	if (is_cpml_yn)
	{
		if (is_TEz)
			update_electric_fields_on_cpu_CPML_TEz_yn();

		if (is_TMz)
			update_electric_fields_on_cpu_CPML_TMz_yn();
	}

	if (is_cpml_yp)
	{
		if (is_TEz)
			update_electric_fields_on_cpu_CPML_TEz_yp();

		if (is_TMz)
			update_electric_fields_on_cpu_CPML_TMz_yp();
	}
}

bool fdtdIterationOnCpu()
{
	if (is_TEz) update_magnetic_fields_on_cpu_TEz();

	if (is_TMz)	update_magnetic_fields_on_cpu_TMz();

	update_impressed_magnetic_currents_on_cpu(time_step);
		
	update_magnetic_fields_for_CPML_on_cpu();

	capture_sampled_magnetic_fields_on_cpu(time_step);

	if (is_TEz)	update_electric_fields_on_cpu_TEz();

	if (is_TMz)	update_electric_fields_on_cpu_TMz();

	update_impressed_electric_currents_on_cpu(time_step);

	update_electric_fields_for_CPML_on_cpu();

	capture_sampled_electric_fields_on_cpu(time_step);

	if (time_step%100 == 0)
		printf("timestep: %d \n", time_step);
	
	time_step++;

	return true;
}
