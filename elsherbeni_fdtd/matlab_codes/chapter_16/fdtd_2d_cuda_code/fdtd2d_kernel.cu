#pragma once
#ifndef _fdtd2d_KERNEL_H_
#define _fdtd2d_KERNEL_H_

#define TILE_SIZE 16

__global__ void
update_magnetic_fields_on_gpu_TEz(float* Chzh, float* Chzex, float* Chzey, float* Hz,  float* Ex, float* Ey, int nxx)
{
	__shared__ float sEy[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;

	sEy[ty][tx] = Ey[ci];
	sEy[ty][tx+TILE_SIZE] = Ey[ci+TILE_SIZE];

	__syncthreads();
		Hz[ci] = Chzh[ci] * Hz[ci] + Chzex[ci] * (Ex[ci+nxx]-Ex[ci])  
					+ Chzey[ci] * (sEy[ty][tx+1]-sEy[ty][tx]); 
}

__global__ void
update_magnetic_fields_on_gpu_TMz(float* Chxh, float* Chxez, float* Chyh, float* Chyez, float* Hx,  float* Hy, float* Ez, int nxx)
{
	__shared__ float sEz[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;

	sEz[ty][tx] = Ez[ci];
	sEz[ty][tx+TILE_SIZE] = Ez[ci+TILE_SIZE];

	__syncthreads();

	Hx[ci] = Chxh[ci] * Hx[ci] + Chxez[ci] * (Ez[ci+nxx]-sEz[ty][tx]); 
	Hy[ci] = Chyh[ci] * Hy[ci] + Chyez[ci] * (sEz[ty][tx+1] - sEz[ty][tx]);
}

__global__ void
update_impressed_M_on_gpu(float* H, int* M_indices,  float* M_Chm, int ind, float M_value)
{
	int field_index = M_indices[ind+threadIdx.x];

	H[field_index] = H[field_index] + M_value * M_Chm[ind + threadIdx.x];
}

__global__ void
update_impressed_J_on_gpu(float* E, int* J_indices,  float* J_Cej, int ind, float J_value)
{
	int field_index = J_indices[ind+threadIdx.x];

	E[field_index] = E[field_index] + J_value * J_Cej[ind + threadIdx.x];
}

__global__ void
update_electric_fields_on_gpu_TEz(float* Cexe, float* Cexhz, float* Ceye, float* Ceyhz,  float* Ex, float* Ey, float* Hz, int nxx)
{
	__shared__ float sHz[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;
	
	sHz[ty][tx] = Hz[ci-TILE_SIZE];
	sHz[ty][tx+TILE_SIZE] = Hz[ci];

	__syncthreads();

	Ex[ci] = Cexe[ci] * Ex[ci] + Cexhz[ci] * (Hz[ci]-Hz[ci-nxx]);
	Ey[ci] = Ceye[ci] * Ey[ci] + Ceyhz[ci] * (sHz[ty][tx+TILE_SIZE]-sHz[ty][tx+TILE_SIZE-1]);
}

__global__ void
update_electric_fields_on_gpu_TMz(float* Ceze, float* Cezhy, float* Cezhx, float* Hx,  float* Hy, float* Ez, int nxx)
{
	__shared__ float sHy[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;
	
	sHy[ty][tx+TILE_SIZE] = Hy[ci];
	sHy[ty][tx] = Hy[ci-TILE_SIZE];

	__syncthreads();
	Ez[ci] = Ceze[ci] * Ez[ci] + Cezhy[ci] * (sHy[ty][tx+TILE_SIZE]-sHy[ty][tx+TILE_SIZE-1]) + Cezhx[ci] * (Hx[ci]-Hx[ci-nxx]);
}

__global__ void
capture_sampled_electric_fields_on_gpu(float* Ex, float* Ey, float* Ez, char* component, int* is, 
											int* js, float* sampled_value, int time_step, int number_of_time_steps, int nxx)
{
	int tx = threadIdx.x;
	int ci = js[tx]*nxx+is[tx]-1;
	float value;

	switch (component[tx])
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
	
	sampled_value[number_of_time_steps*tx+time_step] = value;
}


__global__ void
capture_sampled_magnetic_fields_on_gpu(float* Hx, float* Hy, float* Hz, char* component, int* is, 
											int* js, float* sampled_value, int time_step, int number_of_time_steps, int nxx)
{
	int tx = threadIdx.x;
	int ci = js[tx]*nxx+is[tx]-1;
	float value;

	switch (component[tx])
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
	
	sampled_value[number_of_time_steps*tx+time_step] = value;
}

__constant__ float dvcpml_a_ex_xn[TILE_SIZE];
__constant__ float dvcpml_b_ex_xn[TILE_SIZE];
__constant__ float dvcpml_a_mx_xn[TILE_SIZE];
__constant__ float dvcpml_b_mx_xn[TILE_SIZE];

__constant__ float dvcpml_a_ex_xp[2*TILE_SIZE];
__constant__ float dvcpml_b_ex_xp[2*TILE_SIZE];
__constant__ float dvcpml_a_mx_xp[2*TILE_SIZE];
__constant__ float dvcpml_b_mx_xp[2*TILE_SIZE];

__constant__ float dvcpml_a_ey_yn[TILE_SIZE];
__constant__ float dvcpml_b_ey_yn[TILE_SIZE];
__constant__ float dvcpml_a_my_yn[TILE_SIZE];
__constant__ float dvcpml_b_my_yn[TILE_SIZE];

__constant__ float dvcpml_a_ey_yp[2*TILE_SIZE];
__constant__ float dvcpml_b_ey_yp[2*TILE_SIZE];
__constant__ float dvcpml_a_my_yp[2*TILE_SIZE];
__constant__ float dvcpml_b_my_yp[2*TILE_SIZE];

__global__ void
update_magnetic_fields_on_gpu_CPML_TEz_xn(float* Psi_hzx_xn, float* CPsi_hzx_xn, float* Hz, float* Ey, int nxx)
{
	__shared__ float sEy[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;
	
	int ci = (j+1)*nxx+i;
	int ti = (j+1)*TILE_SIZE+i;

	sEy[ty][tx] = Ey[ci];
	sEy[ty][tx+TILE_SIZE] = Ey[ci+TILE_SIZE];
	__syncthreads();

	Psi_hzx_xn[ti] = dvcpml_b_mx_xn[i] * Psi_hzx_xn[ti] + dvcpml_a_mx_xn[i] * (sEy[ty][tx+1]-sEy[ty][tx]); 	

	Hz[ci] = Hz[ci] + CPsi_hzx_xn[ti] * Psi_hzx_xn[ti];
}

__global__ void
update_magnetic_fields_on_gpu_CPML_TEz_xp(float* Psi_hzx_xp, float* CPsi_hzx_xp, float* Hz, float* Ey, int nxx, int cpml_shift_xp)
{
	__shared__ float sEy[TILE_SIZE][3*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;
	
	int ci = (j+1)*nxx+i+cpml_shift_xp;
	int ti = (j+1)*TILE_SIZE*2+i;

	sEy[ty][tx] = Ey[ci];
	sEy[ty][tx+TILE_SIZE] = Ey[ci+TILE_SIZE];
	__syncthreads();

	Psi_hzx_xp[ti] = dvcpml_b_mx_xp[i] * Psi_hzx_xp[ti] + dvcpml_a_mx_xp[i] * (sEy[ty][tx+1]-sEy[ty][tx]); 	

	Hz[ci] = Hz[ci] + CPsi_hzx_xp[ti] * Psi_hzx_xp[ti];
}

__global__ void
update_magnetic_fields_on_gpu_CPML_TEz_yn(float* Psi_hzy_yn, float* CPsi_hzy_yn, float* Hz, float* Ex, int nxx)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;
	
	int ci = (j+1)*nxx+i;
	int ti = j*nxx+i;

	Psi_hzy_yn[ti] = dvcpml_b_my_yn[j] * Psi_hzy_yn[ti] + dvcpml_a_my_yn[j] * (Ex[ci+nxx] - Ex[ci]); 	

	Hz[ci] = Hz[ci] + CPsi_hzy_yn[ti] * Psi_hzy_yn[ti];
}

__global__ void
update_magnetic_fields_on_gpu_CPML_TEz_yp(float* Psi_hzy_yp, float* CPsi_hzy_yp, float* Hz, float* Ex, int nxx, int cpml_shift_yp)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;
	
	int ci = (j+1+cpml_shift_yp)*nxx+i;
	int ti = j*nxx+i;

	Psi_hzy_yp[ti] = dvcpml_b_my_yp[j] * Psi_hzy_yp[ti] + dvcpml_a_my_yp[j] * (Ex[ci+nxx] - Ex[ci]); 	

	Hz[ci] = Hz[ci] + CPsi_hzy_yp[ti] * Psi_hzy_yp[ti];
}

///////////////////////////////////////////////////////////////////////

__global__ void
update_magnetic_fields_on_gpu_CPML_TMz_xn(float* Psi_hyx_xn, float* CPsi_hyx_xn, float* Hy, float* Ez, int nxx)
{
	__shared__ float sEz[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;
	int ti = (j+1)*TILE_SIZE+i;

	sEz[ty][tx] = Ez[ci];
	sEz[ty][tx+TILE_SIZE] = Ez[ci+TILE_SIZE];

	__syncthreads();

	Psi_hyx_xn[ti] = dvcpml_b_mx_xn[i] * Psi_hyx_xn[ti] + dvcpml_a_mx_xn[i]*(sEz[ty][tx+1] - sEz[ty][tx]); 

	Hy[ci] = Hy[ci] + CPsi_hyx_xn[ti] * Psi_hyx_xn[ti];
}

__global__ void
update_magnetic_fields_on_gpu_CPML_TMz_xp(float* Psi_hyx_xp, float* CPsi_hyx_xp, float* Hy, float* Ez, int nxx, int cpml_shift_xp)
{
	__shared__ float sEz[TILE_SIZE][3*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i+cpml_shift_xp;
	int ti = (j+1)*TILE_SIZE*2+i;

	sEz[ty][tx] = Ez[ci];
	sEz[ty][tx+TILE_SIZE] = Ez[ci+TILE_SIZE];

	__syncthreads();

	Psi_hyx_xp[ti] = dvcpml_b_mx_xp[i] * Psi_hyx_xp[ti] + dvcpml_a_mx_xp[i]*(sEz[ty][tx+1] - sEz[ty][tx]); 

	Hy[ci] = Hy[ci] + CPsi_hyx_xp[ti] * Psi_hyx_xp[ti];
}

__global__ void
update_magnetic_fields_on_gpu_CPML_TMz_yn(float* Psi_hxy_yn, float* CPsi_hxy_yn, float* Hx, float* Ez, int nxx)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;
	
	int ci = (j+1)*nxx+i;
	int ti = j*nxx+i;

	Psi_hxy_yn[ti] = dvcpml_b_my_yn[j] * Psi_hxy_yn[ti] + dvcpml_a_my_yn[j] * (Ez[ci+nxx] - Ez[ci]); 	

	Hx[ci] = Hx[ci] + CPsi_hxy_yn[ti] * Psi_hxy_yn[ti];
}


__global__ void
update_magnetic_fields_on_gpu_CPML_TMz_yp(float* Psi_hxy_yp, float* CPsi_hxy_yp, float* Hx, float* Ez, int nxx, int cpml_shift_yp)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;
	
	int ci = (j+1+cpml_shift_yp)*nxx+i;
	int ti = j*nxx+i;

	Psi_hxy_yp[ti] = dvcpml_b_my_yp[j] * Psi_hxy_yp[ti] + dvcpml_a_my_yp[j] * (Ez[ci+nxx] - Ez[ci]); 	

	Hx[ci] = Hx[ci] + CPsi_hxy_yp[ti] * Psi_hxy_yp[ti];
}
//------------------------------------------------------------------------

__global__ void
update_electric_fields_on_gpu_CPML_TEz_xn(float* Psi_eyx_xn, float* CPsi_eyx_xn, float* Ey, float* Hz, int nxx)
{
	__shared__ float sHz[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;
	int ti = (j+1)*TILE_SIZE+i;

	sHz[ty][tx] = Hz[ci-TILE_SIZE];
	sHz[ty][tx+TILE_SIZE] = Hz[ci];

	__syncthreads();

    Psi_eyx_xn[ti] = dvcpml_b_ex_xn[i] * Psi_eyx_xn[ti] + dvcpml_a_ex_xn[i]*(sHz[ty][tx+TILE_SIZE]-sHz[ty][tx+TILE_SIZE-1]); 

	Ey[ci] = Ey[ci] + CPsi_eyx_xn[ti] * Psi_eyx_xn[ti];

}

__global__ void
update_electric_fields_on_gpu_CPML_TEz_xp(float* Psi_eyx_xp, float* CPsi_eyx_xp, float* Ey, float* Hz, int nxx, int cpml_shift_xp)
{
	__shared__ float sHz[TILE_SIZE][3*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i+cpml_shift_xp;
	int ti = (j+1)*TILE_SIZE*2+i;

	sHz[ty][tx] = Hz[ci-TILE_SIZE];
	sHz[ty][tx+TILE_SIZE] = Hz[ci];

	__syncthreads();

    Psi_eyx_xp[ti] = dvcpml_b_ex_xp[i] * Psi_eyx_xp[ti] + dvcpml_a_ex_xp[i]*(sHz[ty][tx+TILE_SIZE]-sHz[ty][tx+TILE_SIZE-1]); 

	Ey[ci] = Ey[ci] + CPsi_eyx_xp[ti] * Psi_eyx_xp[ti];

}

__global__ void
update_electric_fields_on_gpu_CPML_TEz_yn(float* Psi_exy_yn, float* CPsi_exy_yn, float* Ex, float* Hz, int nxx)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;
	int ti = j*nxx+i;

    Psi_exy_yn[ti] = dvcpml_b_ey_yn[j] * Psi_exy_yn[ti] + dvcpml_a_ey_yn[j]*(Hz[ci]-Hz[ci-nxx]); 

	Ex[ci] = Ex[ci] + CPsi_exy_yn[ti] * Psi_exy_yn[ti];
}

__global__ void
update_electric_fields_on_gpu_CPML_TEz_yp(float* Psi_exy_yp, float* CPsi_exy_yp, float* Ex, float* Hz, int nxx, int cpml_shift_yp)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1+cpml_shift_yp)*nxx+i;
	int ti = j*nxx+i;

    Psi_exy_yp[ti] = dvcpml_b_ey_yp[j] * Psi_exy_yp[ti] + dvcpml_a_ey_yp[j]*(Hz[ci]-Hz[ci-nxx]); 

	Ex[ci] = Ex[ci] + CPsi_exy_yp[ti] * Psi_exy_yp[ti];
}
////////////////////////////////////////////////////////
__global__ void
update_electric_fields_on_gpu_CPML_TMz_xn(float* Psi_ezx_xn, float* CPsi_ezx_xn, float* Ez, float* Hy, int nxx)
{
	__shared__ float sHy[TILE_SIZE][2*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;
	int ti = (j+1)*TILE_SIZE+i;

	sHy[ty][tx] = Hy[ci-TILE_SIZE];
	sHy[ty][tx+TILE_SIZE] = Hy[ci];

	__syncthreads();

    Psi_ezx_xn[ti] = dvcpml_b_ex_xn[i] * Psi_ezx_xn[ti] + dvcpml_a_ex_xn[i]*(sHy[ty][tx+TILE_SIZE]-sHy[ty][tx+TILE_SIZE-1]); 

	Ez[ci] = Ez[ci] + CPsi_ezx_xn[ti] * Psi_ezx_xn[ti];
}

__global__ void
update_electric_fields_on_gpu_CPML_TMz_xp(float* Psi_ezx_xp, float* CPsi_ezx_xp, float* Ez, float* Hy, int nxx, int cpml_shift_xp)
{
	__shared__ float sHy[TILE_SIZE][3*TILE_SIZE+1];

	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i+cpml_shift_xp;
	int ti = (j+1)*TILE_SIZE*2+i;

	sHy[ty][tx] = Hy[ci-TILE_SIZE];
	sHy[ty][tx+TILE_SIZE] = Hy[ci];

	__syncthreads();

    Psi_ezx_xp[ti] = dvcpml_b_ex_xp[i] * Psi_ezx_xp[ti] + dvcpml_a_ex_xp[i]*(sHy[ty][tx+TILE_SIZE]-sHy[ty][tx+TILE_SIZE-1]); 

	Ez[ci] = Ez[ci] + CPsi_ezx_xp[ti] * Psi_ezx_xp[ti];

}

__global__ void
update_electric_fields_on_gpu_CPML_TMz_yn(float* Psi_ezy_yn, float* CPsi_ezy_yn, float* Ez, float* Hx, int nxx)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1)*nxx+i;
	int ti = j*nxx+i;

    Psi_ezy_yn[ti] = dvcpml_b_ey_yn[j] * Psi_ezy_yn[ti] + dvcpml_a_ey_yn[j]*(Hx[ci]-Hx[ci-nxx]); 

	Ez[ci] = Ez[ci] + CPsi_ezy_yn[ti] * Psi_ezy_yn[ti];
}

__global__ void
update_electric_fields_on_gpu_CPML_TMz_yp(float* Psi_ezy_yp, float* CPsi_ezy_yp, float* Ez, float* Hx, int nxx, int cpml_shift_yp)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int i = blockIdx.x * blockDim.x + tx;
	int j = blockIdx.y * blockDim.y + ty;

	int ci = (j+1+cpml_shift_yp)*nxx+i;
	int ti = j*nxx+i;

    Psi_ezy_yp[ti] = dvcpml_b_ey_yp[j] * Psi_ezy_yp[ti] + dvcpml_a_ey_yp[j]*(Hx[ci]-Hx[ci-nxx]); 

	Ez[ci] = Ez[ci] + CPsi_ezy_yp[ti] * Psi_ezy_yp[ti];
}

#endif
