 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fdtd2d_kernel.cu"
#include "global.h"

texture<float4, 2, cudaReadModeElementType> inTex;

__constant__ unsigned int  dvrgb[256];

__global__ void 
find_min_and_max_on_gpu(int nblocks, float* field, float* minimum_field_value, float* maximum_field_value)
{
	__shared__ float minarr[512];
	__shared__ float maxarr[512];
	
	int i = threadIdx.x;
	int nTotalThreads = blockDim.x;

	minarr[i] = field[i];
	maxarr[i] = minarr[i];
	for (int j=1;j<nblocks;j++)
	{
		minarr[i+nTotalThreads] = field[i+nTotalThreads*j];
		if (minarr[i] > minarr[i+nTotalThreads])
				minarr[i] = minarr[i+nTotalThreads];

		if (maxarr[i] < minarr[i+nTotalThreads])
				maxarr[i] = minarr[i+nTotalThreads];
		__syncthreads();
	}
		__syncthreads();

	while(nTotalThreads > 1)
	{
		int halfPoint = (nTotalThreads >> 1);	// divide by two
		if (threadIdx.x < halfPoint)
		{
			float temp = minarr[i + halfPoint];

			if (temp < minarr[i]) minarr[i] = temp;

			temp = maxarr[i + halfPoint];
			if (temp > maxarr[i]) maxarr[i] = temp;
		}
		__syncthreads();
		nTotalThreads = (nTotalThreads >> 1);	
	}
	if (i == 0)
	{
		minimum_field_value[0] = minarr[0];
		maximum_field_value[0] = maxarr[0];
	}
}

void createColormapOnGpu()
{
	cudaError_t et;	

unsigned int rgb[] = {4286775296,4287037440,4287299584,4287561728,4287823872,4288086016,4288348160,
4288610304,4288872448,4289134592,4289396736,4289658880,4289921024,4290183168,
4290445312,4290707456,4290969600,4291231744,4291493888,4291756032,4292018176,
4292280320,4292542464,4292804608,4293066752,4293328896,4293591040,4293853184,
4294115328,4294377472,4294639616,4294901760,4294902784,4294903808,4294904832,
4294905856,4294906880,4294907904,4294908928,4294909952,4294910976,4294912000,
4294913024,4294914048,4294915072,4294916096,4294917120,4294918144,4294919168,
4294920192,4294921216,4294922240,4294923264,4294924288,4294925312,4294926336,
4294927360,4294928384,4294929408,4294930432,4294931456,4294932480,4294933504,
4294934528,4294935296,4294936320,4294937344,4294938368,4294939392,4294940416,
4294941440,4294942464,4294943488,4294944512,4294945536,4294946560,4294947584,
4294948608,4294949632,4294950656,4294951680,4294952704,4294953728,4294954752,
4294955776,4294956800,4294957824,4294958848,4294959872,4294960896,4294961920,
4294962944,4294963968,4294964992,4294966016,4294967040,4294704900,4294442760,
4294180620,4293918480,4293656340,4293394200,4293132060,4292869920,4292607780,
4292345640,4292083500,4291821360,4291559220,4291297080,4291034940,4290772800,
4290510660,4290248520,4289986380,4289724240,4289462100,4289199960,4288937820,
4288675680,4288413540,4288151400,4287889260,4287627120,4287364980,4287102840,
4286840700,4286644096,4286381955,4286119815,4285857675,4285595535,4285333395,
4285071255,4284809115,4284546975,4284284835,4284022695,4283760555,4283498415,
4283236275,4282974135,4282711995,4282449855,4282187715,4281925575,4281663435,
4281401295,4281139155,4280877015,4280614875,4280352735,4280090595,4279828455,
4279566315,4279304175,4279042035,4278779895,4278517755,4278255615,4278254591,
4278253567,4278252543,4278251519,4278250495,4278249471,4278248447,4278247423,
4278246399,4278245375,4278244351,4278243327,4278242303,4278241279,4278240255,
4278239231,4278238207,4278237183,4278236159,4278235135,4278234111,4278233087,
4278232063,4278231039,4278230015,4278228991,4278227967,4278226943,4278225919,
4278224895,4278223871,4278223103,4278222079,4278221055,4278220031,4278219007,
4278217983,4278216959,4278215935,4278214911,4278213887,4278212863,4278211839,
4278210815,4278209791,4278208767,4278207743,4278206719,4278205695,4278204671,
4278203647,4278202623,4278201599,4278200575,4278199551,4278198527,4278197503,
4278196479,4278195455,4278194431,4278193407,4278192383,4278191359,4278190335,
4278190331,4278190327,4278190323,4278190319,4278190315,4278190311,4278190307,
4278190303,4278190299,4278190295,4278190291,4278190287,4278190283,4278190279,
4278190275,4278190271,4278190267,4278190263,4278190259,4278190255,4278190251,
4278190247,4278190243,4278190239,4278190235,4278190231,4278190227,4278190223,
4278190219,4278190215,4278190211,4278190208};

	et = cudaMemcpyToSymbol( dvrgb, rgb, 256*sizeof(int), 0, cudaMemcpyHostToDevice); 
}

__global__ void 
create_image_on_gpu(unsigned int* g_odata, float* Ez, int nxx, float minval, float maxval)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int cind;
	float F;

	int ci = j*nxx+i;
	int ti = (j+1)*nxx+i;
	F = Ez[ti] - minval;
	cind = floor(255 * F/(maxval-minval));
	if (cind > 255) cind = 255;
	g_odata [ci] = dvrgb[cind];
}

extern "C" void
createImageOnGpu(unsigned int* g_odata)
{
    dim3 block(TILE_SIZE, TILE_SIZE, 1);
    dim3 grid(nxx/block.x, nyy/block.y, 1);
	dim3 gridm  = dim3(1,1,1);
	dim3 blockm = dim3(TILE_SIZE*TILE_SIZE,1,1);
	int  nblocks = grid.x * grid.y;
	float minval;
	float maxval;
	float *dvF;

	if (show_Ez) dvF=dvEz; else dvF=dvHz;

	find_min_and_max_on_gpu<<< gridm, blockm>>>(nblocks, dvF, dvminimum_field_value, dvmaximum_field_value);

	cudaMemcpy(&minval, dvminimum_field_value, sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&maxval, dvmaximum_field_value, sizeof(float), cudaMemcpyDeviceToHost);
	
	if (minval>0.0) minval = 0.0;
	if (maxval<0.0) maxval = 0.0;
	if (abs(minval)>maxval) maxval = -minval; else minval = -maxval;
	if (minval<global_min_field) global_min_field = minval;
	if (maxval>global_max_field) global_max_field = maxval;

	create_image_on_gpu<<< grid, block>>>(g_odata, dvF, nxx, global_min_field, global_max_field);
}

bool check_device()
{
	int currentDevice;
    cudaDeviceProp deviceProp;
	int multiProcessorCount;
	int compute_capability_major;
	int compute_capability_minor;

    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if( deviceCount == 0 )
    {
        printf("No gpu devices found!\n");
        return false;
    } 

	if (cudaGetDevice (&currentDevice) != cudaSuccess) 
	{
        printf("No gpu devices available!\n");
        return false;
	}

	if (cudaGetDeviceProperties(&deviceProp, currentDevice) != cudaSuccess) 
	{
        printf("No gpu devices available!\n");
        return false;
	}

	maximum_threads_per_block  = deviceProp.maxThreadsPerBlock;
	multiProcessorCount = deviceProp.multiProcessorCount;
	compute_capability_major = deviceProp.major;
	compute_capability_minor = deviceProp.minor;

	printf("Current Device %d: %s\n", currentDevice, deviceProp.name);
	printf("    Multi Processor Count: %d\n", multiProcessorCount);
	printf("    Maximum Threads per Block: %d\n", maximum_threads_per_block);
	printf("    Compute Capability: %d.%d \n", compute_capability_major, compute_capability_minor);

	return true;
}


int copyFdtdArraysToGpuMemory()
{
	int size_int   = sizeof(int);
	int size_char  = sizeof(char);
	int size_float = sizeof(float);
	int array_size  = number_of_cells * size_float; 
	int array_size_with_pads  = number_of_cells_with_pads * size_float; 

	cudaError_t et;	

	et = cudaMalloc((void**)&dvEx, array_size_with_pads);     if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvEy, array_size_with_pads);     if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvEz, array_size_with_pads);     if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvHx, array_size_with_pads);     if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvHy, array_size_with_pads);     if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvHz, array_size_with_pads);     if (et == cudaErrorMemoryAllocation) return 1;

	et = cudaMalloc((void**)&dvCexe,  array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvCexhz, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvCeye,  array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvCeyhz, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvCeze,  array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvCezhy, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvCezhx, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;

	et = cudaMalloc((void**)&dvChxh,  array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvChxez, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvChyh,  array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvChyez, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvChzh,  array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvChzex, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;
	et = cudaMalloc((void**)&dvChzey, array_size_with_pads);    if (et == cudaErrorMemoryAllocation) return 1;

	cudaMemcpy(dvEx, Ex, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvEy, Ey, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvEz, Ez, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvHx, Hx, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvHy, Hy, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvHz, Hz, array_size_with_pads, cudaMemcpyHostToDevice);

	cudaMemcpy(dvCexe,  Cexe,  array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvCexhz, Cexhz, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvCeye,  Ceye,  array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvCeyhz, Ceyhz, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvCeze,  Ceze,  array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvCezhy, Cezhy, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvCezhx, Cezhx, array_size_with_pads, cudaMemcpyHostToDevice);

	cudaMemcpy(dvChxh,  Chxh,  array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvChxez, Chxez, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvChyh,  Chyh,  array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvChyez, Chyez, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvChzh,  Chzh,  array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvChzex, Chzex, array_size_with_pads, cudaMemcpyHostToDevice);
	cudaMemcpy(dvChzey, Chzey, array_size_with_pads, cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvJ_indices, total_number_of_cej_components*sizeof(int));     if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvJ_indices,  J_indices, total_number_of_cej_components*sizeof(int), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvimpressed_J_Cej, total_number_of_cej_components*sizeof(int));     if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvimpressed_J_Cej,  impressed_J_Cej, total_number_of_cej_components*sizeof(int), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvM_indices, total_number_of_chm_components*sizeof(int));     if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvM_indices,  M_indices, total_number_of_chm_components*sizeof(int), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvimpressed_M_Chm, total_number_of_chm_components*sizeof(int));     if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvimpressed_M_Chm,  impressed_M_Chm, total_number_of_chm_components*sizeof(int), cudaMemcpyHostToDevice);


	et = cudaMalloc((void**)&dvsampled_magnetic_fields_component, number_of_sampled_magnetic_fields*sizeof(char));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_magnetic_fields_component,  sampled_magnetic_fields_component, number_of_sampled_magnetic_fields*sizeof(char), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvsampled_magnetic_fields_is, number_of_sampled_magnetic_fields*sizeof(int));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_magnetic_fields_is,  sampled_magnetic_fields_is, number_of_sampled_magnetic_fields*sizeof(int), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvsampled_magnetic_fields_js, number_of_sampled_magnetic_fields*sizeof(int));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_magnetic_fields_js,  sampled_magnetic_fields_js, number_of_sampled_magnetic_fields*sizeof(int), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvsampled_magnetic_fields_sampled_value, number_of_sampled_magnetic_fields*number_of_time_steps*sizeof(float));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_magnetic_fields_sampled_value, sampled_magnetic_fields_sampled_value, number_of_sampled_magnetic_fields*number_of_time_steps*sizeof(float), cudaMemcpyHostToDevice);


	et = cudaMalloc((void**)&dvsampled_electric_fields_component, number_of_sampled_electric_fields*sizeof(char));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_electric_fields_component,  sampled_electric_fields_component, number_of_sampled_electric_fields*sizeof(char), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvsampled_electric_fields_is, number_of_sampled_electric_fields*sizeof(int));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_electric_fields_is,  sampled_electric_fields_is, number_of_sampled_electric_fields*sizeof(int), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvsampled_electric_fields_js, number_of_sampled_electric_fields*sizeof(int));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_electric_fields_js,  sampled_electric_fields_js, number_of_sampled_electric_fields*sizeof(int), cudaMemcpyHostToDevice);

	et = cudaMalloc((void**)&dvsampled_electric_fields_sampled_value, number_of_sampled_electric_fields*number_of_time_steps*sizeof(float));     
	if (et == cudaErrorMemoryAllocation) return 1;
	cudaMemcpy(dvsampled_electric_fields_sampled_value, sampled_electric_fields_sampled_value, number_of_sampled_electric_fields*number_of_time_steps*sizeof(float), cudaMemcpyHostToDevice);

	if (is_cpml_xn)
	{
		et = cudaMemcpyToSymbol(dvcpml_a_ex_xn, cpml_a_ex_xn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_ex_xn, cpml_b_ex_xn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_a_mx_xn, cpml_a_mx_xn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_mx_xn, cpml_b_mx_xn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 

		if (is_TEz)
		{
			et = cudaMalloc((void**)&dvPsi_hzx_xn,  cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hzx_xn , Psi_hzx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hzx_xn, cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hzx_xn , CPsi_hzx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_eyx_xn,  cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_eyx_xn , Psi_eyx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_eyx_xn, cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_eyx_xn , CPsi_eyx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
		}
		if (is_TMz)
		{
			et = cudaMalloc((void**)&dvPsi_ezx_xn,  cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_ezx_xn , Psi_ezx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_ezx_xn, cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_ezx_xn , CPsi_ezx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_hyx_xn,  cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hyx_xn , Psi_hyx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hyx_xn, cpml_2d_array_size_xn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hyx_xn , CPsi_hyx_xn, cpml_2d_array_size_xn, cudaMemcpyHostToDevice);
		}
	}

	if (is_cpml_xp)
	{
		et = cudaMemcpyToSymbol(dvcpml_a_ex_xp, cpml_a_ex_xp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_ex_xp, cpml_b_ex_xp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_a_mx_xp, cpml_a_mx_xp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_mx_xp, cpml_b_mx_xp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 

		if (is_TEz)
		{
			et = cudaMalloc((void**)&dvPsi_hzx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hzx_xp , Psi_hzx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hzx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hzx_xp , CPsi_hzx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_eyx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_eyx_xp , Psi_eyx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_eyx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_eyx_xp , CPsi_eyx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
		}
		if (is_TMz)
		{
			et = cudaMalloc((void**)&dvPsi_ezx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_ezx_xp , Psi_ezx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_ezx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_ezx_xp , CPsi_ezx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_hyx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hyx_xp , Psi_hyx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hyx_xp, cpml_2d_array_size_xp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hyx_xp , CPsi_hyx_xp, cpml_2d_array_size_xp, cudaMemcpyHostToDevice);
		}
	}

	if (is_cpml_yn)
	{
		et = cudaMemcpyToSymbol(dvcpml_a_ey_yn, cpml_a_ey_yn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_ey_yn, cpml_b_ey_yn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_a_my_yn, cpml_a_my_yn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_my_yn, cpml_b_my_yn, cpml_1d_array_size_xnyn, 0, cudaMemcpyHostToDevice); 
		if (is_TEz)
		{
			et = cudaMalloc((void**)&dvPsi_hzy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hzy_yn , Psi_hzy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hzy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hzy_yn , CPsi_hzy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_exy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_exy_yn , Psi_exy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_exy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_exy_yn , CPsi_exy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
		}
		if (is_TMz)
		{
			et = cudaMalloc((void**)&dvPsi_ezy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_ezy_yn , Psi_ezy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_ezy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_ezy_yn , CPsi_ezy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_hxy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hxy_yn , Psi_hxy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hxy_yn, cpml_2d_array_size_yn);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hxy_yn , CPsi_hxy_yn, cpml_2d_array_size_yn, cudaMemcpyHostToDevice);
		}
	}

	if (is_cpml_yp)
	{
		et = cudaMemcpyToSymbol(dvcpml_a_ey_yp, cpml_a_ey_yp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_ey_yp, cpml_b_ey_yp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_a_my_yp, cpml_a_my_yp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 
		et = cudaMemcpyToSymbol(dvcpml_b_my_yp, cpml_b_my_yp, cpml_1d_array_size_xpyp, 0, cudaMemcpyHostToDevice); 
		if (is_TEz)
		{
			et = cudaMalloc((void**)&dvPsi_hzy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hzy_yp , Psi_hzy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hzy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hzy_yp , CPsi_hzy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_exy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_exy_yp , Psi_exy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_exy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_exy_yp , CPsi_exy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
		}
		if (is_TMz)
		{
			et = cudaMalloc((void**)&dvPsi_ezy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_ezy_yp , Psi_ezy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_ezy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_ezy_yp , CPsi_ezy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvPsi_hxy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_hxy_yp , Psi_hxy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
			et = cudaMalloc((void**)&dvCPsi_hxy_yp, cpml_2d_array_size_yp);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvCPsi_hxy_yp , CPsi_hxy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
		}
	}

	et = cudaMalloc((void**)&dvminimum_field_value, sizeof(float)*TILE_SIZE);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_ezy_yp , Psi_ezy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);
	et = cudaMalloc((void**)&dvmaximum_field_value, sizeof(float)*TILE_SIZE);	if (et == cudaErrorMemoryAllocation) return 1;	cudaMemcpy(dvPsi_ezy_yp , Psi_ezy_yp, cpml_2d_array_size_yp, cudaMemcpyHostToDevice);

	return 0;
}

void update_impressed_magnetic_currents_on_gpu(int time_step)
{
	dim3 threads_M, grid_M;
	float *dvH;

	int chm_ind = 0;
	for (int i=0;i<number_of_impressed_M;i++)
	{
		threads_M = dim3(number_of_chm_components[i], 1, 1);
		grid_M    = dim3(1, 1, 1);
		
		if (impressed_M_direction[i]=='x') dvH = dvHx;
		else if (impressed_M_direction[i]=='y')	dvH = dvHy; else dvH = dvHz;

		update_impressed_M_on_gpu<<<grid_M, threads_M>>>
			(dvH,  dvM_indices, dvimpressed_M_Chm, chm_ind, impressed_M_waveform[number_of_time_steps*i+time_step]);

		chm_ind = chm_ind + number_of_chm_components[i];
	}
}

void update_impressed_electric_currents_on_gpu(int time_step)
{
	dim3 threads_J, grid_J;
	float *dvE;

	int cej_ind = 0;
	for (int i=0;i<number_of_impressed_J;i++)
	{
		threads_J = dim3(number_of_cej_components[i], 1, 1);
		grid_J    = dim3(1, 1, 1);
		
		if (impressed_J_direction[i]=='x') dvE = dvEx;
		else if (impressed_J_direction[i]=='y')	dvE = dvEy; else dvE = dvEz;

		update_impressed_J_on_gpu<<<grid_J, threads_J>>>
			(dvE, dvJ_indices, dvimpressed_J_Cej, cej_ind, impressed_J_waveform[number_of_time_steps*i+time_step]);

		cej_ind = cej_ind + number_of_cej_components[i];
	}
}

void update_magnetic_fields_for_CPML_on_gpu()
{
	int n_bx = (nxx/TILE_SIZE) + (nxx%TILE_SIZE == 0 ? 0 : 1);
	int n_by = (nyy/TILE_SIZE) + (nyy%TILE_SIZE == 0 ? 0 : 1);
	dim3 threads_cpml = dim3(TILE_SIZE, TILE_SIZE, 1);
	dim3 grid_cpml_xn = dim3(1, n_by, 1);
	dim3 grid_cpml_xp = dim3(2, n_by, 1);
	dim3 grid_cpml_yn = dim3(n_bx, 1, 1);
	dim3 grid_cpml_yp = dim3(n_bx, 2, 1);

	if (is_cpml_xn)
	{
		if (is_TEz)
			update_magnetic_fields_on_gpu_CPML_TEz_xn<<<grid_cpml_xn, threads_cpml>>>(dvPsi_hzx_xn, dvCPsi_hzx_xn, dvHz, dvEy, nxx);

		if (is_TMz)
			update_magnetic_fields_on_gpu_CPML_TMz_xn<<<grid_cpml_xn, threads_cpml>>>(dvPsi_hyx_xn, dvCPsi_hyx_xn, dvHy, dvEz, nxx);
	}

	if (is_cpml_xp)
	{
		if (is_TEz)
			update_magnetic_fields_on_gpu_CPML_TEz_xp<<<grid_cpml_xp, threads_cpml>>>(dvPsi_hzx_xp, dvCPsi_hzx_xp, dvHz, dvEy, nxx, cpml_shift_xp);

		if (is_TMz)
			update_magnetic_fields_on_gpu_CPML_TMz_xp<<<grid_cpml_xp, threads_cpml>>>(dvPsi_hyx_xp, dvCPsi_hyx_xp, dvHy, dvEz, nxx, cpml_shift_xp);
	}

	if (is_cpml_yn)
	{
		if (is_TEz)
			update_magnetic_fields_on_gpu_CPML_TEz_yn<<<grid_cpml_yn, threads_cpml>>>(dvPsi_hzy_yn, dvCPsi_hzy_yn, dvHz, dvEx, nxx);

		if (is_TMz)
			update_magnetic_fields_on_gpu_CPML_TMz_yn<<<grid_cpml_yn, threads_cpml>>>(dvPsi_hxy_yn, dvCPsi_hxy_yn, dvHx, dvEz, nxx);
	}

	if (is_cpml_yp)
	{
		if (is_TEz)
			update_magnetic_fields_on_gpu_CPML_TEz_yp<<<grid_cpml_yp, threads_cpml>>>(dvPsi_hzy_yp, dvCPsi_hzy_yp, dvHz, dvEx, nxx, cpml_shift_yp);

		if (is_TMz)
			update_magnetic_fields_on_gpu_CPML_TMz_yp<<<grid_cpml_yp, threads_cpml>>>(dvPsi_hxy_yp, dvCPsi_hxy_yp, dvHx, dvEz, nxx, cpml_shift_yp);
	}
}

void update_electric_fields_for_CPML_on_gpu()
{
	int n_bx = (nxx/TILE_SIZE) + (nxx%TILE_SIZE == 0 ? 0 : 1);
	int n_by = (nyy/TILE_SIZE) + (nyy%TILE_SIZE == 0 ? 0 : 1);
	dim3 threads_cpml = dim3(TILE_SIZE, TILE_SIZE, 1);
	dim3 grid_cpml_xn = dim3(1, n_by, 1);
	dim3 grid_cpml_xp = dim3(2, n_by, 1);
	dim3 grid_cpml_yn = dim3(n_bx, 1, 1);
	dim3 grid_cpml_yp = dim3(n_bx, 2, 1);

	if (is_cpml_xn)
	{
		if (is_TEz)
			update_electric_fields_on_gpu_CPML_TEz_xn<<<grid_cpml_xn, threads_cpml>>>(dvPsi_eyx_xn, dvCPsi_eyx_xn, dvEy, dvHz, nxx);

		if (is_TMz)
			update_electric_fields_on_gpu_CPML_TMz_xn<<<grid_cpml_xn, threads_cpml>>>(dvPsi_ezx_xn, dvCPsi_ezx_xn, dvEz, dvHy, nxx);
	}

	if (is_cpml_xp)
	{
		if (is_TEz)
			update_electric_fields_on_gpu_CPML_TEz_xp<<<grid_cpml_xp, threads_cpml>>>(dvPsi_eyx_xp, dvCPsi_eyx_xp, dvEy, dvHz, nxx, cpml_shift_xp);

		if (is_TMz)
			update_electric_fields_on_gpu_CPML_TMz_xp<<<grid_cpml_xp, threads_cpml>>>(dvPsi_ezx_xp, dvCPsi_ezx_xp, dvEz, dvHy, nxx, cpml_shift_xp);
	}

	if (is_cpml_yn)
	{
		if (is_TEz)
			update_electric_fields_on_gpu_CPML_TEz_yn<<<grid_cpml_yn, threads_cpml>>>(dvPsi_exy_yn, dvCPsi_exy_yn, dvEx, dvHz, nxx);

		if (is_TMz)
			update_electric_fields_on_gpu_CPML_TMz_yn<<<grid_cpml_yn, threads_cpml>>>(dvPsi_ezy_yn, dvCPsi_ezy_yn, dvEz, dvHx, nxx);
	}

	if (is_cpml_yp)
	{
		if (is_TEz)
			update_electric_fields_on_gpu_CPML_TEz_yp<<<grid_cpml_yp, threads_cpml>>>(dvPsi_exy_yp, dvCPsi_exy_yp, dvEx, dvHz, nxx, cpml_shift_yp);

		if (is_TMz)
			update_electric_fields_on_gpu_CPML_TMz_yp<<<grid_cpml_yp, threads_cpml>>>(dvPsi_ezy_yp, dvCPsi_ezy_yp, dvEz, dvHx, nxx, cpml_shift_yp);

	}
}

bool fdtd_time_marching_loop_on_gpu()
{
	int n_bx = (nxx/TILE_SIZE) + (nxx%TILE_SIZE == 0 ? 0 : 1);
	int n_by = (nyy/TILE_SIZE) + (nyy%TILE_SIZE == 0 ? 0 : 1);
	dim3 threads = dim3(TILE_SIZE, TILE_SIZE, 1);
	dim3 grid = dim3(n_bx, n_by, 1);

	dim3 threads_sef = dim3(number_of_sampled_electric_fields, 1, 1);
	dim3 threads_smf = dim3(number_of_sampled_magnetic_fields, 1, 1);
	dim3 grid_sef = dim3(1, 1, 1);
	dim3 grid_smf = dim3(1, 1, 1);
	
	for (int time_step = 0; time_step<number_of_time_steps;time_step++) 
	{
		if (is_TEz)
			update_magnetic_fields_on_gpu_TEz<<<grid, threads>>>(dvChzh, dvChzex, dvChzey, dvHz, dvEx, dvEy, nxx);

		if (is_TMz)
			update_magnetic_fields_on_gpu_TMz<<<grid, threads>>>(dvChxh, dvChxez, dvChyh, dvChyez, dvHx,  dvHy, dvEz, nxx);

		update_impressed_magnetic_currents_on_gpu(time_step);
		
		update_magnetic_fields_for_CPML_on_gpu();

		capture_sampled_magnetic_fields_on_gpu<<<grid_smf, threads_smf>>>
			(dvHx, dvHy, dvHz, dvsampled_magnetic_fields_component, dvsampled_magnetic_fields_is, 
			dvsampled_magnetic_fields_js, dvsampled_magnetic_fields_sampled_value, time_step, number_of_time_steps, nxx);


		if (is_TEz)
		update_electric_fields_on_gpu_TEz<<<grid, threads>>>(dvCexe, dvCexhz, dvCeye, dvCeyhz,  dvEx, dvEy, dvHz, nxx);

		if (is_TMz)
		update_electric_fields_on_gpu_TMz<<<grid, threads>>>(dvCeze, dvCezhy, dvCezhx, dvHx,  dvHy, dvEz, nxx);
		
		update_impressed_electric_currents_on_gpu(time_step);
		
		update_electric_fields_for_CPML_on_gpu();

		capture_sampled_electric_fields_on_gpu<<<grid_sef, threads_sef>>>
			(dvEz, dvEy, dvEz, dvsampled_electric_fields_component, dvsampled_electric_fields_is, 
			dvsampled_electric_fields_js, dvsampled_electric_fields_sampled_value, time_step, number_of_time_steps, nxx);

		runIterationAndDisplay();
		
		if (time_step%100 == 0)
			printf("timestep: %d \n", time_step);
	}

	return true;
}
bool fdtdIterationOnGpu()
{
	int n_bx = (nxx/TILE_SIZE) + (nxx%TILE_SIZE == 0 ? 0 : 1);
	int n_by = (nyy/TILE_SIZE) + (nyy%TILE_SIZE == 0 ? 0 : 1);
	dim3 threads = dim3(TILE_SIZE, TILE_SIZE, 1);
	dim3 grid = dim3(n_bx, n_by, 1);

	dim3 threads_sef = dim3(number_of_sampled_electric_fields, 1, 1);
	dim3 threads_smf = dim3(number_of_sampled_magnetic_fields, 1, 1);
	dim3 grid_sef = dim3(1, 1, 1);
	dim3 grid_smf = dim3(1, 1, 1);
	
	if (is_TEz)
		update_magnetic_fields_on_gpu_TEz<<<grid, threads>>>(dvChzh, dvChzex, dvChzey, dvHz, dvEx, dvEy, nxx);

	if (is_TMz)
		update_magnetic_fields_on_gpu_TMz<<<grid, threads>>>(dvChxh, dvChxez, dvChyh, dvChyez, dvHx,  dvHy, dvEz, nxx);

	update_impressed_magnetic_currents_on_gpu(time_step);
	
	update_magnetic_fields_for_CPML_on_gpu();

	capture_sampled_magnetic_fields_on_gpu<<<grid_smf, threads_smf>>>
		(dvHx, dvHy, dvHz, dvsampled_magnetic_fields_component, dvsampled_magnetic_fields_is, 
		dvsampled_magnetic_fields_js, dvsampled_magnetic_fields_sampled_value, time_step, number_of_time_steps, nxx);

	if (is_TEz)
	update_electric_fields_on_gpu_TEz<<<grid, threads>>>(dvCexe, dvCexhz, dvCeye, dvCeyhz,  dvEx, dvEy, dvHz, nxx);

	if (is_TMz)
	update_electric_fields_on_gpu_TMz<<<grid, threads>>>(dvCeze, dvCezhy, dvCezhx, dvHx,  dvHy, dvEz, nxx);
	
	update_impressed_electric_currents_on_gpu(time_step);
	
	update_electric_fields_for_CPML_on_gpu();

	capture_sampled_electric_fields_on_gpu<<<grid_sef, threads_sef>>>
		(dvEz, dvEy, dvEz, dvsampled_electric_fields_component, dvsampled_electric_fields_is, 
		dvsampled_electric_fields_js, dvsampled_electric_fields_sampled_value, time_step, number_of_time_steps, nxx);

	if (time_step%100 == 0)
		printf("timestep: %d \n", time_step);
	
	time_step++;

	return true;
}

bool fetchResultsFromGpuMemory()
{
	cudaMemcpy(sampled_electric_fields_sampled_value, dvsampled_electric_fields_sampled_value, number_of_sampled_electric_fields*number_of_time_steps*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(sampled_magnetic_fields_sampled_value, dvsampled_magnetic_fields_sampled_value, number_of_sampled_magnetic_fields*number_of_time_steps*sizeof(float), cudaMemcpyDeviceToHost);

	return true;
}

bool deallocateCudaArrays()
{
	cudaFree(dvEx);
	cudaFree(dvEy);
	cudaFree(dvEz);
	cudaFree(dvHx);
	cudaFree(dvHy);
	cudaFree(dvHz);

	cudaFree(dvCexe);
	cudaFree(dvCexhz);
	cudaFree(dvCeye);
	cudaFree(dvCeyhz);
	cudaFree(dvCeze);
	cudaFree(dvCezhy);
	cudaFree(dvCezhx);

	cudaFree(dvChxh);
	cudaFree(dvChxez);
	cudaFree(dvChyh);
	cudaFree(dvChyez);
	cudaFree(dvChzh);
	cudaFree(dvChzex);
	cudaFree(dvChzey);

	cudaFree(dvimpressed_J_Cej);
	cudaFree(dvJ_indices);

	cudaFree(dvimpressed_M_Chm);
	cudaFree(dvM_indices);

	cudaFree(dvsampled_electric_fields_component);
	cudaFree(dvsampled_electric_fields_is);
	cudaFree(dvsampled_electric_fields_js);
	cudaFree(dvsampled_electric_fields_sampled_value);

	cudaFree(dvsampled_magnetic_fields_component);
	cudaFree(dvsampled_magnetic_fields_is);
	cudaFree(dvsampled_magnetic_fields_js);
	cudaFree(dvsampled_magnetic_fields_sampled_value);

	cudaFree(dvminimum_field_value);
	cudaFree(dvmaximum_field_value);

	return true;
}