#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

int main(void)
{
	char size_str[] = "DOMAIN_SIZE";

	max_file_t *maxfile = GOcean_init();
	max_engine_t *engine = max_load(maxfile, "*");

    long size = max_get_constant_uint64t(maxfile, size_str);
    size *= size;
	//int size = 256*256;

	// lmem has to be aligned to the burst size
	const int lmem_align = 48; // *8 = 384 bytes
	int lmem_size = size;
	if(size%lmem_align != 0){
		lmem_size = lmem_align*(size/lmem_align + 1);
	}

	int sizeBytes = lmem_size * sizeof(double);
	double *u = malloc(sizeBytes);
	double *uold = malloc(sizeBytes);
	double *v = malloc(sizeBytes);
	double *vold = malloc(sizeBytes);
	double *p = malloc(sizeBytes);
	double *pold = malloc(sizeBytes);
	// Output buffers
	double *unew = malloc(sizeBytes);
	double *vnew = malloc(sizeBytes);
    double *pnew = malloc(sizeBytes);
	double *uold_new = malloc(sizeBytes);
	double *vold_new = malloc(sizeBytes);
    double *pold_new = malloc(sizeBytes);

	// TODO Generate input data using stream function as in original shallow
	for(int i = 0; i < lmem_size; ++i) {
		u[i] = (double)(random())/(double)RAND_MAX;
		uold[i] = u[i];
		v[i] = (double)(random())/(double)RAND_MAX;
		vold[i] = v[i];
		p[i] = (double)(random())/(double)RAND_MAX;
		pold[i] = p[i];

		unew[i] = (double)i;
	}


	GOcean_actions_t run_scalar;

	run_scalar.param_N = size;
    run_scalar.param_NLmem = lmem_size;

	// Input buffers - stream in field{,_old}

	// Data we stream straight into kernels
	run_scalar.instream_p = p;
	run_scalar.instream_pold = pold;
	run_scalar.instream_u = u;
	run_scalar.instream_uold = uold;
	run_scalar.instream_v = v;
	run_scalar.instream_vold = vold;

	// Output buffers - stream out field{,_old}
	run_scalar.outstream_pout = pnew;
	run_scalar.outstream_pold_out = pold_new;
	run_scalar.outstream_vout = vnew;
	run_scalar.outstream_vold_out = vold_new;
	run_scalar.outstream_uout = unew;
	run_scalar.outstream_uold_out = uold_new;

	// Routing string - now set in the default interface because it seems I need
	// to set 'ignoreall' (in the interface) otherwise my streams have an undefined size.
//	run_scalar.routing_string = "PressureFan -> p1, PressureFan -> p2, "
//			"PressureFan -> p3, PressureFan -> p4, PressureFan -> p5, "
//			"uFan -> u1, uFan -> u2, uFan -> u3, uFan -> u4, "
//			"vFan -> v1, vFan -> v2, vFan -> v3, vFan -> v4, "
//			"cvFan -> cvout1, cvFan -> cvout2, "
//			"cuFan -> cuout1, cuFan -> cuout2, "
//			"zFan -> zout1, zFan -> zout2, "
//			"hFan -> hout1, hFan -> hout2, "
//			"pnewFan -> pnewout, unewFan -> unewout, vnewFan -> vnewout, "
//			"uoldFan -> uold1, uoldFan -> uold2, "
//	        "voldFan -> vold1, voldFan -> vold2, "
//	        "poldFan -> pold1 ";

	// Write to LMEM
	printf("Writing %d bytes to LMEM\n", sizeBytes);
	GOcean_writeLMem_actions_t write_s;
	write_s.param_start = 0;
	write_s.param_size  = lmem_size;
	write_s.instream_cpu2lmem = unew;

	GOcean_writeLMem_run(engine, &write_s);

	printf("Running on DFE (interface: GOcean_run)...\n");

	GOcean_run(engine, &run_scalar);

	// Read back from LMEM
	printf("Reading %d bytes back from LMEM, sizeBytes\n", sizeBytes);
	GOcean_readLMem_actions_t read_s;
	GOcean_readLMem(lmem_size, 0, unew);

	max_unload(engine);
	printf("Done.\n");
	return 0;
}
