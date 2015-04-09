/**\file */
#ifndef SLIC_DECLARATIONS_GOcean_H
#define SLIC_DECLARATIONS_GOcean_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define GOcean_PCIE_ALIGNMENT (16)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface readLMem ----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'readLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_lmem2cpu The stream should be of size (param_size * 8) bytes.
 */
void GOcean_readLMem(
	int32_t param_size,
	int32_t param_start,
	double *outstream_lmem2cpu);

/**
 * \brief Basic static non-blocking function for the interface 'readLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_lmem2cpu The stream should be of size (param_size * 8) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *GOcean_readLMem_nonblock(
	int32_t param_size,
	int32_t param_start,
	double *outstream_lmem2cpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'readLMem'
 * 
 */
typedef struct { 
	int32_t param_size; /**<  [in] Interface Parameter "size". */
	int32_t param_start; /**<  [in] Interface Parameter "start". */
	double *outstream_lmem2cpu; /**<  [out] The stream should be of size (param_size * 8) bytes. */
} GOcean_readLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void GOcean_readLMem_run(
	max_engine_t *engine,
	GOcean_readLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'readLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_readLMem_run_nonblock(
	max_engine_t *engine,
	GOcean_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void GOcean_readLMem_run_group(max_group_t *group, GOcean_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_readLMem_run_group_nonblock(max_group_t *group, GOcean_readLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void GOcean_readLMem_run_array(max_engarray_t *engarray, GOcean_readLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_readLMem_run_array_nonblock(max_engarray_t *engarray, GOcean_readLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* GOcean_readLMem_convert(max_file_t *maxfile, GOcean_readLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*--------------------------- Interface writeLMem ----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'writeLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_cpu2lmem The stream should be of size (param_size * 8) bytes.
 */
void GOcean_writeLMem(
	int32_t param_size,
	int32_t param_start,
	const double *instream_cpu2lmem);

/**
 * \brief Basic static non-blocking function for the interface 'writeLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_cpu2lmem The stream should be of size (param_size * 8) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *GOcean_writeLMem_nonblock(
	int32_t param_size,
	int32_t param_start,
	const double *instream_cpu2lmem);

/**
 * \brief Advanced static interface, structure for the engine interface 'writeLMem'
 * 
 */
typedef struct { 
	int32_t param_size; /**<  [in] Interface Parameter "size". */
	int32_t param_start; /**<  [in] Interface Parameter "start". */
	const double *instream_cpu2lmem; /**<  [in] The stream should be of size (param_size * 8) bytes. */
} GOcean_writeLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void GOcean_writeLMem_run(
	max_engine_t *engine,
	GOcean_writeLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'writeLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_writeLMem_run_nonblock(
	max_engine_t *engine,
	GOcean_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void GOcean_writeLMem_run_group(max_group_t *group, GOcean_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_writeLMem_run_group_nonblock(max_group_t *group, GOcean_writeLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void GOcean_writeLMem_run_array(max_engarray_t *engarray, GOcean_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_writeLMem_run_array_nonblock(max_engarray_t *engarray, GOcean_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* GOcean_writeLMem_convert(max_file_t *maxfile, GOcean_writeLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_N Interface Parameter "N".
 * \param [in] instream_p The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_pold The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_u The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_uold The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_v The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_vold The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_pold_out The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_pout The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_uold_out The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_uout The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_vold_out The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_vout The stream should be of size (param_N * 8) bytes.
 */
void GOcean(
	int64_t param_N,
	const double *instream_p,
	const double *instream_pold,
	const double *instream_u,
	const double *instream_uold,
	const double *instream_v,
	const double *instream_vold,
	double *outstream_pold_out,
	double *outstream_pout,
	double *outstream_uold_out,
	double *outstream_uout,
	double *outstream_vold_out,
	double *outstream_vout);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_N Interface Parameter "N".
 * \param [in] instream_p The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_pold The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_u The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_uold The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_v The stream should be of size (param_N * 8) bytes.
 * \param [in] instream_vold The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_pold_out The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_pout The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_uold_out The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_uout The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_vold_out The stream should be of size (param_N * 8) bytes.
 * \param [out] outstream_vout The stream should be of size (param_N * 8) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *GOcean_nonblock(
	int64_t param_N,
	const double *instream_p,
	const double *instream_pold,
	const double *instream_u,
	const double *instream_uold,
	const double *instream_v,
	const double *instream_vold,
	double *outstream_pold_out,
	double *outstream_pout,
	double *outstream_uold_out,
	double *outstream_uout,
	double *outstream_vold_out,
	double *outstream_vout);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	int64_t param_N; /**<  [in] Interface Parameter "N". */
	const double *instream_p; /**<  [in] The stream should be of size (param_N * 8) bytes. */
	const double *instream_pold; /**<  [in] The stream should be of size (param_N * 8) bytes. */
	const double *instream_u; /**<  [in] The stream should be of size (param_N * 8) bytes. */
	const double *instream_uold; /**<  [in] The stream should be of size (param_N * 8) bytes. */
	const double *instream_v; /**<  [in] The stream should be of size (param_N * 8) bytes. */
	const double *instream_vold; /**<  [in] The stream should be of size (param_N * 8) bytes. */
	double *outstream_pold_out; /**<  [out] The stream should be of size (param_N * 8) bytes. */
	double *outstream_pout; /**<  [out] The stream should be of size (param_N * 8) bytes. */
	double *outstream_uold_out; /**<  [out] The stream should be of size (param_N * 8) bytes. */
	double *outstream_uout; /**<  [out] The stream should be of size (param_N * 8) bytes. */
	double *outstream_vold_out; /**<  [out] The stream should be of size (param_N * 8) bytes. */
	double *outstream_vout; /**<  [out] The stream should be of size (param_N * 8) bytes. */
} GOcean_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void GOcean_run(
	max_engine_t *engine,
	GOcean_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_run_nonblock(
	max_engine_t *engine,
	GOcean_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void GOcean_run_group(max_group_t *group, GOcean_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_run_group_nonblock(max_group_t *group, GOcean_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void GOcean_run_array(max_engarray_t *engarray, GOcean_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *GOcean_run_array_nonblock(max_engarray_t *engarray, GOcean_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* GOcean_convert(max_file_t *maxfile, GOcean_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* GOcean_init(void);

/* Error handling functions */
int GOcean_has_errors(void);
const char* GOcean_get_errors(void);
void GOcean_clear_errors(void);
/* Free statically allocated maxfile data */
void GOcean_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int GOcean_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int GOcean_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_GOcean_H */

