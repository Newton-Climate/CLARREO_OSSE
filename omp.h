/**********************************************************************/
/*                                                                    */
/* COMPONENT_NAME:  smprt (SMP Runtime Library)                       */
/* COMPONENT_ID:    5765-F71-00, 5765-F57-00, 5765-F56-00             */
/*                                                                    */
/* Licensed Materials - Property of IBM                               */
/* (C) Copyright IBM Corp. 1999, 2002. All Rights Reserved.           */
/*                                                                    */
/* US Government Users Restricted Rights - Use, duplication or        */
/* disclosure restricted by GSA ADP Schedule Contract with IBM Corp.  */
/*                                                                    */
/**********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif


#if !defined(_OMP_H_)
#define	_OMP_H_
/* The #endif for this is at the end of the file. */


#include <pthread.h>


typedef int omp_lock_t;

typedef struct {
  int       lock;
  int       count;
  pthread_t thread;
} omp_nest_lock_t;


#if defined(_OPENMP) || defined(_IBMSMP) 

void omp_set_num_threads(int num);
int  omp_get_num_threads(void);
int  omp_get_max_threads(void);
int  omp_get_thread_num(void);
int  omp_get_num_procs(void);
int  omp_in_parallel(void);

void omp_set_dynamic(int flag);
int  omp_get_dynamic(void);

void omp_set_nested(int flag);
int  omp_get_nested(void);

void omp_init_lock(omp_lock_t *lock);
void omp_init_nest_lock(omp_nest_lock_t *lock);

void omp_destroy_lock(omp_lock_t *lock);
void omp_destroy_nest_lock(omp_nest_lock_t *lock);

void omp_set_lock(omp_lock_t *lock);
void omp_set_nest_lock(omp_nest_lock_t *lock);

void omp_unset_lock(omp_lock_t *lock);
void omp_unset_nest_lock(omp_nest_lock_t *lock);

int  omp_test_lock(omp_lock_t *lock);
int  omp_test_nest_lock(omp_nest_lock_t *lock);

#else 

#if __cplusplus
#define __omp_inline  static inline
#else
#define __omp_inline  static __inline
#endif

__omp_inline void omp_set_num_threads(int num) {
}

__omp_inline int omp_get_num_threads(void) {
  return 1;
}

__omp_inline int omp_get_max_threads(void) {
  return 1;
}

__omp_inline int omp_get_thread_num(void) {
  return 0;
}

__omp_inline int omp_get_num_procs(void) {
  return 1;
}

__omp_inline int omp_in_parallel(void) {
  return 0;
}

__omp_inline void omp_set_dynamic(int flag) {
}

__omp_inline int omp_get_dynamic(void) {
  return 0;
}

__omp_inline void omp_set_nested(int flag) {
}

__omp_inline int omp_get_nested(void) {
  return 0;
}

__omp_inline void omp_init_lock(omp_lock_t  *lock) {
}

__omp_inline void omp_init_nest_lock(omp_nest_lock_t *lock) {
  lock->count = 0;
}

__omp_inline void omp_destroy_lock(omp_lock_t *lock) {
}

__omp_inline void omp_destroy_nest_lock(omp_nest_lock_t *lock) {
}

__omp_inline void omp_set_lock(omp_lock_t *lock) {
}

__omp_inline void omp_set_nest_lock(omp_nest_lock_t *lock) {
  lock->count++;
}

__omp_inline void omp_unset_lock(omp_lock_t *lock) {
}

__omp_inline void omp_unset_nest_lock(omp_nest_lock_t *lock) {
  lock->count--;
}

__omp_inline int omp_test_lock(omp_lock_t *lock) {
  return 1;
}

__omp_inline int omp_test_nest_lock(omp_nest_lock_t *lock) {
  lock->count++;
  return lock->count;
}

#undef __omp_inline
#endif  /* (_OPENMP) || defined(_IBMSMP) */

#endif  /* !defined(_OMP_H_) */


#ifdef __cplusplus
}
#endif

