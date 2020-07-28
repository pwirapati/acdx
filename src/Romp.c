/*
 * OpenMP controls callable from R functions
 *
 */

#include <omp.h>

void
Romp_get_num_procs(int *num_procs_)
{
  *num_procs_ = omp_get_num_procs();
}

void
Romp_get_max_threads(int *num_threads_)
{
  *num_threads_ = omp_get_max_threads();
}

void
Romp_set_num_threads(int *num_threads_)
{
  omp_set_num_threads(*num_threads_);
  *num_threads_ = omp_get_max_threads();
}

void
Romp_get_schedule(int *kind_, int *chunk_size_ )
{
  omp_sched_t kind;
  omp_get_schedule( &kind, chunk_size_ );
  switch(kind)
    {
    case omp_sched_static:
      *kind_ = 1;
      break;
    case omp_sched_dynamic:
      *kind_ = 2;
      break;
    case omp_sched_guided:
      *kind_ = 3;
      break;
    case omp_sched_auto:
      *kind_ = 4;
      break;
    default:
      *kind_ = 0; // we don't know
      break;
    }
}

void
Romp_set_schedule(int *kind_, int *chunk_size_)
{
  omp_sched_t kind;
  switch(*kind_)
    {
    case 1:
      kind = omp_sched_static;
      break;
    case 2:
      kind = omp_sched_dynamic;
      break;
    case 3:
      kind = omp_sched_guided;
      break;
    case 4:
      kind = omp_sched_auto;
      break;
    default:
      kind = omp_sched_auto;
    }

  omp_set_schedule(kind, *chunk_size_);
  Romp_get_schedule( kind_, chunk_size_ );
}
