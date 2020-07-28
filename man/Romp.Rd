\name{Romp}
\alias{Romp}
\alias{Romp_get_num_procs}
\alias{Romp_get_max_threads}
\alias{Romp_set_num_threads}
\alias{Romp_get_schedule}
\alias{Romp_set_schedule}

\title{R OpenMP utilities}
\description{Interface to some OpenMP functions for controlling
runtime parameters}
\usage{
  Romp_get_num_procs()
  Romp_get_max_threads()
  Romp_set_num_threads(n)
  Romp_get_schedule()
  Romp_set_schedule(kind,chunk_size)
}
\details{
  \code{Romp_get_num_procs} return the number of processors. On Intel
  under linux, the number includes virtual processors provided by
  hyperthreading.

  \code{Romp_get_max_threads} return the number of threads that will
  be used by the next parallel section. This number may be the same
  as \code{Romp_get_num_procs} by default, but modified by environment
  variable \code{OMP_NUM_THREADS} set before running R. It is also
  modified by \code{Romp_set_num_threads}.

  \code{Romp_set_num_threads} set the number of threads used by
  subsequent parallel sections.

  \code{Romp_get_schedule} returns a list of the \code{kind} of schedule
  (one of strings \code{static}, \code{dynamic}, \code{guided} or \code{auto}),
  and the \code{chunk_size} parameter associated with the kind
  of scheduling. See OpenMP manual for details.

  \code{Romp_set_schedule} sets the schedule.
}
\value{
  \code{Romp_get_num_procs}, \code{Romp_get_max_threads}, and
  \code{Romp_set_num_threads} return an integer of the number
  of threads.
 
  \code{Romp_get_schedule} and \code{Romp_set_schedule}
  return a list of \code{kind}
  (a string of schedule name) and \code{chunk_size} (an integer).
}
