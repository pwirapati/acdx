.onLoad <- function(libname,pkgname)
{
	Romp_set_num_threads( Romp_get_num_procs()/2 )
}
