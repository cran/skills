.onLoad = function(libname,pkgname){
	sets_options("quote",FALSE)
}
.onUnload = function(libname,pkgname){
	sets_options("quote",TRUE)	
}
