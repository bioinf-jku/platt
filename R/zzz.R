# Copyright (C) 2014 Klambauer Guenter 


.onLoad <- function(libname, pkgname) {
	library.dynam("platt", pkgname, libname)
}

.onUnload <- function(libpath)
{
	library.dynam.unload("platt", libpath)
}

