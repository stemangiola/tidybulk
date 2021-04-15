#' @importFrom utils packageDescription
.onAttach = function(libname, pkgname) {
	version = packageDescription(pkgname, fields = "Version")
	
	msg = paste0("========================================
", pkgname, " version ", version, "
If you use TIDYBULK in published research, please cite:

Mangiola et al. tidybulk: an R tidy framework for modular 
transcriptomic data analysis. Genome Biology 2021.

This message can be suppressed by:
  suppressPackageStartupMessages(library(tidybulk))
========================================
")	
	
	packageStartupMessage(msg)
}

rv = R.Version()

if(getRversion() >= "4.0.0" && as.numeric(rv$`svn rev`) >= 77889) {
	unitType = get("unitType", envir = asNamespace("grid"))
} else {
	unitType = function(x, recurse = TRUE) attr(x, "unit")
}