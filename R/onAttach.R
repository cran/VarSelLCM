.onAttach <-function(lib,pkg){
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),
  "Version")
  packageStartupMessage("\nPackage 'VarSelLCM' version ",
as.character(ver),
" loaded\nPlease report bugs at: matthieu.marbac_at_inserm.fr or mohammed.sedki_at_inserm.fr \n",
"Type 'citation(\"VarSelLCM\")' for citing this package.")
}