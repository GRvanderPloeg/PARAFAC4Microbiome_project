function dump(Models, Cons, VarExps, Boots, BootVarExps, Tuckers, path_start, prefix)

writecell(Models, path_start+prefix+"Models.csv");
writecell(Cons, path_start+prefix+"Cons.csv");
writecell(VarExps, path_start+prefix+"VarExps.csv");
writecell(Boots, path_start+prefix+"Boots.csv");
writecell(BootVarExps, path_start+prefix+"BootVarExps.csv");
writematrix(Tuckers{1}, path_start+prefix+"Tuckers1.csv");
writecell(Tuckers{2}, path_start+prefix+"Tuckers2.csv");
writecell(Tuckers{3}, path_start+prefix+"Tuckers3.csv");