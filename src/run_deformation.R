library(bio3d)

args<-commandArgs(TRUE)
pdb_file = args[1]
res_number = as.numeric(args[2])
chain = args[3]
force_field = args[4]

pdb <- read.pdb(pdb_file,rm.alt=FALSE)
pdb_parsed = trim.pdb(pdb,atom.select(pdb,chain=chain))

wt_ind <- atom.select(pdb_parsed,resno=res_number,chain=chain,elety="CA")
all_inds <- atom.select(pdb_parsed,elety="CA")
modes_position = which(all_inds$atom == wt_ind$atom)

if(length(modes_position) != 1){
    if(length(wt_ind$atom) > 1){
        wt_position <- wt_ind$atom[1]
        modes_position = which(all_inds$atom == wt_position)
    }
}

modes <- nma(pdb_parsed,ff=force_field)

# Deformation energies
defe <- deformation.nma(modes)
defsums <- rowSums(defe$ei[,7:12])

cat(defsums[modes_position])
