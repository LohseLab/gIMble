parallel -j1 'bedtools sort -i {} | bedtools merge -i - -c 4 -d -1 -o distinct | grep -vP "[,|\|]" > {.}.4D.bed' ::: *degeneracy.bed
