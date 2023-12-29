from Bio import Entrez

Entrez.email = "pg53501@alunos.uminho.pt"

genes=["ABCB11","COG7","EMCN", "ITIH5L"] #genes em estudo
nucleotide_id=["NG_007374.2","NG_021287.1","NA","NG_013240.1"] # informação a nível do genoma sobre o respetivo gene
"""
Na lista nucleotide_id existe um valor igual a NA e isto acontece, pois na base de dados nucleotide não aparece informação
"""
mRNA_id=["NM_003742.4","NM_153603.4","NM_001159694.2","NM_198510.3"] # informação sobre o mRNA para os respetivos genes
protein_id=["NP_003733.2","NP_705831.1","NP_001153166.1","NP_940912.1"] #informação sobre a proteina codificada pelo gene

for I in range(len(genes)):
    if nucleotide_id[I]=="NA":
        pass
    else:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id[I], rettype="gb", retmode="text")
        doc=open(f'genes_information\{genes[I]}.gb',"w")
        doc.writelines(handle)
        doc.close()
        handle.close()

    handle = Entrez.efetch(db="nucleotide", id=mRNA_id[I], rettype="gb", retmode="text")
    doc=open(f'genes_information\{genes[I]}_mRNA.gb',"w")
    doc.writelines(handle)
    doc.close()
    handle.close()

    handle = Entrez.efetch(db="Protein", id=protein_id[I], rettype="gb", retmode="text")
    doc=open(f'genes_information\{genes[I]}genes[I]_protein.gb',"w")
    doc.writelines(handle)
    doc.close()
    handle.close()