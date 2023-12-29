from Bio import SeqIO
from Bio import SeqFeature

genes=["ABCB11","COG7","EMCN","ITIH5L"] # ITIH6 é sinónimo do ITIH5l

for gene in genes: # filtragem de informação do ficheiro de mRNA
    record=SeqIO.read(f"genes_information\{gene}_mRNA.gb",'genbank')
    name=record.name
    seq=record.seq
    tamanho=len(record.seq)
    des=record.description
    organismo=record.annotations['source']
    tipo_molecula=record.annotations["molecule_type"]
    features=record.features
    comentarios=record.annotations["comment"]
    feat_gene = [] 
    for I,feature in enumerate(features):
        if feature.type=='CDS':
            gene=(feature.qualifiers['gene'])
            notas=(feature.qualifiers['note'][0].split(";"))
            produto=feature.qualifiers['product']
           
            with open(f"{gene[0]}.txt","w") as f:
                f.writelines(f"Gene: {gene[0]}\nFonte do gene: {organismo}\n\n***INFORMACAO RELATIVA AO mRNA***\n\nTipo de molecula: {tipo_molecula}\nAccession number mRNA: {name}\n")
                f.writelines(f'\nTamanho do mRNA: {tamanho} bp\nDescricao do gene: {des}\nNumero de features: {len(features)}\n\nNotas:\n')
                for nota in notas:
                    f.writelines(f"- {nota.strip()}\n")
                f.writelines(f"\nProduto resultante: {produto[0]}\n")
                f.writelines(f"\nComentarios:{comentarios}\n")
        
        if feature.type=='exon':
           feat_gene.append(I)


    f=open (f"{gene[0]}.txt","a")
    f.writelines(f"\nNumero de exoes: {len(feat_gene)}\n\n")

    for I,indice in enumerate(feat_gene):
        localizacao=features[indice].location
        f.writelines(f"Exao {I+1}: {localizacao}\n") # formato do output [nucleotido inicial:nucleotido final] (strand)   

    f.writelines(f'\nSequencia mRNA:\n')

    for I in range(0,len(seq),100):
        seq_=seq[I:I+100]
        f.writelines(f"{seq_}\n")

    f.close()          

genes=["ABCB11","COG7","EMCN","ITIH6"] 

for gene in genes: #filtragem informação proteína
    f=open(f'{gene}.txt','a')
    if gene=="ITIH6":
        gene="ITIH5L"
    f.writelines('\n***INFORMACAO RELATIVA A PROTEINA***\n\n')
    record_prot=SeqIO.read(f"genes_information\{gene}_protein.gb",'genbank')
    id=record_prot.name
    seq=record_prot.seq
    tamanho=len(record_prot.seq)
    des=record_prot.description
    tipo_molecula=record_prot.annotations["molecule_type"]
    features=record_prot.features
    f.writelines(f'Tipo de molecula: {tipo_molecula}\nAccession number proteina: {id}\nTamanho proteina: {len(seq)} aa\n')
    f.writelines(f'Descricao da proteina: {des}\n')
    feat_site=[]
    for I,feature in enumerate(features):
        if feature.type=='Protein':
            peso_molecualr=feature.qualifiers["calculated_mol_wt"][0]
            f.writelines(f'Peso molecular: {peso_molecualr} Dalton\nNumero de features da proteina: {len(features)}\n')
        if feature.type=='Site':    
            feat_site.append(I)

    if len(feat_site)!=0:
        f.writelines('\nLocais de interesse da proteina:\n\n')
        for indice in feat_site:
            feature=features[indice]        
            localizacao=feature.location
            notas_site=feature.qualifiers["note"]
            site_type=feature.qualifiers['site_type']
            f.writelines(f'Site {indice+1}:\n \t-Localizacao: {localizacao}\n\t-Notas: {notas_site[0]}\n\t-Tipo: {site_type[0]}\n')
    
    f.writelines(f'\nSequencia proteina:\n')
    for I in range(0,len(seq),100):
        seq_=seq[I:I+100]
        f.writelines(f"{seq_}\n")

    f.close()
