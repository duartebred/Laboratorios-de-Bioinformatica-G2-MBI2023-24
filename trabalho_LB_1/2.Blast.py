from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import SeqIO

NCBIWWW.email = "PG53481@alunos.uminho.pt"
NCBIXML.email = "PG53481@alunos.uminho.pt"

genes=["ABCB11","COG7","EMCN", "ITIH5L"] # lista de genes em estudo

#Realizacao do blastp
# Este ciclo for realiza o blast para as proteínas codificadas pelos genes de interesse
#Optamos por usar a base de dados swissprot para restringir a procura apenas a proteinas conhecidas e com as suas caracteristicas estudadas

for gene in genes: 
    record=SeqIO.read(f'genes_information\{gene}_protein.gb','genbank')
    sequencia=record.seq
    result=NCBIWWW.qblast('blastp','swissprot',sequencia,hitlist_size=15,expect=5)
    save_result=open(f"{gene}_blast.xml","w")
    save_result.writelines(result.read())
    save_result.close()
    result.close()

#realizacao do blast do gene EMCN alterando a base de dados 
#na swissport não existe muita informacao para este gene o que vai causar problemas no alinhamento multiplo
record=SeqIO.read(f'genes_information\EMCN_protein.gb','genbank')
sequencia=record.seq
result=NCBIWWW.qblast('blastp','nr',sequencia,hitlist_size=10,expect=5)
save_result=open(f"EMCN_2_blast.xml","w")
save_result.writelines(result.read())
save_result.close()
result.close()

genes=["ABCB11","COG7","EMCN","EMCN_2","ITIH5L"] # EMCN_2 nao e um gene novo, mas sim o resultado do blast para a base de dados nr para o gene EMCN_2

#Escolha dos 5 melhores resultados
for gene in genes:
    if gene=="EMCN_2":
        gene_=genes[2]
        record=SeqIO.read(f'genes_information/{gene_}_protein.gb','genbank')
        
    else:
        record=SeqIO.read(f'genes_information/{gene}_protein.gb','genbank')
    
    id=record.id
    
    f=open(f'Blast_results\{gene}_blast_results.txt','w')

    result_handle=open(f'{gene}_blast.xml','r')
    record_blast=NCBIXML.read(result_handle)
    result_handle.close()

    matriz=record_blast.matrix
    gap_pen=record_blast.gap_penalties
    database=record_blast.database
    n_hits=len(record_blast.alignments)
    f.writelines('***RESULTADOS DO BLAST PARA A PROTEINA CODIFICADA PELO GENE {}***\n\n'.format(gene))
    f.write(f"Accession number (NCBI) da proteina codificada pelo gene {gene}: {id}\nNumero de hits: {n_hits}\n")
    f.writelines(f'Matriz usada no alinhamento: {matriz}\nPenalidades de espacamentos (abertura,propagacao): {gap_pen}\nBase de proteinas usada: {database}\n') 

    alinhamentos=record_blast.alignments

    if len(alinhamentos)>=5:
        f.write('\n---MELHORES ALINHAMENTOS---\n\n')
        for I in range(0,6):
                if I==0:
                    alinhamento=alinhamentos[I]
                    accession_match=alinhamento.accession 
                    hit_def=alinhamento.hit_def
                    f.write(f'Alinhameto {I+1}: Proteina codifica pelo gene {gene}\n')
                    f.writelines(f'\tAccession number (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\n')

                else:
                    alinhamento=alinhamentos[I]
                    accession_match=alinhamento.accession 
                    hit_def=alinhamento.hit_def 
                    hit_len=alinhamento.length
                    f.write(f'Alinhamento {I+1}\n')
                    f.writelines(f'\tAccession number hit (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\tTamanho do hit: {hit_len} aa\n')
                    hsps=alinhamento.hsps
                    f.write(f'\tNumero de hsps: {len(hsps)}\n')
                    for N,hsp in enumerate(hsps):
                        E_value=hsp.expect
                        Score=hsp.score
                        Length= hsp.align_length
                        identities= hsp.identities
                        positives= hsp.positives
                        gaps=hsp.gaps
                        f.writelines(f'\t\tNumero hsp: {N+1}\n\t\t\tE-value: {E_value}\n\t\t\tScore: {Score}\n\t\t\tTamanho hsp: {Length} aa\n')
                        f.writelines(f'\t\t\tNumero de aa identicos: {identities}\n\t\t\tNumero de aa identicos e positivos: {positives}\n\t\t\tNumero de gaps: {gaps}\n')           
                    f.write('\n')

    else:
        f.write('\n---MELHORES ALINHAMENTOS---\n\n')
        for I,alinhamento in enumerate(alinhamentos):
            if I==0:
                accession_match=alinhamento.accession 
                hit_def=alinhamento.hit_def
                f.write(f'Alinhameto {I+1}: Proteina codifica pelo gene {gene}\n')
                f.writelines(f'\tAccession number (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\n')

            else:
                accession_match=alinhamento.accession 
                hit_def=alinhamento.hit_def 
                hit_len=alinhamento.length
                f.write(f'Alinhamento {I+1}\n')
                f.writelines(f'\tAccession number hit (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\tTamanho do hit: {hit_len} aa\n')
                hsps=alinhamento.hsps
                f.write(f'\tNumero de hsps: {len(hsps)}\n')
                for N,hsp in enumerate(hsps):
                    E_value=hsp.expect
                    Score=hsp.score
                    Length= hsp.align_length
                    identities= hsp.identities
                    positives= hsp.positives
                    gaps=hsp.gaps
                    f.writelines(f'\t\tNumero hsp: {N+1}\n\t\t\tE-value: {E_value}\n\t\t\tScore: {Score}\n\t\t\tTamanho hsp: {Length} aa\n')
                    f.writelines(f'\t\t\tNumero de aa identicos: {identities}\n\t\t\tNumero de aa identicos e positivos: {positives}\n\t\t\tNumero de gaps: {gaps}\n')           
                f.write('\n')    
    
    f.close()

# Leitura dos resultados do Blast
for gene in genes:
    
    f=open(f'Blast_results\{gene}_blast_results.txt','a')

    result_handle=open(f'{gene}_blast.xml','r')
    record_blast=NCBIXML.read(result_handle)
    result_handle.close()
    
    alinhamentos=record_blast.alignments
    f.write('\n---ALINHAMENTOS---\n\n')

    for I,alinhamento in enumerate(alinhamentos):
        f.write(f'Numero alinhamento: {I+1}\n')
        accession_match=alinhamento.accession #accession number da sequencia que faz match com a query
        hit_def=alinhamento.hit_def #definicao da sequencia match
        hit_len=alinhamento.length
        n_hsps=len(alinhamento.hsps)
        hsps=alinhamento.hsps
        
        if I==0: 
            f.writelines(f'\tProteina codificada pelo gene accession number (relativo a base de dados usado no blast): {accession_match}\n')
            f.writelines(f'\tDefinicao: {hit_def}\n\tTamanho da proteina: {hit_len} aa\n\tNumero de hsps: {n_hsps}\n')
            
        else:
            f.writelines(f'\tHit accession number (relativo a base de dados usada no blast): {accession_match}\n')
            f.writelines(f'\tDefinicao: {hit_def}\n\tTamanho do hit: {hit_len} aa\n\tNumero de hsps: {n_hsps}\n')

        for N,hsp in enumerate(hsps):
            E_value=hsp.expect
            Score=hsp.score
            Length= hsp.align_length
            identities= hsp.identities
            positives= hsp.positives
            gaps=hsp.gaps
            query_seq=hsp.query
            hit_seq=hsp.sbjct
            match_seq=hsp.match
            f.writelines(f'\t\tNumero do hsp:{N+1}\n\t\t\tE-value: {E_value}\n\t\t\tScore: {Score}\n\t\t\tTamanho do hsp: {Length}\n')
            f.writelines(f'\t\t\tNumero de aa identicos: {identities}\n\t\t\tNumero de aa identicos e positivos: {positives}\n\t\t\tNumero de gaps: {gaps}\n')
            if I==0:
                f.write(f'\t\t\tSequencia da proteina:\n')
                for I in range(0,len(query_seq),100):
                    seq_=query_seq[I:I+100]
                    f.writelines(f"\t\t\t\t{seq_}\n")
                f.write('\n')
            else:
                f.write(f'\t\t\tSequencia do hit:\n')
                for I in range(0,len(hit_seq),100):
                    seq_=hit_seq[I:I+100]
                    f.writelines(f"\t\t\t\t{seq_}\n")
                f.write(f'\t\t\tSequencia do alinhamento:\n')
                for I in range(0,len(match_seq),100):
                    seq_=match_seq[I:I+100]
                    f.writelines(f"\t\t\t\t{seq_}\n")
                f.write('\n')
    f.close()

genes=["ABCB11","COG7","EMCN_2","ITIH5L"] 

for gene in genes:
    if gene=="EMCN_2":
        gene_='EMCN'
        record=SeqIO.read(f"genes_information/{gene_}_protein.gb",'genbank')
        seq=record.seq
        name=record.name
    else:
        record=SeqIO.read(f"genes_information/{gene}_protein.gb",'genbank')
        seq=record.seq
        name=record.name

    file = open(f'{gene}_blast.xml','r')
    record_blast = NCBIXML.read(file)
    file.close()
    homo = open(f'Blast_results/{gene}_homology.fa', 'w')
    homo.write(f'>{name}\n{seq}\n')
    alinhamentos=record_blast.alignments
    for I,alinhamento in enumerate(alinhamentos):
       if I!=0:
            for hsp in alinhamento.hsps:
                print(hsp)
                homo.writelines('>' + alinhamento.hit_id + '\n')
                homo.writelines(hsp.sbjct + '\n')
    
    homo.close()
