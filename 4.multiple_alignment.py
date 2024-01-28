from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

genes = ["ABCB11", "COG7", "EMCN", "ITIH5L"]
blast_files = [f"{gene}_blast.xml" for gene in genes]

for gene, blast_file in zip(genes, blast_files):
    # Load BLAST results from XML file
    blast_results = NCBIXML.parse(open(blast_file))

    sequences = []
    id_counts = {}

    for result in blast_results:
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                # Modify the sequence identifier
                sequence_id = alignment.hit_id
                if sequence_id in id_counts:
                    id_counts[sequence_id] += 1
                    sequence_id = f"{sequence_id}_hsps{id_counts[sequence_id]}"
                else:
                    id_counts[sequence_id] = 1

                # Include the modified sequence name in the aligned file
                sequence = f">{sequence_id}\n{hsp.sbjct}\n"
                sequences.append(sequence)

    # Write sequences to a temporary file
    fasta_file_path = f"{gene}_aligned.fasta"
    with open(fasta_file_path, "w") as fasta_file:
        fasta_file.writelines(sequences)

    # Perform multiple sequence alignment using Clustal Omega
    clustalomega_cline = ClustalOmegaCommandline("C:/Users/joaop/Desktop/clustal-omega-1.2.2-win64/clustalo.exe", infile=fasta_file_path, outfile=f"{gene}_aligned.fasta", force=True)
    clustalomega_cline()

    print(f"Alignment for {gene} saved to {gene}_aligned.fasta")

    # Read the aligned sequences
    aligned_file_path = AlignIO.read(f"{gene}_aligned.fasta", "fasta")

    # Calculate distances for the current gene
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(aligned_file_path)

    # Build the phylogenetic tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Save the phylogenetic tree
    tree_file_path = f"{gene}_phylogenetic_tree.newick"
    Phylo.write(tree, tree_file_path, "newick")

    print(f"Phylogenetic tree for {gene} saved to {tree_file_path}")