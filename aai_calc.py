
import os, subprocess, sys, time
from itertools import combinations
from pathlib import Path
from collections import defaultdict



"""Chunk below compares all combinations through best reciprocoal 'BLAST'-hits using
DIAMOND BLASTP. Also, requires an identity > 50% and alignment length >= 60% of
BOTH sequences."""
# TO DO: Put into class ... RBBH with OrthoFinder2 calls on orthology using
# high inflation parameter. Permits tree-walking based comparisons, may be
# attractive to try...

# Perform the best reciprocal BLAST hits for each pair of fasta-files.
def AAI_RBBHs(aa_folder):
    Path('RBBH_Files/').mkdir(exist_ok=True)

    aai_rbbh = {}
    aa_fastas = sorted(i for i in os.listdir(aa_folder) if '.CDS.AA.' in i)

    # Get all unique combinations of the FASTA files [Taxon-Name.CDS.AA.fas].
    taxon_comps = [tuple(map(str, comb)) for comb in combinations(aa_fastas,2)]

    for comp in taxon_comps:
        comp_tsv_1 = run_diamond(f'{aa_folder}{comp[0]}',f'{aa_folder}{comp[1]}')
        while not os.path.exists(comp_tsv_1):
            time.sleep(0.1)
        comp_tsv_2 = run_diamond(f'{aa_folder}{comp[1]}',f'{aa_folder}{comp[0]}')
        while not os.path.exists(comp_tsv_2):
            time.sleep(0.1)

        # Tries to calc AAI before DIAMOND finishes writing file.
        # Put parsing the RBBH files into new function to force slowdown. - XXMA
        aai_rbbh[comp] = get_rbbh(comp_tsv_1, comp_tsv_2)

    # Update the output TSV to possess a unique "run" name. -XXMA
    # TSV file that summarizes the AAI between each pair of taxa.
    with open('AAI.RBBH.Summary.tsv','w+') as w:
        w.write('Taxon-1\tTaxon-2\tNum-Comparisons\tAAI\t100-AAI\n')

        for k, v in aai_rbbh.items():

            overall_aai = sum(v)/len(v)

            taxon_names = '\t'.join([i.split('/')[-1].split('.CDS')[0] for i in k])

            w.write(f'{taxon_names}\t{len(v)}\t{overall_aai:.3f}\t{100-overall_aai:.3f}\n')

    # TSV file with EACH comparison
    with open('AAI.RBBH.AllComparisons.tsv','w+') as w:
        w.write('Taxon-1\tTaxon-2\tNum-Comparisons\tAA-Identity\t100-AA-Identity\n')

        for k, v in aai_rbbh.items():

            taxon_names = '\t'.join([i.split('/')[-1].split('.CDS')[0] for i in k])
            for aa_identity in v:
                w.write(f'{taxon_names}\t{len(v)}\t{aa_identity}\t{100-aa_identity:.3f}\n')

# Performs the "BLASTP" comparative search, with mostly default parameters.
# min_cov is used for both the QUERY and SUBJECT coverage in the alignment.
# No minimum identity is required here, but a weak e-value is used (1e-3).
# Returns the output of DIAMOND for easy access.
def run_diamond(fas1, fas2, min_cov=60):

    f1 = fas1.split('/')[-1].split('.CDS')[0]
    f2 = fas2.split('/')[-1].split('.CDS')[0]

    dmnd_cmd = f'diamond blastp -e 1e-3 ' \
        f'-q {fas1} -d {fas2} ' \
        f'--query-cover {min_cov} ' \
        f'--subject-cover {min_cov} ' \
        f'-o RBBH_Files/{f1}_{f2}.DIAMOND_Comp.tsv ' \
        f'-f 6 qseqid sseqid qlen slen length pident evalue bitscore'

    dmnd_call = subprocess.Popen(dmnd_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        shell=True)

    return f'RBBH_Files/{f1}_{f2}.DIAMOND_Comp.tsv'

# Parse the reciprocal blast hits to grab the "best". Returns list with
# the %AA identity based on the RBBHs for a given taxonomic comparison.
def get_rbbh(tsv_1, tsv_2):
    rbbhs = []
    t1_hits, t2_hits = defaultdict(list), defaultdict(list)

    for line in open(tsv_1).readlines():
        t1_hits[line.split('\t')[0]].append(line)

    for line in open(tsv_2).readlines():
        t2_hits[line.split('\t')[0]].append(line)

    for k, v in t1_hits.items():
        h = v[0].split('\t')[1]

        if h in t2_hits.keys():
            rh = t2_hits[h][0].split('\t')[1]

            # Update to connect the sequences compared with their AA identity. - XXMA
            if rh == k:
                rbbhs.append(float(v[0].split('\t')[-3]))

    return rbbhs

if __name__ == '__main__':
    print('In-Progress.... Will not run -XXMA 2022-07-26')
    sys.exit()
    if len(sys.argv) != 2:
        print('Usage:\n   python3 comp_genera.py [Folder-with-AA-FASTAs]')
        sys.exit(1)
    AA_folder = sys.argv[1]
    AAI_RBBHs(AA_folder)

"""Below is a work in progress to use OrthoFinder2 to cluster putative orthologous
gene families (OGs), which will then be used as the basis for calculating %AAI.
Note that this is NOT finished and may be deleted."""

# def run_of2(aa_folder, of2_results_folder, inflation=10, threads=4):
#     of2_outf = f'{of2_results_folder}_I{inflation}'
#     of2_cmd = f'orthofinder.py -I {inflation} -f {aa_folder} -o {of2_outf}'
#     of2_call = subprocess.Popen(of2_cmd,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.DEVNULL,
#         shell=True)
#     return of2_outf
#
# def mafft_align(msa_file):
#     eins_align_fas = f'{msa_file.split(".fa")[0]}.EINSI.fas'
#     einsi_cmd = f'einsi {msa_file} > {eins_align_fas}'
#     einsi_call = subprocess.Popen(einsi_cmd,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.DEVNULL,
#         shell=True)
#
# def calc_aai(msa_file):
#     aseq, bseq = [i for i in SeqIO.parse(msa_file,'fasta')]
#     num_shared_aa = sum(i==j for i,j in zip(f'{aseq.seq}',f'{bseq.seq}'))
#     num_comp_aa = min((len(aseq) - aseq.seq.count('-')), (len(bseq) - bseq.seq.count('-')))
#     comp_aai = f'{100*num_shared_aa/num_comp_aa:.3f}'
#     comp_names = tuple(sorted([aseq.id.split('_XX')[0], bseq.id.split('_XX')[0]]))
#     return comp_aai, comp_names
#
# def pause_for_file(file):
#     while not os.path.exists(file):
#         time.sleep(0.1)
#     return 0
#
# def comp_taxa(sco_folder, taxon_list):
# taxa = sorted(i.rstrip() for i in open(taxon_list).readlines())
# taxon_comps = {tuple(map(str, comb)):[] for comb in combinations(taxa,2)}
# for f in os.listdir(sco_folder):
#     seqs = [f'>{i.id}\n{i.seq}\n' for i in SeqIO.parse(f'{sco_folder}{f}','fasta')]
#     comps = [tuple(map(str, comb)) for comb in combinations(seqs,2)]
#     for comp_seqs in comps:
#         print(comp_seqs)
#         temp_name = '_oo_'.join([i.split('\n')[0].strip('>') for i in comp_seqs])
#         temp_file_path = f'Temp/{temp_name}'
#         with open(f'{temp_file_path}.fas','w+') as w:
#             w.write(''.join(comp_seqs))
#         pause_for_file(f'{temp_file_path}.EINSI.fas')
#         mafft_align(f'{temp_file_path}.fas')
#         pause_for_file(f'{temp_file_path}.EINSI.fas')
#         aai, comp_pair = calc_aai(f'{temp_file_path}.EINSI.fas')
#         taxon_comps[comp_pair].append(aai)
