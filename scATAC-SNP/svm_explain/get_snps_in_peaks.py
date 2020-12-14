import sys
import subprocess
import os
import pandas as pd
import pybedtools
import tabix
import pysam
import argparse

#snp_file = '/home/ye/Work/BioAligment/SNP/Shi/scripts/snp_prepare_explain/expaned_master_SNP.txt'
#snp_file="/home/ye/Work/BioAligment/SNP/Shi/duplicated_expand_master_SNP.txt"
dbsnp_file = '/home/ye/Data/SNP_List/00-All.vcf.gz'
dbsnp = tabix.open(dbsnp_file)

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main(args):
    print("Reading SNP File")
    snps = pd.read_csv(args.snp_file, sep='\t')
    snps = prep_snps(snps)
    snps = get_alleles(snps)
    #snps_in_peaks(snps,"/home/ye/Work/BioAligment/SNP/Shi/reproduciblePeaks/vkh/Baso.bed")
    snps_in_peaks(snps,args.overlap_peaks,args.outdir)


def prep_snps(snps):
    print("Preparing SNP File")
    #snps.rename(columns = {'hg38_Position': 'end', 'Effect_Allele': 'effect',
    #                       'Noneffect_Allele': 'noneffect', 'hg38_Chromosome': 'chr',
    #                       'SNP_rsID': 'rsid', 'Direction_of_Association': 'direction',
    #                       'GWAS_pvalue': 'pvalue', 'SNP_Type': 'source_gwas',
    #                       'Locus_Number': 'locus_num', 'LD_Tag_chr': 'ld_tag_chr',
    #                       'LD_Tag_pos_hg38': 'ld_tag_pos', 'r2_with_LD_Tag': 'r2_with_ld_tag'},
    #                       inplace=True)

    #snps.rename(columns={"LD_snp_chrom":"ld_tag_chr","LD_snp_0ind_pos":"ld_tag_pos0",
    #    "LD_snp_1ind_pos":"ld_tag_pos","LD_rs":"LD_rs","LD_val":"r2_with_ld_tag",
    #    "chrom_hg19":"chr","snp_pos_hg19_0":"start","snp_pos_hg19":"end","snp_id":"rsid",
    #    "indexSNP_riskallele":"riskallele","Disease":"disease"},inplace=True)
    #snps.rename(columns={"chrom_hg19":"chr","snp_pos_hg19_0":"start",
    #    "snp_pos_hg19":"end","snp_id":"rsid","indexSNP_riskallele":"riskallele","Disease":"disease"},inplace=True)
    #snps['chr'] = snps['chr'].apply(lambda x: 'chr' + str(x).strip('chr'))
    snps['chr'] = snps['chr'].astype(str)
    #snps['start'] = snps['end'] - 1
    snps['start'] = snps['start'].astype(int)
    snps['end'] = snps['end'].astype(int)
    #snps['effect'] = snps['effect'].apply(lambda x: str(x).upper())
    #snps['noneffect'] = snps['noneffect'].apply(lambda x: str(x).upper())
    #snps = snps[['chr', 'start', 'end', 'rsid', 'effect', 'noneffect', 'direction', 'pvalue', 'source_gwas', 'locus_num', 'ld_tag_chr', 'ld_tag_pos', 'r2_with_ld_tag']]
    #snps.sort_values(by=['chr', 'start', 'end', 'rsid', 'effect', 'noneffect'], inplace=True)
    snps.sort_values(by=['chr', 'start', 'end', 'rsid'], inplace=True)
    return snps


def get_alleles(snps):
    print("Getting dbSNP Alleles")
    ref_list = []
    alt_list = []
    major_list = []
    minor_list = []
    counter = 0
    for index,row in snps.iterrows():
        counter += 1
        if counter % 100 == 0:
            print(counter)
        chrom = row['chr']
        start = row['start']
        end = row['end']
        rsid = row['rsid']
        #effect = row['effect']
        #noneffect = row['noneffect']
        matches = dbsnp.query(chrom.strip('chr'), start, end)
        #matches = dbsnp.query(chrom, start, end)
        gotmatches = False
        for match in matches:
            if match[2] == rsid:
                ref = match[3]
                if ',' in match[4]:
                    alleles = [match[3]] + match[4].split(',')
                    alt = match[4]
                else:
                    alleles = [match[3], match[4]]
                    alt = match[4]
                gotmatches = True
                break
        if gotmatches:
            if 'TOPMED' in match[7] and 'CAF' in match[7]:
                topmed_freqs = match[7].split('TOPMED=')[1].split(',')
                topmed_max_freq = max(topmed_freqs)
                topmed_max_ind = topmed_freqs.index(topmed_max_freq)
                caf_freqs = match[7].split('CAF=')[1].split(';')[0].split(',')
                caf_max_freq = max(caf_freqs)
                caf_max_ind = caf_freqs.index(caf_max_freq)
                major = alleles[topmed_max_ind]
                minor = alleles
                minor.remove(major)
            elif 'TOPMED' in match[7]:
                topmed_freqs = match[7].split('TOPMED=')[1].split(',')
                topmed_max_freq = max(topmed_freqs)
                topmed_max_ind = topmed_freqs.index(topmed_max_freq)
                major = alleles[topmed_max_ind]
                minor = alleles
                minor.remove(major)
            elif 'CAF' in match[7]:
                caf_freqs = match[7].split('CAF=')[1].split(';')[0].split(',')
                caf_max_freq = max(caf_freqs)
                caf_max_ind = caf_freqs.index(caf_max_freq)
                major = alleles[caf_max_ind]
                minor = alleles
                minor.remove(major)
            else:
                major = '.'
                minor = ['.']
            minor = ','.join(minor)
        else:
            ref = '.'
            alt = '.'
            major = '.'
            minor = '.'
        ref_list.append(ref)
        alt_list.append(alt)
        major_list.append(major)
        minor_list.append(minor)
    snps['ref'] = ref_list
    snps['alt'] = alt_list
    snps['major'] = major_list
    snps['minor'] = minor_list
    #snps = snps[['chr', 'start', 'end', 'rsid', 'effect', 'noneffect', 'ref', 'alt', 'major', 'minor', 'direction', 'pvalue', 'source_gwas', 'locus_num', 'ld_tag_chr', 'ld_tag_pos', 'r2_with_ld_tag']]
    #snps.sort_values(by=['chr', 'start', 'end', 'rsid', 'effect', 'noneffect'], inplace=True)
    #snps.drop_duplicates(subset=['chr', 'start', 'end', 'rsid', 'effect', 'noneffect'], inplace=True)
    snps.sort_values(by=['chr', 'start', 'end', 'rsid'],inplace=True)
    snps.drop_duplicates(subset=['chr', 'start', 'end', 'rsid',"disease"],inplace=True)
    return snps


def snps_in_peaks(snps,overlap_peaks,output):
    print("Intersecting SNPs with Peaks")
    snps_bed = pybedtools.BedTool.from_dataframe(snps)
    
    overlap_peaks=coord_hg38Tohg19(overlap_peaks)
    overlap_peak_bed = pybedtools.BedTool(overlap_peaks)
    overlap_intersect_bed = snps_bed.intersect(overlap_peak_bed, u=True, wa=True)
    
    if len(overlap_intersect_bed) > 0:
            overlap_intersect_df = pybedtools.BedTool.to_dataframe(overlap_intersect_bed, header=None)
    else:
            overlap_intersect_df = snps[0:0].copy()
    overlap_intersect_df.columns = ['chr', 'start', 'end', 'rsid',"disease",'ref', 'alt', 'major', 'minor']
    #overlap_intersect_df.sort_values(by=['chr', 'start', 'end', 'rsid', 'effect', 'noneffect'], inplace=True)
    overlap_intersect_df.sort_values(by=['chr', 'start', 'end', 'rsid'],inplace=True)
    overlap_intersect_df.drop_duplicates(subset=['chr', 'start', 'end', 'rsid'], inplace=True)

    makedir(output)
    filename=os.path.join(output,"overlap.expanded.snps.hg19.bed")
    overlap_intersect_df.to_csv(filename, sep='\t', index=False)

def coord_hg38Tohg19(peak_bed):
    liftover="/home/ye/anaconda3/envs/scatac/bin/liftOver"
    chain="/home/ye/Work/BioAligment/SNP/Shi/chainTO/hg38ToHg19.over.chain.gz"
    output=os.path.join(os.path.dirname(peak_bed),"hg19_"+os.path.basename(peak_bed))
    unmap=os.path.join(os.path.dirname(peak_bed),"unmap_"+os.path.basename(peak_bed))
    cmd="{} {} {} {} {}".format(liftover,peak_bed,chain,output,unmap)
    submitter(cmd)
    return output

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A program to get Allele between ld-expand bed file and overlap peaks file')
    parser.add_argument('--snp_file',type=str,default="LD_expand/duplicated_expand_master_SNP.txt",help='ld expand snp file')
    parser.add_argument('--overlap_peaks',type=str, help='scATAC cluster peak file to be overlapped')
    parser.add_argument('--outdir',type=str,default="./LD_exapnd",help='path to save result')
    args = parser.parse_args()
    main(args)
