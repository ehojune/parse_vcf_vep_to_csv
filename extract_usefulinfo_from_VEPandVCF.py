#------------------------------------------------------------------------------
# Title     : extract_usefulinfo_from_VEPoutput.py
# Usage     : python3 extract_usefulinfo_from_VEPoutput.py
# Function     : run_provean_massively: the main function to run provean massively
#		 make_qsub_task_file: the function to make *.sh file with provean's command
#		 run_task_file: the function to run *.sh file.
# Example     : 
# Options     : 
# Author     : Hojune Lee ( ehojune@unist.ac.kr )
# Category     :
# Version     : 2021.12.10
#------------------------------------------------------------------------------

## assign_denovo is hard - coded : dependes on the location where sample names are located in jointcalled vcf.


### variables ###

variant_list = []



### parameters ###
Sample_name = "KU10K-05709"
InterVartxt_path = '/BiO/Research/Project1/KOGIC-KU10K_Rare-Genome-2020-11/Analysis/TrioAnalysis/Results/snv/InterVar/KU10K-05709.chr17_2021Dec22.vcf.hg38_multianno.txt'
output_tsv_path = '/BiO/Access/ehojune/extracted_usefulinfos.chr17.NF1patient.fromVEP.tsv5'
vcf_path = '/BiO/Research/Project1/KOGIC-KU10K_Rare-Genome-2020-11/Results/SNVandIndel.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.AltAware.by.GATK.HaplotypeCaller/KU10K-05709/KU10K-05709.chr17.variant.vcf'
jointcall_vcf_path = '/BiO/Research/Project1/KOGIC-KU10K_Rare-Genome-2020-11/Results/JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller/Korea10K-Rare-1/chr17.recal.vcf'


### modules ###

### classes ###

class Variant:
    CHROM = 'na'
    POS = 'na'
    REF = 'na'
    ALT = 'na'
    GT = 'na'
    AD = 'na'
    GQ = 'na'
    DP = 'na'
    AB = 'na'
    VQSR = 'na'
    AF_1K = 'na'
    gene_name = 'na'
    rsID = 'na'
    variant_type = 'na'
    variant_classification = 'na'
    denovo = "False"
    other_sibiling_overlap = "False"
    clinival_significance = 'na'
    dbname = 'na'
    sequence_nearby = 'na'
    repeat = 'na'
    temp_meta = "False"

    def assign_variant_type(self):
        if self.REF != "." and self.ALT != "." and len(self.REF) == len(self.ALT) == 1:   self.variant_type = "SNV"
        elif self.REF == "." or len(self.REF) < len(self.ALT):    self.variant_type = "DEL"
        else: self.variant_type = "INS"

    def assign_AlleleBalance(self):
        # Homo
        if self.GT[0] ==  self.GT[-1]:    pass
        # na or Unread/Unmapped
        elif self.DP == "na" or self.DP == "0" or self.AD == "na":   pass
        # Hetero
        else:
            AB = int(self.AD.split(',')[0])/int(self.DP)
            self.AB = str(min(AB, 1-AB))[:5]

    def assign_denovo(self, line, formatline):
        GT_idx = formatline.index('GT')
        father_GT = line[9].split(":")[GT_idx][0], line[9].split(":")[GT_idx][2]
        mother_GT = line[10].split(":")[GT_idx][0],  line[10].split(":")[GT_idx][2]
        daughter_GT = line[11].split(":")[GT_idx][0], line[11].split(":")[GT_idx][2]
        son_GT = line[12].split(":")[GT_idx][0], line[12].split(":")[GT_idx][2]
        if son_GT[0] not in father_GT and son_GT[0] not in mother_GT and son_GT[0] not in daughter_GT:
            self.denovo = "TRUE"
        if son_GT[0] not in father_GT and son_GT[0] not in mother_GT and son_GT[0] in daughter_GT:
            self.denovo = "TRUE"
            self.other_sibiling_overlap = "TRUE"
        if son_GT[1] not in father_GT and son_GT[1] not in mother_GT and son_GT[1] not in daughter_GT:
            self.denovo = "TRUE"
        if son_GT[1] not in father_GT and son_GT[1] not in mother_GT and son_GT[1] in daughter_GT:
            self.denovo = "TRUE"
            self.other_sibiling_overlap = "TRUE"

    def replace_variable(self, tempVariant):
        self = tempVariant




### functions ###


def extract_info_from_InterVarOutput(variant_list, InterVartxt_path):

    with open(InterVartxt_path, 'r') as fr1:
        
        line = fr1.readline()        
        while True:
            line = fr1.readline()[:-1].split('\t')
            if len(line) < 2:   break
            #if line[4] == '-':  continue
            tempVariant = Variant()
            tempVariant.CHROM = line[0]
            tempVariant.POS = line[1]
            tempVariant.REF = line[3]
            tempVariant.ALT = line[4]
            tempVariant.gene_name = line[6]
            tempVariant.rsID = line[12]
            tempVariant.variant_classification = line[-12] if line[-12] != "." else "na"
            tempVariant.clinival_significance = line[-24] if line[-24] != "." else "na"
            tempVariant.assign_variant_type()
            variant_list.append(tempVariant)



def update_info_with_vcf(variant_list, vcf_path, Sample_name):
    with open(vcf_path, 'r') as fr2:
        
        line = fr2.readline()
        while line:
            if not line.startswith('#CHROM'):
                line = fr2.readline()
            else:   break

        assert Sample_name in line, f"{line} doesn't have {Sample_name}."
        sample_idx = line.split().index(Sample_name)

        prev_tempVariant_POS = ''
        idx = 0

        while True:

            tempVariant = variant_list[idx]
            if tempVariant.POS != prev_tempVariant_POS:
                
                line = fr2.readline()[:-1].split('\t')
                POS = line[1]
                print("BBB2", POS, tempVariant.POS)
                #print(POS, tempVariant.POS)
            if len(line) < 2:   break
            print("CCC", POS, tempVariant.POS)
            if int(tempVariant.POS) - int(POS) == 1:
                tempVariant.POS = POS
                tempVariant.REF = line[3]
                tempVariant.ALT = line[4]
                #print(POS, tempVariant.POS)
                print("DDD", POS, tempVariant.POS)

            print("ZZZ", POS, tempVariant.POS) 
            while POS != tempVariant.POS:
                line = fr2.readline()[:-1].split('\t')
                POS = line[1]

                #print("AAA", POS, tempVariant.POS)

            print("BBB", POS, tempVariant.POS)              
            #if  POS != tempVariant.POS: print(POS, tempVariant.POS)
            assert POS == tempVariant.POS

            formatline = line[8].split(':')
            infoline = line[sample_idx].split(':')
            tempVariant.VQSR = line[6]
            if 'GT' in formatline:
                tempVariant.GT = infoline[formatline.index('GT')].replace('|', '/')
            if 'GQ' in formatline:
               tempVariant.GQ = infoline[formatline.index('GQ')]
            if 'DP' in formatline:
                tempVariant.DP = infoline[formatline.index('DP')]
            if 'AD' in  formatline:
                tempVariant.AD = infoline[formatline.index('AD')]

            variant_list[idx].replace_variable(tempVariant)

            idx += 1
            prev_tempVariant_POS = POS


def update_info_with_jointcalledvcf(variant_list, jointcall_vcf_path, Sample_name):
    with open(jointcall_vcf_path, 'r') as fr2:
        
        line = fr2.readline()
        while line:
            if not line.startswith('#CHROM'):
                line = fr2.readline()
            else:   break
            
        assert Sample_name in line, f"{line} doesn't have {Sample_name}."
        sample_idx = line.split().index(Sample_name)

        prev_tempVariant_POS = ''
        idx = 0
        while True:

            tempVariant = variant_list[idx]

            if tempVariant.POS != prev_tempVariant_POS:
                line = fr2.readline()[:-1].split('\t')
                POS = line[1]

            if len(line) < 2:   break


            while (POS != tempVariant.POS):
                line = fr2.readline()[:-1].split('\t')
                if len(line) < 2 and POS != tempVariant.POS: return 0
                POS = line[1]
                print(POS, tempVariant.POS)

            assert POS == tempVariant.POS


            formatline = line[8].split(':')
            infoline = line[sample_idx].split(':')

            tempVariant.assign_AlleleBalance()
            tempVariant.assign_denovo(line, formatline)
            variant_list[idx].replace_variable(tempVariant)

            idx += 1
            prev_tempVariant_POS = POS



def write_variant_to_tsv(variant_list, output_tsv_path):
    with open(output_tsv_path, 'w') as fw:
        for variant in variant_list:
            if variant.GT == './.' or variant.GT == '0/0': continue
            #if variant.gene_name != 'NF1': continue
            fw.write(variant.CHROM+'\t'+variant.POS+'\t'+variant.REF+'\t'+variant.ALT+'\t'+variant.GT+'\t'+variant.GQ+'\t'+variant.DP+'\t'+variant.AD+'\t'+variant.AB+'\t'+variant.VQSR+'\t'+'*'+'\t'+variant.gene_name+'\t'+variant.rsID+'\t'+variant.variant_type+'\t'+variant.variant_classification+'\t'+variant.denovo+'\t'+variant.other_sibiling_overlap+'\t'+variant.clinival_significance+'\t'+'*'+'\n')


def main():
    extract_info_from_InterVarOutput(variant_list, InterVartxt_path)
    update_info_with_vcf(variant_list, vcf_path, Sample_name)
    update_info_with_jointcalledvcf(variant_list, jointcall_vcf_path, Sample_name)
    write_variant_to_tsv(variant_list, output_tsv_path)


if __name__ == '__main__':
    main()
