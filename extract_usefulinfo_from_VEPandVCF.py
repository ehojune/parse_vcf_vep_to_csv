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
patient_name = "KU10K-05709"
father_name = "KU10K-05706"
mother_name = "KU10K-05707"
sibiling_name = "KU10K-05708"
Samples = [father_name, mother_name, patient_name, sibiling_name]



tempfile_path = 'BiO/Access/ehojune'
InterVartxt_path = '/BiO/Research/Project1/KOGIC-KU10K_Rare-Genome-2020-11/Analysis/TrioAnalysis/Results/snv/InterVar/KU10K-05709.chr17_2021Dec22.vcf.hg38_multianno.txt'
output_tempfile_path = 'BiO/Access/ehojune/temp.jointcalledvcf.chr17.NF1patient.fromVEP.tsv5'
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
        sibiling_GT = line[11].split(":")[GT_idx][0], line[11].split(":")[GT_idx][2]
        son_GT = line[12].split(":")[GT_idx][0], line[12].split(":")[GT_idx][2]
        if son_GT[0] not in father_GT and son_GT[0] not in mother_GT and son_GT[0] not in sibiling_GT:
            self.denovo = "TRUE"
        if son_GT[0] not in father_GT and son_GT[0] not in mother_GT and son_GT[0] in sibiling_GT:
            self.denovo = "TRUE"
            self.other_sibiling_overlap = "TRUE"
        if son_GT[1] not in father_GT and son_GT[1] not in mother_GT and son_GT[1] not in sibiling_GT:
            self.denovo = "TRUE"
        if son_GT[1] not in father_GT and son_GT[1] not in mother_GT and son_GT[1] in sibiling_GT:
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

def extract_jointcallVCF_of_family_from_bulkVCF(jointcall_vcf_path, output_tempfile_path, patient_name, father_name, mother_name, sibiling_name):
    with open(jointcall_vcf_path, 'r') as fr, open(output_tempfile_path, 'w') as fw:
        with open(jointcall_vcf_path, 'r') as fr2:
        
            line = fr.readline()
            while line:
                if not line.startswith('#CHROM'):
                    line = fr.readline()
                else:   break

            vcf_column = line.split()
            for sample in patient_name, father_name, mother_name, sibiling_name:
                assert sample in vcf_column, f"{line} doesn't have {sample}."

            patient_idx = vcf_column.index(patient_name)
            father_idx = vcf_column.index(father_name)
            mother_idx = vcf_column.index(mother_name)
            sibiling_idx = vcf_column.index(sibiling_name)
            fw.write('\t'.join(line.split()[:8]+[patient_name, father_name, mother_name, sibiling_name]))

            while line:
                line = fr.readline()[:-1].split('\t')
                patient_info = line[patient_idx]
                father_info = line[father_idx]
                mother_info = line[mother_idx]
                sibiling_info = line[sibiling_idx]

                patient_GT = patient_info[:2]
                father_GT = father_info[:2]
                mother_GT = mother_info[:2]
                sibiling_GT = sibiling_info[:2]

                #가족 중 한명이 '0','.' 이외의 Genotype을 가질 때
                if len({patient_GT[0], patient_GT[-1], father_GT[0], father_GT[-1], mother_GT[0], mother_GT[-1], sibiling_GT[0], sibiling_GT[-1]}.diff({'0', '.'})):
                    
                    #line = line[:8]
                    #line.append(patient_info)
                    #line.append(father_info)
                    #line.append(mother_info)
                    #line.append(sibiling_info)

                #InterVar 방식에 맞게 index 수정
                    temp_ALT = line[4].split(',')
                    num_of_ALT = len(temp_ALT)
                    for i in range(num_of_ALT):
                        temp_REF = line[3]
                        idx_move = 0
                        while temp_REF[0] == temp_ALT[i][0]:
                            idx_move += 1
                            try:
                                temp_REF = temp_REF[1:]
                            except:
                                temp_REF = '-'
                                break
                            try:
                                temp_ALT[i] = temp_ALT[i][1:]
                            except:
                                temp_ALT = '-'
                                break
                        line_with_modified_idx = line[0] + line
                        fw.write()
                










def assert_all_POS_in_InterVarOutput_are_in_JointCalled_VCF(variant_list, jointcall_vcf_path):

    JointcallVCF_POS_list = []

    with open(vcf_path, 'r') as fr:
        line = fr.readline()
        while line:
            if not line.startswith('#CHROM'):
                line = fr.readline()
            else:   break

        while True:
            line = fr.readline()[:-1].split('\t')
            if len(line) < 2: break
            JointcallVCF_POS_list.append(line[1])

        for variant in variant_list:
            try:    JointcallVCF_POS_list.remove(variant.POS)
            except:
                print(variant.POS)
                JointcallVCF_POS_list.remove(str(int(variant.POS) - 1))
            finally:    print("doesn't exist")


def update_info_with_jointcalledvcf(variant_list, jointcall_vcf_path, patient_name, father_name, mother_name, sibiling_name):
    with open(jointcall_vcf_path, 'r') as fr2:
        
        line = fr2.readline()
        while line:
            if not line.startswith('#CHROM'):
                line = fr2.readline()
            else:   break

        vcf_column = line.split()
        for sample in patient_name, father_name, mother_name, sibiling_name:
            assert sample in vcf_column, f"{line} doesn't have {sample}."

        patient_idx = vcf_column.index(patient_name)
        father_idx = vcf_column.index(father_name)
        mother_idx = vcf_column.index(mother_name)
        sibiling_idx = vcf_column.index(sibiling_name)


        prev_tempVariant_POS = ''
        idx = 0
        POS = '0'
        while True:

            tempVariant = variant_list[idx]

            if tempVariant.ALT == '-':
                tempVariant.POS = [str(int(tempVariant.POS) - 1), str(int(tempVariant.POS) - 2), str(int(tempVariant.POS) - 3)]
            print(type(POS))
            stack = 0
            while POS != tempVariant.POS:
                if type(tempVariant.POS) == list and POS in tempVariant.POS:
                    break
                # if int(tempVariant.POS) - int(POS) == 1:
                #     tempVariant.POS = POS
                #     print('aaa', tempVariant.ALT)
                #     tempVariant.REF = line[3][0] + tempVariant.REF
                #     tempVariant.ALT = line[3][0] + tempVariant.ALT
                #     print('AAA', tempVariant.ALT)
                line = fr2.readline()[:-1].split('\t')

                if len(line) < 2:
                    print("EOF")
                    return 0

                POS = line[1]

                stack += 1
                print(tempVariant.POS, POS)
                if stack == 50:    return 0

            print("out of the loop")

            tempVariant.POS = POS

            assert len(line) > 2, "end of file"

            if '-' in (tempVariant.REF, tempVariant.ALT):
                REF_in_vcf = line[3][0]
                tempVariant.REF = REF_in_vcf
                tempVariant.ALT = REF_in_vcf + tempVariant.ALT.replace('-', '')


            tempVariant.VQSR = line[6]
            formatline = line[8].split(':')

            patient_line = line[patient_idx].split(':')
 
            if 'GT' in formatline:
               # GT_array = [line[]]
                GT_idx = formatline.index('GT')
                tempVariant.GT =  patient_line[GT_idx][0]
            if 'GQ' in formatline:
               tempVariant.GQ =  patient_line[formatline.index('GQ')]
            if 'DP' in formatline:
                tempVariant.DP = patient_line[formatline.index('DP')]
            if 'AD' in  formatline:
                tempVariant.AD = patient_line[formatline.index('AD')]

            variant_list[idx].replace_variable(tempVariant)

            idx += 1
            print('processed:', POS)


'''
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
'''



def write_variant_to_tsv(variant_list, output_tsv_path):
    with open(output_tsv_path, 'w') as fw:
        for variant in variant_list:
            if variant.GT == './.' or variant.GT == '0/0': continue
            #if variant.gene_name != 'NF1': continue
            fw.write(variant.CHROM+'\t'+variant.POS+'\t'+variant.REF+'\t'+variant.ALT+'\t'+variant.GT+'\t'+variant.GQ+'\t'+variant.DP+'\t'+variant.AD+'\t'+variant.AB+'\t'+variant.VQSR+'\t'+'*'+'\t'+variant.gene_name+'\t'+variant.rsID+'\t'+variant.variant_type+'\t'+variant.variant_classification+'\t'+variant.denovo+'\t'+variant.other_sibiling_overlap+'\t'+variant.clinival_significance+'\t'+'*'+'\n')


def main():
    extract_info_from_InterVarOutput(variant_list, InterVartxt_path)
    #assert_all_POS_in_InterVarOutput_are_in_JointCalled_VCF(variant_list, jointcall_vcf_path)

    update_info_with_jointcalledvcf(variant_list, jointcall_vcf_path, patient_name, father_name, mother_name, sibiling_name)
    write_variant_to_tsv(variant_list, output_tsv_path)


if __name__ == '__main__':
    main()
