#------------------------------------------------------------------------------
# Title     : match_af.py
# Usage     : provean program gets an input of protein variant file (HGVC) and reference protein fasta file, to estimate a predicted damage.
# Function     : run_provean_massively: the main function to run provean massively
#		 make_qsub_task_file: the function to make *.sh file with provean's command
#		 run_task_file: the function to run *.sh file.
# Example     : 
# Options     : 
# Author     : Hojune Lee ( ehojune@unist.ac.kr )
# Category     :
# Version     : 2021.11.26
#------------------------------------------------------------------------------


### parameters ###
patient_vcf_path = '/BiO/Access/ehojune/test.chr22.patient.vcf'
trio_vcf_path = '/BiO/Access/ehojune/test.chr22.trio.vcf'
output_vcf_path = '/BiO/Access/ehojune/test.chr22.patient.annotatedwithAFvcf'


### modules ###


### functions ###
def annotate_AF(patient_vcf_path, trio_vcf_path, output_vcf_path):

    patient_VCF_dict = {}                                                        
    with open(patient_vcf_path, 'r') as fr1:
        
        line = fr1.readline()
        while line:
            if line.startswith('#CHROM'):
                break
            else:   line = fr1.readline()
        
        while line:
            if len(line) > 10:
                line = fr1.readline().split()
                POS = int(line[1])
                ALT_genotype = line[4].split()
                AF = [0] * len(ALT_genotype)
                patient_VCF_dict[POS] = ALT_genotype, AF
                print("hi")
            else: break
            
    print("1/3")

    with open(trio_vcf_path, 'r') as fr2:

        line = fr2.readline()
        while line:
            if line.startswith('#CHROM'):
                break
            else:   line = fr2.readline()
        
        patient_POS_list = patient_VCF_dict.keys()

        curr_idx = 0
        while line:
            line = fr2.readline().split()
            trio_vcf_POS = int(line[1])
            patient_vcf_POS = patient_POS_list[curr_idx]
                  
            # Novel Loci in 1K vcf
            if patient_vcf_POS < trio_vcf_POS:
                
                # patient_vcf_POS 가 trio_vcf_POS에 없는 경우, 처음에 patient_vcf 정보를 dict로 불러올 때 AF값을 0으로 초기화했으므로 그냥 넘어가거나,
                # 아래의 코드를 통해 AF를 적는 list를 ['novel']이라는 list로 바꾸는 두 가지 경우의 수를 고려함. (추후 1K dataset에서의 novel한 loci를 검색하기가 용이할 것)
                # patient_POS_list[curr_idx][1] = ['novel']

                curr_idx += 1
                continue
            

            # existing Loci in 1K
            elif patient_vcf_POS == trio_vcf_POS:
                trio_vcf_line_ALT = line[4].split(',')  # [A, G]
                trio_vcf_line_AF =  [ i[2:].split() for i in line[7].split(';') if i.startswith("AF") ] #[1.00, 0.70]
                
                for idx, patient_ALT in enumerate(patient_VCF_dict[patient_vcf_POS][0]):

                    # if patient_ALT is in ALT from VCF, write AF info to patient VCF dict file.
                    if patient_ALT in trio_vcf_line_ALT:
                        patient_POS_list[curr_idx][1][idx] = trio_vcf_line_AF[trio_vcf_line_ALT.index(patient_ALT)]
                    
                   # if patient_ALT is not in ALT from VCF, remain patient_AF value 0.
                    else:
                        # patient_POS_list[curr_idx][1][idx] = 'novel'
                        curr_idx += 1
        
            # read next line in 1K dataset
            else:
                continue

    print("2/3")

    with open(trio_vcf_path, 'r') as fr, open(output_vcf_path, 'w') as fw:

        line = fr.readline()
        while line:
            fw.write(line)
            if line.startswith('#CHROM'):
                break
            else:
                line = fr.readline()

        while line:
            # 수정하면서 적기
            line = fr.readline().split()
            trio_vcf_POS = int(line[2])

            line[7] += f";1K_AF={','.join(patient_VCF_dict[trio_vcf_POS][1])}"
            fw.write("\t".join(line) + "\n")


def main():
	annotate_AF(patient_vcf_path, trio_vcf_path, output_vcf_path)
	
if __name__ == "__main__":
	main()
