#!/usr/bin/python

###########################################################################
## This script takes as input:
# 1) the detailed report file from Biomarks SNPTrace for multiple samples (so that both sample swap and contamination can be tested)
# 2) rs\ genotypes
# and output a formatted vcf file as an input for VerifyBAMID, which requires that it
# a) contains only SNPs
# b) contains both genotype and allele frequency
# c) contains autosomal only (no X, Y)
###########################################################################

import sys
import argparse
import re

#parser = argparse.ArgumentParser()
#parser.add_argument()
input_rs = sys.argv[1] # rs\ genotypes
input_biomarks = sys.argv[2] # Biomarks detailed report file


# hu name vs rs no conversion
conversion = {"hu1":"rs1471939","hu2":"rs4666200","hu3":"rs7554936","hu4":"rs9530435","hu5":"rs6104567","hu7":"rs2272998","hu9":"rs560681","hu11":"rs6591147","hu12":"rs321198","hu14":"rs870347","hu15":"rs2946788","hu16":"rs4891825","hu17":"rs10108270","hu18":"rs2397060","hu20":"rs7229946","hu21":"rs13182883","hu22":"rs1876482","hu23":"rs315791","hu24":"rs7205345","hu25":"rs798443","hu26":"rs4717865","hu27":"rs2416791","hu28":"rs2125345","hu29":"rs4746136","hu31":"rs13218440","hu32":"rs1523537","hu33":"rs1058083","hu34":"rs1344870","hu35":"rs7704770","hu36":"rs1410059","hu37":"rs5768007","hu38":"rs260690","hu39":"rs13400937","hu41":"rs4918842","hu42":"rs9809104","hu44":"rs1821380","hu45":"rs279844","hu46":"rs952718","hu47":"rs447818","hu48":"rs13134862","hu49":"rs4463276","hu51":"rs3943253","hu52":"rs6548616","hu53":"rs731257","hu54":"rs9319336","hu55":"rs1019029","hu56":"rs1358856","hu58":"rs1823718","hu59":"rs2503107","hu61":"rs10236187","hu62":"rs1513181","hu63":"rs7657799","hu64":"rs2504853","hu65":"rs772262","hu66":"rs3737576","hu68":"rs445251","hu69":"rs10488710","hu70":"rs722869","hu71":"rs1109037","hu72":"rs3780962","hu73":"rs7997709","hu74":"rs4670767","hu75":"rs9522149","hu76":"rs4908343","hu78":"rs12629908","hu79":"rs1336071","hu80":"rs740598","hu81":"rs12997453","hu82":"rs2352476","hu83":"rs1554472","hu84":"rs10007810","hu85":"rs1760921","hu86":"rs1040045","hu87":"rs10496971","hu88":"rs7803075","hu89":"rs987640","hu90":"rs6444724","hu91":"rs10092491","hu92":"rs735612","hu93":"rs985492","hu95":"rs9951171","hu96":"rs3907047","hu98Y":"rs1865680","hu103X":"rs525869","hu107X":"rs2040962","hu109X":"rs530501","hu111Y":"rs2032624","hu200":"rs1296819","hu201":"rs316598","hu202":"rs722290","hu203":"rs1872575","hu204":"rs18579","hu205":"rs891700","hu207":"rs8113143","hu208":"rs1008730","hu209Y":"rs17307398"}

# process rs\ genotypes file to get the AF information
aDict = {}
with open(input_rs) as rs_in:
    for line in rs_in:
        line = line.strip()
        if re.match(r'^[0-9]',line): # extract all the autosomal SNPs, not the X or Y.
            #print line
            arr = line.split("\t")
            #print arr[-1]
            matchObj = re.search(r"GMAF=(0.[0-9]*);",arr[-1])
            # chrom pos ref alt AF
            aDict[arr[2]] =  [arr[0],arr[1],arr[3],arr[4],matchObj.group(1)]


#print VCF title
header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
content = ""
# declare a snp based dictionary with the element like "'rs88374':'0|0\t1|0\t1|1.....'"
snp_dict = {}
sample_list = []
with open(input_biomarks) as biomarks_in:
    for line in biomarks_in:
        line = line.strip()
        # start with sample lines
        if re.search(r'^S.*Unknown',line):
            arr = line.split(",")
            # Populate sample list
            sample_ID = arr[0][0:4]+arr[4]
            if sample_ID not in sample_list:
                sample_list.append(sample_ID)
            
            # Populate the non-sex oriented genotype information for al samples in Biomarks
            if not re.search(r'.*["X","Y"]',arr[1]):
                ID = conversion[arr[1]]
                #print ID
                if ID not in snp_dict:
                    if arr[7] == "XX":
                        snp_dict[ID] = "0|0\t"
                    elif arr[7] == "XY":
                        snp_dict[ID] = "0|1\t"
                    elif arr[7] == "YY":
                        snp_dict[ID] = "1|1\t"
                else:
                    if arr[7] == "XX":
                        snp_dict[ID] += "0|0\t"
                    elif arr[7] == "XY":
                        snp_dict[ID] += "0|1\t"
                    elif arr[7] == "YY":
                        snp_dict[ID] += "1|1\t"
                    #print snp_dict
                
# Print out the multiple sample VCF
header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(sample_list)
print header
#print snp_dict
for ID in snp_dict:
    print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(aDict[ID][0],aDict[ID][1],ID,arr[2],arr[3],".","PASS","AF="+aDict[ID][4],"GT",snp_dict[ID])



#            content = content+'{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(aDict[ID][0],aDict[ID][1],ID,arr[2],arr[3],".","PASS","AF="+aDict[ID][4],"GT",geno)
           
    
