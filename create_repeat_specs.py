#!/usr/bin/env python3

import argparse
import subprocess
import re
import os
import sys
import json
import operator 
from subprocess import Popen, PIPE

#User inputs region, delimited by ":|-" (eg. 8:100651852-100652541), and motif (eg. AAG)
parser = argparse.ArgumentParser(description='Identify reference STR tracks to generate a repeat-spec file as input for  Expansion Hunter Original.')
parser.add_argument('-r', '--region', type = str, nargs = 1, required = True, help = 'Region to retrieve reference genomic sequence')
parser.add_argument('-m', '--motif', type = str, nargs = 1, required = True, help = 'STR motif of interest')

args = parser.parse_args()
region = args.region[0]
chrom, start, end = re.split(':|-', region)
start, end = int(start), int(end)
EH_motif = args.motif[0]
chr_chrom = "chr"+chrom


#Runs samtools faidx to pull genomic sequence defined by user (eg. 8:100651852-100652541) and remove line breaks
#Output is stored as a Python string
#Pad the search region by 500bp 
output = subprocess.Popen("samtools faidx /hpf/largeprojects/tcagstor/users/zwang/ref/ngs/hs37d5/human_g1k_v37_decoy.fasta {} | tail -n +2 |  tr -d '\n'".format(region),
stdout = subprocess.PIPE, shell = True).communicate()[0]
ref_DNA = str(output, 'utf-8')

#Finds all repeat tracks found in input DNA sequence
#Returns list of motifs and associated track length [('A', 3.0), ('AT', 4.0),...,('CCG', 3.0)]
def create_repeat_specs(DNA):
    STR_list = []
    r = re.compile(r"(.+?)\1+")
    for match in r.finditer(DNA):
        STR_list.append([match.group(1), len(match.group(0))/len(match.group(1))])
    return STR_list

#Finds all iterations of input_motif; ie. given EHDN motif, find all corresponding motifs
#Returns list of similar motifs which consists of all permutations of input and reverse complement motif
#Example: input_motif = AAG will return ['GAA', 'AGA', 'AAG', 'TTC', 'TCT', 'CTT']
def generate_all_motifs(input_motif):
    base_to_number_dictionary = {"A": 1, "C" :2, "T": 3, "G": 4}
    numeric_motif = [base_to_number_dictionary[letter] for letter in input_motif]
    original_motifs = []
    rev_comp_motifs = []
    offset_vector = [offset for offset in range(1,(len(numeric_motif)+1))]
    #k represents offset value
    k = 0
    while k < len(numeric_motif):
        #Generate offset vector for given offset
        offset_vector = [((x-1 + (len(numeric_motif)-1)) % (len(numeric_motif))+1)
                         for x in offset_vector]
        #Add new motif to list of motifs
        original_motifs.append([numeric_motif[val_off-1] for val_off in offset_vector])
        k += 1
    comp_motifs = [[x + 2 if x<3 else x - 2 for x in m] for m in original_motifs]
    rev_comp_motifs = [motif_to_flip[::-1] for motif_to_flip in comp_motifs]
    original_motifs.extend(rev_comp_motifs)
    trans_motif_back = [[list(base_to_number_dictionary.keys())[list(base_to_number_dictionary.values()).index(x)]
                         for x in motif] for motif in original_motifs]
    motif_strings = [''.join(motif) for motif in trans_motif_back]
    return motif_strings
#Find the number of matches between two strings
def distance(a, b): 
    return sum(x == y for x, y in zip(a, b))

#Given all repeat tracks in user defined region, find those with a similar motif to input motif
#Returns list of repeat tracks and repeat lengths that match the user defined motif 
#Example repeat_motif = AAG, returns [['GAA', 7.0], ['GAA', 3.0]]
def find_repeat_tracks(genomic_DNA, repeat_motif):
    DNA_repeats_all = create_repeat_specs(genomic_DNA)
    motifs, num_repeats = zip(*DNA_repeats_all)
    possible_motifs = generate_all_motifs(repeat_motif)
    matched_repeats = [[motifs[i],num_repeats[i]] for i,x in enumerate(motifs) if x in possible_motifs]
    matched_repeats_with_mismatch = [[motifs[i],num_repeats[i]] for i,x in enumerate(motifs) for pos_motif in possible_motifs if len(x) == len(repeat_motif) and distance(pos_motif, x) >= len(repeat_motif)*(2/3)]
    return DNA_repeats_all, matched_repeats, matched_repeats_with_mismatch

#Show user list of repeat tracks (eg. [['GAA', 7.0], ['GAA', 3.0]])
#that correspond to input motif (eg. AAG) found in user defined region
#Prompt user to pick a repeat track that will be used to create a repeat-specification file
#Options include:
# - [max] (the largest repeat track, eg. ['GAA', 7.0])
#- manual specification (eg. type [GAA,3])
#- [more] - this will return all repeat tracks found in region, not just matches to the input motif
#Following [more], user can pick a repeat track, or type [none] if no appropriate repeat track
def choose_motif(repeats_of_DNA, matched_repeat_full):
    print("Choose repeat motif and track size by inputting [motif],[size]:"+
    "\n"+"- Type [max] for largest track"+"\n"+"- Type [more] if list is empty or for more choices"+
    "\n" + "List of matched EHDN repeat tracks:" + "\n", matched_repeat_full)
    chosen_motif = input()
    if chosen_motif == "max":
        repeat_motif = max(matched_repeat_full, key=operator.itemgetter(1))
    elif chosen_motif == "more":
        print("Input [motif],[size] from all possible repeat tracks of selected region"+"\n"+
    "- Type [none] if no suitable options"+"\n", repeats_of_DNA)
        chosen_motif_more = input()
        if chosen_motif_more == "none":
            repeat_motif=[]
        else:
            repeat_motif = [x for x in chosen_motif_more.split(",")]
    else:
        repeat_motif = [x for x in chosen_motif.split(",")]

    return repeat_motif

#Given the chosen repeat track defined by user, return the DNA sequence associated with track
#If [none] chosen, returns empty string
#If region specified, returns the initial region defined by user with highlighted repeat track "***>>>[repeat]<<<***"
#Asks if user would like to "Proceed with creating repeat-specification file? [y] or [n]"
#Example: User region: 8_100651852_100652541, User motif: AAG, chosen repeat: [max]
#TGTTATATATATATAAGTTATTGCATATATATATATTATTATATATAAAAATATGTTTTTTTTTA
#CCAATCAAATTAAGGAGTTTTTAAAATAATAAGCAGCTTTGGGAGGCTGAGGCAGGAGGATTGCTTGAGGCCA
#GGAGTTTGAGACCAGGTTGGGCTACACTATGAGACCTCATTTCTACAAAAAATGTAAAAATAAGCTAGCCATGG
#TGGTGCATGCCTGTAATCCTAGCTCCTCCGAAGGCTGAAATGGGAGGCTTCCTTGAGCCACCACTGAGGGGTTG
#TAGGCTGCAGTGAGCTATGATCACACCACTGCACTCCAACCTGGGCAACAGAACAAGACTCTCTCTTTAAAAAT
#AGTAATAATAATAATAATAATAATAATAATAATAATAATAATAA***>>>GAAGAAGAAGAAGAAGAAGAA<<<***
#GAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAGAAGAAGCAATA
#CCCAGTTACCCAGTACTTTGGGAGACATTTTTTATCGCTGATAGTAATGTAACTTAGTGCAGCCTTTTTGGATAGAT
#GTTTGACAATATTCAGTAACTTTAAACTTTTTCGTACCATTTGCCTAATAATTCCACTTTTTAAAAATATTTTTGTT
#TTTGTTTTTGTTTTTGAGATAGGGTCTCACTCTG
def chosen_repeat_track(chosen_repeat_motif):
    start_track = 0
    end_track = 0
    repeat_track = chosen_repeat_motif[0]*int(chosen_repeat_motif[1])
    start_track = ref_DNA.index(repeat_track)+start
    end_track = start+ref_DNA.index(repeat_track)+len(repeat_track)-1
    return start_track, end_track

#Create a repeat-specification file for the repeat track selected
def create_repeat_spec(chrom_rt, start_rt, end_rt, motif_rt):
    repeat_id = '{}_{}_{}'.format(chrom_rt, start_rt, end_rt)
    rec = {'RepeatId' : repeat_id,
           'RepeatUnit': motif_rt,
           'TargetRegion' : '{}:{}-{}'.format(chrom_rt, start_rt, end_rt),
           'CommonUnit' : 'true'}
    with open('repeat-specs-test_padding/{}.json'.format(repeat_id), 'w') as spec_file:
        print(json.dumps(rec, sort_keys=True, indent=4), file=spec_file)
    return rec
#print(ref_DNA, EH_motif, "printing all this bullshit")

repeat_tracks = find_repeat_tracks(ref_DNA, EH_motif)
#print(repeat_tracks, "printing repeat_tracks")
#If there is a STR in the region that matches the input motif a repeat-spec file with the matching STR will be generated
if repeat_tracks[1]:
    input_motif = max(repeat_tracks[1], key=operator.itemgetter(1))
#If there is no STR in the region that matches the input motif a repeat-spec file with a mismatched  STR will be generated
elif repeat_tracks[2]:
    input_motif = max(repeat_tracks[2], key=operator.itemgetter(1))
#If there no STR that is a match or mismatch, no repeat-spec file will be generated
else:
    input_motif = []
#If there is a STR with similar input motif, the longest STR track will be identified by incrementing the region
if input_motif:
    increment_STR_length = 1
    while increment_STR_length >= 0: 
        if input_motif[0]*(int(input_motif[1])+increment_STR_length) in ref_DNA: 
            increment_STR_length += 1 
        else: 
            increment_STR_length -= 1
            break
    input_motif_final = [input_motif[0], int(input_motif[1])+increment_STR_length]
    region_create_spec = chosen_repeat_track(input_motif_final)

#Create a directory to hold repeat-spec files
if not os.path.exists('repeat-specs-test_padding/'):
    os.makedirs('repeat-specs-test_padding/')

if input_motif != []:
    create_repeat_spec(chrom, region_create_spec[0], region_create_spec[1], input_motif_final[0])
else:
    curr_file = open('repeat-specs-test_padding/STRs_not_found.txt','a')
    curr_file.write(region+'\t'+EH_motif+'\n')
    curr_file.close()

