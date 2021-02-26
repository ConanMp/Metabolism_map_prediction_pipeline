#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 09:55:10 2019

@author: mconan
"""

import os
import argparse
import time
import itertools
from decimal import Decimal
from ast import literal_eval
import json
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor

def read_SMIRKS_label_dic(SMIRKS_label_dic_path):
    nats_Accerlys_names=[]
    cyps_Accerlys_names=[]
    ugts_Accerlys_names=[]
    sults_Accerlys_names=[]
    gsts_Accerlys_names=[]

    SMIRKS_label_dic_file=open(SMIRKS_label_dic_path,'r')
    lines=SMIRKS_label_dic_file.readlines()
    for line in lines:
        line=line.replace('\n','')
        splitted_line=line.split('\t')
        rule_name=splitted_line[0]
        enzyme_family=splitted_line[1]
        if enzyme_family=="CYPs":
            cyps_Accerlys_names.append(rule_name)
        elif enzyme_family=="NATs":
            nats_Accerlys_names.append(rule_name)
        elif enzyme_family=="SULTs":
            sults_Accerlys_names.append(rule_name)
        elif enzyme_family=="UGTs":
            ugts_Accerlys_names.append(rule_name)
        elif enzyme_family=="GSTs":
            gsts_Accerlys_names.append(rule_name)
    return(cyps_Accerlys_names,ugts_Accerlys_names,nats_Accerlys_names,gsts_Accerlys_names,sults_Accerlys_names)

def get_all_enzyme(json_dic, toolname_list, SMIRKS_label_dic_path):
    enzymes=set([])
    cyp_enzymes=set([])
    reactions=json_dic["links"]
    cyps_Accerlys_names, ugts_Accerlys_names, nats_Accerlys_names, gsts_Accerlys_names, sults_Accerlys_names = read_SMIRKS_label_dic(SMIRKS_label_dic_path)
    for reaction in reactions:
        reaction_name=reaction["reaction_name"]
        if reaction_name in nats_Accerlys_names:
            enzymes.add("NATs")
        elif reaction_name in ugts_Accerlys_names:
            enzymes.add("UGTs")
        elif reaction_name in gsts_Accerlys_names:
            enzymes.add("GSTs")
        elif reaction_name in sults_Accerlys_names:
            enzymes.add("SULTs")
        elif reaction_name in cyps_Accerlys_names:
            enzymes_tool=reaction["enzymes"]
            for tool in enzymes_tool.keys():
                if tool in toolname_list:
                    tool_predicted_enzymes=enzymes_tool[tool]
                    for an_enzyme in tool_predicted_enzymes.keys():
                        an_enzyme=an_enzyme.replace("CYP","")
                        enzymes.add(an_enzyme)
                        cyp_enzymes.add(an_enzyme)
    return (list(enzymes), list(cyp_enzymes))

def define_letters_for_nodes(Json_obj):
    alphabet=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    nodes=Json_obj["nodes"]
    allready_used_letters=set([])
    for node in nodes:
        node_id=int(node["id"])
        if node_id>=26:
            while node_id>=26:
                node_id=node_id-26
        letter_corresponding=alphabet[node_id]
        while letter_corresponding in allready_used_letters:
            letter_corresponding=letter_corresponding+letter_corresponding
        node["Letter"]=str(node["id"])
        allready_used_letters.add(letter_corresponding)
    Json_obj["nodes"]=nodes
    return Json_obj

def generation_des_kmers(une_liste_alphabet):
    k_mer=[]
    for i in range(0, len(une_liste_alphabet)):
        letter=une_liste_alphabet[i]
        if len(k_mer)==0:
            k_mer.append(letter)
            k_mer.append("~"+letter)
        else:
            new_kmer_list=[]
            for j in range(0, len(k_mer)):
                k_mer_str=k_mer[j]
                k_mer_list=k_mer_str.split("\t")
                if letter in k_mer_list:
                    continue
                else:
                    new_kmer_string=k_mer_str+"\t"+letter
                    new_kmer_string2=k_mer_str+"\t"+"~"+letter
                    new_kmer_list.append(new_kmer_string)
                    new_kmer_list.append(new_kmer_string2)
            k_mer=new_kmer_list
    return k_mer

def clean_json(json_dic, rank1_tool_name, rank2_tool_name):
    rooted_nodes=set([0])
    for links in json_dic["links"]:
        source=links["source"]
        if source==0:
            enzyme=links["enzymes"]
            if rank1_tool_name in list(enzyme.keys()):
                rooted_nodes.add(int(links["target"]))
    for links in json_dic["links"]:
        source=links["source"]
        if source!=0 and source in rooted_nodes:
            enzyme=links["enzymes"]
            if rank2_tool_name in list(enzyme.keys()):
                rooted_nodes.add(int(links["target"]))
    to_del_links=[]
    to_del_nodes=[]
    for link_tuple in enumerate(json_dic["links"]):
        link_index=link_tuple[0]
        link=link_tuple[1]
        source = link["source"]
        target = link["target"]
        if source not in rooted_nodes or target not in rooted_nodes:
            to_del_links.append(link_index)

    for node_tuple in enumerate(json_dic["nodes"]):
        node_index=node_tuple[0]
        node=node_tuple[1]
        node_id=node["id"]
        if node_id not in rooted_nodes:
            to_del_nodes.append(node_index)

    to_del_links=sorted(to_del_links, reverse=True)
    to_del_nodes=sorted(to_del_nodes, reverse=True)
    for node_id in to_del_nodes:
        del json_dic["nodes"][node_id]
    for link_id in to_del_links:
        del json_dic["links"][link_id]

def get_nodes_dependencies(json_dic, rank1_tool_name, rank2_tool_name, cyto_list, SMIRKS_label_dic_path):
    node_depencies={}
    cyps_Accerlys_names, ugts_Accerlys_names, nats_Accerlys_names, gsts_Accerlys_names, sults_Accerlys_names = read_SMIRKS_label_dic(SMIRKS_label_dic_path)
    nodes=json_dic["nodes"]
    links=json_dic["links"]
    for node in nodes:
        node_letter=node["Letter"]
        node_id=node["id"]
        if node_id==0:
            node_depencies[node_letter]={}
            continue
        for reaction in links:
            target=reaction["target"]
            if target == node_id:
                source=reaction["source"]
                parent_letter=None
                for node2 in nodes:
                    if node2["id"]==source:
                        parent_letter=node2["Letter"]
                        break
                enzymes=reaction["enzymes"]
                for tool in enzymes.keys():
                    tool_enzymes_score=enzymes[tool]
                    if tool == rank2_tool_name and source!=0:
                        reaction_name=reaction["reaction_name"]
                        if reaction_name in nats_Accerlys_names and rank2_tool_name=="FAME3":
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                a_node_dependencie[(parent_letter, "NATs")]=tool_enzymes_score
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                node_depencies[node_letter]={(parent_letter, "NATs"):tool_enzymes_score}

                        elif reaction_name in cyps_Accerlys_names and rank2_tool_name=="FAME3":
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                for cyto in cyto_list:
                                    a_node_dependencie[(parent_letter, cyto)]=tool_enzymes_score
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                node_depencies[node_letter]={}
                                a_node_dependencie=node_depencies[node_letter]
                                for cyto in cyto_list:
                                    a_node_dependencie[(parent_letter, cyto)]=tool_enzymes_score
                                node_depencies[node_letter]=a_node_dependencie
                                
                        elif reaction_name in cyps_Accerlys_names and rank2_tool_name!="FAME3":
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                for cytochrome in tool_enzymes_score.keys():
                                    score_for_isoform=tool_enzymes_score[cytochrome]
                                    a_node_dependencie[(parent_letter, cytochrome)]=score_for_isoform
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                node_depencies[node_letter]={}
                                a_node_dependencie=node_depencies[node_letter]
                                for cytochrome in tool_enzymes_score.keys():
                                    score_for_isoform=tool_enzymes_score[cytochrome]
                                    a_node_dependencie[(parent_letter, cytochrome)]=score_for_isoform
                                node_depencies[node_letter]=a_node_dependencie
                        

                        elif reaction_name in ugts_Accerlys_names and rank2_tool_name=="FAME3":
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                a_node_dependencie[(parent_letter, "UGTs")]=tool_enzymes_score
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                node_depencies[node_letter]={(parent_letter, "UGTs"):tool_enzymes_score}
                        
                        elif reaction_name in ugts_Accerlys_names and rank2_tool_name!="FAME3":
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                score_for_isoform=tool_enzymes_score["UGT"]
                                a_node_dependencie[(parent_letter, "UGTs")]=score_for_isoform
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                score_for_isoform=tool_enzymes_score["UGT"]
                                node_depencies[node_letter]={(parent_letter, "UGTs"):score_for_isoform}
                        

                        elif reaction_name in gsts_Accerlys_names and rank2_tool_name=="FAME3":
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                a_node_dependencie[(parent_letter, "GSTs")]=tool_enzymes_score
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                node_depencies[node_letter]={(parent_letter, "GSTs"):tool_enzymes_score}

                        elif reaction_name in sults_Accerlys_names and rank2_tool_name=="FAME3":
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                a_node_dependencie[(parent_letter, "SULTs")]=tool_enzymes_score
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                node_depencies[node_letter]={(parent_letter, "SULTs"):tool_enzymes_score}

                    elif tool==rank1_tool_name and source==0:
                        for enzyme_name in tool_enzymes_score.keys():
                            enzyme_score=tool_enzymes_score[enzyme_name]
                            if node_letter in node_depencies.keys():
                                a_node_dependencie=node_depencies[node_letter]
                                a_node_dependencie[(parent_letter, enzyme_name.replace("CYP","").replace("UGT", "UGTs"), tool)]=enzyme_score
                                node_depencies[node_letter]=a_node_dependencie
                            else:
                                node_depencies[node_letter]={(parent_letter, enzyme_name.replace("CYP","").replace("UGT", "UGTs"), tool):enzyme_score}
    return node_depencies

def get_phaseI_node(json_obj):
    node_id_set=set([])
    for link in json_obj["links"]:
        source=link["source"]
        if source==0:
            node_id_set.add(link["target"])
    return(sorted(list(node_id_set)))

def make_cluster(json_obj, phaseI_nodes):
    clusters_dictionnary={}
    dictionnary_nodes_to_cluster_index={}
    for i in range(0,len(phaseI_nodes)):
        node_id=phaseI_nodes[i]
        dictionnary_nodes_to_cluster_index[node_id]=set([i])
        clusters_dictionnary[i]=set([node_id])

    for a_reaction in json_obj["links"]:
        source=a_reaction["source"]
        target=a_reaction["target"]
        if source in dictionnary_nodes_to_cluster_index.keys():
            all_cluster_index_list=list(dictionnary_nodes_to_cluster_index[source])
            for a_cluster_index in all_cluster_index_list:
                clusters_dictionnary[a_cluster_index].add(target)
                if target in dictionnary_nodes_to_cluster_index:
                    dictionnary_nodes_to_cluster_index[target].add(a_cluster_index)
                else:
                    dictionnary_nodes_to_cluster_index[target]=set([a_cluster_index])
    return clusters_dictionnary

def merge_cluster(clusters_dictionnary):
    clusters_to_merge={}
    for a_cluster_index1 in clusters_dictionnary.keys():
        for a_cluster_index2 in clusters_dictionnary.keys():
            if a_cluster_index1==a_cluster_index2:
                continue
            cluster_nodes_set_1=clusters_dictionnary[a_cluster_index1]
            cluster_nodes_set_2=clusters_dictionnary[a_cluster_index2]
            cluster_intersection=cluster_nodes_set_1.intersection(cluster_nodes_set_2)
            not_empty=bool(cluster_intersection)
            if not_empty:
                if a_cluster_index1 not in clusters_to_merge.keys():
                    clusters_to_merge[a_cluster_index1]=set([a_cluster_index2])
                if a_cluster_index2 not in clusters_to_merge.keys():
                    clusters_to_merge[a_cluster_index2]=set([a_cluster_index1])
                for cluster_indexes_also_to_merge in clusters_to_merge[a_cluster_index1]:
                    if cluster_indexes_also_to_merge!=a_cluster_index2:
                        clusters_to_merge[cluster_indexes_also_to_merge].add(a_cluster_index2)
                for cluster_indexes_also_to_merge in clusters_to_merge[a_cluster_index2]:
                    if cluster_indexes_also_to_merge!=a_cluster_index1:
                        clusters_to_merge[cluster_indexes_also_to_merge].add(a_cluster_index1)
                clusters_to_merge[a_cluster_index1].add(a_cluster_index2)
                clusters_to_merge[a_cluster_index2].add(a_cluster_index1)
    keeped_cluster_indexes=set([])

    for key in clusters_to_merge:
        a_merge_group=set([key])
        a_merge_group=a_merge_group.union(clusters_to_merge[key])
        to_keep_cluster=min(list(a_merge_group))
        keeped_cluster_indexes.add(to_keep_cluster)

    for cluster_index in keeped_cluster_indexes:
        cluster_to_keep_set=clusters_dictionnary[cluster_index]
        to_merge_indexes=clusters_to_merge[cluster_index]
        for index in to_merge_indexes:
            to_merge_cluster_set=clusters_dictionnary[index]
            cluster_to_keep_set=cluster_to_keep_set.union(to_merge_cluster_set)
            del clusters_dictionnary[index]
        clusters_dictionnary[cluster_index]=cluster_to_keep_set
    output=[]
    for key in clusters_dictionnary.keys():
        to_str_convert=list(sorted(clusters_dictionnary[key]))
        str_converted=[]
        for elem in to_str_convert:
            new_elem=str(elem)
            str_converted.append(new_elem)
        output.append(str_converted)
    return(output)

def get_new_alphabets(json_dic):
    phaseI=get_phaseI_node(json_dic)
    clusters=make_cluster(json_dic, phaseI)
    new_alphabet=merge_cluster(clusters)
    return(new_alphabet)

def generate_all_cases_nodes2(Json_obj):
    all_cases_list=[]
    clustered_alphabet=get_new_alphabets(Json_obj)
    for alphabet_cluster in clustered_alphabet:
        cas=generation_des_kmers(alphabet_cluster)
        all_cases_list.append(cas)
    return(all_cases_list)

def get_present_parent_and_enzyme(enzyme_case, node_case, letter_dependencie_dic):
    parent_enzymes=[]
    parent_nodes=[]
    for key in letter_dependencie_dic.keys():
        parent=key[0]
        enzyme=key[1]
        parent_enzymes.append(enzyme)
        parent_nodes.append(parent)
    present_parent=[]
    prensent_enzyme=[]
    enzyme_case=enzyme_case.split("\t")
    node_case=node_case.split("\t")
    for an_enz in enzyme_case:
        if "~" in an_enz:
            continue
        else:
            if an_enz in parent_enzymes:
                prensent_enzyme.append(an_enz)
    for a_parent in node_case:
        if "~" in a_parent:
            continue
        else:
            if a_parent in parent_nodes:
                present_parent.append(a_parent)
    to_delete_parent=[]
    for parent in present_parent:
        del_parent=True
        for key in letter_dependencie_dic.keys():
            if key[0]==parent:
                enz=key[1]
                if enz in prensent_enzyme:
                    del_parent=False
        if del_parent:
            to_delete_parent.append(parent)
    for parent in to_delete_parent:
        while parent in present_parent:
            present_parent.remove(parent)
    return (present_parent, prensent_enzyme)

def sort_cases(all_cases, node_dependencies_dic):
    new_all_cases=[]
    for a_familly_cases in all_cases:
        new_familly_cases=[]
        for a_case in a_familly_cases:
            keep_case=True
            present_node=set([])
            a_case_splited=a_case.split('\t')
            for node in a_case_splited:
                if '~' not in node:
                    present_node.add(node)
            for node in present_node:
                has_parent=False
                node_dependencies=node_dependencies_dic[node]
                for parent_tuple in node_dependencies:
                    parent=parent_tuple[0]
                    if parent=='0':
                        has_parent=True
                        break
                    if parent in present_node:
                        has_parent=True
                        break
                if not has_parent:
                    keep_case=False
                    break
            if keep_case:
                new_familly_cases.append(a_case)
        new_all_cases.append(new_familly_cases)
    return new_all_cases

def sort_cases_for_enzymes(enzymes_cases, sorted_cases, node_dependencies_dic):
    output={}
    for an_enzyme_case in enzymes_cases:
        splited_enzyme_case=an_enzyme_case.split('\t')
        present_enzyme=set([])
        for enzyme in splited_enzyme_case:
            if '~' not in enzyme:
                present_enzyme.add(enzyme)
        new_sorted_case=[]
        for case_familly in sorted_cases:
            new_case_familly=[]
            for case in case_familly:
                keep_case=True
                splitted_case=case.split('\t')
                present_node=set([])
                for node in splitted_case:
                    if '~' not in node:
                        present_node.add(node)
                for node in present_node:
                    has_parent=False
                    node_dependencies=node_dependencies_dic[node]
                    for parent_tuple in node_dependencies:
                        parent=parent_tuple[0]
                        enzyme=parent_tuple[1]
                        if parent=='0':
                            if enzyme in present_enzyme:
                                has_parent=True
                                break
                        else:
                            if enzyme in present_enzyme and parent in present_node:
                                has_parent=True
                                break
                    if not has_parent:
                        keep_case=False
                        break
                if keep_case:
                    new_case_familly.append(case)
            if new_case_familly!=[]:
                new_sorted_case.append(new_case_familly)
        output[an_enzyme_case]=new_sorted_case
    return output

def check_json(json_obj):
    reactions=json_obj["links"]
    nodes=json_obj["nodes"]
    is_empty=False
    if reactions==[]:
        if len(nodes)==1:
            principal_node=nodes[0]
            if principal_node["id"]==0:
                is_empty=True
    return(is_empty)

def write_general_probability(directory, general_probability_dic, rank1_tool_name, rank2_tool_name):
    if directory[len(directory)-1]!="/":
        directory=directory+"/"
    directory=directory+"Probabilities_results/"+str(rank1_tool_name)+"_as_rank1_and_"+rank2_tool_name+"_as_rank2/"
    for letter in general_probability_dic.keys():
        if '~' in letter:
            continue
        letter_dic=general_probability_dic[letter]
        if not os.path.exists(directory):
            os.makedirs(directory)
        if not os.path.exists(directory+letter+"/"):
            os.makedirs(directory+letter+"/")
        new_path=directory+letter+"/"+"general_probability_results.tsv"
        file_output=open(new_path,'w')
        for enzyme_case in letter_dic.keys():
            probability=letter_dic[enzyme_case]
            to_print_probability=str(round(probability,4))
            file_output.write(enzyme_case.replace('\t','|')+'\t'+to_print_probability+'\n')
        file_output.close()

def duplicate_dependencie_dic(present_parent, present_enzymes, letter_dependencie_dic):
    new_depend={}
    for key in letter_dependencie_dic.keys():
        parent=key[0]
        enzyme=key[1]
        if parent in present_parent and enzyme in present_enzymes:
            new_depend[key]=letter_dependencie_dic[key]
    return new_depend

def reduce_combinaisons(list_of_combinaison, present_enzymes, present_parent):
    duplicate_tuple={}
    for a_combinaison_tuple1 in enumerate(list_of_combinaison):
        for a_combinaison_tuple2 in enumerate(list_of_combinaison):
            index1=a_combinaison_tuple1[0]
            index2=a_combinaison_tuple2[0]
            if index1 not in duplicate_tuple.keys():
                duplicate_tuple[index1]=set([])
            if index2 not in duplicate_tuple.keys():
                duplicate_tuple[index2]=set([])
            if index1!=index2:
                tuple_list1=a_combinaison_tuple1[1]
                tuple_list2=a_combinaison_tuple2[1]
                tuple1_set=set()
                tuple2_set=set()
                for a_tuple1_str in tuple_list1:
                    a_tuple1=literal_eval(a_tuple1_str)
                    a_tuple1=list(a_tuple1)
                    tuple1_set.add('-'.join(a_tuple1))
                for a_tuple2_str in tuple_list2:
                    a_tuple2=literal_eval(a_tuple2_str)
                    a_tuple2=list(a_tuple2)
                    tuple2_set.add('-'.join(a_tuple2))
                tuple_difference=bool(tuple1_set.difference(tuple2_set))
                if not tuple_difference:
                    if index1 in duplicate_tuple.keys():
                        to_complete=duplicate_tuple[index1]
                        to_complete.add(index2)
                        duplicate_tuple[index1]=to_complete
                    if index1 not in duplicate_tuple.keys():
                        duplicate_tuple[index1]=set([index2])
                    if index2 in duplicate_tuple.keys():
                        to_complete=duplicate_tuple[index2]
                        to_complete.add(index1)
                        duplicate_tuple[index2]=to_complete
                    if index2 not in duplicate_tuple.keys():
                        duplicate_tuple[index2]=set([index1])

    new_combinaison=[]
    used_index=set([])
    for a_combinaison_tuple in enumerate(list_of_combinaison):
        index=a_combinaison_tuple[0]
        combinaison_case=a_combinaison_tuple[1]
        if index not in used_index:
            new_combinaison.append(combinaison_case)
            other_indexes=duplicate_tuple[index]
            used_index.add(index)
            for another_index in other_indexes:
                used_index.add(another_index)
    return(new_combinaison)

def generation_des_combinaisons(dependencie_letter_keys, combinaison_length, present_enzymes, present_parent):
    str_dependencie_letter_keys=[]
    for a_tuple_key in dependencie_letter_keys:
        str_dependencie_letter_keys.append(str(a_tuple_key))
    k_mer=['-'.join(p) for p in itertools.product(str_dependencie_letter_keys, repeat=combinaison_length)]
    combinaison=[]
    for a_combinaison in k_mer:
        combinaison_break=False
        a_combinaison=a_combinaison.split('-')
        for i in range(0,2):
            to_check_set=set()
            for a_tuple in a_combinaison:
                a_tuple=literal_eval(a_tuple)
                to_check_set.add(a_tuple[i])
            to_check_set=list(to_check_set)
            if len(to_check_set)!=combinaison_length:
                combinaison_break=True
        if combinaison_break:
            continue
        else:
            combinaison.append(a_combinaison)
    combinaison=reduce_combinaisons(combinaison, present_enzymes, present_parent)
    return (combinaison)

def calculate_combinaison_probability(a_combinaison, letter_dependencie_dic, node_exist):
    not_exist_probability=1.00
    try:
        for element in a_combinaison:
            element_tuple=literal_eval(element)
            element_probability = letter_dependencie_dic[element_tuple]
            not_exist_probability=not_exist_probability*(1-element_probability)
        if node_exist:
            probability=1.00-not_exist_probability
            return probability
        else:
            return not_exist_probability
    except TypeError:
        print(a_combinaison)

def calculate_node_case_probability(present_parents, present_enzymes, letter_dependencies, node_state):
    case_dependencie_dic=duplicate_dependencie_dic(present_parents, present_enzymes, letter_dependencies)
    if case_dependencie_dic.keys():
        # Il y a une source parent + enzyme minimum:
        combinaisons_length=[len(present_parents), len(present_enzymes)]
        combinaisons_length=min(combinaisons_length)
        combinaisons_parents_enzymes=generation_des_combinaisons(list(case_dependencie_dic.keys()), combinaisons_length, present_enzymes, present_parents)
        maximal_score=0.00
        final_combinaison=combinaisons_parents_enzymes[0]
        for a_combinaison in combinaisons_parents_enzymes:
            combinaison_score=0
            for elem in a_combinaison:
                elem_tuple=literal_eval(elem)
                score=letter_dependencies[elem_tuple]
                combinaison_score=combinaison_score+score
            if combinaison_score>maximal_score:
                maximal_score=combinaison_score
                final_combinaison=a_combinaison
        if "~" in node_state:
            probability=calculate_combinaison_probability(final_combinaison, letter_dependencies, False)
            return probability
        else:
            probability=calculate_combinaison_probability(final_combinaison, letter_dependencies, True)
            return probability
    else:
        if "~" in node_state:
            probability=1.00
            return probability
        else:
            probability=0.00
            return probability

def fonction_unitaire(cases_list, an_enzyme_case, node_dependencie_dic, node_general_prob_dic):
    for a_case in cases_list:
        a_case_with_source="0\t"+a_case
        probability=Decimal(1.0)               
        splitted_case=a_case_with_source.split('\t')
        for letter in splitted_case:
            if letter=='0':
                continue
            elif letter=='~0':
                probability=0
                break
            else:
                present_parent, prensent_enzyme=get_present_parent_and_enzyme(an_enzyme_case, a_case_with_source, node_dependencie_dic[letter.replace("~","")])
                letter_proba=calculate_node_case_probability(present_parent, prensent_enzyme, node_dependencie_dic[letter.replace("~","")], letter)
                probability=probability*Decimal(letter_proba)
        for letter in splitted_case:                                    
            if letter not in list(node_general_prob_dic.keys()):
                node_general_prob_dic[letter]={}
            letter_dic=node_general_prob_dic[letter]
            if an_enzyme_case not in list(letter_dic.keys()):
                letter_dic[an_enzyme_case]=probability
            elif an_enzyme_case in list(letter_dic.keys()):
                letter_dic[an_enzyme_case]=letter_dic[an_enzyme_case]+probability    

def get_parents(node, node_dependencies_dic):
    parents=set([])
    node_dependencies=node_dependencies_dic[node]
    for parent_tuple in node_dependencies.keys():
        parent=parent_tuple[0]
        parents.add(parent)
    parents=list(parents)
    if '0' in parents:
        parents=set(parents)
        return(parents)
    else:
        for a_parent in parents:
            parent_of_parents=get_parents(a_parent, node_dependencies_dic)
            parent_of_parents=set(parent_of_parents)
            parents=set(parents)
            parents=parents.union(parent_of_parents)
        return(parents)

def new_sort_case(node_parents_cases, final_node):
    new_cases=[]
    if node_parents_cases==[]:
        new_cases=[final_node]
        return new_cases
    for case in node_parents_cases:
        at_least_a_present_parent=False
        a_case=case.split('\t')
        for node in a_case:
            if '~' not in node:
                at_least_a_present_parent=True
                break
        if at_least_a_present_parent:
            a_case.append(final_node)
            to_appened='\t'.join(a_case)
            new_cases.append(to_appened)
    return(new_cases)

def new_sort_cases_for_enzymes(enzymes_cases, sorted_cases, node_dependencies_dic):
    output={}
    for an_enzyme_case in enzymes_cases:
        splited_enzyme_case=an_enzyme_case.split('\t')
        present_enzyme=set([])
        for enzyme in splited_enzyme_case:
            if '~' not in enzyme:
                present_enzyme.add(enzyme)
        new_sorted_case=[]
        for case in sorted_cases:
            keep_case=True
            splitted_case=case.split('\t')
            present_node=set([])
            for node in splitted_case:
                if '~' not in node:
                    present_node.add(node)
            for node in present_node:
                has_parent=False
                node_dependencies=node_dependencies_dic[node]
                for parent_tuple in node_dependencies:
                    parent=parent_tuple[0]
                    enzyme=parent_tuple[1]
                    if parent=='0':
                        if enzyme in present_enzyme:
                            has_parent=True
                            break
                    else:
                        if enzyme in present_enzyme and parent in present_node:
                            has_parent=True
                            break
                if not has_parent:
                    keep_case=False
                    break
            if keep_case:
                new_sorted_case.append(case)
        output[an_enzyme_case]=new_sorted_case
    return output

def calculate_node_probability(node, node_dependencies, enzyme_cases):
    parents=get_parents(node, node_dependencies)
    parents.discard('0')
    parents=sorted(list(parents))
    all_parents_cases=generation_des_kmers(parents)
    sorted_parents_cases=new_sort_case(all_parents_cases, node)
    enzyme_sorted_cases=new_sort_cases_for_enzymes(enzyme_cases, sorted_parents_cases, node_dependencies)
    node_probability_dictionnary={}
    for an_enzyme_case2 in enzyme_sorted_cases.keys(): 
            node_cases_list=enzyme_sorted_cases[an_enzyme_case2]
            letter_proba_f=Decimal(0.0)
            for a_case in node_cases_list:
                a_case="0\t"+a_case
                case_probability=Decimal(1.0)
                a_splitted_case=a_case.split('\t')
                for letter in a_splitted_case:
                    if letter=='0':
                        continue
                    present_parent, prensent_enzyme=get_present_parent_and_enzyme(an_enzyme_case2, a_case, node_dependencies[letter.replace("~","")])
                    letter_proba=calculate_node_case_probability(present_parent, prensent_enzyme, node_dependencies[letter.replace("~","")], letter)
                    case_probability=case_probability*Decimal(letter_proba)
                letter_proba_f=letter_proba_f+case_probability
            node_probability_dictionnary[an_enzyme_case2]=letter_proba_f
    return(node_probability_dictionnary)

def get_output_directory_from_map(metabolism_map_path):
    output_dir=None
    splitted_metabolism_map_path=metabolism_map_path.split('/')
    output_dir=splitted_metabolism_map_path[:len(splitted_metabolism_map_path)-1]
    output_dir='/'.join(output_dir)
    output_dir=output_dir+"/Metabolites_probabilities/"
    return(output_dir)

#stime=time.time()

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("map_file", help="the .json file path where the metabolic map is saved with score of FAME3, W2D and Xenosite on reactions")
    parser.add_argument("SMIRKS_label_file", help="the path to file where SMIRKS labels are associated to enzymes")
    parser.add_argument("--output", help="the output directory where probabilities will be written", default=False)
    parser.add_argument("--rank1_tool", help="the name of the tool you want to use for reactions of rank 1. The name has to be \'Xenosite\' OR \'Way2Drug\' with this exact spelling (lower or upper case has no importance)", default="Way2Drug")
    parser.add_argument("--rank2_tool", help="the name of the tool you want to use for reactions of rank 2. The name has to be \'Xenosite\' OR \'Way2Drug\' OR \'FAME3\' with this exact spelling (lower or upper case has no importance)", default="FAME3")

    args = parser.parse_args()
    rank1_tool_name=args.rank1_tool
    rank2_tool_name=args.rank2_tool
    SMIRKS_label_dic_path=args.SMIRKS_label_file
    path_to_metabolism_map=args.map_file
    output_directory=args.output
    if not output_directory:
        output_directory=get_output_directory_from_map(path_to_metabolism_map)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    metabolism_map_file=open(path_to_metabolism_map, 'r')
    metabolism_map=json.loads(metabolism_map_file.read())
    metabolism_map_file.close()

    is_metabolism_map_empty=check_json(metabolism_map)
    if is_metabolism_map_empty:
        general_prob_test={"0":{"empty graph":1.0}}
        write_general_probability(output_directory, general_prob_test, rank1_tool_name, rank2_tool_name)
    else:
        enzyme_list, cytochrome_list=get_all_enzyme(metabolism_map, [rank1_tool_name, rank2_tool_name], SMIRKS_label_dic_path)
        metabolism_map=define_letters_for_nodes(metabolism_map)
        all_enzyme_cases=generation_des_kmers(enzyme_list)
        clean_json(metabolism_map, rank1_tool_name, rank2_tool_name)
        all_node_dependencies=get_nodes_dependencies(metabolism_map, rank1_tool_name, rank2_tool_name, cytochrome_list, SMIRKS_label_dic_path)
        general_probability_dictionnary={}
        for node in all_node_dependencies.keys():
            if node=="0":
                general_probability_dictionnary["0"]={}
                for enzyme_case in all_enzyme_cases:
                    general_probability_dictionnary["0"][enzyme_case]=Decimal(1.0)
            else:
                node_probabilities=calculate_node_probability(node, all_node_dependencies, all_enzyme_cases)
                general_probability_dictionnary[node]=node_probabilities
        write_general_probability(output_directory, general_probability_dictionnary, rank1_tool_name, rank2_tool_name)

main()
