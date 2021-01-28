#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 14:08:41 2020

@author: mconan
"""
import os
import json
import sygma
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# Each step in a scenario lists the ruleset and the number of reaction cycles to be applied
scenario2 = sygma.Scenario([
    [sygma.ruleset['phase1'], 1],
    [sygma.ruleset['phase2'], 1]])

scenario = sygma.Scenario([
    [sygma.ruleset['phase1'], 1]])

scenario3 = sygma.Scenario([
    [sygma.ruleset['phase1'], 2]])

# Sygma utilise des réactions SMARTS en utilisant RDKIT et une classe de RDKIT avec la fonction runreactants.
# Celle ci retourne un tuple de tuples mais SyGMa ne retient que le premier élément de ce tuple ( [0] )
# Il n'y a donc jamais l'information de l'atome qui réagit il faut donc bien le retrouver autrement ou alors retravailler SyGMa

def predict_metabolite(smile, scenario, image_output_dir):
    dict_to_return={"nodes":[],"links":[]}
    parent = Chem.MolFromSmiles(smile)
    metabolic_tree = scenario.run(parent)
    nodes=metabolic_tree.nodes
    source_parent_key=metabolic_tree.parentkey
    source_parent_node=nodes[source_parent_key]
    source_parent_mol=source_parent_node.mol
    source_parent_smile=Chem.MolToSmiles(source_parent_mol)
    image_output=image_output_dir+str(0)+".svg"
    Draw.MolToFile(source_parent_mol, image_output, size=(170, 170), includeAtomNumbers=False)
    Draw.MolToFile(source_parent_mol, image_output.replace(".svg","_atom_num.svg"), size=(170, 170), includeAtomNumbers=True)
    parent_node_dic={"id":0, "Original_smile":smile,"smile_in_sygma_tree":source_parent_smile, "node_key":source_parent_key, "img_svg":"0"}
    dict_to_return["nodes"].append(parent_node_dic)
    id_index=0
    reaction_id=0
    for key_node in nodes.keys():
        if key_node!=source_parent_key:
            id_index+=1
            node=nodes[key_node]
            node_chem_obj=node.mol
            node_smile=Chem.MolToSmiles(node_chem_obj)
            parents=node.parents
            node_image_output=image_output_dir+str(id_index)+".svg"
            Draw.MolToFile(node_chem_obj, node_image_output, size=(170, 170), includeAtomNumbers=False)
            Draw.MolToFile(node_chem_obj, node_image_output.replace(".svg","_atom_num.svg"), size=(170, 170), includeAtomNumbers=True)
            node_dic={"id":id_index,"smile_in_sygma_tree":node_smile, "node_key":key_node, "img_svg":str(id_index)}
            dict_to_return["nodes"].append(node_dic)
            for parent_id in parents:
                reaction_id+=1
                rule_used=parents[parent_id]
                node_parent=nodes[parent_id]
                parent_mol=node_parent.mol
                parent_smile=Chem.MolToSmiles(parent_mol)
                rule_name=rule_used.rulename
                reaction_dic={"id":reaction_id, "source_key":parent_id, "target_key":key_node, "reaction_name":rule_name, "source_smile":parent_smile, "target_smile":node_smile, "atom_number_of_reaction":0}
                for node_json in dict_to_return["nodes"]:
                    if node_json["node_key"]==parent_id:
                        reaction_dic["source"]=node_json["id"]
                    if node_json["node_key"]==key_node:
                        reaction_dic["target"]=node_json["id"]
                dict_to_return["links"].append(reaction_dic)
    return(dict_to_return)

def write_new_json(json_object, output):
    output_file=open(output, 'w')
    output_file.write(json.dumps(json_object, indent=4))
    output_file.close()
    
scenario_AHA = sygma.Scenario([
    [sygma.ruleset['phase1'], 2],
    [sygma.ruleset['phase2'], 2]])
# On pourrait refaire les cartes comme la caffeine mais alors la stratégie serait différente de celle pourtant déjà adopté par l'équipe en 2019.
# Pour être cohérent et ne pas à avoir à justifier la stratégie, et pouvoir se comparer à la fin, on définis un scénario avec 2 cycles de phase I e 2 cycles de phase II en même temps.
# La caffeine justifiait une nouvelle stratégie pour être au plus proche de ce que l'on connait de son métabolisme, ici on ne cherchait pas à découvrir de nouvelles choses mais à valider en suivant à la lettre son métabolisme.
# On se mettait dans la meilleur des situations car c'est la pire situation pour le modèle bayesien car peut de surprédiction et des molécule avec une forte confiance.
# Ce n'est pas le cas des AHAs

AHAs_smiles={
    "1,5,6-TMIP":"N=1C=2C=C(C(=CC2N(C1N)C)C)C",
    "1,6-DMIP":"N=1C=2C=CC(=CC2N(C1N)C)C",
    "3,5,6-TMIP":"N=1C=2C=C(C(=NC2N(C1N)C)C)C",
    "4,7,8-TriMeIQx":"N=1C=2C=C(C3=C(N=C(N)N3C)C2N=C(C1C)C)C",
    "4,8-DiMeIQx":"N1=CC(=NC2=C1C=C(C3=C2N=C(N)N3C)C)C",
    "4-CH2OH-8-MeIQx":"OCC1=CC=2N=CC(=NC2C=3N=C(N)N(C31)C)C",
    "4'-OH-PhIP":"OC=1C=CC(=CC1)C=2C=NC=3N=C(N)N(C3C2)C",
    "6,7-DiMeIgQx":"N=1C=2C=C3N=C(C(=NC3=CC2N(C1N)C)C)C",
    "7,8-DiMeIQx":"N=1C=2C=CC3=C(N=C(N)N3C)C2N=C(C1C)C",
    "7,9-DiMeIgQx":"N1=CC(=NC=2C1=CC=3N=C(N)N(C3C2C)C)C",
    "7-MeIgQx":"N1=CC(=NC2=CC3=C(N=C(N)N3C)C=C12)C",
    "AaC":"N=1C(N)=CC=C2C1NC=3C=CC=CC32",
    "AMPNH":"N=1C=CC=2C=3C=CC=CC3N(C4=CC=C(N)C(=C4)C)C2C1",
    "APNH":"N=1C=CC=2C=3C=CC=CC3N(C4=CC=C(N)C=C4)C2C1",
    "GluP1":"N1=C(N)C=CC=2N=C3C(=CC=CN3C12)C",
    "GluP2":"N1=C(N)C=CC=2N=C3C=CC=CN3C12",
    "Harman":"N=1C=CC=2C=3C=CC=CC3NC2C1C",
    "IFP":"N1=C2N=C(N)N(C2=CC=3OC(=CC13)C)C",
    "IgQx":"N1=CC=NC2=CC3=C(N=C(N)N3C)C=C12",
    "IQ":"N1=CC=CC2=C1C=CC3=C2N=C(N)N3C",
    "IQ[4,5-b]":"N1=C2N=C(N)N(C2=CC3=CC=CC=C13)C",
    "IQx":"N1=CC=NC2=C1C=CC3=C2N=C(N)N3C",
    "MeAaC":"N=1C(N)=C(C=C2C1NC=3C=CC=CC32)C",
    "MeIQ":"N1=CC=CC2=C1C=C(C3=C2N=C(N)N3C)C",
    "MeIQx":"N1=CC(=NC2=C1C=CC3=C2N=C(N)N3C)C",
    "NorHarman":"N=1C=CC2=C(C1)NC=3C=CC=CC32",
    "PheP1":"N1=CC(=CC=C1N)C2=CC=CC=C2",
    "PhIP":"N1=CC(=CC2=C1N=C(N)N2C)C=3C=CC=CC3",
    "TrP1":"N=1C(N)=C(C=2NC=3C=CC=CC3C2C1C)C",
    "TrP2":"N1=C(N)C=C2NC=3C=CC=CC3C2=C1C"
}

Aha_result_dir="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/All_AHAs/"
for AHA in AHAs_smiles.keys():
    if not os.path.exists(Aha_result_dir+AHA):
        os.makedirs(Aha_result_dir+AHA+"/Mol_img/")
    SMILES=AHAs_smiles[AHA]
    AHA_output=Aha_result_dir+AHA+"/"
    AHA_img_output=Aha_result_dir+AHA+"/Mol_img/"
    AHA_predicted_map=predict_metabolite(SMILES, scenario_AHA, AHA_img_output)
    AHA_map_output=AHA_output+AHA+".json"
    write_new_json(AHA_predicted_map, AHA_map_output)

#AaC_predicted_map=predict_metabolite("N=1C(N)=CC=C2C1NC=3C=CC=CC32", scenario_AHA)
#MeIQx_predicted_map=predict_metabolite("N1=CC(=NC2=C1C=CC3=C2N=C(N)N3C)C", scenario_AHA)
#PhIP_predicted_map=predict_metabolite("N1=CC(=CC2=C1N=C(N)N2C)C=3C=CC=CC3", scenario_AHA)

#AaC_predicted_map_output="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/AaC_predicted_map_14_10_20.json"
#MeIQx_predicted_map_output="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/MeIQx_predicted_map_14_10_20.json"
#PhIP_predicted_map_output="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/PhIP_predicted_map_14_10_20.json"

# test1=Chem.MolFromSmiles("Cn1c(N(O)C2OC(C(=O)O)C(O)C(O)C2O)[n+](C2OC(C(=O)O)C(O)C(O)C2O)c2ncc(-c3ccccc3)c(O)c21")
# Draw.MolToFile(test1, "/home/mconan/Documents/Bayesian_model/AHA_with_sygma/"+"test_name.svg", size=(200, 200), includeAtomNumbers=False)


#write_new_json(AaC_predicted_map, AaC_predicted_map_output)
#write_new_json(MeIQx_predicted_map, MeIQx_predicted_map_output)
#write_new_json(PhIP_predicted_map, PhIP_predicted_map_output)
