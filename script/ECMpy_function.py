# -*- coding: utf-8 -*-
# This code is used to introduce enzyme concentration constraint in GEMs
# by COBRApy and to calculate the parameters that need to be entered
# during the construction of the enzyme-constrained model.
#from warnings import warn

import cobra
import math
import random
import statistics
import sys
import pandas as pd
import json
import re
import copy
from copy import deepcopy
import numpy as np
from typing import Any, Dict, List
from cobra.core import Reaction
from cobra.util.solver import set_objective
from cobra.io.dict import model_to_dict
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from optlang.symbolics import Zero, add
import plotly.graph_objects as go

def json_load(path: str) -> Dict[Any, Any]:
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder

def convert_to_irreversible(model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.

    """
    #warn("deprecated, not applicable for optlang solvers", DeprecationWarning)
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in reaction._metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)
    
def get_genes_and_gpr(model,gene_outfile,gpr_outfile):
    """Retrieving genes and gene_reaction_rule from GEM.

    Arguments
    ----------
    * model: cobra.Model ~ A genome scale metabolic network model for
        constructing the enzyme-constrained model.

    :return: all genes and gpr in model.
    """
    model_dict = model_to_dict(model, sort=False)
    genes = pd.DataFrame(model_dict['genes']).set_index(['id'])
    genes.to_csv(gene_outfile)
    all_gpr = pd.DataFrame(model_dict['reactions']).set_index(['id'])
    all_gpr.to_csv(gpr_outfile)
    return [genes, all_gpr]

def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reaction(r_add)
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model

def get_reaction_mw(sbml_path,project_folder,project_name,json_output_file):
        model = cobra.io.read_sbml_model(sbml_path)
        basepath: str = project_folder + project_name
        # READ REACTIONS<->KEGG ID XLSX
        protein_id_mass_mapping: Dict[str, float] = json_load(
            basepath + "_protein_id_mass_mapping.json")
        #subunit_num	1 and 1 and 1 
        reaction_mw={}
        for r in model.reactions:
                if re.search('_num',r.id):
                    r_id=r.id.split('_num')[0]
                else:
                    r_id=r.id            
            #if r_id in reactions_kcat_mapping_database.keys():
                #print(r.id,r.gene_reaction_rule)
                mass_sum = .0
                if re.search(' and ',r.gene_reaction_rule):
                    genelist=r.gene_reaction_rule.split(' and ')
                    for eachgene in genelist:
                        enzyme_unit_number=1
                        if eachgene in protein_id_mass_mapping.keys():
                            mass_sum += protein_id_mass_mapping[eachgene] * enzyme_unit_number
                    #print(mass_sum)
                    reaction_mw[r.id]=mass_sum
                else:  # Single enzyme
                    eachgene=r.gene_reaction_rule
                    enzyme_unit_number = 1
                    if eachgene in protein_id_mass_mapping.keys():
                        #print(protein_id_mass_mapping[eachgene] * enzyme_unit_number)
                        reaction_mw[r.id]=protein_id_mass_mapping[eachgene] * enzyme_unit_number
        json_write(json_output_file, reaction_mw) 
                        
# PUBLIC FUNCTIONS
def get_reaction_kcat_mw(model_file: str,project_folder: str, project_name: str,enzyme_unit_number_file: str,
                                       type_of_default_kcat_selection: str = "median") -> None:
    """Adds proteomic constraints according to sMOMENT to the given stoichiometric model and stores it as SBML.

    Arguments
    ----------

    * model: cobra.Model ~ A cobra Model representation of the metabolic network. This model will
      be changed using cobrapy functions in order to add the proteomic constraints.
    * project_folder: str ~ The folder in which the spreadsheets and JSONs with the model's supplemental
      data can be found.
    * project_name: str ~ The sMOMENTed model creation's name, which will be added at the beginning
      of the created SBML's name.
    * type_of_default_kcat_selection: str ~ The type of selection of default kcat values. Can be "mean",
      "median" or "random". Is "median" by default.

    Output
    ----------
    An SBML in the given folder with the given name, which describes the given stoichiometric model
    enhanced by the protein constraint introduction with this function.
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # Set folder path for newly created SBML and name for the reaction ID addition (added at the end,
    # and used in order to have a programatically convinient way to separate additions such as 'reverse'
    # from the 'actual' reaction ID).
    basepath: str = project_folder + project_name
    id_addition: str = "_num"
    # READ REACTIONS<->KEGG ID XLSX
    protein_id_mass_mapping: Dict[str, float] = json_load(
        basepath + "_protein_id_mass_mapping.json")
    if enzyme_unit_number_file != 'none':
        enzyme_unit_number=json_load(enzyme_unit_number_file) 
    # Make model irreversible, separating all reversible reactions to which a gene rule is given
    # in order to save some reactions.
    model = cobra.io.read_sbml_model(model_file)
    convert_to_irreversible(model)
    #split isoenzyme
    model = isoenzyme_split(model)
    
    # Read reaction <-> kcat mapping :-)
    reactions_kcat_mapping_database = json_load(
        basepath + "_reactions_kcat_mapping_combined.json")

    # sMOMENT :D
    # Get all kcats which are not math.nan and calculate the median of them, which will be used as default kcat
    all_kcats = [x["forward"] for x in reactions_kcat_mapping_database.values()] + \
                [x["reverse"] for x in reactions_kcat_mapping_database.values()]
    all_kcats = [x for x in all_kcats if not math.isnan(x)]

    if type_of_default_kcat_selection == "median":
        default_kcat = statistics.median(all_kcats)
    elif type_of_default_kcat_selection == "mean":
        default_kcat = statistics.mean(all_kcats)
    elif type_of_default_kcat_selection == "max":
        default_kcat = np.max(all_kcats)
    elif type_of_default_kcat_selection == "random":
        default_kcat = random.choice(all_kcats)
    else:
        default_kcat = 'Null'

    print(f"Default kcat is: {default_kcat}")

    # Get all reaction IDs of the given model
    model_reaction_ids = [x.id for x in model.reactions]

    # Main loop :D, add enzyme constraints to reactions \o/
    reaction_kcat_mw={}
    for model_reaction_id in model_reaction_ids:
        # Get the reaction and split the ID at the ID addition
        reaction = model.reactions.get_by_id(model_reaction_id)
        splitted_id = reaction.id.split(id_addition)

        # If the reaction has no name, ignore it
        if splitted_id[0] == "":
            continue
        # Take the reaction ID from the first part of the split
        reaction_id = splitted_id[0]
        # Remove GPRSPLIT name addition from reactions with measured protein concentrations
        if "_GPRSPLIT_" in reaction_id:
            reaction_id = reaction_id.split("_GPRSPLIT_")[0]

        # Retrieve the reaction's forward and reverse kcats from the given reaction<->kcat database
        if re.search('_reverse',reaction_id):
            reaction_id=reaction_id.split('_reverse')[0]
        if reaction_id not in reactions_kcat_mapping_database.keys():
            continue
        forward_kcat = reactions_kcat_mapping_database[reaction_id]["forward"]
        reverse_kcat = reactions_kcat_mapping_database[reaction_id]["reverse"]

        # If the given reaction<->kcat database contains math.nan as the reaction's kcat,
        # set the default kcat as math.nan means that no kcat could be found.
        if math.isnan(forward_kcat):
            forward_kcat = default_kcat
        if math.isnan(reverse_kcat):
            reverse_kcat = default_kcat

        # Add the given forward or reverse kcat is the reaction was
        # splitted due to its reversibility.
        # If the reaction is not splitted, add the forward kcat (this
        # is the only possible direction for non-splitted=non-reversible
        # reactions)
        if model_reaction_id.endswith(id_addition + "forward"):
            reaction_kcat = forward_kcat
        elif model_reaction_id.endswith(id_addition + "reverse"):
            reaction_kcat = reverse_kcat
        else:
            reaction_kcat = forward_kcat

        reaction_kcat_mw[model_reaction_id]={}
        if reaction_kcat=='Null':
            continue
        reaction_kcat_mw[model_reaction_id]['kcat']=reaction_kcat
        
        #MW
        #subunit_num	1 and 1 and 1 
        reaction_mw={}
        for r in model.reactions:
                #print(r.id)
                if re.search('_num',r.id):
                    r_id=r.id.split('_num')[0]
                else:
                    r_id=r.id            
            #if r_id in reactions_kcat_mapping_database.keys():
                #print(r.id,r.gene_reaction_rule)
                mass_sum = .0
                if re.search(' and ',r.gene_reaction_rule):
                    genelist=r.gene_reaction_rule.split(' and ')
                    for eachgene in genelist:
                        if eachgene in protein_id_mass_mapping.keys():
                            if enzyme_unit_number_file != 'none':
                                if r.id in enzyme_unit_number.keys():
                                    if eachgene in enzyme_unit_number[r.id].keys():
                                        mass_sum += protein_id_mass_mapping[eachgene] * int(enzyme_unit_number[r.id][eachgene])
                            else:
                                #print(eachgene)
                                mass_sum += protein_id_mass_mapping[eachgene]
                    #print(mass_sum)
                    if mass_sum>0:
                        reaction_mw[r.id]=mass_sum
                else:  # Single enzyme
                    eachgene=r.gene_reaction_rule
                    #enzyme_unit_number = 1
                    if eachgene in protein_id_mass_mapping.keys():
                        #print(protein_id_mass_mapping[eachgene] * enzyme_unit_number)
                        #reaction_mw[r.id]=protein_id_mass_mapping[eachgene]['mw'] * enzyme_unit_number
                        if enzyme_unit_number_file != 'none':
                            if r.id in enzyme_unit_number.keys():
                                if eachgene in enzyme_unit_number[r.id].keys(): 
                                    reaction_mw[r.id]=protein_id_mass_mapping[eachgene] * int(enzyme_unit_number[r.id][eachgene])
                        else:
                            reaction_mw[r.id]=protein_id_mass_mapping[eachgene]
                # isoform enzyme's mw should be modified
                if re.search('ASPK',r.id):
                    reaction_mw[r.id] = 126584
                if re.search('ACGS',r.id):
                    reaction_mw[r.id] = 79428
                if re.search('ASP1DC',r.id):
                    reaction_mw[r.id] = 56588   
        #print(model_reaction_id,reaction_mw.keys())
        
        if model_reaction_id in reaction_mw.keys():
            # print(reaction_mw.keys())
            print(model_reaction_id,reaction_mw[model_reaction_id])
            reaction_kcat_mw[model_reaction_id]['MW']=reaction_mw[model_reaction_id]
            
            reaction_kcat_mw[model_reaction_id]['kcat_MW']=reaction_kcat_mw[model_reaction_id]['kcat']*3600000/reaction_mw[model_reaction_id]
        
    reaction_kcat_mw_df = pd.DataFrame(reaction_kcat_mw)
    reaction_kcat_mw_df_T=reaction_kcat_mw_df.T
    reaction_kcat_mw_df_T_select=reaction_kcat_mw_df_T[abs(reaction_kcat_mw_df_T['kcat_MW'])>0]
    reaction_kcat_mw_df_T_select.to_csv(project_folder + 'reaction_kcat_MW.csv')
    
def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reaction(r_add)
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model

def trans_model2enz_json_model_split_isoenzyme(model_file, reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Tansform cobra model to json mode with  
    enzyme concentration constraintat.

    Arguments
    ----------
    * model_file:   The path of sbml model
    * reaction_kcat_mw_file: The path of storing kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint. 
    将眼镜蛇模型转换为json模式，其中包括  酶浓度的约束。
    参数
    ----------
    * model_file。  sbml模型的路径
    * reaction_kcat_mw_file。存储GEM模型中催化每个反应的酶的kcat/MW值的路径。
     反应的kcat/MW值。
    * f: 酶的质量分数 
    * ptot: 细胞中的总蛋白质部分。 
    * sigma: 酶的近似平均饱和度。
    * lowerbound:  酶浓度约束的下限。
    * upperbound:  酶浓度限制的上限。
    * 输出一个酶约束模型


    """
    model = cobra.io.read_sbml_model(model_file)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        #if re.search('_num',reaction_id):
        #    reaction_id=reaction_id.split('_num')[0]
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    json_write(json_output_file, dictionary_model)
    
def get_enzyme_constraint_model(json_model_file):#将模型读取成酶约束模型
    """using enzyme concentration constraint
    json model to create a COBRApy model.

    Arguments
    ----------
    * json_model_file: json Model file.

    :return: Construct an enzyme-constrained model.
    """

    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)

    coefficients = dict()
    for rxn in model.reactions:
        for eachr in dictionary_model['reactions']:
            if rxn.id == eachr['id']:
                if eachr['kcat_MW']:
                    coefficients[rxn.forward_variable] = 1 / float(eachr['kcat_MW'])
                break

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']
    print(upperbound)
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    model.reactions.get_by_id("R1975_num1").bounds=(2.8,2.8)
    model.reactions.get_by_id("R1975_num2").bounds=(0,0)
    model.reactions.get_by_id("R1975_num3").bounds=(0,0)
    model.reactions.get_by_id("R1975_num4").bounds=(0,0)
    model.reactions.get_by_id("R1975_num5").bounds=(0,0)
    return model

def get_enzyme_constraint_model_percent(json_model_file,percent):
    """using enzyme concentration constraint
    json model to create a COBRApy model.

    Arguments
    ----------
    * json_model_file: json Model file.
    使用酶浓度约束
    json模型来创建一个COBRApy模型。

    论据
    ----------
    * json_model_file: json模型文件。

    :返回。构建一个受酶约束的模型。
    :return: Construct an enzyme-constrained model.
    """

    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)

    coefficients = dict()
    for rxn in model.reactions:
        for eachr in dictionary_model['reactions']:
            if rxn.id == eachr['id']:
                if eachr['kcat_MW']:
                    coefficients[rxn.forward_variable] = 1 / float(eachr['kcat_MW'])
                break

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']*percent
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model

def get_fluxes_detail_in_model(model,model_pfba_solution,fluxes_outfile,json_model_file):
    """Get the detailed information of each reaction

    Arguments
    ----------
    * model: cobra.Model.
    * fluxes_outfile: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.

    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    ""获得每个反应的详细信息

    论据
    ----------
    * model: cobra.Model。
    * fluxes_outfile：反应通量文件。
    * reaction_kcat_mw_file：反应kcat/mw文件。

    :返回：数据帧中的通量、kcat、MW和kcat_MW。
    """
    dictionary_model = json_load(json_model_file)
    model_pfba_solution = model_pfba_solution.to_frame()
    model_pfba_solution_detail = pd.DataFrame()
    for index, row in model_pfba_solution.iterrows():
        reaction_detail = model.reactions.get_by_id(index)
        model_pfba_solution_detail.loc[index, 'fluxes'] = row['fluxes']
        for eachreaction in dictionary_model['reactions']:
            if index ==eachreaction['id']:
                if 'annotation' in eachreaction.keys():
                    if 'ec-code' in eachreaction['annotation'].keys():
                        if isinstance (eachreaction['annotation']['ec-code'],list):
                            model_pfba_solution_detail.loc[index, 'ec-code'] = (',').join(eachreaction['annotation']['ec-code'])
                        else:
                            model_pfba_solution_detail.loc[index, 'ec-code'] = eachreaction['annotation']['ec-code']    
                if 'kcat_MW' in eachreaction.keys():
                    if eachreaction['kcat_MW']:
                        model_pfba_solution_detail.loc[index, 'kcat_MW'] = eachreaction['kcat_MW']
                        model_pfba_solution_detail.loc[index, 'E'] = float(row['fluxes'])/float(eachreaction['kcat_MW'])
                break
        model_pfba_solution_detail.loc[index, 'equ'] = reaction_detail.reaction
    print('Enzyme cost total is:'+str(np.sum(model_pfba_solution_detail['E'])))
    model_pfba_solution_detail.to_csv(fluxes_outfile)
    return model_pfba_solution_detail


def json_load(path):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary


def json_write(path, dictionary):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)
        
def GENENAME_2_ACC_from_uniprot(query,outfile):
    #print(' '.join(query).replace('511145.',''))
    url = 'https://legacy.uniprot.org/uploadlists/'
    params = {
        'from': 'GENENAME',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(query),
        'columns':'id,entry name,protein names,genes,organism,ec,mass,database(PDB)'
    }
    data = urlencode(params).encode()
    request = Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urlopen(request)
    page = response.read()
    # print(page)
    outFile = open(outfile,'w') 
    namesRegex = re.compile(r'yourlist:(.*)\n')
    outFile.write(namesRegex.sub('Gene ID\n',page.decode('utf-8')))
    #print(namesRegex.sub('Protein AC\t',page.decode('utf-8')))
    outFile.close()
    
def GENENAME_2_ACC_from_uniprot_byID(query,outfile):
    #print(' '.join(query).replace('511145.',''))
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from': 'ID',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(query),
        'columns':'id,entry name,protein names,genes,organism,ec,mass,database(PDB)'
    }
    data = urlencode(params).encode()
    request = Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urlopen(request)
    page = response.read()
    outFile = open(outfile,'w') 
    namesRegex = re.compile(r'yourlist:(.*)\n')
    outFile.write(namesRegex.sub('Gene ID\n',page.decode('utf-8')))
    #print(namesRegex.sub('Protein AC\t',page.decode('utf-8')))
    outFile.close()    
    
def calculate_f(uni_model_gene_list, gene_abundance_file, gene_mw_file, gene_mw_colname,gene_abundance_colname):
    """Calculating f (the mass fraction of enzymes that are accounted
    in the model out of all proteins) based on the protein abundance
    which can be obtained from PAXdb database.

    Arguments
    ----------
    * genes: All the genes in the model.
    * gene_abundance_file: The protein abundance of each gene
     in the E. coli genome.
    * gene_mw_file: The molecular weight of the protein expressed by each gene.

    :return: The enzyme mass fraction f.
    """
    gene_abundance = pd.read_csv(gene_abundance_file, index_col=0)
    gene_mw = pd.read_csv(gene_mw_file, sep='\t', index_col=gene_mw_colname)
    enzy_abundance = 0
    pro_abundance = 0
    for gene_i in gene_abundance.index:
        if gene_i in gene_mw.index:
            if isinstance(gene_mw.loc[gene_i,'Mass'],str):
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * int(gene_mw.loc[gene_i,'Mass'].replace(',',''))
            else:
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * int(gene_mw.loc[gene_i,'Mass'][0].replace(',',''))
            pro_abundance += abundance
            if gene_i in uni_model_gene_list:
                enzy_abundance += abundance
    f = enzy_abundance/pro_abundance
    return f

def calculate_f_special(gene_abundance_file,modelgene2uniprot_file,paxdbgene2uniprot_file,gene_abundance_colname):
    """Calculating f (the mass fraction of enzymes that are accounted
    in the model out of all proteins) based on the protein abundance
    which can be obtained from PAXdb database.
    根据蛋白质丰度，计算f（在模型中占所有蛋白质的酶的质量分数
    在模型中占所有蛋白质的质量分数），根据蛋白质丰度
    可以从PAXdb数据库获得。
    Arguments
    ----------
    * genes: All the genes in the model.
    * gene_abundance_file: The protein abundance of each gene
     in the E. coli genome.
    * gene_mw_file: The molecular weight of the protein expressed by each gene.

    :return: The enzyme mass fraction f.
    * 基因。模型中的所有基因。
    * gene_abundance_file。每个基因的蛋白质丰度
     在大肠杆菌基因组中的蛋白质丰度。
    * gene_mw_file。每个基因所表达的蛋白质的分子量。

    :返回。酶的质量分数f。
    """
    gene_abundance = pd.read_csv(gene_abundance_file, index_col='gene_id',sep='\t')
    model_gene = pd.read_csv(modelgene2uniprot_file, sep='\t', index_col='Gene ID')
    paxdb_gene = pd.read_csv(paxdbgene2uniprot_file, sep='\t', index_col='Gene ID')
    paxdb_gene=paxdb_gene[paxdb_gene['Organism']=='Myceliophthora thermophila (strain ATCC 42464 / BCRC 31852 / DSM 1799) (Sporotrichum thermophile)']
    model_gene=model_gene[model_gene['Organism']=='Myceliophthora thermophila (strain ATCC 42464 / BCRC 31852 / DSM 1799) (Sporotrichum thermophile)']
    enzy_abundance = 0
    pro_abundance = 0
    for gene_i in gene_abundance.index:
        # print(gene_i)
        if gene_i in list(paxdb_gene.index):
            # print(gene_i)
            if str(paxdb_gene.loc[gene_i,'Mass']) !='nan':
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * float(paxdb_gene.loc[gene_i,'Mass'].replace(',',''))
                # print(abundance)
                pro_abundance += abundance
                if gene_i in list(model_gene.index):
                    # print(gene_i)
                    enzy_abundance += abundance
    # print(enzy_abundance,pro_abundance)
    f = enzy_abundance/pro_abundance
    return f

def calculate_f_special_new(model,gene_abundance_file):
    gene_abundance = pd.read_csv(gene_abundance_file)
    enzy_abundance = 0
    pro_abundance = 0
    model_genes = [gene.id for gene in model.genes]
    model_genes
    for gene_i in range(len(gene_abundance['gene_id'].index)):
            # print(gene_i)
                if str(gene_abundance.loc[gene_i,'Mass']) !='nan':
                    abundance = gene_abundance.loc[gene_i, 'normalized_counts'] * float(gene_abundance.loc[gene_i,'Mass'])
                    # print(abundance)
                    pro_abundance += abundance
#                     print(pro_abundance)
                if gene_abundance.loc[gene_i, 'gene_id'] in model_genes:
                    # print(gene_i)
                    enzy_abundance += abundance
    # print(enzy_abundance,pro_abundance)
    f = enzy_abundance/pro_abundance
    return f

def get_model_substrate_obj(use_model):
    ATPM='No' 
    substrate_list=[]
    concentration_list=[]
    #这些应当是底物的信息
    EX_exclude_reaction_list=['R2147','R2150','R2152','R2155','R2391','R2392','R2393','R2394','R2395','R2396','R2397','R2398','R2400','R2401','R2402','R2403','R2404','R2405','R2406','R2407','R2408','R2409','R2410','R2411','R2412','R2413','R2414','R2415','R2416','R2417','R2418','R2419','R2420','R2421','R2422','R2423','R2424','R2425','R2426','R2427','R2428','R2429','R2430','R2431','R2432','R2433','R2434','R2435','R2436','R2437','R2438','R2439','R2440','R2441','R2442','R2443','R2444','R2445','R2446','R2447','R2448','R2449','R2450','R2451','R2452','R2453','R2454','R2455','R2456','R2457','R2458','R2459','R2460','R2461','R2462','R2463','R2464','R2465','R2466','R2467','R2468','R2469','R2470','R2471','R2472','R2473','R2474','R2475','R2476','R2477','R2478','R2479','R2480','R2481','R2482','R2483','R2484','R2485','R2486','R2487','R2488','R2489','R2490','R2491','R2492','R2493','R2494','R2495','R2496','R2497','R2498','R2499','R2500','R2501','R2502','R2503','R2504','R2505','R2506','R2507','R2508','R2509','R2510','R2511','R2512','R2513','R2514','R2515','R2516','R2517','R2518','R2519','R2520','R2521','R2522','R2523','R2524','R2525','R2526','R2539','R2542','R2545','R2548','R2551','R2558','R2568','R2583','R2593','R2597','R2599']#不要的反应
    EX_exclude_reaction_list=EX_exclude_reaction_list+[i+'_reverse' for i in EX_exclude_reaction_list]
    for r in use_model.reactions:
        if r.objective_coefficient == 1:
            obj=r.id #Product name
        #elif not r.lower_bound==0 and not r.lower_bound==-1000 and not r.lower_bound==-999999 and abs(r.lower_bound)>0.1:#排除很小的值
        elif not r.upper_bound==0 and not r.upper_bound==1000 and not r.upper_bound==999999 and abs(r.upper_bound)>0.1:#排除很小的值
            #print(r.id,r.upper_bound,r.lower_bound)
            #判断反应中的ATPM是否复合要求
            if r.id=='R1975_num1'or r.id=='R1975_num2'or r.id=='R1975_num3'or r.id=='R1975_num4'or r.id=='R1975_num5':
                if r.upper_bound>0:
                    ATPM='Yes' #ATP maintenance requirement
            #获取反应中的底物信息
            elif r.id not in EX_exclude_reaction_list:
                #print(r.id,r.upper_bound,r.lower_bound)
                #substrate=r.id #Substrate name
                substrate_list.append(r.id)
                #concentration=r.upper_bound #Substrate uptake rate  
                concentration_list.append(r.upper_bound)
    return(obj,substrate_list,concentration_list,ATPM)

def parse_sabio_rk_for_eclist(ec_numbers_list: List[str], json_output_path: str, bigg_id_name_mapping_path: str) -> None:
    """Retrieves kcats from SABIO-RK for the given model and stores it in a JSON for the given model in the given path.

    Algorithm
    ----------
    Using the SABIO-RK REST API (as of 2019/30/04, it is explained under
    http://sabiork.h-its.org/layouts/content/docuRESTfulWeb/RESTWebserviceIntro.gsp),


    Arguments
    ----------
    * eclist: List[str] ~ eclist.
    * json_output_path: str ~ The path of the JSON that shall be created

    Output
    ----------
    * A JSON in the given project folder with the following structure:
    <pre>
        {
            "$EC_NUMBER_OR_KEGG_REACTION_ID": {
                "$SUBSTRATE_WITH_BIGG_ID_1": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                },
                (...),
                "REST": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                }
            }
            (...),
        }
    </pre>
    'REST' stands for a substrate without found BIGG ID.
    """
    # GET KCATS FOR EC NUMBERS
    ec_number_kcat_mapping = get_ec_number_kcats_wildcard_search(
        ec_numbers_list, bigg_id_name_mapping_path)

    json_write(json_output_path, ec_number_kcat_mapping)
    
    
def get_protein_mass_mapping_from_local(sbml_path: str, project_folder: str, project_name: str,uniprot_data_file: str) -> None:
    """Returns a JSON with a mapping of protein IDs as keys, and as values the protein mass in kDa.

    The protein masses are calculated using the amino acid sequence from UniProt (retrieved using
    UniProt's REST API).

    Arguments
    ----------
    * model: cobra.Model ~ The model in the cobrapy format
    * project_folder: str ~ The folder in which the JSON shall be created
    * project_name: str ~ The beginning of the JSON's file name
    * uniprot_data_file: str ~ The gene information obtained from uniprot
    Output
    ----------
    A JSON file with the path project_folder+project_name+'_protein_id_mass_mapping.json'
    and the following structure:
    <pre>
    {
        "$PROTEIN_ID": $PROTEIN_MASS_IN_KDA,
        (...),
    }
    </pre>
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # The beginning of the created JSON's path :D
    basepath: str = project_folder + project_name
    model= cobra.io.read_sbml_model(sbml_path)
    # GET UNIPROT ID - PROTEIN MAPPING
    uniprot_id_protein_id_mapping: Dict[str, List[str]] = {}
    for gene in model.genes:
        # Without a UniProt ID, no mass mapping can be found
        if "uniprot" not in gene.annotation:
            continue
        uniprot_id = gene.annotation["uniprot"]
        if uniprot_id in uniprot_id_protein_id_mapping.keys():
            uniprot_id_protein_id_mapping[uniprot_id].append(gene.id)
        else:
            uniprot_id_protein_id_mapping[uniprot_id] = [gene.id]

    # GET UNIPROT ID<->PROTEIN MASS MAPPING
    uniprot_id_protein_mass_mapping = json_load(uniprot_data_file)
    
    # Create the final protein ID <-> mass mapping
    protein_id_mass_mapping: Dict[str, float] = {}
    for uniprot_id in list(uniprot_id_protein_mass_mapping.keys()):
        try:
            protein_ids = uniprot_id_protein_id_mapping[uniprot_id]
        except Exception:
            #print(f"No mass found for {uniprot_id}!")
            continue
        for protein_id in protein_ids:
            protein_id_mass_mapping[protein_id] = uniprot_id_protein_mass_mapping[uniprot_id]

    # Write protein mass list JSON :D
    #print("Protein ID<->Mass mapping done!")
    json_write(basepath+"_protein_id_mass_mapping.json", protein_id_mass_mapping)
    
def change_enz_model_by_enz_usage(enz_ratio,json_model_path, reaction_flux_file, EC_max_file, reaction_kcat_mw, need_change_reaction_list, changed_reaction_list,f, ptot, sigma, lowerbound, upperbound, json_output_file):

    """Get new enzyme model using enzyme usage to calibration

    Arguments
    ----------
    * enz_ratio: enzyme ratio which needed change.
    * json_model_path: The file storing json model.
    * reaction_flux_file: reaction-flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file： enzyme usage of each reaction.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_mw_outfile: changed file stored reaction kcat/mw.

    :return: new enzyme model
    使用酶的使用情况获得新的酶模型以进行校准

    参数
    ----------
    * enz_ratio: 需要改变的酵素比率。
    * json_model_path: 存储json模型的文件。
    * reaction_flux_file: 反应流量的文件。
    * reaction_kcat_mw_file: 反应kcat/mw文件。
    * reaction_enz_usage_file：每个反应的酶用量。
    * kcat_database_combined_file：从autoPACMEN得到的合并kcat文件。
    * model_file：眼镜蛇模型。
    * f: 酶的质量分数 
    * ptot: 细胞中的总蛋白质部分。 
    * sigma: 酶的近似平均饱和度。
    * lowerbound:  酶浓度约束的下限。
    * upperbound:  酶浓度约束的上限值。 
    * json_output_file: json文件存储json模型。
    * reaction_mw_outfile: 更改后的文件，存储反应kcat/mw。

    :返回：新的酶模型
    """ 
    reaction_fluxes = pd.read_csv(reaction_flux_file, index_col=0)
    reaction_fluxes['enz ratio'] = reaction_fluxes['E']/np.sum(reaction_fluxes['E'])
    reaction_fluxes=reaction_fluxes.sort_values(by="enz ratio", axis=0, ascending=False)
    #每次只选第一个
    i=0
    select_reaction = reaction_fluxes.index[0]
    # if select_reaction == 'PPNDH':
    #     select_reaction = reaction_fluxes.index[1]
    while (select_reaction in need_change_reaction_list):
        i=i+1
        #print(i)
        select_reaction = reaction_fluxes.index[i+1]
        
    print('Need changing reaction: ')
    print(select_reaction)
    [need_change_reaction_list, changed_reaction_list,reaction_kcat_mw] = adj_reaction_kcat_by_database(json_model_path,[select_reaction], need_change_reaction_list, changed_reaction_list, EC_max_file, reaction_kcat_mw)
    print('Changed reaction: ')
    print(changed_reaction_list)

    adj_trans_model2enz_model(json_model_path, reaction_kcat_mw, f, ptot, sigma, lowerbound, upperbound, json_output_file)

    enz_model = get_enzyme_constraint_model(json_output_file)
    print('Enzyme cost total is:'+str(np.sum(reaction_fluxes['E'])))
    return (enz_model,reaction_kcat_mw,need_change_reaction_list, changed_reaction_list)


def adj_reaction_kcat_by_database_old(json_model_path,select_reactionlist, need_change_reaction_list, changed_reaction_list,EC_max_file, reaction_kcat_mw):
    """Use the kcat in database to change reaction kcat in model

    Arguments
    ----------
    * json_model_path: The file storing json model.
    * select_reactionlist: reaction list need to change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_kapp_change_file: changed file stored reaction kcat/mw.

    :return: a dataframe stored new reaction kcat/mw .
    """
    Brenda_sabio_combined_select = json_load(EC_max_file)
    
    json_model=cobra.io.load_json_model(json_model_path)
    for eachreaction in select_reactionlist:
        need_change_reaction_list.append(eachreaction)
        select_reaction = json_model.reactions.get_by_id(eachreaction)
        if "ec-code" in select_reaction.annotation.keys():
            ec_number = select_reaction.annotation["ec-code"]
            kcat_max_list = []
            if isinstance(ec_number, str):
                if ec_number in Brenda_sabio_combined_select.keys():
                    reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat_max']
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max * 3600:
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max * 3600#h_1
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                        changed_reaction_list.append(eachreaction) 
            else:
                for eachec in ec_number:
                    if eachec in Brenda_sabio_combined_select.keys():
                        kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat_max'])
                reaction_kcat_max = np.max(kcat_max_list)     
                if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max * 3600:
                    reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max * 3600
                    reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                    changed_reaction_list.append(eachreaction)    
                
    return(need_change_reaction_list,changed_reaction_list,reaction_kcat_mw)

def adj_reaction_kcat_by_database(json_model_path,select_reactionlist, need_change_reaction_list, changed_reaction_list,EC_max_file, reaction_kcat_mw):
    """Use the kcat in database to change reaction kcat in model

    Arguments
    ----------
    * json_model_path: The file storing json model.
    * select_reactionlist: reaction list need to change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_kapp_change_file: changed file stored reaction kcat/mw.

    :return: a dataframe stored new reaction kcat/mw .
    """
    Brenda_sabio_combined_select = json_load(EC_max_file)
    
    json_model=cobra.io.load_json_model(json_model_path)
    for eachreaction in select_reactionlist:
        need_change_reaction_list.append(eachreaction)
        select_reaction = json_model.reactions.get_by_id(eachreaction)
        if "ec-code" in select_reaction.annotation.keys():
            ec_number = select_reaction.annotation["ec-code"]
            kcat_max_list = []
            if isinstance(ec_number, str):
                if ec_number in Brenda_sabio_combined_select.keys():
                    reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat_max']
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max:
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max#h_1
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                        changed_reaction_list.append(eachreaction) 
            else:
                for eachec in ec_number:
                    if eachec in Brenda_sabio_combined_select.keys():
                        kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat_max'])
                reaction_kcat_max = np.max(kcat_max_list)     
                if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max :
                    reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max 
                    reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                    changed_reaction_list.append(eachreaction)    
                
    return(need_change_reaction_list,changed_reaction_list,reaction_kcat_mw)

def adj_trans_model2enz_model(model_file, reaction_kcat_mw, f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Tansform cobra model to json mode with  
    enzyme concentration constraintat.

    Arguments
    ----------
    * model_file:   The path of sbml model
    * reaction_kcat_mw_file: The path of storing kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint. 

    """
    if re.search('\.xml',model_file):
        model = cobra.io.read_sbml_model(model_file)
    elif re.search('\.json',model_file):
        model = cobra.io.json.load_json_model(model_file)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        #if re.search('_num',reaction_id):
        #    reaction_id=reaction_id.split('_num')[0]
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    json_write(json_output_file, dictionary_model)

    
    
def draw_cdf_fig(data_cdf_data,output_file,x_name,y_name,y_index,nticks):
    trace0 = go.Scatter(x=data_cdf_data,y=y_index,mode='lines',marker={'color': 'blue'})
    trace1 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines',xaxis='x2',yaxis="y2")
    data1 = [trace0, trace1]
    layout = go.Layout(plot_bgcolor='lightgrey',
            xaxis=dict(title=dict(text=x_name,font=dict(size=20, family='Times New Roman')),
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   linecolor='black',ticks='inside',tickcolor='black',zeroline=False,
                   showexponent = 'all',exponentformat =  "power"),
            xaxis2=dict(linecolor='black',showticklabels=False,type="log",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   rangemode="tozero",overlaying='x', side='top',nticks=nticks, zeroline=False,
                   showexponent = 'all', exponentformat =  "power",gridcolor="white"),
            yaxis=dict(title=dict(text=y_name,font=dict(size=20, family='Times New Roman')),range=[0, 1],
                   showgrid=False,zeroline=False,rangemode="tozero", tickfont=dict(color='black', size=20, family='Times New Roman'),
                   ticks='inside',tickcolor='black',linecolor='black'),
            yaxis2=dict(range=[0, 1],linecolor='black',showgrid=False,zeroline=False,tickfont=dict(color='black', size=20, family='Times New Roman'),
                    showticklabels=False,overlaying='y',side='right'),
            showlegend=False,height=450,width=750,margin=go.layout.Margin(l=10, r=10, b=10, t=10))

    fig = go.Figure(data1, layout=layout)
    fig.add_hline(y=0.5,line_width=2,line_color="orange")
    fig.write_image(output_file)
    return fig

def draw_cdf_fig_kcat(data_cdf_data,output_file,x_name,y_name,y_index,nticks):
    trace0 = go.Scatter(x=data_cdf_data,y=y_index,mode='lines',marker={'color': 'blue'})
    trace1 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines',line={'color': 'blue', 'width': 3},xaxis='x2',yaxis="y2")
    data1 = [trace0, trace1]
    layout = go.Layout(plot_bgcolor='lightgrey',xaxis=dict(title=dict(text=x_name,font=dict(size=20, family='Times New Roman')),
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   linecolor='black',ticks='inside',tickcolor='black',zeroline=False,
                   showexponent = 'all',exponentformat =  "power", gridcolor="yellow"),
            xaxis2=dict(linecolor='black',showticklabels=False,type="log",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   rangemode="tozero",overlaying='x', side='top',nticks=nticks, zeroline=False,
                   showexponent = 'all', exponentformat =  "power",gridcolor="white"),
            yaxis=dict(title=dict(text=y_name,font=dict(size=20, family='Times New Roman')),range=[0, 1],
                   showgrid=False,zeroline=False,rangemode="tozero", tickfont=dict(color='black', size=20, family='Times New Roman'),
                   ticks='inside',tickcolor='black',linecolor='black'),
            yaxis2=dict(range=[0, 1],linecolor='black',showgrid=False,zeroline=False,tickfont=dict(color='black', size=20, family='Times New Roman'),
                    showticklabels=False,overlaying='y',side='right'),
            showlegend=False,height=450,width=750,margin=go.layout.Margin(l=10, r=10, b=10, t=10))

    fig = go.Figure(data1, layout=layout)
    fig.add_hline(y=0.5,line_width=2,line_color="orange")
    fig.write_image(output_file)
    return fig
def get_PhPP_dataold(model_file,model_type,obj,number,outputfile):
    if model_type=='GEM':
        cg_model = cobra.io.json.load_json_model(model_file)
    else:
        cg_model=get_enzyme_constraint_model(model_file)
    objlist = []
    glclist = []
    o2list = []
    exlist = list(np.linspace(0,30,number))
    print(exlist)
    exlistn = []
    for i in exlist:
        i = format(i,'.4f')
        exlistn.append(float(i))
    exlistm=deepcopy(exlistn)
    exlistn.insert(0,0.0)
    ectest = pd.DataFrame(exlistn)
    o2df = pd.DataFrame(exlistn).T
    df = pd.concat([o2df,ectest],axis=0)
    df = df.reset_index(drop = True)
    df1 = df.drop(1,axis=0)
    df1 = df1.reset_index(drop = True)
    df1
    k=1
    v=1
    for i in exlistm:
        condi = -i
        for j in exlistm:
            condj = -j

            cg_model.reactions.get_by_id('EX_o2_e').bounds=(condj,0)
            cg_model.reactions.get_by_id('EX_o2_e_reverse').bounds=(0,0)
            cg_model.reactions.get_by_id('EX_glc_e').bounds=(condi,0)
            cg_model.reactions.get_by_id('EX_glc_e_reverse').bounds=(0,0)

            cg_model.objective=obj
            enz_model_pfba_solution = cobra.flux_analysis.pfba(cg_model)
            # print(k,v)
            df1.iloc[k,v] = enz_model_pfba_solution.fluxes[obj]
            k = k+1
            if k == number+1:
                k=1
                v = v+1
        if v == number+1:
            break
#     print(df1)
    df1.to_csv(outputfile)
    return(df1)

def get_PhPP_data(model_file,model_type,obj,number,outputfile):
    if model_type=='GEM':
        cg_model = cobra.io.json.load_json_model(model_file)
    else:
        cg_model=get_enzyme_constraint_model(model_file)
    objlist = []
    glclist = []
    o2list = []
    exlisto2 = []
    exlist = list(np.linspace(0,10,number))
    exlisto2 = list(np.linspace(0,20,number))
    print(exlist)

    exlistn = []
    exlistno2 = []
    for i in exlist:
        i = format(i,'.4f')
        exlistn.append(float(i))
    exlistm=deepcopy(exlistn)
    for i in exlisto2:
        i = format(i,'.4f')
        exlistno2.append(float(i))
    exlistmo2=deepcopy(exlistno2)    
    
    exlistn.insert(0,0.0)
    exlistno2.insert(0,0.0)

    ectest = pd.DataFrame(exlistn)
    o2df = pd.DataFrame(exlistno2).T
    df = pd.concat([o2df,ectest],axis=0)
    df = df.reset_index(drop = True)
    df1 = df.drop(1,axis=0)
    df1 = df1.reset_index(drop = True)
    df1
    k=1
    v=1
    for i in exlistmo2:
        condi = -i
        for j in exlistm:
            condj = -j

            cg_model.reactions.get_by_id('R2396').bounds=(condi,0)#氧气
            cg_model.reactions.get_by_id('R2396_reverse').bounds=(0,0)
            cg_model.reactions.get_by_id('R2399').bounds=(condj,0)#葡萄糖
            cg_model.reactions.get_by_id('R2399_reverse').bounds=(0,0)

            cg_model.objective=obj
            enz_model_pfba_solution = cobra.flux_analysis.pfba(cg_model)
            print(k,v)
            df1.iloc[k,v] = enz_model_pfba_solution.fluxes[obj]
            k = k+1
            if k == number+1:
                k=1
                v = v+1
        if v == number+1:
            break
#     print(df1)
    df1.to_csv(outputfile)
    return(df1)

def draw_cdf_fig_mw(data_cdf_data,output_file,x_name,y_name,y_index,nticks):
    trace0 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines',line={'color': 'blue', 'width': 3},xaxis='x2',yaxis="y2")
    trace1 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines',line={'color': 'blue', 'width': 3},xaxis='x2',yaxis="y2")
    data1 = [trace0,trace1]
    layout = go.Layout(plot_bgcolor='lightgrey',
            xaxis=dict(title=dict(text=x_name,font=dict(size=20, family='Times New Roman')),
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   linecolor='black',ticks='inside',tickcolor='black',zeroline=False,showgrid=False,
                   showexponent = 'all',exponentformat =  "power"),
            xaxis2=dict(showticklabels=False,
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   linecolor='black',overlaying='x', side='top',tickcolor='black', zeroline=False,nticks = nticks,
                   showexponent = 'all', exponentformat =  "power"),
            yaxis=dict(title=dict(text=y_name,font=dict(size=20, family='Times New Roman')),range=[0, 1],
                   showgrid=False,zeroline=False,rangemode="tozero", tickfont=dict(color='black', size=20, family='Times New Roman'),
                   ticks='inside',tickcolor='black',linecolor='black'),
            yaxis2=dict(range=[0, 1],linecolor='black',showgrid=False,zeroline=False,tickfont=dict(color='black', size=20, family='Times New Roman'),
                    showticklabels=False,overlaying='y',side='right'),
            showlegend=False,height=450,width=750,margin=go.layout.Margin(l=10, r=10, b=10, t=10))

    fig = go.Figure(data1, layout=layout)
    fig.add_hline(y=0.5,line_width=2,line_color="orange")
    fig.write_image(output_file)
    return fig

def get_robustness_data(model_file,model_type,obj,number,outputfile):
    if model_type=='GEM':
        cg_model = cobra.io.json.load_json_model(model_file)
    else:
        cg_model=get_enzyme_constraint_model(model_file)
    objlist = []
    glclist = []
    o2list = []
    exlist = list(np.linspace(0,30,number))
    print(exlist)
    exlistn = []
    for i in exlist:
        i = format(i,'.4f')
        exlistn.append(float(i))
    exlistm=deepcopy(exlistn)
    exlistn.insert(0,0.0)
    ectest = pd.DataFrame(exlistn)
    o2df = pd.DataFrame(exlistn).T
    df = pd.concat([o2df,ectest],axis=0)
    df = df.reset_index(drop = True)
    df1 = df.drop(1,axis=0)
    df1 = df1.reset_index(drop = True)
    df1
    k=1
    v=1
    for i in exlistm:
        condi = -i
        for j in exlistm:
            condj = -j

            cg_model.reactions.get_by_id('EX_o2_e').bounds=(condj,0)
            cg_model.reactions.get_by_id('EX_o2_e_reverse').bounds=(0,0)
            cg_model.reactions.get_by_id('EX_glc_e').bounds=(condi,0)
            cg_model.reactions.get_by_id('EX_glc_e_reverse').bounds=(0,0)

            cg_model.objective=obj
            enz_model_pfba_solution = cobra.flux_analysis.pfba(cg_model)
            # print(k,v)
            df1.iloc[k,v] = enz_model_pfba_solution.fluxes[obj]
            k = k+1
            if k == number+1:
                k=1
                v = v+1
        if v == number+1:
            break
#     print(df1)
    df1.to_csv(outputfile)
    return(df1)

def get_min_enzyme_cost(model, dict_coeff):
    """Get model flux using Minimum enzyme cost algorithm

    Arguments
    ----------
    * model: cobra model.
    * dict_coeff: {reaction ID: coeffient}.
    
    :return: cobra solution.
    """
    with model:
        bounds = (model.slim_optimize(), model.slim_optimize())
        cons_obj = model.problem.Constraint(
            model.objective.expression,
            lb=min(bounds), ub=max(bounds))
        model.add_cons_vars(cons_obj)

        dict_obj = dict()
        for r in model.reactions:
            if r.id in list(dict_coeff.index):
                #print(dict_coeff.loc[r.id,'kcat_MW'])
                dict_obj[r.forward_variable] = 1 / dict_coeff.loc[r.id,'kcat_MW']

        model_obj = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective = model_obj
        model.objective.set_linear_coefficients(dict_obj)
        solution = model.optimize()
    return solution

def draw_3d_rbas(z_data,out_fig_file):
    layout = go.Layout(template="none",plot_bgcolor='lightgrey')
    exlist = list(np.linspace(0,30,31))
    fig = go.Figure(data=[go.Surface(y=exlist,x=exlist,z=z_data.values)],layout=layout)
    fig.update_layout(scene = dict(
        xaxis = dict(backgroundcolor = "lightgrey",title=dict(text="Glucose uptake rates<br>(mmol/gDW/h)",font=dict(size=15, family='Times New Roman'),)),
        yaxis = dict(backgroundcolor = "lightgrey",title=dict(text="O2 uptake rates<br>mmol/gDW/h",font=dict(size=15, family='Times New Roman'))),
        zaxis = dict(range=[0,1], backgroundcolor = "grey", gridcolor = "white", title=dict(text="Growth rates (1/h)",font=dict(size=15, family='Times New Roman')))))

    fig.update_traces(contours_z=dict(usecolormap=True, highlightcolor="mistyrose", project_z=True))
    fig.update_xaxes(type='linear', side='top',showgrid=False, 
                     title={'font': {'size': 18}, 'text': 'Sepal Length', 'standoff': 10},automargin=True)
    # fig.update_traces(hovertemplate="none")
    fig.update_layout(title='Robustness analysis', autosize=False,scene_camera_eye=dict(x=-1.8, y=-1.5, z=0.5),
        width=800, height=800,margin=dict(l=65, r=50, b=65, t=90))
    fig.update_scenes(yaxis_tickangle=0)
    fig.update_scenes(xaxis_tickangle=0)

    fig.write_image(out_fig_file) 
    return fig
def get_enz_model_use_c13(reaction_kcat_mw_file,json_model_path, c13reaction_file, percentage, kcat_database_combined_file, model_file, f, ptot, sigma, lowerbound, upperbound, json_output_file, reaction_mw_outfile):
    """Get new enzyme model using C13 reaction to calibration
    Arguments
    ----------
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * json_model_path: The file storing json model.
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * percentage:  percentage which needed change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.  
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_mw_outfile: changed file stored reaction kcat/mw.
    :return: new enzyme model.
    reaction_kcat_mw_file: 反应kcat/mw文件。
    * json_model_path: 存储json模型的文件。
    * c13reaction_file: 该文件包含有c13通量的反应列表。
    * percentage: 需要改变的百分比。
    * kcat_database_combined_file: 从autoPACMEN得到的kcat合并文件。 
    * model_file: 眼镜蛇模型。
    * f: 酶的质量分数。
    * ptot: 细胞中的总蛋白质部分。 
    * sigma: 酶的近似平均饱和度。
    * lowerbound:  酶浓度约束的下限。
    * upperbound:  酶浓度约束的上限值。 ptot*f
    * json_output_file: json文件存储json模型。
    * reaction_mw_outfile: 更改后的文件，存储反应kcat/mw。
    :返回：新的酶模型。
    """
    c13reaction_selecet = select_calibration_reaction_by_c13(
        reaction_kcat_mw_file, c13reaction_file, upperbound, percentage, sigma,json_model_path)
    print('need changing reaction: ')
    print(c13reaction_selecet)

    #if isinstance(df_reaction_select, pd.DataFrame):
    #    reaction_kcat_mw_file = "./analysis/reaction_change_by_biomass.csv"

    reaction_kapp_change_file = reaction_mw_outfile
    #c13reaction_selecet=['CS','ACONTa','ACONTb','ICDHyr','MALS', 'MDH', 'ICL', 'SUCOAS_reverse', 'SUCDi', 'AKGDH']
    change_reaction_list_round1 = change_reaction_kcat_by_database(json_model_path,
        c13reaction_selecet, kcat_database_combined_file, reaction_kcat_mw_file, reaction_kapp_change_file)
    print('changed reaction: ')
    print(change_reaction_list_round1)
    reaction_kcat_mw_file = reaction_mw_outfile
    trans_model2enz_json_model_split_isoenzyme(
        model_file, reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file)
    enz_model = get_enzyme_constraint_model(json_output_file)
    return enz_model
def select_calibration_reaction_by_c13(reaction_kcat_mw_file, c13reaction_file, enzyme_amount, percentage, sigma,json_model_path):
    """Get reaction list need change kcat using c13 data
    Arguments
    ----------
    * reaction_kcat_mw_file: original file stored kcat/mw of reaction.
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * enzyme_amount:  total enzyme amount(e.g. 0.227).
    * percentage:  percentage which needed change.
    * sigma: The approximated average saturation of enzyme. 
    * reaction_kcat_mw_file：存储反应的kcat/mw的原始文件。
    * c13reaction_file。该文件包含有c13通量的反应列表。
    * enzyme_amount: 酶的总量（如0.227）。
    * percentage：需要改变的百分比。
    * sigma: 酶的近似平均饱和度。
    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """    
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    c13reaction = pd.read_csv(c13reaction_file, index_col=0)
#     print(c13reaction)
    c13reaction_selecet = []
    for index, row in c13reaction.iterrows():
        if index in reaction_kcat_mw.index:
            enz_model=get_enzyme_constraint_model(json_model_path)
            enz_model_solution = cobra.flux_analysis.pfba(enz_model)
            print('通量为：',enz_model_solution[index])
            ECMpy_c13_reaction_flux = (reaction_kcat_mw.loc[index,'kcat_MW']/1000)*enzyme_amount*percentage*sigma
#             print(enzyme_amount,percentage,sigma)
#             print('ECMpy_c13_reaction_flux:',ECMpy_c13_reaction_flux)
            print(index,'ECMpy_c13_reaction_flux:',ECMpy_c13_reaction_flux,'flux:',row['Flux norm'])
            if enz_model_solution[index]< row['Flux norm']:
                if ECMpy_c13_reaction_flux < row['Flux norm']:
                      c13reaction_selecet.append(index)
    return(c13reaction_selecet)
    
def change_reaction_kcat_by_database(json_model_path,select_reactionlist, kcat_database_combined_file, reaction_kcat_mw_file, reaction_kapp_change_file):
    """Use the kcat in database to change reaction kcat in model
    Arguments
    ----------
    * json_model_path: The file storing json model.
    * select_reactionlist: reaction list need to change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_kapp_change_file: changed file stored reaction kcat/mw.
    :return: a dataframe stored new reaction kcat/mw .
    """
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    Brenda_sabio_combined_select = json_load(kcat_database_combined_file)
    json_model=cobra.io.load_json_model(json_model_path)
    reaction_change_accord_database = []
    for eachreaction in select_reactionlist:
        select_reaction = json_model.reactions.get_by_id(eachreaction)
        if "ec-code" in select_reaction.annotation.keys():
            ec_number = select_reaction.annotation["ec-code"]
            kcat_max_list = []
            if isinstance(ec_number, str):
                if ec_number in Brenda_sabio_combined_select.keys():
                    reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat_max']
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max :
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max 
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                        reaction_change_accord_database.append(eachreaction) 
            else:
                for eachec in ec_number:
                    if eachec in Brenda_sabio_combined_select.keys():
                        kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat_max'])
                reaction_kcat_max = np.max(kcat_max_list)     
                if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max :
                    reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max
                    reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                    reaction_change_accord_database.append(eachreaction)              
    reaction_kcat_mw.to_csv(reaction_kapp_change_file)
    return(reaction_change_accord_database)