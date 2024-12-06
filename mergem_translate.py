f ='/ibex/user/niuk0a/funcarve/cobra/iAF987.xml'


import cobra
import mergem
from mergem import translate, load_model, save_model, map_localization, map_metabolite_univ_id, map_reaction_univ_id, get_metabolite_properties, get_reaction_properties, update_id_mapper
# def translate_model(input):
#     if type(input) == str:
#         model = cobra.io.read_sbml_model(input)
#     elif type(input) == cobra.Model:
#         model = input
    # results = mergem.merge(model, set_objective='merge', exact_sto=False, use_prot=False, extend_annot=False, trans_to_db=None)
    # merged_model = results['merged_model']
    # jacc_matrix = results['jacc_matrix']
    # num_met_merged = results['num_met_merged']
    # num_reac_merged = results['num_reac_merged']
    # met_sources = results['met_sources']
    # reac_sources = results['reac_sources']

mm = mergem.load_model(f)
mrs = [r.id for r in mm.reactions]
print('bf',mrs[:100])
tm = mergem.translate(mm,'seed')
rs = [r.id for r in tm.reactions]
print('aft',rs[:100])
#save the translated model
save_model(tm, 'iAF987_seed.xml')