# generate biggr2ec
reactionf = '/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_models_reactions.txt'
reactions = pd.read_csv(reactionf, sep='\t')
ecs = []

for link in list(reactions.database_links.values):
    link = database_links_reformat(link)
    ids, rhea_ids, mnxs, seeds, biocyc, ec, keggr = links_to_id(link)
    ecs.append(ec)

biggr2ec = {}
biggec2r = {}
for index, row in reactions.iterrows():
    id = row['bigg_id']
    oldid = row['old_bigg_ids']
    
    ec = ecs[index]
    for eitem in ec : 
        try:
            biggec2r[eitem].append(id)
            if '; ' in oldid:
                oids = oldid.split('; ')
                for oid in oids:
                    biggec2r[eitem].append(oid) 
            else:
                biggec2r[eitem].append(oldid)
        except:
            biggec2r[eitem] = [id]
            if '; ' in oldid:
                oids = oldid.split('; ')
                biggec2r[eitem] = oids
            else:
                biggec2r[eitem] = [oldid]
    try:
        biggr2ec[id] += ec
        if '; ' in oldid:
            oids = oldid.split('; ')
            for oid in oids:
                biggr2ec[oid] += ec
        else:
            biggr2ec[oldid] += ec
    except:
        biggr2ec[id] = ec
        if '; ' in oldid:
            oids = oldid.split('; ')
            for oid in oids:
                biggr2ec[oid] = ec
        else:
            biggr2ec[oldid] = ec
       
#save biggr2ec


# Save biggr2ec dictionary to a file
with open('biggr2ec.pkl', 'wb') as f:
    pickle.dump(biggr2ec, f)
with open('biggec2r.pkl', 'wb') as f:
    pickle.dump(biggec2r, f)
print("biggr2ec dictionary saved successfully.")