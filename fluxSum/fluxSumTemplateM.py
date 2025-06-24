import cobra

def isAuto(rxn, meta):
    pro=[x.id for x in rxn.products]
    whole=[x.id for x in rxn.metabolites]
    con=list(set(whole)-set(pro))

    proNc=[x[:-2] for x in pro]
    conNc=[x[:-2] for x in con]

    both=set(proNc)&set(conNc)
    whole=set(proNc)|set(conNc)

    if not meta in whole:
        raise ValueError('input metabolite is not related with input rxn')

    if meta in both:
        return True
    else:
        return False

def run(templateModelPath, autoFlag=False):
    result={}
    model=cobra.io.read_sbml_model(templateModelPath)
    for m in model.metabolites:
        temp={'pre':[], 'next':[]}
        for r in m.reactions:
            if autoFlag:
                if isAuto(r, m.id[:-2]):
                    continue

            pro=[x.id for x in r.products]
            if m.id in pro:
                temp['pre'].append(r.id)
            else:
                temp['next'].append(r.id)
        result[m.id]=temp
    Template=result
    if autoFlag:
        Template={x[:-2]:{'pre':[], 'next':[]} for x in result}
        for mc in result:
            m=mc[:-2]
            Template[m]['pre']+=result[mc]['pre']
            Template[m]['next']+=result[mc]['next']

        for m in Template:
            Template[m]['pre']=list(set(Template[m]['pre']))
            Template[m]['next']=list(set(Template[m]['next']))
    return Template
