import copy


def isEnd(R_i, excepts):
    nextCandi=[]
    nextExcept=[]
    for nextMeta in R_i.metabolites:
        if nextMeta.id == excepts:
            continue
        R_i2=[x for x in nextMeta.reactions if x.id!=R_i.id]
        if len(R_i2)==1:
            nextCandi.append(R_i2[0])
            nextExcept.append(nextMeta.id)
    if len(set([x.id for x in nextCandi]))==1:
        return nextCandi[0], nextExcept[0]
    else:
        return 'end'
    
    
def getDirection(R_i, M_i):
    P_metas=[x.id for x in R_i.products]
    C_metas=[x.id for x in R_i.reactants]
    if R_i.reversibility:
        return 'B' #both
    elif M_i in P_metas:
        return 'P' #producing
    else:
        return 'C' #consuming


def getLinearLinkage(CtxGem, targets):
    '''
    CtxGems : GEMs for paitients
    targets : list of targets ([metabolite.id])
    
    It return dictionary which keys are target metabolites and values are dictionary
    The dictionarys for values has the keys for each reactions which directly linked to the target metabolites(dist1)
    and the values consists linearly linked reactions to dist1 reaction (starts from dist1 rxns)
    '''
    
    Result={}
    
    for eachTarget in targets:
        Result[eachTarget]={}
        metabolite=CtxGem.metabolites.get_by_id(eachTarget)
        cnt=0
        MetaKey=copy.deepcopy(eachTarget)
        for dist1Rxn in metabolite.reactions:
            eachTarget=metabolite.id
            RxnKey=copy.deepcopy(dist1Rxn.id)
            Result[MetaKey][RxnKey]=[]
            Rhist=[dist1Rxn.id]
            Mhist=[eachTarget]
            
            Direction=[getDirection(dist1Rxn, metabolite.id)]
            while True:
                P=isEnd(dist1Rxn, eachTarget)
                if type(P)==tuple:
                    dist1Rxn=P[0]
                    eachTarget=P[1]
                    if dist1Rxn in Rhist:
                        break
                    if eachTarget in Mhist:
                        break
                    cnt+=1
                    Rhist.append(dist1Rxn.id)
                    Mhist.append(eachTarget)
                    Direction.append(getDirection(dist1Rxn, eachTarget))
                    if 'P' in Direction and 'C' in Direction:
                        Rhist=Rhist[:-1]
                        Mhist=Mhist[:-1]
                        Direction=Direction[:-1]
                        break
                else:
                    break
                if cnt>100:
                    raise
#             Result[MetaKey][RxnKey]={'R':Rhist,'M':Mhist,'D':Direction} #for debug
            Result[MetaKey][RxnKey]=Rhist
    
    return Result


def isRelated(Rxntarget, linkedMap):
    RelMetas=[]
    for Meta in linkedMap:
        for D1Rxn in linkedMap[Meta]:
            if Rxntarget in linkedMap[Meta][D1Rxn]:
                RelMetas.append(Meta)
                
    return list(set(RelMetas))


def evaluate(defaultF, fluxProf, target, weight, Topology, riskScoreMetas):
    
#     Topology=getLinearLinkage(CTXGEM, riskScoreMetas)
    adjMetas=isRelated(target, Topology)
    
    Score=0
    AdjScore=0
    for factor in riskScoreMetas: #for metabolites consists the risk score
        W=weight[factor]
        F=fluxProf[factor]
        
        if factor in adjMetas:
            adjF=defaultF[factor]
        else:
            adjF=fluxProf[factor]
        
        Score+=W*F
        AdjScore+=W*adjF
    
    return AdjScore, Score
    