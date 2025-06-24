import cobra
import pandas as pd
import time
from simulation.fluxSum import fluxSumTempleteM

class ContributionFluxSum():
    
    def __init__(self, fluxPath, templateModel):
        self.Templete=fluxSumTempleteM.run(templateModel, False)
        self.rawflux=pd.read_csv(fluxPath, index_col=0)
        self.model=cobra.io.read_sbml_model(templateModel)
        self.netFlux, self.direction = self.preprocessing()
        self.pathDict={}
        
#     def getPathInfo(self):
#         pathDict={}
#         newPath=pd.read_excel('/data/grlee/NewRecon2M.2Subsystems_210228.xlsx')
#         convert={newPath['Recon2M.2'][x]:newPath['NewReconSubsystem'][x] for x in newPath.index if newPath['Recon2M.2'][x]==newPath['Recon2M.2'][x]}
#         for r in self.model.reactions:
#             pathDict[r.id]=r.subsystem
#         getpras=[x for x in pathDict if pathDict[x]=='']
#         for r in getpras:
#             pathDict[r]='Unassigned'
            
            
#         for r in pathDict:
#             pathOld=pathDict[r]
#             pathNew=convert[pathOld]
#             pathDict[r]=pathNew
#         self.pathDict=pathDict
        
    def preprocessing(self):
        tdf=self.rawflux.fillna(0)
        drops=[]
        for i in tdf.index:
            data=list(tdf.loc[i])
            if max(data)==min(data) and min(data)==0:
                drops.append(i)
        print(len(drops))
        tdf=tdf.drop(drops)
        directionDF=(tdf>=0)*1+(tdf<0)*-1
        
        return tdf, directionDF
    
    def getPutMeta(self):
        meta=[]
        for i in self.netFlux.index:
            r=self.model.reactions.get_by_id(i)
            for m in r.metabolites:
#                 meta.append(m.id[:-2])
                meta.append(m.id)
        meta=list(set(meta))
        
        return meta
    
    def calculateCFluxSum(self, mode='producingRxns'):
        print('mode is',mode, 'among consumingRxns AND producingRxns')
        st=time.time()
        RECORD={}
        targetMetas=self.getPutMeta()
#         self.getPathInfo()
        cnt=0
        for metabol in targetMetas:
            temp1=getDirectionMetaFluxDf(metabol, mode, self.netFlux, self.Templete, self.direction)
            self.checkRxnCoefficient(metabol, temp1)
            temp1=resolvePathwayResolutionSumPath(temp1, self.pathDict)
            temp1=temp1.rename({x:metabol+'^'+x for x in temp1.index})
            temp1=temp1.transpose()
            temp1=temp1.to_dict()
            RECORD.update(temp1)
            
            cnt+=1
            if cnt%100==0:
                print(cnt, 'metabolites done', time.time()-st)
        RESULT=pd.DataFrame.from_dict(RECORD)
        RESULT=RESULT.transpose()
        return RESULT
    
    def checkRxnCoefficient(self, metabolite, df):
        for rxn in df.index:
            calledRxn=self.model.reactions.get_by_id(rxn)
            for meta in calledRxn.metabolites:
#                 if meta.id[:-2]==metabolite:
                if meta.id==metabolite:
                    c_meta=meta.id
                    
            coef=calledRxn.get_coefficient(c_meta)
            coef=abs(coef)
            if coef!=1:
                df.loc[rxn]=df.loc[rxn]*coef

def getSourceUseInfo(qMeta, templete, direction):
    ## For previouis edge
    
    metaTems=templete[qMeta]['pre']
    call=set(metaTems)&set(direction.index)
    preDF=direction.loc[call]
    
    result={}
    for sam in preDF:
        temp={'producingRxns':[], 'consumingRxns':[]}
        rights=preDF[preDF[sam]==1].index
        antis=set(preDF.index)-set(rights)

        temp['producingRxns']+=list(rights)
        temp['consumingRxns']+=list(antis)

        result[sam]=temp
    ## For next edge
    metaTems=templete[qMeta]['next']
    call=set(metaTems)&set(direction.index)
    nextDF=direction.loc[call]
    
    for sam in nextDF:
        rights=nextDF[nextDF[sam]==1].index
        antis=set(nextDF.index)-set(rights)

        result[sam]['producingRxns']+=list(antis)
        result[sam]['consumingRxns']+=list(rights)
        
    return result

def getDirectionMetaFluxDf(testMeta, mode,flxDf, tplete, direc):
    temp=pd.DataFrame(getSourceUseInfo(testMeta, tplete, direc))
    
    tensorMulTem={}
    for sam in temp:
        _={}
        for r in temp[sam][mode]:
            _[r]=1
        tensorMulTem[sam]=_
    tensorMulTem=pd.DataFrame.from_dict(tensorMulTem)
    tensorMulTem=tensorMulTem.fillna(0)
    
    
#     result=flxDf.loc[tensorMulTem.index]*tensorMulTem
    
    result=flxDf.loc[tensorMulTem.index]
    result=result.sort_index(axis=1)
    result=result.sort_index(axis=0)
    tensorMulTem=tensorMulTem.sort_index(axis=1)
    tensorMulTem=tensorMulTem.sort_index(axis=0)
    
    result=result*tensorMulTem
    
    result=abs(result)
    ## opposite direction flux of input mode is zero
    return result

def resolvePathwayResolutionSumPath(data, pathInfo={}):
    pathDict={}
    path='justSum'
    for rxn in data.index:
        if not path in pathDict:
            pathDict[path]=[]
        pathDict[path].append(rxn)
    
    result={}
    for path in pathDict:
        call=pathDict[path]
        temp=data.loc[call].sum().to_dict()
        result[path]=temp
    result=pd.DataFrame(result)
    result=result.transpose()
    
    return result

        
    

    