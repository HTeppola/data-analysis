__author__ = 'dlonardoni'

import matplotlib as m
m.use('Qt5Agg')

import pylab as plt
import cellCultureNet as Net
from tools import tools
import pickle
import numpy as np
G=pickle.load(open("savedGraph/4fgauss_0.005_1024_1_.graph"))
G=G.subgraph(list(range(20)))
ifNMDA=0
magnitude=0.08
mAMPA=1
G.name="1024AMPA_test_%s_%s"%(mAMPA,magnitude)
_peso=37.5
_rpeso=13.9
_tot=1
groupname="IG_final_control2"
Net.stimPar['active']=1
Net.stimPar['positions']=np.array([[.1,.9],[.15,.10],[.9,.1],[.1,.65]])
np.random.shuffle(Net.stimPar['positions'])
HEADER=Net.NetworkMaker(mAMPA=mAMPA,ifNMDA=ifNMDA,replicate=False,graph=G,stp=100,stpinh=100,N=len(G.nodes()),peso=_peso,rpeso=_rpeso,ts=1000*40,tipo=0,sigma=0.0005,percInh=0.2,tot=_tot,r=0,group=groupname,magnitude=magnitude,vblock=-35,kb=6,U=.5)

sname,gname=Net.formatNames(G.name,'tau_w',_tot,0,0,_peso,_rpeso)
gname += "_%d_NMDA_%0.2d_NOBIC" % (ifNMDA, magnitude)
sname += "sfinal_%d_NMDA_%0.2f_NOBIC" % (ifNMDA, magnitude)
Net.SetAttribute(G=G)
tools.SaveGraph(G,fname=gname,folder="Results/"+groupname+"/Graph/")
Net.go()

tools.SaveSpikeTrains(Net.CellList,HEADER,fname=sname,folder="Results/"+groupname+"/Spike/")
spikelist=[]
for cell in Net.CellList:
    try:
        formatData={'type':cell.ID_syn,"spikes":np.array(cell.record["Stimspk"]),"ID":np.array([cell.num]),"x":np.array([cell.pos["x"]]),"y":np.array([cell.pos["y"]])}
        spikelist.append(formatData)
    except:
        pass
np.save(open('Results/'+groupname+'/Spike/'+sname+'_STIM.npy', 'w'),spikelist)



clist=np.load('Results/'+groupname+'/Spike/'+sname+'.npy', allow_pickle=True)

for cell in clist[1:]:
    plt.plot(cell['spikes'],np.ones_like(cell['spikes'])*cell['ID'],'ok')
plt.ion()
plt.show()