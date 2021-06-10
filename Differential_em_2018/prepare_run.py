import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sample', '-s', default=None, help='Output name')

args = parser.parse_args()


if __name__ == "__main__":
    place=""
    datatype=""
    name=[]

    if args.sample=="data":
       place="/hdfs/store/user/caillol/SMHTT_2018_27jul_data/"
       datatype="data"
       name=[
"data_MuonEG_Run2018A-17Sep2018","DataA","0",
"data_MuonEG_Run2018B-17Sep2018","DataB","0",
"data_MuonEG_Run2018C-17Sep2018","DataC","0",
"data_MuonEG_Run2018D-PromptReco","DataD","0"
]

    if args.sample=="V":
       place="/hdfs/store/user/caillol/SMHTT_2018_27jul_data/"
       datatype="mc"
       name=[
"DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","DY1_v1","Z",
"DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","DY1_v2","Z",
"DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","DY2_v1","Z",
"DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","DY2_v2","Z",
"DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","DY3_v1","Z",
"DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","DY3_v2","Z",
"DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","DY4_v1","Z",
"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","DY_v1","Z",
"EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1","EWKWminus","W",
"EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1","EWKWplus","W",
"EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1","EWKZLL","Z",
"EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1","EWKZNuNu","Z",
"W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","W1_v1","W",
"W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","W2_v1","W",
"W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","W2_v2","W",
"W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","W3_v1","W",
"W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","W3_v2","W",
"W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","W4_v1","W",
"W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","W4_v2","W",
"WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1","WG","0",
"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2","W_v1","W"
]

    if args.sample=="VV":
       place="/hdfs/store/user/caillol/SMHTT_2018_27jul_data/"
       datatype="mc"
       name=[
"VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1","VV2L2Nu","0",
"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1","WW1L1Nu2Q","0",
"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1","WZ1L1Nu2Q","0",
"WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1","WZ1L3Nu","0",
"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1","WZ2L2Q","0",
"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1","WZ3L1Nu_v1","0",
"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15_ext1-v2","WZ3L1Nu_v2","0",
"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1","ZZ2L2Q","0",
"ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1","ZZ4L","0",
"ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_-102X_upgrade2018_realistic_v15-v1","ST_t_antitop","0",
"ST_t-channel_top_5f_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1","ST_t_top","0",
"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext1-v1","ST_tW_antitop","0",
"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext1-v1","ST_tW_top","0",
"WW_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v1","WW_v1","0",
"WW_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v2","WW_v2","0",
"WZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v1","WZ_v1","0",
"WZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v3","WZ_v2","0",
"ZZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v2","ZZ_v1","0"
]

    if args.sample=="TT":
       place="/hdfs/store/user/caillol/SMHTT_2018_27jul_data/"
       datatype="mc"
       name=[
"TTToHadronic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1","TTToHadronic_v1","0",
"TTToHadronic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext2-v2","TTToHadronic_v2","0",
"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1","TTTo2L2Nu","0",
"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1","TTToSemiLeptonic_v1","0",
"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext3-v2","TTToSemiLeptonic_v2","0"
]

    if args.sample=="embedded":
       place="/hdfs/store/user/caillol/SMHTT_2018_27jul_embedded/"
       datatype="embedded"
       name=[
"embedded_EmbeddingRun2018A_ElMuFinalState","EmbeddedA","0",
"embedded_EmbeddingRun2018B_ElMuFinalState","EmbeddedB","0",
"embedded_EmbeddingRun2018C_ElMuFinalState","EmbeddedC","0",
"embedded_EmbeddingRun2018D_ElMuFinalState","EmbeddedD","0"
]

    if args.sample=="signal":
       place="/hdfs/store/user/caillol/SMHTT_2018_27jul_data/"
       datatype="mc"
       name=[
"ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v1","GGZHLLTT","0",
"ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v1","GGZHNNTT","0",
"ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v1","GGZHQQTT","0",
"GluGluHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2","GGHTT","Z",
"GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_-102X_upgrade2018_realistic_v15-v1","GGHWW","Z",
"GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5_PSweights_-102X_upgrade2018_realistic_v15-v1","GGZHWW","0",
"HWminusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5_-102X_upgrade2018_realistic_v15-v1","WminusHWW","0",
"HWplusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5_-102X_upgrade2018_realistic_v15-v1","WplusHWW","0",
"HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5_-102X_upgrade2018_realistic_v15-v1","ZHWW","0",
"ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1","ttHTT","0",
"VBFHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15_ext1-v1","VBFHTT","Z",
"VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_-102X_upgrade2018_realistic_v15-v1","VBFHWW","Z",
"WminusHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2","WminusHTT","0",
"WplusHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2","WplusHTT","0",
"ZHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2","ZHTT","0"
]

    datadir="/nfs_scratch/caillol/differentialem2018_3aug/"
    all_File = open("do_submit_"+ args.sample+"_em.sh" , 'w')
    line=""
    for j in range(0,len(name)/3):
       print name[3*j],name[3*j+1],name[3*j+2]
       submit_File = open("Submit_"+name[3*j+1]+"_em.sh" , 'w')
       line+="mkdir -p "+datadir+"Out_"+name[3*j+1]+"\n"
       line+="sh Submit_"+name[3*j+1]+"_em.sh" +"\n"
       f=os.popen("ls -rt " + place+name[3*j] + "/make* ")
       command1=""
       ligne=0
       for i in f.readlines():
           command1=command1+"./skim_em.exe "+datatype + " " +datadir+"Out_"+name[3*j+1]+"/"+name[3*j+1]+str(ligne)+".root " + i[0:-1] + " " + name[3*j+2] +"\n"
           ligne=ligne+1
       submit_File.write(command1)
       submit_File.close()
    all_File.write(line)
    all_File.close()

