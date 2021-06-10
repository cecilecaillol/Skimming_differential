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
      place="/hdfs/store/user/aloeliger/SMHTT_2016_data_7_Jul_2020/"
      datatype="data"
      name=[
"data_SingleElectron_Run2016B_v1","DataBv1","0",
"data_SingleElectron_Run2016B_v2","DataBv2","0",
"data_SingleElectron_Run2016C","DataC","0",
"data_SingleElectron_Run2016D","DataD","0",
"data_SingleElectron_Run2016E","DataE","0",
"data_SingleElectron_Run2016F","DataF","0",
"data_SingleElectron_Run2016G","DataG","0",
"data_SingleElectron_Run2016H","DataH","0"
]


    if args.sample=="embedded":
      place="/hdfs/store/user/aloeliger/SMHTT_2016_embedded_7_Jul_2020/"
      datatype="embedded"
      name=[
"embedded_EmbeddingRun2016B_ElTauFinalState","embeddedB","0",
"embedded_EmbeddingRun2016C_ElTauFinalState","embeddedC","0",
"embedded_EmbeddingRun2016D_ElTauFinalState","embeddedD","0",
"embedded_EmbeddingRun2016E_ElTauFinalState","embeddedE","0",
"embedded_EmbeddingRun2016F_ElTauFinalState","embeddedF","0",
"embedded_EmbeddingRun2016G_ElTauFinalState","embeddedG","0",
"embedded_EmbeddingRun2016H_ElTauFinalState","embeddedH","0"
]


    if args.sample=="V":
      place="/hdfs/store/user/aloeliger/SMHTT_2016_7_Jul_2020/"
      datatype="mc"
      name=[
"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","DY1_v1","Z",
"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","DY2_v1","Z",
"DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","DY3","Z",
"DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","DY4","Z",
"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","DY_v1","Z",
"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2","DY_v2","Z",
"EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","EWKWminus_v1","W",
"EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","EWKWminus_v2","W",
"EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2","EWKWminus_v3","W",
"EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","EWKWplus_v1","W",
"EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","EWKWplus_v2","W",
"EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2","EWKWplus_v3","W",
"EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","EWKZLL_v1","Z",
"EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2","EWKZLL_v2","Z",
"EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","EWKZNuNu_v1","Z",
"EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","EWKZNuNu_v2","Z",
"EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2","EWKZNuNu_v3","Z",
"W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","W1","W",
"W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","W2_v1","W",
"W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","W2_v2","W",
"W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","W3_v1","W",
"W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","W3_v2","W",
"W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","W4_v1","W",
"W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","W4_v2","W",
"W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1","W4_v3","W",
"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1","WG_v1","W",
"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1","WG_v2","W",
"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1","WG_v3","W",
"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","W_v1","W",
"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2","W_v2","W"
]

    if args.sample=="VV":
      place="/hdfs/store/user/aloeliger/SMHTT_2016_7_Jul_2020/"
      datatype="mc"
      name=[
"ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","ST_t_antitop","0",
"ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","ST_t_top","0",
"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1","ST_tW_antitop","0",
"ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1","ST_tW_top","0",
"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","TT","0",
"VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","VV2L2Nu_v1","0",
"VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","VV2L2Nu_v2","0",
"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WW1L1Nu2Q","0",
"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WZ1L1Nu2Q","0",
"WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WZ1L3Nu","0",
"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WZ2L2Q","0",
"WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WZ3L1Nu","0",
"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","ZZ2L2Q","0",
"ZZTo4L_13TeV-amcatnloFXFX-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2","ZZ4L","0"
]

    if args.sample=="signal":
      place="/hdfs/store/user/aloeliger/SMHTT_2016_7_Jul_2020/"
      datatype="mc"
      name=[
"GluGluHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","GGHTT_v1","Z",
"GluGluHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v3","GGHTT_v2","Z",
"GluGluHToWWTo2L2Nu_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","GGHWW","Z",
"GluGluZH_HToWW_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","GGZHWW","0",
"HWminusJ_HToWW_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WminusHWW","0",
"HWplusJ_HToWW_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WplusHWW","0",
"HZJ_HToWW_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","ZHWW","0",
"VBFHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","VBFHTT","Z",
"VBFHToWWTo2L2Nu_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","VBFHWW","Z",
"WminusHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WplusHTT","0",
"WplusHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","WminusHTT","0",
"ZHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","ZHTT","0",
"ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","GGZHLLTT","0",
"ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","GGZHNNTT","0",
"ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1","GGZHQQTT","0",
"ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2","ttHnonbb","0"
]

    datadir="/nfs_scratch/caillol/differentialet2016_3aug/"
    all_File = open("do_submit_"+ args.sample+"_et.sh" , 'w')
    line=""
    for j in range(0,len(name)/3):
       print name[3*j],name[3*j+1],name[3*j+2]
       submit_File = open("Submit_"+name[3*j+1]+"_et.sh" , 'w')
       line+="mkdir -p "+datadir+"Out_"+name[3*j+1]+"\n"
       line+="sh Submit_"+name[3*j+1]+"_et.sh" +"\n"
       f=os.popen("ls -rt " + place+name[3*j] + "/make* ")
       command1=""
       ligne=0
       for i in f.readlines():
	   command1=command1+"./skim_et.exe "+datatype + " " +datadir+"Out_"+name[3*j+1]+"/"+name[3*j+1]+str(ligne)+".root " + i[0:-1] + " " + name[3*j+2] +"\n"
           ligne=ligne+1
       submit_File.write(command1)
       submit_File.close()
    all_File.write(line)
    all_File.close()

