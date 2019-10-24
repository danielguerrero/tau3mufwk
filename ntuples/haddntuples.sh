#samples associated to tau3mu with 0 PU
#hadd  DsTau3muPU0.root  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Ds_Tau3Mu_PU0_7May_Tau3MuSamples/output/ | grep '\.root'`
#hadd Single Muon with 0 PU
#hadd SingleMuPU0.root         -f `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Mu_FlatPt2to100-pythia8-gun_PU0_SampleFullStat/output/ | grep '\.root'`
#hadd SingleMuPU0_EMTF.root    -f `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Mu_FlatPt2to100-pythia8-gun_PU0_SampleFullStat_EMTF/output/ | grep '\.root'`
#samples associated to tau3mu with 200 PU
#hadd  DsTau3muPU200.root  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Tau3Mu_PU200_SampleFullStat/output/ | grep '\.root'`
#sample minBias sample with 200 PU 
#hadd MinBiasPU200.root  -f  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Nu_E10-pythia8-gun_SampleFullStat/output/ | grep '\.root'`
#hadd MinBiasPU200_EMTF.root   -f  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Nu_E10-pythia8-gun_SampleFullStat_EMTF/output/ | grep '\.root'`
#Single Muon with 200 PU
#hadd SingleMuPU200.root -f `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Mu_FlatPt2to100-pythia8-gun_SampleFullStat/output/ | grep '\.root'`
#hadd SingleMuPU200_EMTF.root -f  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Mu_FlatPt2to100-pythia8-gun_SampleFullStat_EMTF/output/ | grep '\.root'`

#samples for tau3mu for TDR
#hadd  DsTau3muPU0_TDR.root  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Tau3Mu_PU0_Config_EMTFmode_TkMuStub1/output/ | grep '\.root'`
#hadd  DsTau3muPU140_TDR.root  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Tau3Mu_PU140_Config_EMTFmode_TkMuStub1/output/ | grep '\.root'`
#hadd  DsTau3muPU200_TDR.root  `xrdfsls -u /store/user/guerrero/L1MuTrks_ntuples/Tau3Mu_PU200_Config_EMTFmode_TkMuStub1/output/ | grep '\.root'`

#Fixed samples 
#hadd DsTau3muPU0_TDR_fixed.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU0_12Sett_Tau3MuSamples/output/ | grep '\.root'`
#hadd MinBias_TDR_fixed.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/NuGun_200PU_11Sett_TDR_MC_EMTFmode_FixCppFlag_SvenTPFix/output | grep '\.root'`

#Lastest samples 28/09/2019
#hadd DsTau3muPU0_TDR.root   `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU0_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`
#hadd DsTau3muPU200_TDR.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU200_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`
#hadd MinBiasPU200_TDR.root  `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/NuGun_200PU_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`
#hadd SingleMuPU0_TDR.root   `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/MuMu_flatPt_0PU_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`
#hadd SingleMuPU200_TDR.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/MuMu_flatPt_200PU_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`

#Lastest samples 22-24/10/2019
#hadd DsTau3muPU0_TDR.root   `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU0_21Ott_TDR_MC_TkMuStubv3/output/ | grep '\.root'`
hadd DsTau3muPU200_TDR.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU200_21Ott_TDR_MC_TkMuStubv3_extDanieListTau3muv1_resub/output/ | grep '\.root'`
#hadd MinBiasPU200_TDR.root  `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/NuGun_200PU_21Ott_TDR_MC_TkMuStubv3/output/ | grep '\.root'`

#hadd DsTau3muPU200_TDR_small.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU200_21Ott_TDR_MC_TkMuStubv3/output/ | grep '\.root'`
#hadd SingleMuPU0_TDR.root   `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/MuMu_flatPt_0PU_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`
#hadd SingleMuPU200_TDR.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/MuMu_flatPt_200PU_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`
