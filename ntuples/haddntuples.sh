#Lastest samples 22-24/10/2019
hadd DsTau3muPU0_TDR.root   `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU0_21Ott_TDR_MC_TkMuStubv3/output/ | grep '\.root'`
hadd DsTau3muPU200_TDR.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU200_21Ott_TDR_MC_TkMuStubv3_extDanieListTau3muv1_resub/output/ | grep '\.root'`
hadd MinBiasPU200_TDR.root  `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/NuGun_200PU_21Ott_TDR_MC_TkMuStubv3/output/ | grep '\.root'`

#hadd DsTau3muPU200_TDR_small.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/TauTo3Mu_TuneCP5_14TeV-pythia8_TDR_PU200_21Ott_TDR_MC_TkMuStubv3/output/ | grep '\.root'`
#hadd SingleMuPU0_TDR.root   `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/MuMu_flatPt_0PU_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`
#hadd SingleMuPU200_TDR.root `xrdfsls -u /store/user/lcadamur/L1MuTrks_ntuples/MuMu_flatPt_200PU_28Sett_TDR_MC_EMTFpp_gentrackinfo_v2/output/ | grep '\.root'`