#Signal studies in fiducial region (PU=0,200)
tau3muAnalysis.exe --config config/DsToTau3Mu_fiducial.cfg --input ntuples/DsTau3muPU0_TDR.root   --tagname Tau3muPU0_fiducial    --is-signal
tau3muAnalysis.exe --config config/DsToTau3Mu_fiducial.cfg --input ntuples/DsTau3muPU200_TDR.root --tagname Tau3muPU200_fiducial  --is-signal
tau3muAnalysis.exe --config config/DsToTau3Mu_fiducialextended.cfg --input ntuples/DsTau3muPU0_TDR.root   --tagname Tau3muPU0_fiducialextended    --is-signal
tau3muAnalysis.exe --config config/DsToTau3Mu_fiducialextended.cfg --input ntuples/DsTau3muPU200_TDR.root --tagname Tau3muPU200_fiducialextended  --is-signal
#Rate studies  (PU 200)
tau3muAnalysis.exe --config config/DsToTau3Mu_fiducialextended.cfg --input ntuples/MinBiasPU200_TDR.root  --tagname MinBiasPU200                  --is-minbias