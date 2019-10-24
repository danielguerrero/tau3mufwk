##Properties default and fiducial cuts
multiplotter.exe --inputFile histograms/Tau3muPU0_fiducial_genproperties.root   --samplename '#tau#rightarrow3#mu' --pileup   0 --configTH1 config/plot_gen_fiducial.cfg
multiplotter.exe --inputFile histograms/Tau3muPU200_fiducial_genproperties.root --samplename '#tau#rightarrow3#mu' --pileup 200 --configTH1 config/plot_gen_fiducial.cfg
multiplotter.exe --inputFile histograms/Tau3muPU0_fiducialextended_genproperties.root   --samplename '#tau#rightarrow3#mu' --pileup   0 --configTH1 config/plot_gen_fiducialextended.cfg
multiplotter.exe --inputFile histograms/Tau3muPU200_fiducialextended_genproperties.root --samplename '#tau#rightarrow3#mu' --pileup 200 --configTH1 config/plot_gen_fiducialextended.cfg
##Signal efficiencies
multiplotter.exe --inputFile histograms/Tau3muPU0_fiducial_muontriggers.root   --samplename '#tau#rightarrow3#mu' --pileup 0   --configEff config/plot_triggers_eff_fiducial.cfg --configTH1 config/plot_triggers_trk_fiducial.cfg
multiplotter.exe --inputFile histograms/Tau3muPU200_fiducial_muontriggers.root --samplename '#tau#rightarrow3#mu' --pileup 200 --configEff config/plot_triggers_eff_fiducial.cfg --configTH1 config/plot_triggers_trk_fiducial.cfg
multiplotter.exe --inputFile histograms/Tau3muPU0_fiducialextended_muontriggers.root   --samplename '#tau#rightarrow3#mu' --pileup 0   --configEff config/plot_triggers_eff_fiducialextended.cfg --configTH1 config/plot_triggers_trk_fiducialextended.cfg
multiplotter.exe --inputFile histograms/Tau3muPU200_fiducialextended_muontriggers.root --samplename '#tau#rightarrow3#mu' --pileup 200 --configEff config/plot_triggers_eff_fiducialextended.cfg --configTH1 config/plot_triggers_trk_fiducialextended.cfg
##Rates
multiplotter.exe --inputFile histograms/MinBiasPU200_muontriggerrates.root    --samplename MinBias --pileup 200 --configTH1 config/plot_triggers_rate_minbias.cfg
multiplotter.exe --inputFile histograms/MinBiasPU200_muontriggerrates.root    --samplename MinBias --pileup 200 --configTH1 config/plot_triggers_trk_minbias.cfg
##Plot ROC curve
source scripts/runROC.sh