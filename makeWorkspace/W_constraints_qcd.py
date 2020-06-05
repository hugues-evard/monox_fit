import ROOT
from counting_experiment import *
from W_constraints import do_stat_unc
# Define how a control region(s) transfer is made by defining cmodel provide, the calling pattern must be unchanged!
# First define simple string which will be used for the datacard 
model = "qcd_wjets"
def cmodel(cid,nam,_f,_fOut, out_ws, diag,year, convention="BU"):
  
  # Some setup
  _fin    = _f.Get("category_%s"%cid)
  _wspace = _fin.Get("wspace_%s"%cid)


  # ############################ USER DEFINED ###########################################################
  # First define the nominal transfer factors (histograms of signal/control, usually MC 
  # note there are many tools available inside include/diagonalize.h for you to make 
  # special datasets/histograms representing these and systematic effects 
  # but for now this is just kept simple 
  processName  = "WJets" # Give a name of the process being modelled
  metname      = 'mjj'    # Observable variable name 
  targetmc     = _fin.Get("signal_qcdwjets")      # define monimal (MC) of which process this config will model
  controlmc    = _fin.Get("Wmn_qcdwjets")  # defines in / out acceptance
  controlmc_e  = _fin.Get("Wen_qcdwjets")  # defines in / out acceptance

  # Create the transfer factors and save them (not here you can also create systematic variations of these 
  # transfer factors (named with extention _sysname_Up/Down
  
  WScales = targetmc.Clone(); WScales.SetName("qcd_wmn_weights_%s"%cid)
  WScales.Divide(controlmc);  _fOut.WriteTObject(WScales)  

  WScales_e = targetmc.Clone(); WScales_e.SetName("qcd_wen_weights_%s"%cid)
  WScales_e.Divide(controlmc_e);  _fOut.WriteTObject(WScales_e)  



  #######################################################################################################

  _bins = []  # take bins from some histogram, can choose anything but this is easy 
  for b in range(targetmc.GetNbinsX()+1):
    _bins.append(targetmc.GetBinLowEdge(b+1))

  # Here is the important bit which "Builds" the control region, make a list of control regions which 
  # are constraining this process, each "Channel" is created with ...
  #   (name,_wspace,out_ws,cid+'_'+model,TRANSFERFACTORS) 
  # the second and third arguments can be left unchanged, the others instead must be set
  # TRANSFERFACTORS are what is created above, eg WScales

  CRs = [
   Channel("qcd_singlemuon",_wspace,out_ws,cid+'_'+model,WScales, convention=convention),
   Channel("qcd_singleelectron",_wspace,out_ws,cid+'_'+model,WScales_e, convention=convention),
  ]


  # See https://docs.google.com/spreadsheets/d/15vq-c2xejGA-Nw6yzZU3mUDftter_l7OOcmJEwuCPyI/edit?usp=sharing
  if year == 2017:
    jes = 0.015
    jer = 0.015
  elif year==2018:
    jes = 0.01
    jer = 0.01
  for c in CRs:
    c.add_nuisance('CMS_VBF_scale_j'.format(YEAR=year), jes)
    c.add_nuisance('CMS_res_j_{YEAR}'.format(YEAR=year), jer)
    c.add_nuisance('CMS_veto{YEAR}_t'.format(YEAR=year),      0.01)
    c.add_nuisance('CMS_veto{YEAR}_m'.format(YEAR=year),      0.015)
    c.add_nuisance('CMS_veto{YEAR}_e'.format(YEAR=year),      0.03)


  # ############################ USER DEFINED ###########################################################
  # Add systematics in the following, for normalisations use name, relative size (0.01 --> 1%)
  # for shapes use add_nuisance_shape with (name,_fOut)
  # note, the code will LOOK for something called NOMINAL_name_Up and NOMINAL_name_Down, where NOMINAL=WScales.GetName()
  # these must be created and writted to the same dirctory as the nominal (fDir)
  do_stat_unc(WScales,proc='qcd_wmn', region='qcd_singlemuon', CR=CRs[0], cid=cid,outfile=_fOut)
  do_stat_unc(WScales_e,proc='qcd_wen', region='qcd_singleelectron', CR=CRs[1], cid=cid,outfile=_fOut)


  #######################################################################################################

  cat = Category(model,cid,nam,_fin,_fOut,_wspace,out_ws,_bins,metname,targetmc.GetName(),CRs,diag, convention=convention)
  cat.setDependant("qcd_zjets","qcd_wjetssignal")  # Can use this to state that the "BASE" of this is already dependant on another process
  # EG if the W->lv in signal is dependant on the Z->vv and then the W->mv is depenant on W->lv, then 
  # give the arguments model,channel name from the config which defines the Z->vv => W->lv map! 
  # Return of course
  return cat

