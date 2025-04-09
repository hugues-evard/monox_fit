- [1. Overview of what is done when running `generate_combine_model.py`:](#1-overview-of-what-is-done-when-running-generate_combine_modelpy)
  - [1.1. Addition of nuisances in `cmodel`](#11-addition-of-nuisances-in-cmodel)
    - [1.1.1. List of nuisances](#111-list-of-nuisances)
    - [1.1.2. Difference with monojet models](#112-difference-with-monojet-models)
  - [1.2. Overview of `init_channels`](#12-overview-of-init_channels)
  - [1.3. Overview of `convert_to_combine_workspace`](#13-overview-of-convert_to_combine_workspace)
- [2. Overview of the old model construction scripts](#2-overview-of-the-old-model-construction-scripts)
  - [2.1. General overview](#21-general-overview)
    - [2.1.1. Beginning of script](#211-beginning-of-script)
    - [2.1.2. inside `my_function` (applying theory variations):](#212-inside-my_function-applying-theory-variations)
    - [2.1.3. back to `cmodel`](#213-back-to-cmodel)
  - [2.2. Recap of `Z_constraints_qcd_withphoton.py`:](#22-recap-of-z_constraints_qcd_withphotonpy)
  - [2.3. Recap of `W_constraints.py`:](#23-recap-of-w_constraintspy)
  - [2.4. Recap of `Z_constraints_ewk_withphoton.py`:](#24-recap-of-z_constraints_ewk_withphotonpy)
  - [2.5. Recap of `W_constraints.py`:](#25-recap-of-w_constraintspy)
- [3. Overview of `Category`, `Channel` and `Bin` classes](#3-overview-of-category-channel-and-bin-classes)
  - [3.2. `Channel`](#32-channel)
  - [3.1. `Category`](#31-category)
  - [3.3. `Bin`](#33-bin)
- [4. Stepping through `init_channels` in more details](#4-stepping-through-init_channels-in-more-details)
- [5. Stepping through `convert_to_combine_workspace` in more details](#5-stepping-through-convert_to_combine_workspace-in-more-details)
- [6. What is a `RooParametricHist` ?](#6-what-is-a-rooparametrichist-)
  - [6.1. Member Variables:](#61-member-variables)
  - [6.2. Constructor Behavior](#62-constructor-behavior)
  - [6.3. Main Methods](#63-main-methods)
    - [6.3.1. `evaluate()`:](#631-evaluate)
    - [6.3.2. `evaluatePartial()`:](#632-evaluatepartial)
      - [6.3.2.1. Code:](#6321-code)
      - [6.3.2.2. Step-by-step:](#6322-step-by-step)
    - [6.3.3. `evaluateFull()`:](#633-evaluatefull)
      - [6.3.3.1. Code:](#6331-code)
      - [6.3.3.2. Step-by-step:](#6332-step-by-step)
    - [6.3.4. `evaluateMorphFunction(int bin)`:](#634-evaluatemorphfunctionint-bin)
      - [6.3.4.1. Code:](#6341-code)
      - [6.3.4.2. Step-by-step:](#6342-step-by-step)

# 1. Overview of what is done when running `generate_combine_model.py`:

This script is where all processes are expressed as a product of transfer factors, nuisances, and parametrized by the same process,
$Z^{\text{QCD}}_{\text{SR}}\to \nu\nu$.

There are three main function calls:
   1. `cmodel`, The model construction scripts
      - This is run for each "model", containing a "target process" 
         (e.g. $Z^{\text{QCD}}_{\text{SR}}\to \nu\nu$ or $W^{\text{EWK}}_{\text{SR}}\to l\nu$) and a list of "control processes".
      - For each control process, it computes a transfer factors, ratio of the yield of the target process and control process.
      - It creates the parameters for the relevant nuisances for each transfer factor 
         (veto, JES/JER, theory, and statistical uncertainties) in the workspace.
   2. `init_channels`:
      - makes the modeled yield of events in each bin as a function of $(Z^{\text{QCD}}_{\text{SR}}\to \nu\nu)$:
         - $(Z^{\text{QCD}}_{\text{SR}}\to \nu\nu) \times \frac{CR}{Z^{\text{QCD}}_{\text{SR}}\to \nu\nu} \times \Pi^{nuis}_{CR}{(1+nuis)}$ 
         for channels in `qcd_zjets` model
         - $(Z^{\text{QCD}}_{\text{SR}}\to \nu\nu) \times \frac{Z^{\text{EWK}}_{\text{SR}}\to \nu\nu}{Z^{\text{QCD}}_{\text{SR}}\to \nu\nu} \times \Pi^{nuis}_{Z^{\text{EWK}}_{\text{SR}}\to \nu\nu}{(1+nuis)} \times \frac{Z^{\text{EWK}}_{\text{diMuon CR}} \to ll}{Z^{\text{EWK}}_{\text{SR}}\to \nu\nu} \times \Pi^{nuis}_{Z^{\text{EWK}}_{\text{diMuon CR}} \to ll}{(1+nuis)}$ 
         for channels in `ewk_zjets` `qcd_wjets` and `ewk_wjets` models
   3. `convert_tocombine_workspace`:
      - stores the distribution of modeled yields above in `RooParametricHist`s and saves them to the workspace

## 1.1. Addition of nuisances in `cmodel`

These are the steps performed when adding a nuisances `nuis` for a given transfer factor:
   - We check if a parameter for this nuisance already exist in the workspace. If not, we create and import the following:
      - The nuisance parameter `nuis`, with its value at 0, and a range of +- 3
      - A constraints named `const_{nuis}`, a gaussian centered at 0 with a sigma of 1
   - A nuisance can either affect all bins uniformely or each bin with a different shape.
   - Each bin of the distribution affected by the nuisance is assigned a function:
      - `sys_function_{nuis}_cat_vbf_2018_{model}_ch_{channel}_bin_{b}`
      - The effect of the function is `nuis * size` is the same for each bin of the process.
         So for a given process, the effect of the nuisance is uniform across all bins. 
      - Different processes can be affect by the same nuisance `nuis`, and with different strenght `size`.
   - Shape nuisances are similar, but need the up and down variations computed beforehand:
      - When adding the shape nuisance, we first fetch the variations called `{channel}_weights_vbf_2018_{nuis}_(Up|Down)`
      - For each bin, the effect of function is 
         $\left(1 + \frac{\Delta(\text{variation up}, \text{variation down})}{2 \times \text{nominal}}\right)^{\pm \text{nuisance}} - 1$
      - The exponant is $+nuisance$ if $\Delta(\text{variation up}, \text{variation down}) > 0$,
         $-nuisance$ if $\Delta(\text{variation up}, \text{variation down}) < 0$,

### 1.1.1. List of nuisances 

The following nuisances are created for each model (for the vbf case):
   - veto:
      - `"CMS_veto2018_e "`
      - `"CMS_veto2018_m "`
      - `"CMS_veto2018_t"`
   - JES/JER (shape):
      - `"jer_2018"`
      - `"jesAbsolute"`
      - `"jesAbsolute_2018"`
      - `"jesBBEC1"`
      - `"jesBBEC1_2018"`
      - `"jesEC2"`
      - `"jesEC2_2018"`
      - `"jesFlavorQCD"`
      - `"jesHF"`
      - `"jesHF_2018"`
      - `"jesRelativeBal"`
      - `"jesRelativeSample_2018"`
   - Theory (shape): 
      - QCD, PDF: `{qcd_label}_{var[1]}_vbf` 
         - `qcd_label` is `"(Photon|ZnunuWJets)_{QCD|EWK}"`
         - `var[1]` is in `["facscale", "pdf", "renscale"]`
         - e.g. `"Photon_EWK_facscale_vbf"`
      - EWK: `{ewk_label}_vbf_bin{b}`
         - `ewk_label` is in `["qcd_ewk", "qcd_photon_ewk", "ewk_ewk", "ewkphoton_ewk"]`
   - Stastistical (shape): `"vbf_2018_stat_error_{region}_bin{b}"`

### 1.1.2. Difference with monojet models

In the monojet script, things are a little different
   - vetos uncertainties are a shape nuisance (read from some systematics file templats)
   - there are trigger shape nuisances
   - there are electron/photon id shape nuisances
   - there are photon scale shape nuisances
   - there are prefiring shape nuisances

## 1.2. Overview of `init_channels`

Once all of theses nuisances are created for all models (qcd z, ewk z, qcd w, ewk w), `init_channels` is ran on each model:

   - We extract one (any) histogram, whose shape will be used as a reference 
      (the mjj bin edges, should be the same for all distributions)
   - We go through every `Channel` in the model
      - For each bin, we create a `Bin` object. This contains:
         - The bin edges
         - some ID to link it to other channels in other models
         - The initial yield for that mjj bin for the target process of that model
         ( so either $Z^{\text{QCD}}_{\text{SR}}\to \nu\nu$, $Z^{\text{EWK}}_{\text{SR}}\to \nu\nu$,
            $W^{\text{QCD}}_{\text{SR}}\to l\nu$ or $W^{\text{EWK}}_{\text{SR}}\to l\nu$)
         - The transfer factor $\frac{\text{control process}}{\text{target process}}$ for that process in that bin.
            - This way, the control process can be expressed as target process $\times$ transfer factor, and parametrized by the target process
            - This is stored as a constant `RooRealVar` `sfactor_cat_vbf_2018_{model}_ch_{channel_name}_bin_{b}`
         - Model the expected number of events.
            - At initialization, two cases:
               - For `qcd_zjets` model, RooRealVar `model_mu_cat_vbf_2018_{model}_bin_{b}` whose value is set to the QCD Znunu (SR) yield, in a range from 0 to 3 times the yield.
               - For other models, it fetches `pmu_cat_vbf_2018_qcd_zjets_bin_{b}`,  (see below)
            - Make the product of all nuisances:
               - Fetch each nuisance `nuis` , and create a `delta` function that has the formula `1+nuis`
               `delta_cat_vbf_2018_{model}_ch_{channel}_bin_{b}_{nuis}`
               - Only nuisances affecting the bin are added (so theory EWK and statistical uncertainties, decorrelated by bin, each contribute one nuisance; contributions for the other bins are skipped).
            - make `pure_mu`: two cases:
               - For `qcd_zjets`, product of `model_mu_cat_vbf_2018_{model}_bin_{b}`, `sfactor_cat_vbf_2018_{model}_ch_{channel_name}_bin_{b}` and all nuisances
               - For all other categories, product of `pmu_cat_vbf_2018_qcd_zjets_bin_{b}`, `sfactor_cat_vbf_2018_{model}_ch_{channel_name}_bin_{b}` and all nuisances

            - In other words, we have, for each bin, 
            $(Z^{\text{QCD}}_{\text{SR}}\to \nu\nu) \times \frac{CR}{Z^{\text{QCD}}_{\text{SR}}\to \nu\nu} \times \Pi^{nuis}_{CR}{(1+nuis)}$ 
            for each CR in the `qcd_zjets` model:
               - $Z^{\text{QCD}}_{\text{diMuon CR}} \to ll$
               - $Z^{\text{QCD}}_{\text{diElectron CR}} \to ll$
               - $W^{\text{QCD}}_{\text{SR}} \to l\nu$
               - $(\gamma + \text{jets})^{\text{QCD}}_{\text{SR}}$
               - $Z^{\text{EWK}}_{\text{SR}} \to \nu\nu$
            - For the other models, we fetch the previous expression from the channel the are linked to, and multiply to
            $\frac{CR}{target} \times \Pi^{nuis}_{CR}{(1+nuis)}$ 
            - For instance, in the `ewk_zjets` model, the target is $Z^{\text{EWK}}_{\text{SR}}\to \nu\nu$, we are linked to the corresponding CR in `qcd_zjets`. If we look for instance at the $Z^{\text{EWK}}_{\text{diMuon CR}} \to ll$ CR, this ends up with
            $(Z^{\text{QCD}}_{\text{SR}}\to \nu\nu) \times \frac{Z^{\text{EWK}}_{\text{SR}}\to \nu\nu}{Z^{\text{QCD}}_{\text{SR}}\to \nu\nu} \times \Pi^{nuis}_{Z^{\text{EWK}}_{\text{SR}}\to \nu\nu}{(1+nuis)} \times \frac{Z^{\text{EWK}}_{\text{diMuon CR}} \to ll}{Z^{\text{EWK}}_{\text{SR}}\to \nu\nu} \times \Pi^{nuis}_{Z^{\text{EWK}}_{\text{diMuon CR}} \to ll}{(1+nuis)}$ 
            - This is save as `pmu_cat_vbf_2018_{model}_bin_{b}`, and also wrapped in the `RooFormulaVar` `mu_cat_vbf_2018_{model}_bin_{b}`
            - The `"observed"` is fetched, and a Poisson PDF is constructed for `observed` using `mu_cat_vbf_2018_{model}_bin_{b}`.
               It is uncleared where `observed` comes from
         - Once this modelling is done for all bins, save all prefit distributions
         - The last step is unclear, maybe a check that all expected distribution exist.
## 1.3. Overview of `convert_to_combine_workspace`
This is where the all distribution used for combine are saved to the workspace.

   - We extract one (any) histogram, whose shape will be used as a reference (`samplehist`)
      (the mjj bin edges, should be the same for all distributions)
   - We extract `mjj`, rename to `mjj_vbf_2018`
   - Convert every histogram `hist` from the input file (`limit_vbf.root`) to a `RooDataHist` called `vbf_2018_{hist}`, as a function of `mjj_vbf_2018` and add it to the workspace
   - Loop over every `Category` (as in "model", qcd zjets, qcd wjets, ewk zjets, ewk wjets):
      - fetch `model_mu_cat_vbf_2018_{model}_bin_{b}` for every bin
      - Only for `qcd_zjets`  model, create the signal distribution:
         - Create a `RooParametricHist` with name `vbf_2018_signal_qcd_zjets_model`, holding the distribution of all
            `model_mu_cat_vbf_2018_{model}_bin_{b}` for every bin as a function of `mjj_vbf_2018` (with shape given by `samplehist`)
         - Create addition of expectations `vbf_2018_signal_qcd_zjets_model_norm` (integral of previous `RooParametricHist` over the whole `mjj` range)
         - import these in the workspace
      - Loop for every `Channel` in the model to get all backgrounds:
         - Create a `RooParametricHist` with name `vbf_2018_{channel}_{model}_model`, holding the distribution of all
            `pmu_cat_vbf_2018_ch_{model}_bin_{b}` for every bin as a function of `mjj_vbf_2018` (with shape given by `samplehist`)
         - Create addition of expectations `vbf_2018_{channel}_{model}_model_norm` (integral of previous `RooParametricHist` over the whole `mjj` range)
         - import these in the workspace
   - Get all parameters in the workspace
   - for all background nuisances: print the line to add in the datacard template, `{param.GetName()} param {param.getVal()} 1`

# 2. Overview of the old model construction scripts

 some notes of what is being done in each model construction script for vbf:
   - Z_constraints_qcd_withphoton.py
   - W_constraints_qcd.py
   - Z_constraints_ewk_withphoton.py
   - W_constraints_ewk.py

## 2.1. General overview

### 2.1.1. Beginning of script

Initializing the inputs:
   - Target process, the one that is used to parametrize the control samples
   - Control MC samples, samples to derive the transfer factor and apply nuisances on.

Compute all transfer factors: for each control sample, compute target divided by control

### 2.1.2. inside `my_function` (applying theory variations):
This is only done for the Z models (EWK and QCD).

   - Getting theory uncertainties from file `sys/vbf_z_w_gjets_theory_unc_ratio_unc.root`

   - Compute $\frac{Z_\text{SR}\to \nu\nu}{W_\text{SR}\to l\nu} \times$ theory uncertainty for QCD and PDF uncertainties

   - For EWK uncertainties, decorrelated among bins:
      - Create one clone of $\frac{Z_\text{SR}\to \nu\nu}{W_\text{SR}\to l\nu}$ for each bin
      - Each clone i get only bin i replace with bin i of $\frac{Z_\text{SR}\to \nu\nu}{W_\text{SR}\to l\nu} \times $ EWK uncertainty

   - Same computations for $\frac{Z_\text{SR}\to \nu\nu}{{\gamma + \text{jets}}_\text{SR}}$:
      - QCD
      - PDF
      - EWK, decorrelated among bins

### 2.1.3. back to `cmodel`

   - extract binning of mjj

   - Create one `Channel` by transfer factor

   - For relevant channels,  add nuisances for vetos: `RooFormulaVar` corresponding to vetoname * value
      - `CMS_veto{YEAR}_t * {value veto t}`
      - `CMS_veto{YEAR}_m * {value veto m}`
      - `CMS_veto{YEAR}_e * {value veto e}`

   - Extract JER and JES uncertainties for transfer factors from `sys/vbf_jes_jer_tf_uncs.root`:
      - for each variation 
      - for all transfer factor (except Znunu QCD / EWK):
      - Add histogram for Transfer factor  * variation uncertainty in output
         (variation uncertainty for transfer factor is a histogram with one bin,
         and the whole transfer factor distribution get Scaled according to the value in this one bin)
      - Add function (quadratic) to model systematic uncertainty on transfer factor

   - Add variation from statistic uncertainty + quadratic model

   - Add variation from theory uncertainties + quadratic model (QCD, PDF, EWK decorellated, both for Z/W and Z/gamma+jet)

   - Return everything as `Category` 

## 2.2. Recap of `Z_constraints_qcd_withphoton.py`:

| Sample name | transfer factor                                                               | Veto | JES/JER | Theory | Stat |
| ----------- | ----------------------------------------------------------------------------- | ---- | ------- | ------ | ---- |
| `Zmm`       | $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{Z^\text{QCD}_{2\mu}\to ll}$          |      | X       |        | X    |
| `Zee`       | $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{Z^\text{QCD}_{2e}\to ll}$            |      | X       |        | X    |
| `WZ`        | $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{W^\text{QCD}_{SR}\to l\nu}$          | X    | X       | X      | X    |
| `EQ`        | $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{Z^\text{EWK}_{SR}\to \nu\nu}$        |      | X       |        | X    |
| `Photon`    | $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{{\gamma+jets}^\text{QCD}_{1\gamma}}$ |      | X       | X      | X    |

Return everything as `Category` called `qcd_zjets`

## 2.3. Recap of `W_constraints.py`:

| Sample name | transfer factor                                                      | Veto | JES/JER | Theory | Stat |
| ----------- | -------------------------------------------------------------------- | ---- | ------- | ------ | ---- |
| `Wmu`       | $\frac{W^\text{QCD}_\text{SR}\to l\nu}{W^\text{QCD}_{1\mu}\to l\nu}$ | X    | X       |        | X    |
| `We`        | $\frac{W^\text{QCD}_\text{SR}\to l\nu}{W^\text{QCD}_{1e}\to l\nu}$   | X    | X       |        | X    |

Return everything as `Category` called `qcd_wjets`,  specifying it is dependant on `WZ` transfer factor of previous `qcd_zjets` category.

Later on, both transfer factor above will be multiplied by $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{W^\text{QCD}_{SR}\to l\nu}$
to parametrize it to $Z^\text{QCD}_\text{SR}\to \nu\nu$

## 2.4. Recap of `Z_constraints_ewk_withphoton.py`:


| Sample name | transfer factor                                                               | Veto | JES/JER | Theory | Stat |
| ----------- | ----------------------------------------------------------------------------- | ---- | ------- | ------ | ---- |
| `Zmm`       | $\frac{Z^\text{EWK}_\text{SR}\to \nu\nu}{Z^\text{EWK}_{2\mu}\to ll}$          |      | X       |        | X    |
| `Zee`       | $\frac{Z^\text{EWK}_\text{SR}\to \nu\nu}{Z^\text{EWK}_{2e}\to ll}$            |      | X       |        | X    |
| `WZ`        | $\frac{Z^\text{EWK}_\text{SR}\to \nu\nu}{W^\text{EWK}_{SR}\to l\nu}$          | X    | X       | X      | X    |
| `Photon`    | $\frac{Z^\text{EWK}_\text{SR}\to \nu\nu}{{\gamma+jets}^\text{EWK}_{1\gamma}}$ |      | X       | X      | X    |

Return everything as `Category` called `ewk_zjets`, specifying it is dependant on `EQ` transfer factor of previous `qcd_zjets` category.

Later on, both transfer factor above will be multiplied by $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{Z^\text{EWK}_{SR}\to \nu\nu}$
to parametrize it to $Z^\text{QCD}_\text{SR}\to \nu\nu$

## 2.5. Recap of `W_constraints.py`:

| Sample name | transfer factor                                                      | Veto | JES/JER | Theory | Stat |
| ----------- | -------------------------------------------------------------------- | ---- | ------- | ------ | ---- |
| `Wmu`       | $\frac{W^\text{EWK}_\text{SR}\to l\nu}{W^\text{EWK}_{1\mu}\to l\nu}$ | X    | X       |        | X    |
| `We`        | $\frac{W^\text{EWK}_\text{SR}\to l\nu}{W^\text{EWK}_{1e}\to l\nu}$   | X    | X       |        | X    |

Return everything as `Category` called `ewk_wjets`, specifying it is dependant on `WZ` transfer factor of previous `ewk_zjets` category.

Later on, both transfer factor above will be multiplied by $\frac{Z^\text{EWK}_\text{SR}\to \nu\nu}{W^\text{EWK}_{SR}\to l\nu}$
which is itself multiplied by $\frac{Z^\text{QCD}_\text{SR}\to \nu\nu}{Z^\text{EWK}_{SR}\to \nu\nu}$
to parametrize it to $Z^\text{QCD}_\text{SR}\to \nu\nu$

# 3. Overview of `Category`, `Channel` and `Bin` classes

These classes are used at the different stages of `generate_combine_model.py`
to store transfer factors, nuisances, and build the distributions
of the different process parametrized by the $Z^\text{QCD}_\text{SR}\to \nu\nu$ yield.

## 3.2. `Channel`

A `Channel` object holds:
   - A transfer factor (ratio of yields from two different processes)
   - an input and an output `RooWorkspace` containting nuisances
   - A name, and some IDs to link to `Category` and `Bin`

And two methods are defined to add:
   - a nuisances that applies on all bins (normalization)
   - a shape nuisance, where the effect 
      of the nuisance in each bin is modeled with a function (lognormal) that depends 
      on the up and down variations the nuisance is derived from.

All channels are created in the different model construction scripts,
which also add the relevent nuisances to them 
(for veto, JES/JER, theory systematics as well as statistical uncertainties).

## 3.1. `Category`

A `Category` object holds:
   - A collection of `Channel`s
   - A name, and some IDs to link to `Category` and `Bin`
   - Optionally, a dependance on `Channel` from another `Category`

It's main method is `init_channels`, 
which is used to build the distributions of all the different process,
parametrized by the $Z^\text{QCD}_\text{SR}\to \nu\nu$ yield and all nuisances affecting each process.
This is done by first creating a collection of `Bin` objects that handle parametrizing the number of events for each bin separately.


## 3.3. `Bin`

A `Bin` object holds:
   - A name, and some IDs to link to `Category` and `Bin`

This is the class that handles, inside the `setup_expect_var` method, fetching the transfer factor, yield of $Z^\text{QCD}_\text{SR}\to \nu\nu$,
and nuisances for a given process and make their product to make a parametrized distribution.

# 4. Stepping through `init_channels` in more details

Loop over categories, do `Category.init_channels()`

Loop over control regions (`Channel`s)

Loop over each bin:
   - extract the edges
   - initialize `Bin` object, link it to control region (using id)
   - set categoryname "bin_number"
   - initialize Y value of bin based on target dataset of this control region:
      ``` 
      INIT Y:  
      mjj>=200.0 && mjj<400.0
      self.rngename = rnge_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0 (a name given to the var in setRange, for later acces)
      self.wspace = Name: wspace_vbf_2018 Title: wspace_vbf_2018
      self.wspace.data(mcdataset) = Name: signal_qcdzjets Title: DataSet - vbf_2018, signal_qcdzjets
      mcdataset = signal_qcdzjets
      ```
   - In other words, value of bin is initialized to mjj yield of category (qcd z, qcd w, ewk z or ewk w)
   - Get "scale factor" (rather transfer factor) for that control region: 
      - given by 1 / cr.scalefactor
      - following exemple in qcd_zjets category, qcd_dimuon channel 
      - transfer factor is qcd Znunu (SR) / qcd Zll (dimuon CR)
      - transfer factor = 20087 / 2097 = 9.577
      - ret_scalefactor gives 1/9.577 = 0.104 (so qcd Zll (dimuon) / qcd Znunu (SR))
   - Set "scale factor":
      - `Bin.sfactor` is a RooRealVar `sfactor_cat_vbf_2018_{model}_ch_{channel_name}_bin_{b}` holding previous value (for qcd zjets, qcd dimuon, 0.104)
      - no range( +-INF), constant
      - TODO: why the use of `ROOT.RooFit.RecycleConflictNodes()`?
   - "Setup expected var", this is where dependences are handled
      - Case of no dependence:
         - initialize `model_mu`, a RooRealVar `model_mu_cat_vbf_2018_{model}_bin_{b}`
         - value: yield of target process (eg qcd Znunu SR)
         - range: from 0 to 3 times yield
      - Case of a dependence:
         - initialize `model_mu` as `pure_mu` from the channel it depends on (number of expected event in the channel it depends on, see bellow)
      - arglist: `model_mu` and `sfactor`
      - fetch nuisances (JES/JER, veto, theory, stat)
      - initialize empty `RooArgList` for nuisances `nuis_args`
      - for each nuisance:
         - get formula of nuisance in that bin, put it in a new `RooArgList`
         - get "`delta`", `RooFormulaVar` giving `1+nuisance`
         - add `delta` to `nuis_args`
         - write `delta` to the workspack
         - in systematic has no effect on the bin (sf up - sf down = 0), in `Channel.add_nuisance_shape`, nuisance is marked as `"temp"` and will not
            be skipped here, so not added to the `nuis_args`. This is likely the case for nuisances that are decorrelated by bins (stat and ewk theory) for bins they don't act on
         - resulting list of nuisances for qcd_dimuon channel in qcd_zjets:
         ```
         RooFormulaVar::delta_cat_vbf_2018_qchttps://root.cern.ch/doc/master/namespaceTMath.html#a9cef8a9eac06254be610bd64b03106d0d_zjets_ch_qcd_dimuon_bin_0_jesBBEC1_2018[ actualVars=(sys_function_jesBBEC1_2018_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesRelativeSample_2018[ actualVars=(sys_function_jesRelativeSample_2018_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesBBEC1[ actualVars=(sys_function_jesBBEC1_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesFlavorQCD[ actualVars=(sys_function_jesFlavorQCD_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesRelativeBal[ actualVars=(sys_function_jesRelativeBal_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jer_2018[ actualVars=(sys_function_jer_2018_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesEC2_2018[ actualVars=(sys_function_jesEC2_2018_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesHF_2018[ actualVars=(sys_function_jesHF_2018_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesAbsolute_2018[ actualVars=(sys_function_jesAbsolute_2018_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesEC2[ actualVars=(sys_function_jesEC2_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesHF[ actualVars=(sys_function_jesHF_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_jesAbsolute[ actualVars=(sys_function_jesAbsolute_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1
         RooFormulaVar::delta_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0_vbf_2018_stat_error_qcd_dimuonCR_bin0[ actualVars=(sys_function_vbf_2018_stat_error_qcd_dimuonCR_bin0_cat_vbf_2018_qcd_zjets_ch_qcd_dimuon_bin_0) formula="1+@0" ] = 1 
         ```
      - make product of deltas
      - add product to `arglist` (so now arglist is model_mu (yield of qcd znunu), sfactor (qcd znunu / qcd zmumu) and product of delta)
      - "`pure_mu`" is a formula for number of expected (signal) events in the bin, given by product of everything in arglist, so yield * transfer factor * effect of nuisances
         - named `pmu_cat_vbf_2018_{model}_ch_{channel}_bin_{b}`
      - "`mu`" is the same ? Number of expected events in the bin, but I don't see any difference with `pure_mu`
      -  add `mu` to the workspace (adds dependencies as well)
      - note about the dependant case: `pure_mu` will can be re-used by other channels. For instance all channels in the `Category` qcd_wjets depend on the qcd_w channel (qcd znunu in SR / qcd wlnu in SR) in qcd_zjets.
         In that case, for instance for the qcd_wmunu channel (qcd wln SR / qcd wln (singlemuon CR)) will initialize its `arglist` with the `pure_mu` of qcd_w. we and up with:
            pure_mu(qcd_w) * sfactor(qcd_wmunu) = (yield (= qcd znunu) * sfactor (= qcd wlnu / qcd znunu) * product of nuisances in qcd_w ) * sfactor (= qcd wmunu / qcd wlnu)
            that finaly can be simplified to yield of qcd wmunu * product of nuisance in qcd_w

            To this we will further add the product of nuisances in qcd_wmunu
      
      - get `observed` (not sure where this is passed)
      - build [poisson PDF](https://root.cern.ch/doc/master/namespaceTMath.html#a9cef8a9eac06254be610bd64b03106d0) from `observed` and `mu`
   - set_initE: returns 0
   - add_to_dataset: does nothing
- Save prefit distributions
- TODO: Last step is unclear, is it a check that all expected distribution exist?

# 5. Stepping through `convert_to_combine_workspace` in more details

   - open the workspace, 
   - For every "category" (as in the year, here `vbf_2018`)
      - open the input directory containing the histograms
      - extract the binning and varname by fetching the first histogram it can find in the directory
         - Varname is taken from histogram.GetXaxis().GetTitle(). This is why all histograms in input file need to have title set to `"mjj"`
      - extract mjj from workspace (RooRealVar `varl`, vary from 200 to 6000, lowest and highest edges in distribution)
      - rename variable to mjj_vbf_2018 
      - Loop over all histograms in directory
         - if integral is 0, set first bin to something small to avoid errors in combine TODO is this needed?
         - Create RooDataHist `vbf_2018_{histogram.GetName()}` holding the distribution of the histogram along one dimension for mjj_vbf_2018
         - Import the RooDataHist in combine workspace
      - for every `Category` (as in "model", qcd zjets, qcd wjets, ewk zjets, ewk wjets):
         - import the model
         - initialize RooArgList of expectations
         - fetch `model_mu_cat_vbf_2018_{model}_bin_{b}` for every bin, add it to list of expectations.
         - Only for `qcd_zjets`  model:
            - **This if for the signal**
            - Create phist, [RooParametricHist](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/interface/RooParametricHist.h)
            (custom class inheriting from RooAbsPdf, see also in `src/HiggsAnalysis/CombineLimits`) with name `vbf_2018_signal_qcd_zjets_model`, RooRealVar `varl` (_x), RooArgList `expectations` (_pars), TH1 `samplehist` (_shape)
               - Histogram that holds the distribution of `expectations` (the expected number of events), with the binning provided by `samplehist` (only used for its shape),
                  and will evaluate the distribution based on the value of `varl` (in this case `mjj_vbf_2018`), it will look on the relevent bin of mjj, and return yield / bin width * nuisances
            - crate addition of expectations `vbf_2018_signal_qcd_zjets_model_norm`
            - import these in the workspace
         - loop again for every `Category` (TODO: is this really needed)
            - check that model name matches the one imported previously, else continue
            - for every `Channel` in the model, create the background distributions:
               - Fetch every `pmu_cat_vbf_2018_{model}_ch_{channel}_bin_{b}`
               - for every bin in the channel:
                  - add expectation of bin to list of expectation of the channel, using the "`pure_mu`" (from `Bin.setup_expect_var`)
               - We obtain list of `pure_mu` expectation for every bin in that channel
               - create p_phist, RooParametricHist with name `vbf_2018_{channel}_{model}_model`, RooRealVar `varl` (`mjj_vbf_2018`, `_x`), RooArgList `cr_expectations` (`_pars`) and `_shape` `samplehist`
               - create addition of expectations `vbf_2018_{channel}_{model}_model_norm`
               - import these in the workspace
   - We are done now out of  every loop
   - Get all parameters in the workspace
   - for all background nuisances: print the line to add in the datacard template, `{param.GetName()} param {param.getVal()} 1`

      
# 6. What is a `RooParametricHist` ?

The `RooParametricHist` class is a **custom RooFit probability density function (PDF)**. It models a **binned histogram-like PDF**, where the bin contents are **free parameters**.

---

The class defines a PDF where:
- The shape is determined by a histogram (`TH1`) with a number of bins.
- Each bin has an associated **parameter** that represents its content (normalized).
- The PDF is evaluated based on which bin the observable `x` falls into.
- Optional: it can **morph** bin values based on nuisance parameters, used to model systematic uncertainties.

---

## 6.1. Member Variables:
- `x`: The observable (in our case, $m_{jj}$).
- `pars`: A list of parameters, each representing the bin content (as a `RooAbsReal`).
- `bins`, `widths`: Vectors holding the bin edges and widths from the input histogram.
- `_coeffList`: Morphing coefficients (nuisance parameters).
- `_sums`, `_diffs`: Data structures for storing the shape variation due to morphing.
- `_hasMorphs`: Flag to enable morphing.
- `_cval`: Cached last value of the PDF.

---

## 6.2. Constructor Behavior
```cpp
RooParametricHist::RooParametricHist(const char *name, const char *title, RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape)
```
- Sets up the histogram-like structure using the shape of `TH1`.
- Validates that number of parameters matches the number of bins.
- Initializes bin edges and widths.

---

## 6.3. Main Methods

### 6.3.1. `evaluate()`:
- Main function called by RooFit to evaluate the PDF.
- Calls either:
  - `evaluateFull()` if morphing is enabled.
  - `evaluatePartial()` otherwise.
- Normalizes by bin width.
- Returns the PDF value at current `x`.

---

### 6.3.2. `evaluatePartial()`:
This function is used when **morphing is not enabled**.
It finds the bin where the current observable value `x` lies, retrieves the corresponding bin parameter (from `pars`),
and returns the normalized value for that bin.

#### 6.3.2.1. Code:
```cpp
double RooParametricHist::evaluatePartial() const
{
  auto it = std::upper_bound(std::begin(bins), std::end(bins), x);
  if ( it == std::begin(bins) ) {
    // underflow
    return 0;
  }
  else if ( it == std::end(bins) ) {
    // overflow
    return 0;
  }
  size_t bin_i = std::distance(std::begin(bins), it) - 1;
  RooAbsReal *retVar = (RooAbsReal*)pars.at(bin_i);

  double ret = retVar->getVal();
  ret /= widths[bin_i];
  return ret;
}
```

#### 6.3.2.2. Step-by-step:
1. **Find bin:** Uses `std::upper_bound` to find the first bin edge that is *greater than* the current value of `x`. This tells us which bin `x` falls into.
2. **Check edge cases:**
   - If `x` is less than the first bin edge → return 0 (underflow).
   - If `x` is greater than the last bin edge → return 0 (overflow).
3. **Get bin index:** Subtract 1 to get the actual bin index.
4. **Get parameter:** Get the corresponding parameter for that bin: `pars.at(bin_i)`.
5. **Normalize by bin width:** Divide by the bin width to get the probability density.
6. **Return result.**

---

### 6.3.3. `evaluateFull()`:
Used when **morphing is enabled**. This is similar to `evaluatePartial`, but it applies a **morphing correction factor** to the bin content.

#### 6.3.3.1. Code:
```cpp
double RooParametricHist::evaluateFull() const
{
  int bin_i;
  if (x < bins[0]) return 0;
  else if (x >= bins[N_bins]) return 0;

  else {
    for (bin_i=0; bin_i<N_bins; bin_i++) {
      if (x >= bins[bin_i] && x < bins[bin_i+1]) break;
    }
  }
  double mVar = evaluateMorphFunction(bin_i);
  RooAbsReal *retVar = (RooAbsReal*)pars.at(bin_i);

  double ret = retVar->getVal() * mVar;
  ret /= widths[bin_i];
  return ret;
}
```

#### 6.3.3.2. Step-by-step:
1. **Check range:**
   - If `x` is before the first bin → return 0.
   - If `x` is past the last bin → return 0.
2. **Find bin index:** Iterate through bins to find the one where `x` falls between `bins[bin_i]` and `bins[bin_i+1]`.
3. **Apply morphing:** Call `evaluateMorphFunction(bin_i)` to get a scale factor for this bin.
4. **Get parameter:** Get the corresponding parameter from `pars`.
5. **Multiply and normalize:** Multiply the parameter by the morph scale, then divide by bin width.
6. **Return result.**

---

### 6.3.4. `evaluateMorphFunction(int bin)`:
- Computes a scale factor to apply to a bin’s parameter based on current values of morph parameters (`_coeffList`).
- Uses `_diffs` and `_sums` to smoothly morph the bin content.

#### 6.3.4.1. Code:
```cpp
double RooParametricHist::evaluateMorphFunction(int j) const
{
    double scale = 1.0;
    if (!_hasMorphs) return scale;

    int ndim = _coeffList.getSize();
    double f0 = static_cast<RooAbsReal*>(pars.at(j))->getVal();

    for (int i = 0; i < ndim; ++i) {
        double x = (dynamic_cast<RooRealVar*>(_coeffList.at(i)))->getVal();
        double a = 0.5 * x;
        double b = smoothStepFunc(x);  // smooth transition
        scale *= 1 + (1.0 / f0) * a * (_diffs[j][i] + b * _sums[j][i]);
    }

    return scale;
}
```

#### 6.3.4.2. Step-by-step:
1. **Skip if no morphing:** If morphing isn’t enabled, return scale = 1 (no effect).
2. **Loop over morphing parameters (`_coeffList`):**
   - Get the current value of each nuisance parameter `x`.
   - Compute a morphing term using `_diffs` and `_sums`, which represent up/down variations in the shape.
   - Multiply `scale` by a factor that adjusts the bin value.
3. **Return scale factor** that can be applied to the bin.
