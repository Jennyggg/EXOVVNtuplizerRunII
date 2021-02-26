#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Fall17" # for MC
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2017B" # for Data from 2017
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2018B" # for Data from 2018
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Autumn18" # for MC from 2018
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Summer16" # for MC from 2016


###### Process initialization ##########

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')


#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.TFileService = cms.Service("TFileService",
#                                    fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/debug.root')
                                    fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/flatTuple_small2_trig.root')
#                                    fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/flatTuple_small1_trig_local.root')
#                                     fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/ZeroBias2018/debug.root')
#                                     fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/MinBiasMC2018_NoPU/flatTuple_small.root')
#                                     fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPHMC/flatTuple.root')
#                                     fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/MinBiasMC/flatTuple.root')
#                                     fileName = cms.string('file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/HardQCD2018/flatTuple_small.root')
                                   )

#from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_data_cfi import config
from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_generic_cfi import config

# change from its original value
#config["DZCUT"] = 0.25
#config["FSIGCUT"] = 3
#config["VPROBCUT"] = 0.1
#config["DNNCUT"] = 0.2
#config["RUNONMC"]=False
####### Config parser ##########

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

options.register( 'RunPeriod',
                  '',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "RunNumber (Default Run2017B)")

options.register( 'runUpToEarlyF',
                  'false',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list                                                                                                                                 
                  VarParsing.VarParsing.varType.bool,          # string, int, or float                                                                                                                        
                  "false")# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile

####
#options.maxEvents = 100
options.maxEvents = -1

#data file

#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v4/230000/322CC3C2-921E-7448-902A-8FCB1F0A2F72.root'
#options.inputFiles = '/store/user/manzoni/BcToJpsiX_TuneCP5_13TeV-pythia8/RunIISummer19UL18_MINIAODSIM_v1/201112_133556/0000/RJpsi-BcToXToJpsiMuMuSelected-RunIISummer19UL18MiniAOD_7.root'

#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/100000/78A9DF86-EAC5-D242-8E7B-171AFD012CC7.root'

#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/OniaAndX_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/01325465-A815-E24E-ABB3-DAB8D4880BDE.root'
#options.inputFiles = 'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/005650EB-AA60-D34C-8100-95959C4B352D.root'
#options.inputFiles = ['file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_16.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_17.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_18.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_19.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_20.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_21.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_22.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_23.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_24.root',
#'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPH/Run2018A_25.root']
#options.inputFiles = 'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/ZeroBias2018/013404B1-9437-A740-9C7A-F3EBE718EFC2.root'
#options.inputFiles = 'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/MinBiasMC/8F79C765-3C14-2641-8737-9AD8980241B8.root'
#options.inputFiles = 'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/BPHMC/0F2BA15C-2CB8-E24C-B58C-AA59AAFA2B07.root'
#options.inputFiles = 'file:///pnfs/psi.ch/cms/trivcat/store/user/jinw/HardQCD2018/0132101A-3696-6D48-A7AC-2B6B1D2E6139.root'
#options.inputFiles = ['/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0109D8A3-FF62-6245-9745-675EE9FD6243.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0187FB10-EA3E-ED4F-8340-8931B659E533.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/01C007B5-2B48-3543-A338-22272E441925.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/022E3F41-8CD9-5246-96F6-5BA07B4180AF.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/028BA34E-416C-E94E-B976-5615D158D99A.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/03A6084D-770F-EA40-8500-126262842B17.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/03CBF487-1DAA-2345-BB25-05EB7244867B.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0447AB7C-9EDC-E046-B41D-F68DE3A1EAB0.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/048E57B6-A15E-D645-8731-4A1F1B27786C.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0923C159-3544-3649-9424-694E69CE855F.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0A4F96E2-8DC1-2240-A606-C8027B7DDF0C.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0AAE8C6F-1197-844C-B10D-34DDD4847F3A.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0B14F6F4-5DD5-B04B-A7EF-7D2B80AD33FA.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0B3FE33C-5B87-974B-A80B-2783E45886AD.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0BC38BDB-A396-FF45-A1B9-EFD9DA98D316.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0C14C0CC-8733-2242-91F3-94C55DBAF8DD.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0DF8969A-A985-4648-808D-CD0A1AD7367F.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0F441A41-FA10-8A4C-B07D-9EC21B82A5E0.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0F544919-1564-2446-81B8-5CF792CAD956.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0F882E56-C6E5-6A49-AA40-FAAEB621F675.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/1189C7B3-4874-704C-A485-F4EBDC21784D.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/11BAF64E-7D44-8B49-BAF3-2AD6AB1CEF79.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/12208844-ED2A-5A4C-B127-4C0444B899A4.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/12280ABA-1297-EA4A-8C59-C429FB0693A9.root']
#options.inputFiles = ['/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/12AC9C66-D582-1D4F-A39F-0BDBFDFCD02A.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/12C659FD-D961-6640-9C3D-48CD12A0D033.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/16216319-7A5E-774C-AAE5-5E36ED1E724E.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/17214445-47CE-B740-837B-01F1CEF607F7.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/18244BBF-8FF4-F847-8E6C-CAC64492F34E.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/18529C8A-BF7E-FB4D-8BE5-481EE498ECFF.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/185BA78A-4197-FC49-9353-845F717FC0C0.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/19F53134-F58E-584F-9471-C290C95F55C1.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/1AFB234E-4D54-CA40-A374-D7D3848E2E7B.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/1B9C68C4-67A2-CB4F-8BCE-9638EC7A2E31.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/1D796FDD-1879-3441-BB24-C0A61C33B5D7.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/1EEBDB7D-3C35-1045-8271-B48535954774.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/200AE1E6-4D9B-F84F-938B-1EEBBC5036E7.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/201AB255-6FB8-6945-B0BC-BC8BD52A04BF.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2129F080-982D-A649-8D20-4944620E99A6.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2184D704-2916-0E4A-913B-DF794BAF22F1.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/22A2D0DA-3725-954A-8A3E-F3286E8CC5D7.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/27098DE0-F01F-EF46-AB9C-E97521F1E994.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/272E571D-8FC9-4D44-B714-801092498538.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/280C18A9-E24D-ED47-B709-F29231A6D378.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/284E9437-65FE-2A44-8FB4-A60E4473AF20.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/28A1A007-F49D-F644-BA2C-85D76DFAA199.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2988A4FE-0604-C847-8C5B-3B575522CDEC.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/29D9A261-AF66-664C-BF86-0D6E34FD3DDE.root',
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2B59D9F7-F59A-9F45-A6C8-992D2D0D0F58.root']
options.inputFiles = ['/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2D654358-78BB-CC47-9BDB-B06F7A921AA4.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2E04833B-1907-D042-AFDB-6C90C998D4CD.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2F1325AF-4A4A-5048-8B66-458A4CD29BCE.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/30D0485A-0BC6-AA45-8248-7F511B0BAA40.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/3173E522-2B70-4447-9F97-86509F860887.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/324CFB6F-04C5-CA43-8FC4-7D71865ACC74.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/32CCB237-4000-0A40-BA5F-7944298C512E.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/33C6D0EF-F07A-9140-AC0D-5BEB78BB5620.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/34A41A95-F0BE-A848-AA40-A34FE8C02957.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/34C995F1-FA2A-C941-A149-6B63A5194E2A.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/357B8A1D-6B20-3846-9FB1-ECD2C0E70281.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/36410479-7E3C-BE40-B625-F0AC4BD555F7.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/3644E20B-EB2B-DD42-B7BE-5D1ADBA3F62D.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/36FCFFBA-9557-CE4C-95A5-A3A3A3B11A85.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/386D8AA5-CD18-D047-9609-AE40C008F3A8.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/3C22B72D-B9A5-3A47-A4BA-EF8D8AA35B31.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/3D6A89D1-76A4-454B-B01B-D2E5D40BAF7B.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/3DC147F2-0677-0B4C-A5F2-4F8BEB140F0F.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/3EDEDE1C-B2CC-DD40-B61F-CE3C81DB8EF8.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/3F1DC056-C085-FA4F-A6CA-9C343A5791BA.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/45575733-8518-7247-8D64-15A887EF57C7.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/46B00DD1-0FF1-134D-8AD2-0AD045BE3159.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/474DF484-D863-B142-BD9D-91541BA36FEA.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/477710C5-6B9B-D348-90F8-6DAB3D30A1D5.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/489A28C0-68F3-6A41-AC4D-2AAF77FC21AD.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/489C65AD-415E-E840-A912-E11E54ECCADB.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/4938462C-59AB-FD4B-8076-80B84EDBB29C.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/4ABD886C-978A-8A46-BFA5-9F8FD5DD615F.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/4B49014B-C005-584B-A0E0-0EC01E8A7851.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/4C3F6DD5-1AD3-0444-80F8-830DBA773403.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/4CBC3552-C55C-424A-904C-EB464F080B89.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/4EED88BB-26E0-4044-977B-956B049AB332.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/4F6D3185-6655-154A-93AE-9A71367B7969.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/502CD439-6A19-1C44-8DBA-75BD6ACCDF3E.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/50889893-7D46-9241-881B-6D8AAFD16328.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/50B48698-91FB-EF46-B78F-31D2D5CF7AB0.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/510A8701-EA1F-E24F-B457-276166E72D04.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/519BCDC5-3D0C-EC41-8EF6-E8C92CDA168D.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/51BA0779-C792-A945-A05E-82DC0B436CD8.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/523C89FE-9BCB-4A4F-B127-F81DC6571EB2.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/53605403-3162-1143-98B5-DB30A0CBEEA3.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/539E9B0C-7383-DB45-B5DD-C1412AFC0EBB.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/54CCBC83-B709-3F47-8698-FD138A8EAFCC.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/553664C7-B7FF-7946-AF23-A8278AC9EF87.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/558D3028-D103-F84B-829F-EAE8897A2D33.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/5610D644-50A6-3245-AD5E-DD2A58EE073B.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/57E1BB3B-E7E0-E840-9038-0E56E8B7A83B.root',
'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/57EF9727-1B58-764A-8A64-94FD70C3C100.root',
]
#options.inputFiles = ['/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/013404B1-9437-A740-9C7A-F3EBE718EFC2.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/0183EE1D-086E-B540-904D-F3B086500E5F.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/01AEC6FB-B54A-8A4B-8053-02C94620FC9D.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/06266F51-9959-AF4B-9D5E-6D82E10C36D9.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/07CB125B-A423-BD47-861E-F09A78B9EECB.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/0861F6F5-8811-3F45-A8A2-D970635E7840.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/08BA9DE0-7264-6A4B-BB7D-24ED3270E32A.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/08D68E7E-CF36-8740-A9E2-D1F413F3C218.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/098EBC6C-4F51-4844-94B3-D6B2E8BB30A1.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/0AED8F0A-EBA9-3F41-8B33-0188A63C23D2.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/0BE90291-6319-154A-BBCF-78430B539542.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/0C329E58-D296-3E4A-86B0-93D9C0E14B12.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/0EE98D40-1118-9E4D-87CB-00A295C979B9.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/115B8AE0-E7E5-B44E-BDA8-D5B0C296F44F.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/11D1ADBF-BD16-434B-81F9-583468A26771.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/122810A5-95A6-9C43-8751-0458C8F17E62.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/12CBD5D8-9900-D94C-940A-5DB8A65B5C62.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1479D0AF-2176-4B47-99F7-F1CD2561C2B6.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/151D4F1C-DF0A-8944-929F-CB21C843A148.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/155A5F55-13B9-3A4E-8C85-9C05E459291C.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/159C2FC7-7129-8947-B757-FB7D7145767A.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/15D079D7-34A9-EF40-9C88-72E73D373F76.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1623EC71-15FB-FE4E-98DC-34E83A7726C9.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/16AB417A-4C2F-134C-98F5-93AB97252C91.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/175D5D85-654B-3042-A3D8-75B8B2606382.root']
#options.inputFiles = ['/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/184836FF-0596-0744-9C8D-A1C79A1C7F29.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1B97BB81-AFC4-1941-B133-EFA7E853CAEA.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1BF8E97A-B319-A84B-A69A-505431938F4D.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1CF8FBB7-8D8A-1348-A414-1D43C2CC70EF.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1E2746AF-627A-AB40-8071-15CDBFF0CD67.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1F23C958-956E-0041-BBF9-885FB16D0934.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1F38074C-8615-244A-8F52-FE2A5C926D34.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1F38826A-78D2-1C45-B673-22FAFD622FF5.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/1FD50E1F-836D-6543-89B9-6E99C7D77045.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/201CC24C-8F0C-DB48-8E42-821B9061FAA1.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/211684D3-FFF8-C749-A39C-AF735255939B.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/21174BE8-3A55-8D44-9BC3-4D2DC9F65D03.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/230E3669-9D31-1746-8C32-2FAA9172A40D.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/2325BE01-4E4F-DD4A-8C76-288BFEBAEA39.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/23B683AD-5CB3-9A40-B48A-CE81E09FC709.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/23B76C7A-C582-EE4A-BCD9-075DA3982338.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/25B75002-195A-7E4C-8FDF-2B1773BCD23B.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/25BD8B19-007B-7C41-848A-2C34F6935C96.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/2600FA2A-2FF3-EC46-9D66-78FADEAEA8AE.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/26DF6479-BD2F-0E49-ACD5-E8ED97FA51F0.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/26E0C5C2-85CD-3447-B669-014E5B8067D6.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/274B7A46-693F-254F-A584-46D753EF745F.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/2B40FF07-C5EC-AB4A-8324-494472347258.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/2D3A445B-FEA9-5F4F-8C65-F963FA143932.root',
#'/store/data/Run2018A/ZeroBias/MINIAOD/17Sep2018-v1/00000/2EE99EAF-6D27-5E40-866F-0C28DED344D6.root']
#options.inputFiles = ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/00272A55-D6D7-014D-AE1C-BA1B4E6B3539.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/009DC2AA-D956-9C4A-9330-8CB339C5BE82.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/0110ED45-ADFF-5046-989A-6DCB04CD03F0.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/04D0F69A-0B70-7E48-828B-13260A84D403.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/05706613-0162-934D-ADB5-F8C81791BE70.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/08FB6778-BCBF-F847-AAB1-D5776044018F.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/095226B1-E4D5-4846-ABCC-266E6E98A444.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/0973B272-5119-3B4F-AAA9-3C2ED5B87A1A.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/0AEBD36D-19C1-3C4B-8EF7-9948B907CF7A.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/0F9EFCE7-65F4-E947-AA1B-F0F3F4C46F61.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/0FA3BE7B-B900-8342-8B53-FA052DDD6316.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/121F402E-C4FB-F047-9EC4-C676C23A30EE.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/127C693A-7F9B-7B4C-9D3B-1E1E9A975DF0.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/138FBDCF-72E5-9747-BC19-2BA9C55449AB.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/15207932-98B4-8C47-AF34-69FD1A3504F2.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/1597FCBC-5E20-8F4F-B183-7D0BEC0609BE.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/178A7726-6EB9-E64B-8B7B-FE6C5D58BA15.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/1843FF8E-8457-1E4F-AFB5-EF722ABEE0FD.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/18B53E0E-8B87-3346-9C24-98650D0D80B9.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/19E5D491-9A3B-AA4E-883E-FFBAD3957B3D.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/1B035962-E443-CA42-96CE-E3FE39AE65B4.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/1BB4D9A7-133E-CE4A-9874-17609817F271.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/1DD3476F-8734-2B43-B6ED-834B65602030.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/237CCDA6-FB92-D446-A6E4-1E8DBC572E20.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/25582CBF-B836-8E41-A349-EC4309B2F4F4.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_102X_upgrade2018_realistic_v15-v1/270000/26F01674-2D5C-464B-9081-259A3A0C2447.root']
#options.inputFiles = ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/00491599-68A3-BA42-838A-2650D6DCBDF9.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/0132101A-3696-6D48-A7AC-2B6B1D2E6139.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/054D0C14-97DE-654B-8FCC-66C8038A7E96.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/058F1445-BF8E-9241-BA76-6F27F6840A75.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/05EA2787-E276-1946-BCD4-24C9497F9AA9.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/060717E3-A3AD-CD41-A852-211C0685A3AE.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/08FB2BDC-8066-BB46-B6CC-0CC245064849.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/0A9A1F36-8D8D-F042-97B6-B701D8EE5A85.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/0ED2F72C-A5C0-8341-8AF9-E9197C276371.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/0FD8F6C2-7510-C94A-B9D8-063E59693AAA.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/10617138-36D0-F540-B102-2CF24A1A11CB.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/12F030DF-CD51-BF42-99EA-898003957411.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/12FBB782-318C-0B43-85E0-A7A79A9B6707.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/13021408-8DC6-EE48-9FD7-AC6E73EB073C.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/141F9009-1A1D-6F41-ACC5-F09BC2E0AC39.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/1B1E914E-2B4B-DC42-9095-ACC847B7B347.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/1CF74F47-F8DA-8D44-8B99-3A913D71E7AC.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/1E0C8FE7-E08B-F146-8DB4-D04648B37C12.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/1F2DF31F-5376-024F-9846-86ABA7C26B8F.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/1FEC3889-C277-AB43-9F97-42DB0FE0ECFF.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/200A1CDB-F11C-A344-9332-C917C95ED196.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/22191761-A204-004C-9268-2D9EC67347AD.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/23A8EE2C-D6B9-A046-AB93-240E0D360B07.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/25DF5A6C-F05F-B04C-9689-8B12FC3309EB.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/29AF66CC-3086-4342-A9E3-E7DF4534D33A.root']
#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v1/110000/02CCAE3F-A631-3D4D-B5A6-6B5B81B4E210.root'
#options.inputFiles = '/store/data/Run2018C/ParkingBPH1/MINIAOD/05May2019-v1/100000/005650EB-AA60-D34C-8100-95959C4B352D.root'
#options.inputFiles = '/store/data/Run2018D/Charmonium/MINIAOD/12Nov2019_UL2018-v1/280000/D7FD376D-30CD-AA48-8D03-E0220043BBDE.root'
#options.inputFiles = '/store/data/Run2018C/ParkingBPH1/MINIAOD/05May2019-v1/100000/005650EB-AA60-D34C-8100-95959C4B352D.root'
#options.inputFiles = '/store/data/Run2017F/Charmonium/MINIAOD/09Aug2019_UL2017-v1/20000/00BACB48-9B0F-8F48-A68B-2F08A3E9E681.root'
#options.inputFiles = '/store/data/Run2016B/Charmonium/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/240000/0333D5C7-28C0-7641-994A-ADE29A1EBAAD.root'
#options.inputFiles = '/store/data/Run2016B/Charmonium/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/240000/00251310-7FD5-BE47-B127-8CFC5B8DFE6E.root'

options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True),
                     )

process.options.numberOfThreads=cms.untracked.uint32(2)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# run, lumi, event ID
#                              eventsToProcess = cms.untracked.VEventRange('1:94:182460'),

if config["RUNONMC"]:
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(options.inputFiles),
                              # Run, lumi, Event
#                              eventsToProcess = cms.untracked.VEventRange('1:35137:14345'),
                              duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),                              
                              ) 
else:                    
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(options.inputFiles),
#                              skipEvents=cms.untracked.uint32(23000)
                              ) 

print " process source filenames %s" %(process.source) 
######## Sequence settings ##########

hltFiltersProcessName = 'RECO'
if config["RUNONMC"] or config["JSONFILE"].find('reMiniAOD') != -1:
  hltFiltersProcessName = 'PAT'

# ####### Logger ##########
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(1)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

####### Define conditions ##########
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

GT = ''


jetcorr_levels=[]
jetcorr_levels_groomed=[]
if config["RUNONMC"]:
  jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
  jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute']) # NO L1 corretion for groomed jets
else:
  jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
  jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])

   
######### read JSON file for data ##########					                                                             
if not(config["RUNONMC"]) and config["USEJSON"]:

  import FWCore.PythonUtilities.LumiList as LumiList
  import FWCore.ParameterSet.Types as CfgTypes
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  myLumis = LumiList.LumiList(filename = config["JSONFILE"]).getCMSSWString().split(',')
  process.source.lumisToProcess.extend(myLumis) 

  

####### Redo Jet clustering sequence ##########
betapar = cms.double(0.0)
fatjet_ptmin = 100.0

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.TrackJetParameters_cfi import *
from PhysicsTools.PatAlgos.tools.helpers import *
from RecoJets.JetProducers.ak4TrackJets_cfi import *
from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.JetIDProducers_cff import *
from RecoJets.Configuration.RecoTrackJets_cff import *
from RecoJets.Configuration.RecoJetAssociations_cff import *
from RecoJets.Configuration.RecoJPTJets_cff import *
#pattask = getPatAlgosToolsTask(process)
 
#process.jetGlobalReco = cms.Sequence(recoJets*recoJetIds*recoTrackJets)
#process.jetHighLevelReco = cms.Sequence(recoJetAssociations*recoPFJets*recoJPTJets)
                                                                                                         
process.chs = cms.EDFilter("CandPtrSelector",
  src = cms.InputTag('packedPFCandidates'),
  cut = cms.string('fromPV')
)

process.ak4PFJetsCHS = ak4PFJetsCHS.clone( src = 'chs' )
process.ak4PFJetsCHS.doAreaFastjet = True

#process.load("RecoJets.JetProducers.ak4TrackJets_cfi")
#process.ak4TrackJets = ak4TrackJets.clone()

pattask = getPatAlgosToolsTask(process)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#dataFormat=DataFormat.MiniAOD
#switchOnVIDElectronIdProducer(process,dataFormat,task=pattask)

#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
     
#my_id_modules = [
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
#                 ]
           
#add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection,task=pattask)

####### Ntuplizer initialization ##########
jetsAK4 = "slimmedJets"


METS = "slimmedMETs"
METS_EGclean = "slimmedMETsEGClean"
METS_MEGclean = "slimmedMETsMuEGClean"
METS_uncorr = "slimmedMETsUncorrected"

if config["USENOHF"]: METS = "slimmedMETsNoHF"  

##___________________ MET significance and covariance matrix ______________________##


##___________________ Jets ______________________##

  

######## JEC ########
jecLevelsAK8chs = []
jecLevelsAK8Groomedchs = []
jecLevelsAK4chs = []
jecLevelsAK4 = []
jecLevelsAK8Puppi = []
jecLevelsForMET = []

print "1. options.RunPeriod ", options.RunPeriod
if options.RunPeriod=="" : options.RunPeriod=options.inputFiles[0]

if  config["RUNONMC"] :
  JECprefix = ""
  if ("Fall17" in options.RunPeriod):
    JECprefix = "Fall17_17Nov2017_V32"
    GT ='94X_mcRun2_asymptotic_v3'
  elif ("Summer16" in options.RunPeriod):
    JECprefix = "Summer16_07Aug2017_V11"
    GT = '94X_mc2017_realistic_v17'
  elif (("Autumn18" in options.RunPeriod) or ("Fall18" in options.RunPeriod)):
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'
  else:
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'

  #jecAK8chsUncFile = "JEC/%s_MC_Uncertainty_AK8PFchs.txt"%(JECprefix)
  jecAK4chsUncFile = "JEC/%s_MC_Uncertainty_AK4PFchs.txt"%(JECprefix)
 



else : #Data
   JECprefix = ""
   JEC_runDependent_suffix= ""

   if ("2017" in options.RunPeriod):
     if ("Run2017B" in  options.RunPeriod): JEC_runDependent_suffix= "B"
     elif ("Run2017C" in  options.RunPeriod): JEC_runDependent_suffix= "C"
     elif ("Run2017D" in  options.RunPeriod): JEC_runDependent_suffix= "D"
     elif ("Run2017E" in  options.RunPeriod): JEC_runDependent_suffix= "E"
     elif ("Run2017F" in  options.RunPeriod): JEC_runDependent_suffix= "F"
     
     JECprefix = "Fall17_17Nov2017"+JEC_runDependent_suffix+"_V32"
     GT = '94X_dataRun2_v11'
     
     
   elif ("2016" in options.RunPeriod):
     if ("Run2016D" in  options.RunPeriod or "Run2016B" in  options.RunPeriod  or "Run2016C" in  options.RunPeriod  ): JEC_runDependent_suffix= "ABC"
     elif ("Run2016E" in  options.RunPeriod): JEC_runDependent_suffix= "EF"
     elif ("Run2016G" in  options.RunPeriod): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and  not options.runUpToEarlyF): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and   options.runUpToEarlyF): JEC_runDependent_suffix= "EF"


     JECprefix = "Summer16_07Aug2017"+JEC_runDependent_suffix+"_V11"
     GT ='94X_dataRun2_v10'

   elif ("2018" in options.RunPeriod):
     if ("Run2018A" in  options.RunPeriod ): 
       JEC_runDependent_suffix= "A"
       GT="102X_dataRun2_Sep2018ABC_v2" 
     elif ("Run2018B" in  options.RunPeriod): 
       JEC_runDependent_suffix= "B"
       GT="102X_dataRun2_Sep2018ABC_v2"
#       GT="106X_dataRun2_v28"
     elif ("Run2018C" in  options.RunPeriod): 
       JEC_runDependent_suffix= "C"
       GT="102X_dataRun2_Sep2018ABC_v2"
     elif ("Run2018D" in  options.RunPeriod): 
       JEC_runDependent_suffix= "D"
       GT = '102X_dataRun2_Prompt_v16' 

     JECprefix = "Autumn18_Run"+JEC_runDependent_suffix+"_V8"
    
   #jecAK8chsUncFile = "JEC/%s_DATA_Uncertainty_AK8PFchs.txt"%(JECprefix)
   jecAK4chsUncFile = "JEC/%s_DATA_Uncertainty_AK4PFchs.txt"%(JECprefix)
 
#   GT = '106X_dataRun2_v27' 
   print "jec JEC_runDependent_suffix %s ,  prefix %s " %(JEC_runDependent_suffix,JECprefix)

print "jec prefix ", JECprefix

print "doing corrections  to met on the fly %s" ,config["CORRMETONTHEFLY"]

print "*************************************** GLOBAL TAG *************************************************" 
print GT
print "****************************************************************************************************" 
process.GlobalTag = GlobalTag(process.GlobalTag, GT)




if config["CORRMETONTHEFLY"]:  
   if config["RUNONMC"]:
     jecLevelsForMET = [				       
     	 'JEC/%s_MC_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_MC_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_MC_L3Absolute_AK4PFchs.txt'%(JECprefix)
       ]
   else:       					       
     jecLevelsForMET = [
     	 'JEC/%s_DATA_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_DATA_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_DATA_L3Absolute_AK4PFchs.txt'%(JECprefix),
         'JEC/%s_DATA_L2L3Residual_AK4PFchs.txt'%(JECprefix)
       ]	
      			    

                                                                       
################## Ntuplizer ###################

process.ntuplizer = cms.EDAnalyzer("Ntuplizer",
    runOnMC	      = cms.bool(config["RUNONMC"]),
    useDNN	      = cms.bool(config["USEDNN"]),
    useHammer	      = cms.bool(config["USEHAMMER"]),
    doGenParticles    = cms.bool(config["DOGENPARTICLES"]),
    doGenEvent	      = cms.bool(config["DOGENEVENT"]),
    doPileUp	      = cms.bool(config["DOPILEUP"]),
    doJpsiMu	      = cms.bool(config["DOJPSIMU"]),
    doJpsiTau	      = cms.bool(config["DOJPSITAU"]),
    doBsTauTau	      = cms.bool(config["DOBSTAUTAU"]),
    doBsTauTauFH      = cms.bool(config["DOBSTAUTAUFH"]),
    doBsTauTauFH_mr   = cms.bool(config["DOBSTAUTAUFH_mr"]),
    doBsDstarTauNu    = cms.bool(config["DOBSDSTARTAUNU"]),
    doVertices	      = cms.bool(config["DOVERTICES"]),
    doMissingEt       = cms.bool(config["DOMISSINGET"]),
    doGenHist         = cms.bool(config["DOGENHIST"]),
    doInstanton       = cms.bool(config["DOINSTANTON"]),
    verbose           = cms.bool(config["VERBOSE"]),
    dzcut             = cms.double(config['DZCUT']),
    fsigcut           = cms.double(config['FSIGCUT']),
    vprobcut          = cms.double(config['VPROBCUT']),
    dnncut            = cms.double(config['DNNCUT']),
    tau_charge        = cms.uint32(config['TAU_CHARGE']),
    dnnfile_old       = cms.string(config['DNNFILE_OLD']),                        
    dnnfile_perPF     = cms.string(config['DNNFILE_PERPF']),                        
    dnnfile_perEVT    = cms.string(config['DNNFILE_PEREVT']),                        
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    taus = cms.InputTag("slimmedTaus"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),

#    eleHEEPId51Map = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
#    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
#    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
#    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
#    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
#    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),

#    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"),
#    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
#    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
#    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),

#    eleHLTIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"), 
#    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
                                   
#    eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
#    eleMVATightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
#    mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
#    mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
    dupCluster          = cms.InputTag("particleFlowEGammaGSFixed:dupECALClusters"),
    hitsNotReplaced     = cms.InputTag("ecalMultiAndGSGlobalRecHitEB:hitsNotReplaced"),
    mets = cms.InputTag(METS),
    mets_EGclean = cms.InputTag(METS_EGclean),
    mets_MEGclean = cms.InputTag(METS_MEGclean),
    mets_uncorr = cms.InputTag(METS_uncorr),
    mets_puppi = cms.InputTag("slimmedMETsPuppi"),
    mets_mva = cms.InputTag("MVAMET","MVAMET"),
    corrMetPx = cms.string("+0.1166 + 0.0200*Nvtx"),
    corrMetPy = cms.string("+0.2764 - 0.1280*Nvtx"),
    jecAK4forMetCorr = cms.vstring( jecLevelsForMET ),
    jetsForMetCorr = cms.InputTag(jetsAK4),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    genparticles = cms.InputTag("prunedGenParticles"),
    packedgenparticles = cms.InputTag("packedGenParticles"),
#    gentaus = cms.InputTag("tauGenJets"),
    PUInfo = cms.InputTag("slimmedAddPileupInfo"),
    genEventInfo = cms.InputTag("generator"),
    externallheProducer = cms.InputTag("externalLHEProducer"),
    HLT = cms.InputTag("TriggerResults","","HLT"),
    triggerobjects = cms.InputTag("slimmedPatTrigger"),
    triggerprescales = cms.InputTag("patTrigger"),
    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
    jecpath = cms.string(''),
   
    
    ## Noise Filters ###################################
    # defined here: https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),   # both data and MC for 2018,
    noiseFilterSelection_HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
    noiseFilterSelection_HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
    noiseFilterSelection_HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),    # both data and MC for 2018,  
    noiseFilterSelection_ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),  # both data and MC for 2018,
    noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHaloFilter'),
    noiseFilterSelection_CSCTightHalo2015Filter = cms.string('Flag_CSCTightHalo2015Filter'),
    noiseFilterSelection_hcalLaserEventFilter = cms.string('Flag_hcalLaserEventFilter'),
    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),  # both data and MC for 2018,
    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),  # both data and MC for 2018,
    noiseFilterSelection_trackingFailureFilter = cms.string('Flag_trackingFailureFilter'),
    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
    noiseFilterSelection_ecalLaserCorrFilter = cms.string('Flag_ecalLaserCorrFilter'),
    noiseFilterSelection_trkPOGFilters = cms.string('Flag_trkPOGFilters'),
    
    #New for ICHEP 2016
    noiseFilterSelection_CSCTightHaloTrkMuUnvetoFilter = cms.string('Flag_CSCTightHaloTrkMuUnvetoFilter'),
    noiseFilterSelection_globalTightHalo2016Filter = cms.string('Flag_globalTightHalo2016Filter'),
    noiseFilterSelection_globalSuperTightHalo2016Filter = cms.string('Flag_globalSuperTightHalo2016Filter'), # both data and MC for 2018,  
    noiseFilterSelection_HcalStripHaloFilter = cms.string('Flag_HcalStripHaloFilter'),
    noiseFilterSelection_chargedHadronTrackResolutionFilter = cms.string('Flag_chargedHadronTrackResolutionFilter'),
    noiseFilterSelection_muonBadTrackFilter = cms.string('Flag_muonBadTrackFilter'),
    
    #New for Moriond
    noiseFilterSelection_badMuonsFilter = cms.string('Flag_BadPFMuonFilter'),    #('Flag_badMuons'),  # both data and MC for 2018, 
    noiseFilterSelection_duplicateMuonsFilter = cms.string('Flag_duplicateMuons'),
    noiseFilterSelection_nobadMuonsFilter = cms.string('Flag_nobadMuons'),

    # and the sub-filters
    noiseFilterSelection_trkPOG_manystripclus53X = cms.string('Flag_trkPOG_manystripclus53X'),
    noiseFilterSelection_trkPOG_toomanystripclus53X = cms.string('Flag_trkPOG_toomanystripclus53X'),
    noiseFilterSelection_trkPOG_logErrorTooManyClusters = cms.string('Flag_trkPOG_logErrorTooManyClusters'),
    # summary
    noiseFilterSelection_metFilters = cms.string('Flag_METFilters'),

    packedpfcandidates = cms.InputTag('packedPFCandidates'),
    SecondaryVertices = cms.InputTag('slimmedSecondaryVertices'),
#    losttrack = cms.InputTag('lostTracks')
)

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(True)#False
    )




process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


####### Final path ##########
process.p = cms.Path()

process.p += process.ecalBadCalibReducedMINIAODFilter

#if config["RUNONMC"]:
#  process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
#  process.tauGenJets.GenParticles = cms.InputTag('prunedGenParticles')
#  process.p += process.tauGenJets


process.p += process.ntuplizer
#process.p += process.ak4TrackJets
process.p.associate(pattask)

print pattask

#  LocalWords:  tauIdMVAIsoDBoldDMwLT
