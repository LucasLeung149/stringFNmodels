(* ::Package:: *)

(* ::Title:: *)
(*stringFNmodels: GA Environment for Froggatt-Nielsen-type string model building*)


(* ::Section:: *)
(*Start package*)


BeginPackage["stringFNmodels`"];


(* ::Section:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*General help*)


(* ::Subsubsection::Closed:: *)
(*Basic*)


stringFNmodels::usage="\[FilledSmallSquare]  This package realises a GA and RL environment for stringy Froggatt-Nielsen (heterotic line bundle) models.\n"<>
"\[FilledSmallSquare]  Basic parameters of the environment are defined in Options[stringFNmodels], including the number of singlets, the number of U(1)s, etc.\n"<>
"\[FilledSmallSquare]  For more information about the options type ?stringFNmodelsOpts.\n"<>
"\[FilledSmallSquare]  The current version of the code allows linkage to the GA code, `Genetic. \n"<>
"\[FilledSmallSquare]  The O1 coefficients are fixed at the start .\n"<>
"\[FilledSmallSquare]  Type ?stringFNmodelsAuxModules for more information on the auxiliary modules.\n"<>
"\[FilledSmallSquare]  CompleteStateQsector[state] is the key function to compute the full state assocation given a state in its state-vector, bitlist or state-association form.\n"<>
"\[FilledSmallSquare]  The basic guide of how to use the enviroment can be found by typing ?stringFNmodelsHowtoUse.\n"<>
"\[FilledSmallSquare]  Details of the main modules can be found by typing ?stringFNmodelsModules.";
(*"\[FilledSmallSquare]  For the purpose of linking to an RL system a state is given in association form <|\"Input\"->statevector,\"Terminal\"->True/False,\"Value\"->value|>.\n"<>
"\[FilledSmallSquare]  For the purpose of linking to an GA system, the state is given in association form <|\"Bits\"->statebitlist,\"Terminal\"->True/False,\"Fitness\"->value|>.\n"<>*)


(* ::Subsubsection::Closed:: *)
(*Options*)


(* Documentation for opts *)
stringFNmodelsOpts::usage="\[FilledSmallSquare]  Options[stringFNmodels] gives the global environment options for the package.\n"<>
"\[FilledSmallSquare]  For global environment settings type ?stringFNmodelsOptsGlobal for details.\n"<>
"\[FilledSmallSquare]  For fitness factors settings type ?stringFNmodelsOptsFitness for details.\n"<>
"\[FilledSmallSquare]  For SM Quantities type ?stringFNmodelsOptsSM for details.";

stringFNmodelsOptsGlobal::usage=" Global environment options:\n"<>
"\[FilledSmallSquare]  Global/Environment Settings: \n"<>
"  \[FilledVerySmallSquare]  \"nU1s\" - Number of anomalous GS-type U(n)s in the model.\n"<>
"  \[FilledVerySmallSquare]  \"nChargeVec\" - The split vector specifying the \!\(\*SubscriptBox[\(n\), \(a\)]\) of the split structure group.\n"<>
"  \[FilledVerySmallSquare]  \"nSinglet\" - Number of singlets/RH neutrinos in the model.\n"<>
"  \[FilledVerySmallSquare]  \"nNonperts\" - Number of non-perturbative fields \[CapitalPhi] in the model.\n"<>
"  \[FilledVerySmallSquare]  \"\[CapitalPhi]ChargeRange\" - Charge range of non-perturbative fields in the form {a,b}. The charge range takes integral values and must have a range that is a power of 2, i.e. b-a = \!\(\*SuperscriptBox[\(2\), \(n\)]\) where n is a number.\n"<>
"  \[FilledVerySmallSquare]  \"DimOp\" - Number of maximum singlet field insertions in each operator.\n"<>
"  \[FilledVerySmallSquare]  \"DimOpof\[CapitalPhi]\" - Number of maximum non-perturbative fields in the operator expansion.\n"<>
"  \[FilledVerySmallSquare]  \"bitsfor\[Phi]\" - Number of bits specified for each of the singlet moduli VEVs. \n"<>
"  \[FilledVerySmallSquare]  \"bitsfor\[CapitalPhi]\" - Number of bits specified for each of the Kahler moduli VEVs. \n"<>
"\[FilledSmallSquare]  Symbol Settings: \n"<>
"  \[FilledVerySmallSquare]  \"BundleModuliFields\" - Symbol used to denote bundle moduli/ singlet fields in the code. Currently set to Global`\[Phi].\n"<>
"  \[FilledVerySmallSquare]  \"KahlerEffFields\" - Symbol used to denote effective Kahler moduli/ non-perturbative fields in the code. Currently set to Global`\[CapitalPhi].\n"<>
"  \[FilledVerySmallSquare]  \"VEVsScale\" - Symbol used to denote the scale of the VEVs. Currently set to Global`\[Epsilon].\n"<>
"\[FilledSmallSquare]  VEV Settings: \n"<>
"  \[FilledVerySmallSquare]  \"\[Phi]VEVmin\" - The minimum value of the singlet VEVs.\n"<>
"  \[FilledVerySmallSquare]  \"\[Phi]VEVlogGen\" - The logarithmic factor to generate singlet VEVs.\n"<>
"  \[FilledVerySmallSquare]  \"\[CapitalPhi]VEVmin\" - The minimum value of the non-perturbative VEVs.\n"<>
"  \[FilledVerySmallSquare]  \"\[CapitalPhi]VEVlogGen\" - The logarithmic factor to generate non-perturbative VEVs.\n"<>
"\[FilledSmallSquare]  Fixed Charge Settings: \n"<>
"  \[FilledVerySmallSquare]  \"fixedCharges\" - A boolean determining whether the 10 and 5bH reps are fixed.\n"<>
"  \[FilledVerySmallSquare]  \"tenQ\" - The fixed charges for the 10 reps.\n"<>
"  \[FilledVerySmallSquare]  \"fiveBarHiggsQ\" - The fixed charges for the 5bH rep.\n"<>
"\[FilledSmallSquare]  O1 Settings: \n"<>
"  \[FilledVerySmallSquare]  \"O1CoeffRange\" - The range of the O1 coefficients, preset to {0.5,3}.\n"<>
"  \[FilledVerySmallSquare]  \"numOfO1Coeff\" - The number of O1 coefficients generated in the Global variable stringFNmodels`O1CLst at the start of the code.\n"<>
"  \[FilledVerySmallSquare]  \"O1CoeffVarDiff\" - The O(1) variation implemented in the O(1) fitness calculation. This is the factor in the denominator of the numerical differentiation.\n"<>
"\[FilledSmallSquare]  Method of Generating Operators: \n"<>
"  \[FilledVerySmallSquare]  \"MethodGenOps\" - This is the method of generating operating operators. The current available options are \"List\" and \"Graph\".";

stringFNmodelsOptsFitness::usage=" Fitness options:\n"<>
"\[FilledSmallSquare]  Fitness Sector Options: \n"<>
"  \[FilledVerySmallSquare]  \"FitRange\" - This is the minimum of the fitness value for the state to be terminal.\n"<>
"  \[FilledVerySmallSquare]  \"factorfitnessH\" - The factor multiplying the fitness contribution of the Higgs sector.\n"<>
"  \[FilledVerySmallSquare]  \"factorFitnessRG\" - The factor multiplying the fitness contribution of the RG sector.\n"<>
"  \[FilledVerySmallSquare]  \"factorFitnessQ\" - The factor multiplying the fitness contribution of the Quark Mass sector.\n"<>
"  \[FilledVerySmallSquare]  \"factorFitnessCKM\" - The factor multiplying the fitness contribution of the Quark Mixing sector.\n"<>
"  \[FilledVerySmallSquare]  \"factorFitnessGS\" - The factor multiplying the fitness contribution of the Green-Schwarz sector.\n"<>
"  \[FilledVerySmallSquare]  \"factorFitnessNAF\" - The factor multiplying the fitness contribution of the non-Abelian sector.\n"<>
"  \[FilledVerySmallSquare]  \"factorFitnessTEX\" - The factor multiplying the fitness contribution of the texture sector.\n"<>
"  \[FilledVerySmallSquare]  \"factorFitnessO1\" - The factor multiplying the fitness contribution of the O(1) sector.\n"<>
"  \[FilledVerySmallSquare]  \"MaxLogPower\" - This is the maximum contribution in logarithmic contributions in the fitness functions.\n"<>
"\[FilledSmallSquare]  Extra Fitness Penalties Options: \n"<>
"  \[FilledVerySmallSquare]  \"indPenal\" - A boolean determining whether individual penalties are turned on.\n"<>
"  \[FilledVerySmallSquare]  \"diffPenal\" - A boolean determining whether the extra fitness penalties are added due to different between sectors.\n"<>
"  \[FilledVerySmallSquare]  \"indFitMax\" - The maximum fitness contribution for individual penalties in the fitness.\n"<>
"\[FilledSmallSquare]  Individual Fitness Factors Options: \n"<>
"  \[FilledVerySmallSquare]  \"fitDownFac\" - A dim-3 factor to be multiplied to the down-quark mass contributions in the fitness.\n"<>
"  \[FilledVerySmallSquare]  \"fitUpFac\" - A dim-3 factor to be multiplied to the up-quark mass contributions in the fitness.\n"<>
"  \[FilledVerySmallSquare]  \"fitCKMFac\" - A dim-9 factor to be multiplied to the CKM matrix entry contributions in the fitness.\n"<>
"  \[FilledVerySmallSquare]  \"fitTextFac\" - A dim-4 factor to be multiplied to texture contributions in the fitness. This is in order: (powers,variance_of_scale,scale_range,o1_top_range).";

stringFNmodelsOptsSM::usage=" SM Quantities (Data):\n"<>
"\[FilledSmallSquare]  \"v\" - Higgs VEV.\n"<>
"\[FilledSmallSquare]  \"vErr\" - Higgs VEV error.\n"<>
"\[FilledSmallSquare]  \"mu\" - up, charm and top quark masses.\n"<>
"\[FilledSmallSquare]  \"md\" - down, strange and bottom quark masses.\n"<>
"\[FilledSmallSquare]  \"CKMmat\" - CKM matrix entries (magnitude).";


(* ::Subsubsection::Closed:: *)
(*Modules*)


(* Module descriptions and list of functions. *)
stringFNmodelsAuxModules::usage="List of Modules (type ?$$Modules for details):\n"<>
"\[FilledSmallSquare]  BitGen - auxiliary modules for bitlist generations and coupling to GA.\n"<>
"\[FilledSmallSquare]  DataAnal - auxiliary modules for data analysis.\n"<>
"\[FilledSmallSquare]  FindOps - auxiliary modules for finding operators.\n"<>
"\[FilledSmallSquare]  Fit - auxiliary modules for fitness computations.\n"<>
"\[FilledSmallSquare]  O1 - auxiliary modules for O(1) coefficients.\n"<>
"\[FilledSmallSquare]  OpsToMat - auxiliary modules for O(1) coefficients.\n"<>
"\[FilledSmallSquare]  RG - auxiliary modules for RG computations.\n"<>
"\[FilledSmallSquare]  SM - auxiliary modules for SM quantity computations.\n"<>
"\[FilledSmallSquare]  StateGen - auxiliary modules for state generation.\n"<>
"\[FilledSmallSquare]  Text - auxiliary modules for texture computations.";

StateGenModules::usage="StateGenModules contains auxiliary modules for state generation.\n"<>
"List of StateGen modules:\n"<>
"\[FilledSmallSquare]  RandomStatewithVEVs[] generates a random state in state-vector form.\n"<>
"\[FilledSmallSquare]  ConvertStatewithVEVs[stateVec] converts a state from state-vector to association form.\n"<>
"\[FilledSmallSquare]  StateVectoMatwithVEVs[stateVec] converts a state from state-vector form to a charge matrix where the rows are the species and the columns are the U(1) factors.\n"<>
"\[FilledSmallSquare]  StateVectoVEVs[stateVec] converts a state from state-vector form to VEV substitution form.\n"<>
"\[FilledSmallSquare]  AssocFormToStateVec[assoc] converts the state in association form into a state-vector form.\n"<>
"\[FilledSmallSquare]  AssocFormToVEVs[assoc] converts the inputted state in association form into VEV substitution rules.\n"<>
"\[FilledSmallSquare]  ConvertStateVecToVEVPowers[statevec] converts the inputted state-vector into an association of powers of the VEVs.";

BitGenModules::usage="BitGenModules contains auxiliary modules for bitlist generation.\n"<>
"This is primarily used with couplings to the GA module.\n"<>
"List of BitGen modules:\n"<>
"\[FilledSmallSquare]  GADetermDimStateVecwithVEVs[] determines the dimension of the binary vector to be generated in GA.\n"<>
"\[FilledSmallSquare]  GARandomStatewithVEVs[] generates a random bitstring for GA.\n"<>
"\[FilledSmallSquare]  GAConvertBitlstStateVecwithVEVs[bitlst] converts a bitlist into state-vector form.\n"<>
"\[FilledSmallSquare]  The following modules are used when the 10 reps and H rep charges are fixed.\n"<>
"\[FilledSmallSquare]  GADetermDimStateVecFix[] determines the dimension of the binary vector to be generated in GA.\n"<>
"\[FilledSmallSquare]  GARandomStateFix[] generates a random bitstring for GA.\n"<>
"\[FilledSmallSquare]  GAConvertBitlstStateVecFix[bitlst] converts a bitlist into state-vector form.\n"<>
"\[FilledSmallSquare]  GAConvertStateVectoBitlst[stateVec] converts a state-vector into bitlist form.";

FindOpsModules::usage="FindOpsModules contains auxiliary modules for computing operators.\n"<>
"List of FindOps modules:\n"<>
"\[FilledSmallSquare]  QMatRepChange[qvec,nvec] converts the upstairs chargevec rep to the downstairs one.\n"<>
"\[FilledSmallSquare]  FromPowertoOps[powervec] converts a power vector to the approriate operator form.\n"<>
"\[FilledSmallSquare]  FindOpsforWIJawithList\[Phi]VEVsQsector[state,O1vec,O1count] find all the possible operators combinations in the Standard Model part of the superpotential from a given state, using the List method.\n"<>
"\[FilledSmallSquare]  UpdateVectorWithMinimum[vec,posList] updates the inputted vec according to the posList which gives the duplicate entries.\n"<>
"\[FilledSmallSquare]  MinReplaceListRule[vec,posList] gives a binary vector on whether there is a repeated Q entry being replaced in UpdateVectorWithMinimum.\n"<>
"\[FilledSmallSquare]  FindOpsFromGraphMethod[\[Phi]Qs,\[Phi]VEVs,QToSolve] aims to find the leading order operator insertion given the charges and VEVs of the singlet moduli fields to the inputted matched charge.\n"<>
"\[FilledSmallSquare]  GenerateOrderedTuples[sum, length] generates combinations of integer-tuples up to the sum of length.\n"<>
"\[FilledSmallSquare]  SumChargesByPowers[kPows,kQs] sum the charges according to the power.\n"<>
"\[FilledSmallSquare]  ProcessModPowers[vec1,vec2] joins the two vectors element-wise and delete the Null ones.\n"<>
"\[FilledSmallSquare]  FindOpsFromGraphAdd\[CapitalPhi][\[Phi]Qs,\[Phi]VEVs,QToSolve,\[CapitalPhi]Qs] computes the operators using the graph method with the non-perturbative insertions.\n"<>
"\[FilledSmallSquare]  FindOpsforWIJawithGraph\[Phi]VEVsQsector[state,O1vec,O1count] find all the possible operators combinations in the Standard Model part of the superpotential from a given state, using the Graph method.\n";

TextModules::usage="TextModules contains auxiliary modules for computing texture and fine-tuning constraints.\n"<>
"List of Text modules:\n"<>
"\[FilledSmallSquare]  FromWOpsGetLeadingYukawaMatrix[WOps,VEVsPowerListAssoc] computes the leading order powers of the Yukawa matrix entries.\n"<>
"\[FilledSmallSquare]  LowestOrderOfMonomialListOneVar[monomialList,var] finds the lowest order monomial from an inputted list for one variable.\n"<>
"\[FilledSmallSquare]  RemoveCoeffFromMonomial[monomial] removes coefficient from the monomial.\n"<>
"\[FilledSmallSquare]  MassHierarchyFromLeadingOrders[YukLeadPowAssoc] computes the mass hierarchy from the leading order powers of the two Yukawa matrices.\n"<>
"\[FilledSmallSquare]  LocateExtractO1Coeffs[WOpsOut,VEVPowersAssoc,MHAssoc,O1coeffList] locates the leading order positons of the O1 coefficients and the simultaneously the positions of the O(1) factors for the top and bottom couplings in O1coeffList.\n"<>
"\[FilledSmallSquare]  CalcScaleAndHiggsO1Coeff[massHierarchyAssoc] computes the scale and the O(1) coefficient from the leading order data.\n"<>
"\[FilledSmallSquare]  GeometricStandardDeviation[list] computes the geometric SD of the inputted list.\n"<>
"\[FilledSmallSquare]  FitnessTexture[massHierarchyAssoc,VEVPowers] computes the fitness part from the Texture.";

OpsToMatModules::usage="OpsToMatModules contains auxiliary modules for converting the Yukawa matrices from operator form into matrices.\n"<>
"List of OpsToMat modules:\n"<>
"\[FilledSmallSquare]  AdjoinO1coeffWOpsTextQSec[WOps,MHAssoc,VEVsPowers,O1Assoc,O1Vec] adjoins the O(1) coefficients to the operators to obtain the numerical Yukawa matrices.";

RGModules::usage="RGModules contains auxiliary modules for computing the RG evolution.\n"<>
"List of RG modules:\n"<>
"\[FilledSmallSquare]  RGGivenYukMZFindYukMGUT[ytMZ,ybMZ] computes the RG evolution from IR to UV.\n"<>
"\[FilledSmallSquare]  RGGivenYukMGUTFindYukMZ[ytMGUT,ybMGUT] computes the RG evolution from UV to IR.\n"<>
"\[FilledSmallSquare]  RGGivenYukMZFindYukRatiosAll[ytMZ,ybMZ] computes the RG evolution from IR to UV.\n"<>
"\[FilledSmallSquare]  RGGivenYukMZFindYukRatiosAllFixedPlane[ytMZ,ybMZ] computes the RG evolution from IR to UV.\n"<>
"\[FilledSmallSquare]  RGGivenYukGUTFindYukRatiosAll[ytMGUT,ybMGUT] computes the RG evolution from UV to IR.\n"<>
"\[FilledSmallSquare]  ConstructRGFactors[y33Factor,y3iFactor,yijFactor,leadPos] computes the RG factor form (to be multiplied to the GUT Yukawa matrix.\n"<>
"\[FilledSmallSquare]  GivenYukMatsFindRGMatricesPlane[yukMZAssoc] computes the RG matrix to be multiplied to the GUT Yukawa matrices.\n"<>
"\[FilledSmallSquare]  O1LogDevPenal[num] test whether the number is O(1) which range is specified in the Options.\n"<>
"\[FilledSmallSquare]  FitnessRG[YukO1Assoc] computes the relevant fitness in the RG sector.";

SMModules::usage="SMModules contains auxiliary modules for computing SM quantities.\n"<>
"List of SM Computation modules:\n"<>
"\[FilledSmallSquare]  ComputeSMQuantities[Yd,Yu] computes the SM quantities for the quark sector given the down- and up-Yukawa matrices.";

FitModules::usage="FitModules contains auxiliary modules for computing fitness functions.\n"<>
"List of Fitness modules:\n"<>
"\[FilledSmallSquare]  NormalLog[num_,denom_] computes the absolute value of log10 of the faction (num/denom).\n"<>
"\[FilledSmallSquare]  FitnessHiggs[massAssoc] computes the fitness in the Higgs sector.\n"<>
"\[FilledSmallSquare]  FitnessQuarksMassMix[massAssoc] computes the fitness in the Quark Mass and Mixing sector.\n"<>
"\[FilledSmallSquare]  FitnessGSAnomaly[Qmat] computes the fitness in the Green-Schwarz sector.\n"<>
"\[FilledSmallSquare]  CheckCharge[vector] takes in a 2 component vector and checks whether the entries are equal.\n"<>
"\[FilledSmallSquare]  FitnessNonAbelianFactors[stateVec] computes the fitness in the non-abelian sector.\n"<>
"\[FilledSmallSquare]  FitnessCalcQsectorOnly[massAssoc,stateVec,YukMZValsO1] computes all of the fitness functions relevant to the quark sector.";

O1Modules::usage="O1Modules contains auxiliary modules for generating and evaluating O(1) coefficients.\n"<>
"List of O1 Modules:\n"<>
"\[FilledSmallSquare]  SetO1Coeff[O1Vec] sets the O1 coefficient to be the inputted vector O1Vec.\n"<>
"\[FilledSmallSquare]  GetO1Coeff[] obtains the current O(1) coefficients in the settings.\n"<>
"\[FilledSmallSquare]  NewO1Coeff[] generates a vector of O(1) coefficients for the operators.\n"<>
"\[FilledSmallSquare]  ChangeO1VecByAmount[O1Vec,n] changes the O1Vec by a fixed amount set by the Options at position n.\n"<>
"\[FilledSmallSquare]  FitnessO1Coeff[varMassAssoc,massAssoc,O1val] computes the fitness contributions of O1 coefficients.\n";

DataAnalModules::usage="DataAnalModules contain auxiliary modules that are used in analying the data produced by the C-code GA runs.\n"<>
"List of Data Analysis Modules:\n"<>
"\[FilledSmallSquare]  LowestOrderTerms[expression,vars] splits out the lowest order terms given variables and expressions.\n"<>
"\[FilledSmallSquare]  VecDuplSymm[vec] determines the duplicates and their positions in a given vector.\n"<>
"\[FilledSmallSquare]  RedunQSymm[vec] determines the redundant symmetry in the state in statevector form given the settings in Options.\n"<>
"\[FilledSmallSquare]  CalcQScoreFix[fiveBarQ,oneQ,kahlQ] determines whether we need to exchange the charges by canonically defining a direction in the redundant \!\(\*SubscriptBox[\(S\), \(2\)]\) symmetry, if it exists.\n"<>
"\[FilledSmallSquare]  StateNormalForm[stateVec] presents the normal form for the inputted statevector.\n"<>
"\[FilledSmallSquare]  EvalStateVecSoloList[inputList] evaluates the inputted list of bitlists and adjoins the \"StateVec\" association to the states.\n"<>
"\[FilledSmallSquare]  AppendSetNoToStates[inputList] tags the list of input states with the appropriate set number.\n"<>
"\[FilledSmallSquare]  WriteUniqueStates[filename] reads the file and converts the content into a list of states that has been processed.\n"<>
"\[FilledSmallSquare]  SelectVecsSameCharge[listOfStates] delete the duplicates of the list of states after being processed by the StateNormalForm function.\n"<>
"\[FilledSmallSquare]  SelectVecsSameModel[listOfStates] delete the duplicates of the list of states after being processed by the StateNormalForm function.\n"<>
"\[FilledSmallSquare]  ProduceReducedList[evalStates] runs through the inputted this and produces the list of states flattened and reduced using StateNormalForm.\n"<>
"\[FilledSmallSquare]  ProduceReducedModelList[evalStates] runs through the inputted this and produces the list of states flattened and reduced using StateNormalForm.\n"<>
"\[FilledSmallSquare]  ComputeNumModuli[state] computes the number of moduli used in the leading order Yukawa textures.\n"<>
"\[FilledSmallSquare]  ProduceFullRedList[inputFile,outputFile] takes in the input file of a GA output and prints out the list of reduced states in the outputFile.\n"<>
"\[FilledSmallSquare]  ProduceFullRedModelList[inputFile,outputFile] takes in the input file of a GA output and prints out the list of reduced states in the outputFile.\n"<>
"\[FilledSmallSquare]  PlotAccumGraph[accumVec,title] produces the GA accumulation file given the accumulation vector and the title.\n"<>
"\[FilledSmallSquare]  FromRedListGetPlot[filePath,num] produces the GA accumulation file given file path for the states.\n"<>
"\[FilledSmallSquare]  GivenModNumFromRedListGetPlot[filePath,num,numSing,numKahl] produces the GA accumulation file given file path for the states specified to the number of moduli specified.\n"<>
"\[FilledSmallSquare]  ChargeStateNormalForm[QstateVec] presents the normal form for the inputted charge part of statevector (without the VEVs).\n"<>
"\[FilledSmallSquare]  FromStateVecListGiveNum[filePath] reduces the list of states and gives the number of reduced states as output.";

stringFNmodelsModules::usage="Below is a list of the main modules.\n"<>
"\[FilledSmallSquare]  CompleteStateGivenWOps[stateVec,WOpsOut,massHierAssoc,VEVsPowers,O1Vec] computes the state given the WOps.\n"<>
"\[FilledSmallSquare]  CompleteStateO1Coeffs[stateVec,O1Vec] computes the full state given the O1 vector O1Vec.\n"<>
"\[FilledSmallSquare]  CompleteStateQsector[state] computes the full state.\n"<>
"\[FilledSmallSquare]  GAFitnessFuncQsector[state] computes the full state for GA.\n"<>
"\[FilledSmallSquare]  GAInitialiseStateQsector[] randomises a state and computes its full state.";


(* ::Subsubsection::Closed:: *)
(*Guides*)


stringFNmodelsHowtoUse::usage=" This HowtoUse Guide outlines the basic elements of the system. \n"<>
"\[FilledSmallSquare]  This package realises the environment for heterotic line bundle models. \n"<>
"\[FilledSmallSquare]  The basic state of the system is encoded in state-vector form or bitlist form.\n"<>
"\[FilledSmallSquare]  The bitlist form is used to couple to the GA system. The map from bitlist to state-vector is surjective but not necessarily injective.\n"<>
"\[FilledSmallSquare]  The statevector form is the basic form to express a state. It is encoded in the following form:\n"<>
"\[FilledSmallSquare]  The first part of the statevector encodes the charges of the SU(5) GUT charges. Each number correspond to the corresponding U(1) factor the charge appears,
and this number must not exceed the \"nU1s\" number specified in Options.
The first three entries are the charges for the three 10 reps.
The following six entries encode the first charge of the three 5 bar reps, followed by the second charge of the three 5 bar reps.
The final next 2 charges encode two from the bar 5 rep of Hd. The Hu charges are determined from the negative of the Hd charges.\n"<>
"\[FilledSmallSquare]  For example, a GUT chrage pattern of 10_(1),10_(2),10_(3),bar(5)_(1,4),bar(5)_(2,4),bar(5)_(3,4),Hd_(2,4) should have a state vector of the form: {1,2,3,1,2,3,4,4,4,2,4} for \"nU1s\"=4.\n"<>
"\[FilledSmallSquare]  The second part of the charge pattern encodes the moduli charges.\n"<>
"\[FilledSmallSquare]  The first section of this part encodes the charges of the singlet/perturbative moduli. This must be a charge vector of length-\"nSinglet\"*2, and each entry must be a positive integer not exceeding \"nU1s\".
The first \"nSinglet\" charges encode the first charge of the 1 reps and the last \"nSinglet\" charges encode the negative charge of the 1 reps.
As an example, if we have two perturbative moduli, with 1_(2,-3) and 1_(4,-1) the charge vector will be {2,4,3,1}. \n"<>
"\[FilledSmallSquare]  The second section of this moduli part encodes the charges of the non-perturbative moduli. 
This is a length-\"nNonperts\"*\"nU1s\" vector and exactly encodes the charges of these perturbative moduli \[CapitalPhi]=e^-t where t are the Kalher moduli.
Each entry must be an integer taking value from 1-2^(k-1) to 2^(k-1) (k an integer here), this range is specified by \"\[CapitalPhi]ChargeRange\" in the options (and should be ideally of this form).
For example, for a non-pert moduli of charge (-4,1,3,2) for 4 U1 factors the part of chargevector reads {-4,1,3,2}.\n"<>
"\[FilledSmallSquare]  The last part of the charge vector encodes the VEVs of the system. This is a integer vector of length-2*(\"nSinglet\"+\"nNonperts\").
The code reads the vector in the following form. For each \[Phi]-moduli, the code reads two digits of the vector and convert that to a complex number.
For each digit the VEV is encoded using the formula (sign)*\"\[Phi]VEVmin\"*\"\[Phi]VEVlogGen\"^(num) where sign is the sign of the digit of the integer, and num is value of the number. 
This procedure is the same for the non-perturbative part except that with the \[CapitalPhi] counterparts of the quantities mentioned above.\n"<>
"\[FilledSmallSquare]  The statevector is the vector formed from joining the MSSM charge, moduli charge and VEVs parts as described above. The length can be checked using the module GADetermDimStateVecwithVEVs.\n"<>
"\[FilledSmallSquare]  CompleteStateQsector takes the state in either bitlist or statevector form to obtain the full association. Use this to compute the properties of the system.\n"<>
"\[FilledSmallSquare]  Currently there is no reverse map built in the model to convert a statevector to a bitlist. This will be added in a later version.\n"<>
"\[FilledSmallSquare]  See ?stringFNmodelsAuxModules for the individual modules used in the code.";


(* ::Subsection::Closed:: *)
(*O1coeff Module help*)


SetO1Coeff::usage="SetO1Coeff[O1vec] sets the O1 coefficient to the inputted vector.\n"<>
"\[FilledSmallSquare]  O1vec -> This must be a vector of real values used as O1 coefficients.";

GetO1Coeff::usage="GetO1Coeff[] obtains the current O(1) coefficients in the settings.\n"<>
"\[FilledSmallSquare]  This gets the O1CLst global variable of the package.";

NewO1Coeff::usage="NewO1Coeff[] generates a vector of O(1) coefficients for the operators.\n"<>
"\[FilledSmallSquare]  This vector is returned as the output.";


(* ::Subsection::Closed:: *)
(*Auxiliary Module help*)


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - State Generation*)


RandomStatewithVEVs::usage="RandomStatewithVEVs[] generates a random state in a vector form.\n"<>
"\[FilledSmallSquare]  The first part of state vector encodes the charge of the system. This is previously the chargevector of the system.\n"<>
"\[FilledSmallSquare]  The second part of the state vector encodes the VEVs of the moduli fields.";

ConvertStatewithVEVs::usage="ConvertStatewithVEVs[stateVec] converts a state from vector to association form.\n"<>
"\[FilledSmallSquare]  stateVec must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  The output is an association form with the details of the state.";

StateVectoMatwithVEVs::usage="StateVectoMatwithVEVs[stateVec] converts a state from statevector form to a charge matrix where the rows are the species and the columns are the U(1) factors.\n"<>
"\[FilledSmallSquare]  stateVec must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  The output is a dimension (17+\"nSinglet\"+\"nNonperts\")*(\"nU1s\") matrix.";

StateVectoVEVs::usage="StateVectoVEVs[stateVec] finds the VEVs of the moduli fields.\n"<>
"\[FilledSmallSquare]  stateVec must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  The output is the moduli VEVs in substitution form, to be used later in the code.";

AssocFormToStateVec::usage="AssocFormToStateVec[assoc] converts the state from association form to state-vector form.\n"<>
"\[FilledSmallSquare]  assoc - this is an association with the following keys.\n"<>
"\[FilledSmallSquare]  \"Qq\" - this must be a (3*\"nU1s\") matrix. This is the charge of the Q in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qd\" - this must be a (3*\"nU1s\") matrix. This is the charge of the d in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qu\" - this must be a (3*\"nU1s\") matrix. This is the charge of the u in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"QL\" - this must be a (3*\"nU1s\") matrix. This is the charge of the L in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qe\" - this must be a (3*\"nU1s\") matrix. This is the charge of the e in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"QHd\" - this must be a (1*\"nU1s\") matrix. This is the charge of the Hd in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Q\[Phi]\" - this must be a (\"nSinglet\"*\"nU1s\") matrix. This is the charge of the bundle moduli in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Q\[CapitalPhi]\" this must be a (\"nNonperts\"*\"nU1s\") matrix. This is the charge of the Kahler moduli in the upstairs representation.\n"<>
"\[FilledSmallSquare]  \"\[Phi]VEVs\" -> this must be a dimension-(\"nSinglet\") vector. This is the VEVs for the bundle moduli. \n"<>
"\[FilledSmallSquare]  \"\[CapitalPhi]VEVs\" -> this must be a dimension-(\"nNonperts\") vector. This is the VEVs for the Kahler moduli. \n"<>
"\[FilledSmallSquare]  Note that \"Qq\", \"Qu\" and \"Qe\" should be equal since they are from the 10 rep in SU(5). Alternatively the package reads \"Q10\" as a 3-vector which just gives the charge-positions of 10s. \n"<>
"\[FilledSmallSquare]  Note that \"QL\" and \"Qd\" should be equal since they are from the 5 bar rep in SU(5). Alternatively the package reads \"Q5b\" as a (3*2)-matrix which gives the charge positions {e_a,e_b} for each 5 bars. \n"<>
"\[FilledSmallSquare]  Instead of \"QHd\" the package can also take in \"Q5Hb\", a 2-vector specifying the direction of the down-Higgs charges.\n"<>
"\[FilledSmallSquare]  Note that the 1 rep has charge pattern e_a-e_b. In this case there is a charge pattern penalty when the bundle moduli has zero charge in all U1s.\n"<>
"\[FilledSmallSquare]  In this case we will have to specify the charges using \"Q1\". This is a (\"nSinglet\"*2)-matrix which gives the charge positions {e_a,e_b} for each 1s. \n"<>
"\[FilledSmallSquare]  If the dimensions of the entered association is not matched the module will print an error message and exit.\n"<>
"\[FilledSmallSquare]  The output is the state in statevector form with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  If the VEVs coordinates are not integer-values they are rounded to the nearest integer.";

AssocFormToVEVs::usage="AssocFormToVEVs[assoc] converts the inputted state in association form into VEV substitution rules.\n"<>
"\[FilledSmallSquare]  assoc - this is an association with the following keys.\n"<>
"\[FilledSmallSquare]  \"\[Phi]VEVs\" - this must be a dimension-(\"nSinglet\") vector. This is the VEVs for the bundle moduli. \n"<>
"\[FilledSmallSquare]  \"\[CapitalPhi]VEVs\" - this must be a dimension-(\"nNonperts\") vector. This is the VEVs for the Kahler moduli. \n"<>
"\[FilledSmallSquare]  The output is the moduli VEVs in substitution form, to be used later in the code.\n"<>
"\[FilledSmallSquare]  There is no discretisation in this convertion, the VEV value is simply taken from the inputted (valid) VEVs in the association.";

ConvertStateVecToVEVPowers::usage="ConvertStateVecToVEVPowers[statevec] converts the inputted state-vector into an association with the powers of the VEVs.\n"<>
"\[FilledSmallSquare]  stateVec must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  The module converts the state-vector into an association with the powers of the VEVs.\n"<>
"\[FilledSmallSquare]  \"\[Phi]VEVPowers\" encode the powers of the bundle moduli VEVs.\n"<>
"\[FilledSmallSquare]  \"\[CapitalPhi]VEVPowers\" encode the powers of the Kahler moduli VEVs.\n";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - State Generation for GAbitlst*)


GADetermDimStateVecwithVEVs::usage="GADetermDimStateVecwithVEVs[] determines the dimension of the binary vector to be generated in GA.\n"<>
"\[FilledSmallSquare]  The output a dimension 3 vector.\n"<>
"\[FilledSmallSquare]  The first entry specifies the full dimension of the binary bitlist.\n"<>
"\[FilledSmallSquare]  The second entry specifies the dimension of the bitlist used to allocate 1 digit of the \[CapitalPhi] charge.\n"<>
"\[FilledSmallSquare]  The last entry specifies the dimension of the bitlist allocated for the VEVs.";

GARandomStatewithVEVs::usage="GARandomStatewithVEVs[] generates a random bitstring for GA module.\n"<>
"\[FilledSmallSquare]  This random bitstring generated is the state used in the genetic algorithm.\n"<>
"\[FilledSmallSquare]  The output is a bitlist with its dimension determined by the module GADetermDimStateVecwithVEVs.";

GAConvertBitlstStateVecwithVEVs::usage="GAConvertBitlstStateVecwithVEVs[bitlist] converts the GA bitlist state into statevector form to be used in the module.\n"<>
"\[FilledSmallSquare]  The bitlist inputted must be in binary form, with dimension as determined by the module GADetermDimStateVecwithVEVs.\n"<>
"\[FilledSmallSquare]  The output is the statevector with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"*\"bitsfor\[Phi]\"+\"nNonperts\"*\"bitsfor\[CapitalPhi]\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[CapitalPhi]\")-1)}.";

GADetermDimStateVecFix::usage="GADetermDimStateVecFix[] determines the dimension of the binary vector to be generated in GA.\n"<>
"\[FilledSmallSquare]  This module is similar to GADetermDimStateVecwithVEVs[] but accounted for the fact that the charges of the 10s and Higgs are fixed.\n"<>
"\[FilledSmallSquare]  The first entry specifies the full dimension of the binary bitlist.\n"<>
"\[FilledSmallSquare]  The second entry specifies the dimension of the bitlist used to allocate 1 digit of the \[CapitalPhi] charge.\n"<>
"\[FilledSmallSquare]  The last entry specifies the dimension of the bitlist allocated for the VEVs.";

GARandomStateFix::usage="GARandomStateFix[] generates a random bit-string for the GA module.\n"<>
"\[FilledSmallSquare]  This module assumes that the 10 and 5barH reps are fiexed by the Options.\n"<>
"\[FilledSmallSquare]  The output is a bitlist with its dimension determined by the module GADetermDimStateVecFix (the first entry of the output).";

GAConvertBitlstStateVecFix::usage="GAConvertBitlstStateVecFix[bitlist] converts the GA bitlist into statevector form to be used in the main module.\n"<>
"\[FilledSmallSquare]  The bitlist inputted must be in binary form, with the dimension determined by the first output of the module GADetermDimStateVecFix[].\n"<>
"\[FilledSmallSquare]  This module assumes that the 10 and 5barH reps are fiexed by the Options.\n"<>
"\[FilledSmallSquare]  The output is the statevector with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"*\"bitsfor\[Phi]\"+\"nNonperts\"*\"bitsfor\[CapitalPhi]\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[CapitalPhi]\")-1)}.";

GAConvertStateVectoBitlst::usage="GAConvertStateVectoBitlst[stateVec] converts a statevector into a GA bitlist.\n"<>
"\[FilledSmallSquare]  The input is the statevector with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"*\"bitsfor\[Phi]\"+\"nNonperts\"*\"bitsfor\[CapitalPhi]\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[CapitalPhi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The module converts the state-vector form into bitlist form, the map is not 1:1 so one of the bitlist is picked.\n"<>
"\[FilledSmallSquare]  This module checks whether the \"fixCharges\" option is turned on and output the bitlist of correct dimension.\n"<>
"\[FilledSmallSquare]  The output bitlist is in binary form, with the dimension determined by the first output of the module GADetermDimStateVec$$[].";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - Find Operators (Quark Sector)*)


QMatRepChange::usage="QMatRepChange[Qchargevec,nchargevec] change the charge vector representation from the upstairs to the downstairs one.\n"<>
"\[FilledSmallSquare]  Here the upstairs representation correspond to the representation with the q~q+n identification, as outlined in the note.\n"<>
"\[FilledSmallSquare]  Here the downstairs representation represents the representation with the redundancy removed.\n"<>
"\[FilledSmallSquare]  The output is the statevector representation one dimension lower than Qchargevec.";

FromPowertoOps::usage="FromPowertoOps[] exponentiates power vector form into operator form.\n"<>
"\[FilledSmallSquare]  This is an auxiliary module to be used in the find operator modules.\n"<>
"\[FilledSmallSquare]  The input is the power vector, a dimension (\"nSinglet\"+\"nNonperts\") integer vector to be exponentiated.\n"<>
"\[FilledSmallSquare]  The output is the required operator in the form \[Phi][i]^n_i and etc.";

FindOpsforWIJaListwith\[Phi]VEVsQsector::usage="FindOpsforWIJaListwith\[Phi]VEVsQsector[stateVec] find all the possible operators combinations in the Standard Model part of the superpotential from a given state.\n"<>
"\[FilledSmallSquare]  This submodule is used with the environment setting where the moduli VEVs are encoded in the state.\n"<>
"\[FilledSmallSquare]  This submodule ultilises the list method, listing out all available powers and checking whether the charges add up to zero.\n"<>
"\[FilledSmallSquare]  This submodule only computes the relevant terms in the quark sector.\n"<>
"\[FilledSmallSquare]  stateVec - This must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\".\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options\n"<>
"\[FilledSmallSquare]  The middle \"nNonperts\"*\"nU1s\" coordinates must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The last 2*(\"nSinglet\"*\"bitsfor\[Phi]\"+\"nNonperts\"*\"bitsfor\[CapitalPhi]\") digits are binary and encodes the VEVs of the moduli fields.\n"<>
"\[FilledSmallSquare]  The output gives the expansion of the operators in a {18,1} matrix as a list in each entry, which gives the respective operator expansion with the corresponding coefficients for each expansion term.\n"<>
"\[FilledSmallSquare]  If there are no operators then instead a 0 is outputted as the entry to the {18,1} matrix.";

UpdateVectorWithMinimum::usage="UpdateVectorWithMinimum[vec,posList] updates the inputted vec according to the posList which gives the duplicate entries.\n"<>
"\[FilledSmallSquare]  Inputs:\n"<>
"  \[FilledVerySmallSquare]  vec  - a dimension-\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) vector which gives the VEV powers of the system.\n"<>
"  \[FilledVerySmallSquare]  posList  - a nx2 matrix which gives the list of duplicate positions. This is if the moduli fields have the same charge patterns.\n"<>
"\[FilledSmallSquare]  The module checks the minimum of every pair of inputted positions and replace the larger VEV with the smaller VEV.\n"<>
"\[FilledSmallSquare]  Output - a dim-\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) vector giving the VEV powers of the system.\n"<>
"\[FilledSmallSquare]  Module is used in FindOpsFromGraphMethod.\n"<>
"\[FilledSmallSquare]  The module is to canonically assign the edge-weights to the solver.";

MinReplaceListRule::usage="MinReplaceListRule[vec,posList] gives a binary vector on whether there is a repeated Q entry being replaced in UpdateVectorWithMinimum.\n"<>
"\[FilledSmallSquare]  Inputs:\n"<>
"  \[FilledVerySmallSquare]  vec  - a dimension-\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) vector which gives the VEV powers of the system.\n"<>
"  \[FilledVerySmallSquare]  posList  - a nx2 matrix which gives the list of duplicate positions. This is if the moduli fields have the same charge patterns.\n"<>
"\[FilledSmallSquare]  The module checks the minimum of every pair of inputted positions and replace the larger VEV with the smaller VEV.\n"<>
"\[FilledSmallSquare]  Output - a dim-\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) vector giving a binary code on whether a field is being replaced.\n";

FindOpsFromGraphMethod::usage="FindOpsFromGraphMethod[\[Phi]Qs,\[Phi]VEVs,QToSolve] aims to find the leading order operator insertion given the charges and VEVs of the singlet moduli fields to the inputted matched charge.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  \[Phi]Qs - the charges of the singlet operators. This is a dim-(2\!\(\*SubscriptBox[\(xn\), \(\[Phi]\)]\)) vector encoded in the same manner as the state-vector: \!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) positive charges followed by the \!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) negative charges. This encodes the flow directions (edges) in the graph.\n"<>
"  \[FilledVerySmallSquare]  \[Phi]VEVs - the VEvs of the singlet operators. This is a dim-\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) vector encoding the VEV powers in the same manner as the state-vector. This encodes the edge weights in the graph.\n"<>
"  \[FilledVerySmallSquare]  QToSolve - a dim-|n| vector encoding the charge to be matched by the singlet insertions.\n"<>
"\[FilledSmallSquare]  Module uses the Minimum-Cost-Maximum-Flow algorithm to find the required operator.\n"<>
"\[FilledSmallSquare]  Output - a dim-\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) vector encoding the insertion powers of the required insertion. If no operators are returned Null is returned.";

FindOpsforWIJaGraphwith\[Phi]VEVsQsector::usage="FindOpsforWIJaGraphwith\[Phi]VEVsQsector[stateVec] find all the possible operators combinations in the Standard Model part of the superpotential from a given state.\n"<>
"\[FilledSmallSquare]  This submodule is used with the environment setting where the moduli VEVs are encoded in the state.\n"<>
"\[FilledSmallSquare]  This submodule ultilises the graph method, following the Minimum-Cost-Maximum-Flow problem.\n"<>
"\[FilledSmallSquare]  This submodule only computes the relevant terms in the quark sector.\n"<>
"\[FilledSmallSquare]  stateVec - This must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\".\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options\n"<>
"\[FilledSmallSquare]  The middle \"nNonperts\"*\"nU1s\" coordinates must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The last 2*(\"nSinglet\"*\"bitsfor\[Phi]\"+\"nNonperts\"*\"bitsfor\[CapitalPhi]\") digits are binary and encodes the VEVs of the moduli fields.\n"<>
"\[FilledSmallSquare]  The output gives the expansion of the operators in a {18,1} matrix as a list in each entry, which gives the respective operator expansion with the corresponding coefficients for each expansion term.\n"<>
"\[FilledSmallSquare]  If there are no operators then instead a 0 is outputted as the entry to the {18,1} matrix.";

GenerateOrderedTuples::usage="GenerateOrderedTuples[sum, length] generates combinations of integer-tuples up to the sum of length.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  sum - integer, the maximum of the L1 norm of the vectors.\n"<>
"  \[FilledVerySmallSquare]  length - integer, the length of the vectors.\n"<>
"\[FilledSmallSquare]  Output: list of vectors ordered by sum_i sum_j=i.";

SumChargesByPowers::usage="SumChargesByPowers[kPows,kQs] sum the charges according to the power.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  kPows - dim-\!\(\*SubscriptBox[\(n\), \(\[CapitalPhi]\)]\) integer vector, the power of the \[CapitalPhi] vector.\n"<>
"  \[FilledVerySmallSquare]  kQs - dim-(\!\(\*SubscriptBox[\(n\), \(\[CapitalPhi]\)]\)x|n|) integer matrix, the charge matrix of \[CapitalPhi] moduli.\n"<>
"\[FilledSmallSquare]  Output: dim-\!\(\*SubscriptBox[\(n\), \(\[CapitalPhi]\)]\) integer vector that is the charge of the \[CapitalPhi] insertion.";

ProcessModPowers::usage="ProcessModPowers[vec1,vec2] joins the two vectors element-wise and delete the Null ones.\n"<>
"\[FilledSmallSquare]  Inputs: vec1 and vec2 must have the same length with vectors as elements.\n"<>
"\[FilledSmallSquare]  Module checks whether vec1 has Null elements. If not then the corresponding entry in vec2 is joined to the element in vec1.\n"<>
"\[FilledSmallSquare]  Inputs: Returns a vector of vectors.";

FindOpsFromGraphAdd\[CapitalPhi]::usage="FindOpsFromGraphAdd\[CapitalPhi][\[Phi]Qs,\[Phi]VEVs,QToSolve,\[CapitalPhi]Qs] computes the operators using the graph method with the non-perturbative insertions.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  \[Phi]Qs - the charges of the singlet operators. This is a dim-(2\!\(\*SubscriptBox[\(xn\), \(\[Phi]\)]\)) vector encoded in the same manner as the state-vector: \!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) positive charges followed by the \!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) negative charges. This encodes the flow directions (edges) in the graph.\n"<>
"  \[FilledVerySmallSquare]  \[Phi]VEVs - the VEvs of the singlet operators. This is a dim-\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\) vector encoding the VEV powers in the same manner as the state-vector. This encodes the edge weights in the graph.\n"<>
"  \[FilledVerySmallSquare]  QToSolve - a dim-|n| vector encoding the charge to be matched by the singlet insertions.\n"<>
"  \[FilledVerySmallSquare]  \[CapitalPhi]Qs - a dim-(\!\(\*SubscriptBox[\(n\), \(\[CapitalPhi]\)]\)x|n|) matrix encoding the charges of the effective kahler operators (non-perturbative) fields. This is the vector encoded in the same manner as the \[CapitalPhi]Qs-part in the state-vector.\n"<>
"\[FilledSmallSquare]  The module computes the leading operator insertion for each non-perturbative effective field combination using the Minimum-Cost-Maximum-Flow algorithm.\n"<>
"\[FilledSmallSquare]  Output: A list of dim-(\!\(\*SubscriptBox[\(n\), \(\[Phi]\)]\)+\!\(\*SubscriptBox[\(n\), \(\[CapitalPhi]\)]\)) vectors of the operator insertions.";

ComputeWOpsGivenStateVec::usage="ComputeWOpsGivenStateVec[statevec] computes the WOps for the Yukawa couplings given the state vector as an input.\n"<>
"\[FilledSmallSquare]  Depending on the settings for the method of generating operators (\"MethodGenOps\"), this module will choose the correct sub-module to compute the operators.\n"<>
"\[FilledSmallSquare]  For example, choosing \"MethodGenOps\"->\"Graph\" will use the module FindOpsforWIJaGraphwith\[Phi]VEVsQsector.\n"<>
"\[FilledSmallSquare]  And choosing \"MethodGenOps\"->\"List\" will use the module FindOpsforWIJaListwith\[Phi]VEVsQsector.\n"<>
"\[FilledSmallSquare]  The output is the expansion of the operators as a list (in each entry) in a {18,1} matrix, which gives the respective operator expansion with the corresponding coefficients for each expansion term.\n"<>
"\[FilledSmallSquare]  If there are no operators then instead a 0 is outputted as the entry to the {18,1} matrix.";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - Texture Analysis*)


FromWOpsGetLeadingYukawaMatrix::usage="FromWOpsGetLeadingYukawaMatrix[WOps,VEVsPowerListAssoc] computes the leading order powers of the Yukawa matrix entries.\n"<>
"\[FilledSmallSquare]  WOps -> This is a dimension {18,1}-matrix with each entry either being a 0 (no allowed operators) or a list of allowed operators.\n"<>
"\[FilledSmallSquare]  This is nicely the output of the submodule ComputeWOpsGivenStateVec[stateVec].\n"<>
"\[FilledSmallSquare]  VEVsPowerListAssoc -> This is an association with the keys \"\[Phi]VEVPowers\" and \"\[CapitalPhi]VEVPowers\" that encodes the powers of the VEVs of the matrix.\n"<>
"\[FilledSmallSquare]  This is the output of the module ConvertStateVecToVEVPowers[statevec].\n"<>
"\[FilledSmallSquare]  This module computes the leading order contribution in each Yukawa entry assuming that the scale of the VEVs of moduli is set by some number \[Epsilon]<1 and the order is determined by the inputted power list.\n"<>
"\[FilledSmallSquare]  It outputs the list of leading order contribution in symbolic form, using the symbol defined by the option \"VEVsScale\".\n"<>
"\[FilledSmallSquare]  The output is an association with the keys \"YukdLead\" and \"YukuLead\" which are {3,3}-matrices determining the leading order term in powers of the scale.\n"<>
"\[FilledSmallSquare]  The symbolic form is for easy visualisation for analysis purposes.";

LowestOrderOfMonomialListOneVar::usage="LowestOrderOfMonomialListOneVar[monomialList,var] finds the lowest order monomial from an inputted list for one variable.\n"<>
"\[FilledSmallSquare]  monomialList -> A list of monomials {... , ...} in the variable var.\n"<>
"\[FilledSmallSquare]  var -> some symbol.\n"<>
"\[FilledSmallSquare]  This is a submodule to aid the calculation is finding the leading order contribution in the Yukawa entries in FromWOpsGetLeadingYukawaMatrix.\n"<>
"\[FilledSmallSquare]  The output is ONE of the leading order terms in the list of monomials.\n"<>
"\[FilledSmallSquare]  i.e. if there are multiple terms of the same order this submodule will output the first one in the list that satisfies the constraints.";

RemoveCoeffFromMonomial::usage="RemoveCoeffFromMonomial[monomial] removes coefficient from the monomial.\n"<>
"\[FilledSmallSquare]  Input is monomial.\n"<>
"\[FilledSmallSquare]  Output of submodule is monomial without the coefficient.";

MassHierarchyFromLeadingOrders::usage="MassHierarchyFromLeadingOrders[YukLeadPowAssoc] computes the mass hierarchy from the leading order powers of the two Yukawa matrices.\n"<>
"\[FilledSmallSquare]  YukLeadPowAssoc -> This is an association with the keys \"\[Phi]VEVPowers\" and \"\[CapitalPhi]VEVPowers\" that encodes the powers of the VEVs of the matrix.\n"<>
"\[FilledSmallSquare]  Each entry must be a (3,3)-matrix with entries of the form scale^a where scale is the symbol defined by \"VEVsScale\" and a is a number.\n"<>
"\[FilledSmallSquare]  The input is naturally the output of the module of FromWOpsGetLeadingYukawaMatrix.\n"<>
"\[FilledSmallSquare]  This module assumes the eigenvalues are distinct in their hierarchy (no two eigenvalues are of the same order) and computes the order of the eigenvalues of the two Yukawa matrices based on this assumption.\n"<>
"\[FilledSmallSquare]  It then outputs the relevant scales in the form of an association.\n"<>
"\[FilledSmallSquare]  The list of the keys for the output association is as follows:\n"<>
"\[FilledSmallSquare]  \"YukUpEigens\" -> The leading-order eigenvalue powers of the up Yukawa matrix.\n"<>
"\[FilledSmallSquare]  \"YukDownEigens\" -> The leading-order eigenvalue powers of the down Yukawa matrix.\n"<>
"\[FilledSmallSquare]  \"mtPow\" -> The leading-order eigenvalue from the up Yukawa matrix.\n"<>
"\[FilledSmallSquare]  \"mbPow\" -> The leading-order eigenvalue from the down Yukawa matrix.\n"<>
"\[FilledSmallSquare]  \"mu/mcPow\" -> The leading-order scale-difference between the up and charm quarks (eigenvalues).\n"<>
"\[FilledSmallSquare]  \"mc/mtPow\" -> The leading-order scale-difference between the charm and top quarks (eigenvalues).\n"<>
"\[FilledSmallSquare]  \"md/msPow\" -> The leading-order scale-difference between the down and strange quarks (eigenvalues).\n"<>
"\[FilledSmallSquare]  \"ms/mbPow\" -> The leading-order scale-difference between the strange and bottom quarks (eigenvalues).\n"<>
"\[FilledSmallSquare]  \"YukUpMatLead\" -> This is a {3,3}-matrix listing the leading-order term of the up Yukawa matrix.\n"<>
"\[FilledSmallSquare]  \"YukDownMatLead\" -> This is a {3,3}-matrix listing the leading-order term of the down Yukawa matrix.\n"<>
"\[FilledSmallSquare]  \"PosOTop\" -> Diagonal entry position that contributes the top eigenvalue in the up-Yukawa matrix.\n"<>
"\[FilledSmallSquare]  \"PosOTop\" -> Diagonal entry position that contributes the bottom eigenvalue in the down-Yukawa matrix.";

IdentifyLowestPowerPos::usage="IdentifyLowestPowerPos[ListExp] identifies the position of the smallest scale in the WijOp entry.\n"<>
"\[FilledSmallSquare]  Input a list {\[Epsilon]^a,\[Epsilon]^b,...} and the output is the location of the lowest power of scale in the list.\n"<>
"\[FilledSmallSquare]  If there are more than positions returned the first is returned canonically.\n"<>
"\[FilledSmallSquare]  If 0 is inputted 1 is returned.";

LocateExtractO1Coeffs::usage="LocateExtractO1Coeffs[WOpsOut,VEVPowersAssoc,MHAssoc,O1coeffList] locates the leading order positons of the O1 coefficients and the simultaneously the positions of the O(1) factors for the top and bottom couplings in O1coeffList.\n"<>
"\[FilledSmallSquare]  Inputs:\n"<>
"  \[FilledVerySmallSquare]  WOpsOut - a dimension {18,1}-matrix with each entry either being a 0 (no allowed operators) or a list of allowed operators.\n"<>
"  \[FilledVerySmallSquare]  It is the output of the module ComputeWOpsGivenStateVec[stateVec].\n"<>
"  \[FilledVerySmallSquare]  VEVPowersAssoc - an association with keys \"\[Phi]VEVPowers\" and \"\[CapitalPhi]VEVPowers\" that gives the powers of the VEVs of the moduli fields.\n"<>
"  \[FilledVerySmallSquare]  It is the output of the module ConvertStateVecToVEVPowers[stateVec].\n"<>
"  \[FilledVerySmallSquare]  MHAssoc - the mass Hierarchy association. It is naturally is the output of MassHierarchyFromLeadingOrders[YukLeadPowAssoc].\n"<>
"  \[FilledVerySmallSquare]  Only two keys are needed \"PosOTop\" and \"PosOBot\" which gives the position of the top and bottom entry along the diagonal of the respective Yukawa matrices (takes values from 1,2,3).\n"<>
"  \[FilledVerySmallSquare]  O1coeffList - a list of O(1) coefficients which must have a larger dimension then the number of operator insertions.\n"<>
"\[FilledSmallSquare]  Module computes the following outputs as an association: \n"<>
"  \[FilledVerySmallSquare]  \"O1StartPosList\" - the starting O(1) positions of each of the Yukawa matrices. \n"<>
"  \[FilledVerySmallSquare]  \"O1LeadPosList\" - the positions of O(1) corresponding to the leading term in each Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadList\" - the O(1) coefficients corresponding to the leading term in each Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadUpMat\" - the O(1) coefficients corresponding to the leading term in each up-Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadDownMat\" - the O(1) coefficients corresponding to the leading term in each down-Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"OTopVecPos\" - the position of the O(1) coefficient corresponding to the term of the top Yukawa coupling. \n"<>
"  \[FilledVerySmallSquare]  \"OBotVecPos\" - the position of the O(1) coefficient corresponding to the term of the bottom Yukawa coupling.";

CalcScaleAndHiggsO1Coeff::usage="CalcScaleAndHiggsO1Coeff[massHierarchyAssoc] computes the scale and the O(1) coefficient from the leading order data.\n"<>
"\[FilledSmallSquare]  Inputs:\n"<>
"  \[FilledVerySmallSquare]  massHierarchyAssoc - This is an association that contains the leading-order data of the Yukawa matrices.\n"<>
"  \[FilledVerySmallSquare]  This is naturally the output of the module MassHierarchyFromLeadingOrders.\n"<>
"  \[FilledVerySmallSquare]  The required list of the keys for the input association is as follows:\n"<>
"  \[FilledVerySmallSquare]  \"mtPow\" - The leading-order eigenvalue from the up Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"mbPow\" - The leading-order eigenvalue from the down Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"mu/mcPow\" - The leading-order scale-difference between the up and charm quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"mc/mtPow\" - The leading-order scale-difference between the charm and top quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"md/msPow\" - The leading-order scale-difference between the down and strange quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"ms/mbPow\" - The leading-order scale-difference between the strange and bottom quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"YukUpMatLead\" - This is a {3,3}-matrix listing the leading-order tmer of the up Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"YukDownMatLead\" - This is a {3,3}-matrix listing the leading-order tmer of the down Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"PosOTop\" - Diagonal entry position that contributes the top eigenvalue in the up-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"PosOBot\" - Diagonal entry position that contributes the bottom eigenvalue in the down-Yukawa matrix.\n"<>
"\[FilledSmallSquare]  The submodule then computes the best scale and the O(1) coefficient associated to the top mass by fixing the scale of Electroweak Breaking.\n"<>
"\[FilledSmallSquare]  The best scale is found by computing the geometric mean of the five quark mass scales as described above.\n"<>
"\[FilledSmallSquare]  The O(1) coefficient associated with the bottom mass is set to 1 (by fixing the scale).\n"<>
"\[FilledSmallSquare]  The O(1) coefficient associated with the top mass is solved from the triangular Higgs relation. This coefficient can be complex.\n"<>
"\[FilledSmallSquare]  The output is the input association with the following extra keys:\n"<>
"  \[FilledVerySmallSquare]  \"ListOfScales\" -> The list of scales between the quark masses solved after comparing to the quark hierarchies. They are canonically in the order: {mu/mc,mc/mt,md/ms,ms/mb}.\n"<>
"  \[FilledVerySmallSquare]  \"BestScale\" -> This is best scale (geometric mean) by fixing the scales to the actual quark hierarchies.\n"<>
"  \[FilledVerySmallSquare]  \"Higgs\" -> This is a {2,1}-matrix encoding <Hu> and <Ht> without the O(1) coefficient.\n"<>
"  \[FilledVerySmallSquare]  \"BestO1Top\" -> This is the best O(1) coefficient associated to the top mass (diagonal entry in Yukawa matrix that gives the top-mass).\n"<>
"  \[FilledVerySmallSquare]  This O(1) coefficient can be complex (purely real or imaginary).\n"<>
"  \[FilledVerySmallSquare]  If the O(1) coefficient is real then the EW breaking scale can be fixed.\n"<>
"  \[FilledVerySmallSquare]  If the O(1) coefficient is imaginary then the EW breaking scale cannot be fixed. The model is not viable.";

GeometricStandardDeviation::usage="GeometricStandardDeviation[list] computes the geometric SD of the inputted list.\n"<>
"\[FilledSmallSquare]  The module returns a number.";

FitnessTexture::usage="FitnessTexture[massHierarchyAssoc,VEVPowers] computes the fitness part from the Texture.\n"<>
"\[FilledSmallSquare]  Any contributions coming from this fitness module signals the deviation from the \"Froggart-Nielsen\" spirit model building.\n"<>
"\[FilledSmallSquare]  Inputs:\n"<>
"  \[FilledVerySmallSquare]  massHierarchyAssoc - This is an association that contains the leading-order data of the Yukawa matrices.\n"<>
"  \[FilledVerySmallSquare]  This is naturally the output of the module CalcScaleAndHiggsO1Coeff[mHAssoc].\n"<>
"  \[FilledVerySmallSquare]  The required list of the keys for the input association is as follows:\n"<>
"  \[FilledVerySmallSquare]  \"mu/mcPow\" - The leading-order scale-difference between the up and charm quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"mc/mtPow\" - The leading-order scale-difference between the charm and top quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"md/msPow\" - The leading-order scale-difference between the down and strange quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"ms/mbPow\" - The leading-order scale-difference between the strange and bottom quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"BestScale\" - This is best scale (geometric mean) by fixing the scales to the actual quark hierarchies.\n"<>
"  \[FilledVerySmallSquare]  \"BestO1Top\" - This is the best O(1) coefficient associated to the top mass (diagonal entry in Yukawa matrix that gives the top-mass).\n"<>
"  \[FilledVerySmallSquare]  \"ListOfScales\" - The list of scales between the quark masses solved after comparing to the quark hierarchies. They are canonically in the order: {mu/mc,mc/mt,md/ms,ms/mb}.\n"<>
"  \[FilledVerySmallSquare]  \"HiggsSquared\" - This is a {2,1}-matrix encoding <Hu>^2 and <Hd> ^2 without the O(1) coefficient.\n"<>
"  \[FilledVerySmallSquare]  VEVPowersAssoc - an association with keys \"\[Phi]VEVPowers\" and \"\[CapitalPhi]VEVPowers\" that gives the powers of the VEVs of the moduli fields.\n"<>
"  \[FilledVerySmallSquare]  It is the output of the module ConvertStateVecToVEVPowers[stateVec].\n"<>
"\[FilledSmallSquare]  There are four contributions to the Fitness.\n"<>
"\[FilledSmallSquare]  1. If the powers of the mass ratios are not positive (eigenvalues are not-correctly ordered).\n"<>
"\[FilledSmallSquare]  2. The variance (standard deviation) of the scale calculated from the five different mass ratios.\n"<>
"\[FilledSmallSquare]  3. Scale penalty - if the scale is not smaller than 1 (super-Planckian masses) then a Log contribution is added.\n"<>
"\[FilledSmallSquare]  4. O(1)-top penalty - the O(1) coefficient must be order-1. This is another log-contribution.\n"<>
"\[FilledSmallSquare]  The output is an association with the related Fitness values.\n"<>
"\[FilledSmallSquare]  The output association has the following keys:\n"<>
"  \[FilledVerySmallSquare]  \"FitnessTextSec\" -> This is the total fitness contribution from texture analysis.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessPowers\" -> Fitness from negative powers in scales.\n"<>
"  \[FilledVerySmallSquare]  \"VarsScale\" -> Fitness from variance of calculated scales.\n"<>
"  \[FilledVerySmallSquare]  \"ScalePenalty\" -> Penalty from super-Planckian scales.\n"<>
"  \[FilledVerySmallSquare]  \"OTopPenalty\" -> Penalty from O(1).\n"<>
"  \[FilledVerySmallSquare]  The total fitness is calculated with the scale factor (multiplicative) as specified by the option \"fitnessTexturefactor\" to the four contributions.\n"<>
"  \[FilledVerySmallSquare]  The O(1) coefficient penalty only takes into account the magnitude of the O(1) coefficent. The EW-breaking scale will be correctly asccounted in the FitnessHiggs module.";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - Operators to Yukawa Matrices*)


AdjoinO1coeffWOpsTextQSec::usage="AdjoinO1coeffWOpsTextQSec[WOps,MHAssoc,VEVsPowers,O1Assoc,O1Vec] adjoins the O(1) coefficients to the operators to obtain the numerical Yukawa matrices.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  WOps -> This is a dimension {18,1}-matrix with each entry either being a 0 (no allowed operators) or a list of allowed operators.\n"<>
"  \[FilledVerySmallSquare]  This is nicely the output of the submodule ComputeWOpsGivenStateVec[stateVec].\n"<>
"  \[FilledVerySmallSquare]  massHierarchyAssoc - This is an association that contains the leading-order data of the Yukawa matrices.\n"<>
"  \[FilledVerySmallSquare]  This is naturally the output of the module CalcScaleAndHiggsO1Coeff[mHAssoc].\n"<>
"  \[FilledVerySmallSquare]  The required list of the keys for this association is as follows:\n"<>
"  \[FilledVerySmallSquare]  \"BestScale\" - This is best scale (geometric mean) by fixing the scales to the actual quark hierarchies.\n"<>
"  \[FilledVerySmallSquare]  \"BestO1Top\" - This is the best O(1) coefficient associated to the top mass (diagonal entry in Yukawa matrix that gives the top-mass).\n"<>
"  \[FilledVerySmallSquare]  \"PosOTop\" -> Diagonal entry position that contributes the top eigenvalue in the up-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"PosOTop\" -> Diagonal entry position that contributes the bottom eigenvalue in the down-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare] VEVsPowers -> This is an association with the keys \"\[Phi]VEVPowers\" and \"\[CapitalPhi]VEVPowers\" that encodes the powers of the VEVs of the matrix.\n"<>
"  \[FilledVerySmallSquare] This is the output of the module ConvertStateVecToVEVPowers[statevec].\n"<>
"  \[FilledVerySmallSquare]  O1Assoc - is the O(1) association for O(1) coefficients.\n"<>
"  \[FilledVerySmallSquare]  This is the output of the module LocateExtractO1Coeffs[WOpsOut,VEVPowersAssoc,MHAssoc,O1coeffList]. \n"<>
"  \[FilledVerySmallSquare]  The keys required is \"O1LeadList\" - the O(1) coefficients corresponding to the leading term in each Yukawa matrix element.\n"<>
"\[FilledSmallSquare]  The module substitutes the O(1) coefficients and VEVs into the operators.\n"<>
"\[FilledSmallSquare]  The top and bottom O(1) coefficients are automatically updated.\n"<>
"\[FilledSmallSquare]  The output is an association with the following keys.\n"<>
"  \[FilledVerySmallSquare]  \"YdMatVal\" - a 3x3 matrix giving the down-Yukawa numerical matrix.\n"<>
"  \[FilledVerySmallSquare]  \"YuMatVal\" - a 3x3 matrix giving the up-Yukawa numerical matrix.\n"<>
"  \[FilledVerySmallSquare]  \"YdMatOps\" - a 3x3 matrix giving the down-Yukawa matrix of operators.\n"<>
"  \[FilledVerySmallSquare]  \"YuMatOps\" - a 3x3 matrix giving the up-Yukawa matrix of operators.\n"<>
"  \[FilledVerySmallSquare]  \"ytMZVal\" - the top Yukawa coupling at \!\(\*SubscriptBox[\(M\), \(Z\)]\)-scale.\n"<>
"  \[FilledVerySmallSquare]  \"ybMZVal\" - the bottom Yukawa coupling at \!\(\*SubscriptBox[\(M\), \(Z\)]\)-scale.\n"<>
"  \[FilledVerySmallSquare]  \"O1StartPosList\" - the starting O(1) positions of each of the Yukawa matrices. \n"<>
"  \[FilledVerySmallSquare]  \"O1LeadPosList\" - the positions of O(1) corresponding to the leading term in each Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadList\" - the updated O(1) coefficients corresponding to the leading term in each Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadUpMat\" - the updated O(1) coefficients corresponding to the leading term in each up-Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadDownMat\" - the updated O(1) coefficients corresponding to the leading term in each down-Yukawa matrix element.\n"<>
"  \[FilledVerySmallSquare]  \"OTopVecPos\" - the position of the O(1) coefficient corresponding to the term of the top Yukawa coupling. \n"<>
"  \[FilledVerySmallSquare]  \"OBotVecPos\" - the position of the O(1) coefficient corresponding to the term of the bottom Yukawa coupling.\n"<>
"  \[FilledVerySmallSquare]  \"PosUpLead\" - Diagonal entry position that contributes the top eigenvalue in the up-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"PosDownLead\" - Diagonal entry position that contributes the bottom eigenvalue in the down-Yukawa matrix.";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - RG Computations*)


RGGivenYukMZFindYukMGUT::usage="RGGivenYukMZFindYukMGUT[ytMZ,ybMZ] computes the RG evolution from IR to UV.\n"<>
"\[FilledSmallSquare]  ytMZ and ybMZ are the values of the yu33 and yd33 entries at low energies.\n"<>
"\[FilledSmallSquare]  The evolution is computed by assuming that the electroweak contribution is negliglble, i.e. only g3 effects are considered.\n"<>
"\[FilledSmallSquare]  Also it is assumed that only the y33 entries matter in both Yukawa matrices, and ye33 is 0.\n"<>
"\[FilledSmallSquare]  Returns {ytMGUT,ybMGUT}, the Yukawa y33 entries at \!\(\*SubscriptBox[\(M\), \(\(GUT\)\(\\\ \)\)]\)scale.";

RGGivenYukMGUTFindYukMZ::usage="RGGivenYukMGUTFindYukMZ[ytMGUT,ybMGUT] computes the RG evolution from UV to IR.\n"<>
"\[FilledSmallSquare]  ytMGUT and ybMGUT are the values of the yu33 and yd33 entries at low energies.\n"<>
"\[FilledSmallSquare]  The evolution is computed by assuming that the electroweak contribution is negliglble, i.e. only g3 effects are considered.\n"<>
"\[FilledSmallSquare]  Also it is assumed that only the y33 entries matter in both Yukawa matrices, and ye33 is 0.\n"<>
"\[FilledSmallSquare]  Returns {ytMZ,ybMZ}, the Yukawa y33 entries at \!\(\*SubscriptBox[\(M\), \(\(Z\)\(\\\ \)\)]\)scale.";

RGGivenYukMZFindYukRatiosAll::usage="RGGivenYukMZFindYukRatiosAll[ytMZ,ybMZ] computes the RG evolution from IR to UV.\n"<>
"\[FilledSmallSquare]  ytMZ and ybMZ are the values of the yu33 and yd33 entries at low energies.\n"<>
"\[FilledSmallSquare]  The evolution is computed by assuming that the electroweak contribution is negliglble, i.e. only g3 effects are considered.\n"<>
"\[FilledSmallSquare]  Also it is assumed that only the y33 entries matter in both Yukawa matrices, and ye33 is 0.\n"<>
"\[FilledSmallSquare]  Returns association describing the ratios between the yukawa couplings, where the ratios are \!\(\*SubscriptBox[\(\[Lambda]\), \(GUT\)]\)/\!\(\*SubscriptBox[\(\[Lambda]\), \(Z\)]\).\n"<>
"\[FilledSmallSquare]  The keys are of the form \"yxnnR\", where x indicates up or down and nn indicates the position.,";

RGGivenYukGUTFindYukRatiosAll::usage="RGGivenYukGUTFindYukRatiosAll[ytMGUT,ybMGUT] computes the RG evolution from UV to IR.\n"<>
"\[FilledSmallSquare]  ytMGUT and ybMGUT are the values of the yu33 and yd33 entries at high energies.\n"<>
"\[FilledSmallSquare]  The evolution is computed by assuming that the electroweak contribution is negliglble, i.e. only g3 effects are considered.\n"<>
"\[FilledSmallSquare]  Also it is assumed that only the y33 entries matter in both Yukawa matrices, and ye33 is 0.\n"<>
"\[FilledSmallSquare]  Returns association describing the ratios between the yukawa couplings, where the ratios are \!\(\*SubscriptBox[\(\[Lambda]\), \(Z\)]\)/\!\(\*SubscriptBox[\(\[Lambda]\), \(GUT\)]\).\n"<>
"\[FilledSmallSquare]  The keys are of the form \"yxnnR\", where x indicates up or down and nn indicates the position.,";

RGGivenYukMZFindYukRatiosAllFixedPlane::usage="RGGivenYukMZFindYukRatiosAllFixedPlane[ytMZ,ybMZ] computes the RG evolution from IR to UV.\n"<>
"\[FilledSmallSquare]  ytMZ and ybMZ are the values of the yu33 and yd33 entries at low energies.\n"<>
"\[FilledSmallSquare]  The evolution is computed by assuming that the electroweak contribution is negliglble, i.e. only g3 effects are considered.\n"<>
"\[FilledSmallSquare]  Also it is assumed that only the y33 entries matter in both Yukawa matrices, and ye33 is 0.\n"<>
"\[FilledSmallSquare]  The ratio is approximated by finding a plane in the \!\(\*SubscriptBox[\(\[Lambda]\), \(t33, MZ\)]\)-\!\(\*SubscriptBox[\(\[Lambda]\), \(b33, MZ\)]\)-\!\(\*SubscriptBox[\(\[Lambda]\), \(MZ\)]\)/\!\(\*SubscriptBox[\(\[Lambda]\), \(MGUT\)]\) Space.\n"<>
"\[FilledSmallSquare]  Returns association describing the ratios between the yukawa couplings, where the ratios are \!\(\*SubscriptBox[\(\[Lambda]\), \(GUT\)]\)/\!\(\*SubscriptBox[\(\[Lambda]\), \(Z\)]\).\n"<>
"\[FilledSmallSquare]  The keys are of the form \"yxnnR\", where x indicates up or down and nn indicates the position.,";

ConstructRGFactors::usage="ConstructRGFactors[y33Factor,y3iFactor,yijFactor,leadPos] computes the RG factor form (to be multiplied to the GUT Yukawa matrix.\n"<>
"\[FilledSmallSquare]  y33Factor,y3iFactor,yijFactor are the \!\(\*SubscriptBox[\(\[Lambda]\), \(MZ\)]\)/\!\(\*SubscriptBox[\(\[Lambda]\), \(MGUT\)]\) factors for the appropriate positions respectively.\n"<>
"\[FilledSmallSquare]  leadPos must be 1,2 or 3 and gives the leading order on the diagonal.\n";
"\[FilledSmallSquare]  Returns a 3x3 matrix.";

GivenYukMatsFindRGMatricesPlane::usage="GivenYukMatsFindRGMatricesPlane[yukMZAssoc] computes the RG matrix to be multiplied to the GUT Yukawa matrices.\n"<>
"\[FilledSmallSquare]  Input: yukMZAssoc is an association that has the keys \"ytMZVal\" and \"ybMZVal\".\n"<>
"\[FilledSmallSquare]  It also must have the keys \"PosDownLead\", \"PosUpLead\" giving the diagonal leading order term in the respective Yukawa matrices.\n";
"\[FilledSmallSquare]  These are the values of the top- and bottom-Yukawa couplings and is computed in module AdjoinO1coeffWOpsTextQSec.\n"<>
"\[FilledSmallSquare]  The module computes the relevant RG Yukawa matrices which consists of entries \!\(\*SubscriptBox[\(\[Lambda]\), \(MZ\)]\)/\!\(\*SubscriptBox[\(\[Lambda]\), \(MGUT\)]\) in the appropriate entries.\n"<>
"\[FilledSmallSquare]  The output is the relevant RG matrices in asssociation form, with the keys \"RGUpMat\" and \"RGDownMat\".";

O1LogDevPenal::usage="O1LogDevPenal[num] test whether the number is O(1) which range is specified in the Options.\n"<>
"\[FilledSmallSquare]  If the num is in range the module returns 0, otherwise it returns the logarithmic deviation away from it.";

FitnessRG::usage="FitnessRG[YukO1Assoc] computes the relevant fitness in the RG sector.\n"<>
"\[FilledSmallSquare]  Input: YukO1Assoc is an association that is the output of AdjoinO1coeffWOpsTextQSec and must contain the following keys:\n"<>
"  \[FilledVerySmallSquare]  \"ytMZVal\" - leading order diagonal Yukawa term in up.\n"<>
"  \[FilledVerySmallSquare]  \"ybMZVal\" - leading order diagonal Yukawa term in down.\n"<>
"  \[FilledVerySmallSquare]  \"PosUpLead\" - position along diagonal identifying the leading order eigenvalue in up.\n"<>
"  \[FilledVerySmallSquare]  \"PosDownLead\" - position along diagonal identifying the leading order eigenvalue in down.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadUpMat\" - 3x3 matrix giving the O(1) coefficients for the leading order terms in the up-yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadDownMat\" - 3x3 matrix giving the O(1) coefficients for the leading order terms in the down-yukawa matrix.\n"<>
"\[FilledSmallSquare]  This is naturally the output of the module MassHierarchyFromLeadingOrders.\n"<>
"\[FilledSmallSquare]  The module then computes the RG running upwards given the MZ-scale Yukawa matrices and relevant quantities.\n"<>
"\[FilledSmallSquare]  The relevant O(1) coefficients are checked to see if they are still O(1) - otherwise they will incur a logarithmic penalty.\n"<>
"\[FilledSmallSquare]  The module returns a number which is the fitness relevant to the RG running.";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - Mass Computation*)


ComputeSMQuantities::usage="ComputeSMQuantities[Yd,Yu] computes the SM quantities for the quark sector given the down- and up-Yukawa matrices.\n"<>
"\[FilledSmallSquare]  Inputs: Yd and Yu are 3x3 numerical Yukawa matrices.\n"<>
"\[FilledSmallSquare]  Ideally it is the output of the AdjoinO1coeffWOpsTextQSec module, with keys \"YdMatVal\" and \"YuMatVal\".\n"<>
"\[FilledSmallSquare]  Output: An association with the following keys:\n"<>
"  \[FilledVerySmallSquare]  \"Hd\" - VEV of down-Higgs (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"Hu\" - VEV of up-Higgs (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"Higgs\" - VEV of Higgs v (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"tan\[Beta]\" - the factor Hu/Hd.\n"<>
"  \[FilledVerySmallSquare]  \"CKMmat\" - the CKM matrix.\n"<>
"  \[FilledVerySmallSquare]  \"md\" - the down, strange and bottom quark masses in GeV.\n"<>
"  \[FilledVerySmallSquare]  \"mu\"- the up, charm and top quark masses in GeV.";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - Fitness Computation*)


NormalLog::usage="NormalLog[num_,denom_] computes the absolute value of log10 of the faction (num/denom).\n"<>
"\[FilledSmallSquare]  The result is bounded at 10, if Abs(Log10(n/d)) >= 10 then the output is 10.\n";

FitnessHiggs::usage="FitnessHiggs[massAssoc] computes the fitness in the Higgs sector.\n"<>
"\[FilledSmallSquare]  Input: massAssoc an association with the keys:\n"<>
"  \[FilledVerySmallSquare]  \"Higgs\" - VEV of Higgs in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"Hu\" - VEV of up-Higgs in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"Hd\" - VEV of down-Higgs in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"tan\[Beta]\" - tan\[Beta] factor. \n"<>
"\[FilledSmallSquare]  Output: an association with the keys:\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHSec\" - total fitness from Higgs sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHiggs\" - fitness contributions from Higgs VEVs.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessTan\[Beta]\" - fitness contributions from tan\[Beta].\n"<>
"\[FilledSmallSquare]  Since the \!\(\*SubscriptBox[\(o\), \(t\)]\) should fix the Higgs VEV this is currently turned off, i.e. all output values are set to 1,.";

FitnessQuarksMassMix::usage="FitnessQuarksMassMix[massAssoc] computes the fitness in the Quark Mass and Mixing sector.\n"<>
"\[FilledSmallSquare]  Input: massAssoc an association with the keys:\n"<>
"  \[FilledVerySmallSquare]  \"mu\" - masses of up, charm and top quarks in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"md\" - masses of down, strange and bottom quarks in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"CKM\" - CKM matrix. \n"<>
"\[FilledSmallSquare]  Output: an association with the keys:\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMassSec\" - total fitness from Quark Mass sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMixSec\" - total fitness from Quark Mix sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQmass\" - individual fitness contributions from quark masses, down families first.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessCKM\" - individual fitness contributions from CKM entries.";

FitnessGSAnomaly::usage="FitnessGSAnomaly[Qmat] computes the fitness in the Green-Schwarz sector.\n"<>
"\[FilledSmallSquare]  It detects whether the U(1) anomaly can be cancelled via a Green-Schwarz mechanism.\n"<>
"\[FilledSmallSquare]  Input: Qmat ->This must be a (17+\"nSinglet\"+\"nNonperts\")*(\"nU1s\") integer-valued matrix, specifically the charge matrix of the state.\n"<>
"\[FilledSmallSquare]  This is the output of the module StateVectoMatwithVEVs[stateVec].\n"<>
"\[FilledSmallSquare]  If an integer solution can be found then the contribution is 0.\n"<>
"\[FilledSmallSquare]  If an integer solution cannot be found but a real solution exists then the contribution is L1-rounding error to the nearest integer solution.\n"<>
"\[FilledSmallSquare]  If no solution can be found the contribution is the error in the least square problem.\n"<>
"\[FilledSmallSquare]  The output is an association with one key \"FitnessGSSec\" which gives the fitness contribution in this sector.";

CheckCharge::usage="CheckCharge[vector] takes in a 2 component vector and checks whether the entries are equal.\n"<>
"\[FilledSmallSquare]  It outputs a zero if the entries are non-equal and gives the value of the entry if the entries are equal.";

FitnessNonAbelianFactors::usage="FitnessNonAbelianFactors[stateVec] computes the fitness in the non-abelian sector.\n"<>
"\[FilledSmallSquare]  The input stateVec is the state-vector, a vector of dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  The module checks whether the charges of the 5, 5b and 1 reps charge directions lie in the same direction as the non-abelian deformations of the split bundle are specified by the vector n.\n"<>
"\[FilledSmallSquare]  If the condition is not satisfied an extra 1 is added to the contribution for each instance of violation.\n"<>
"\[FilledSmallSquare]  The output is a real number to be added in the total fitness.";

FitnessCalcQsectorOnly::usage="FitnessCalcQsectorOnly[massAssoc,stateVec,YukMZValsO1] computes all of the fitness functions relevant to the quark sector.\n"<>
"\[FilledSmallSquare] Inputs:\n"<>
"  \[FilledVerySmallSquare]  massAssoc - this is the output of ComputeSMQuantities and contains the computed SM quantities.\n"<>
"  \[FilledVerySmallSquare]  \"Higgs\" - VEV of Higgs in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"Hu\" - VEV of up-Higgs in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"Hd\" - VEV of down-Higgs in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"tan\[Beta]\" - tan\[Beta] factor. \n"<>
"  \[FilledVerySmallSquare]  \"mu\" - masses of up, charm and top quarks in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"md\" - masses of down, strange and bottom quarks in GeV. \n"<>
"  \[FilledVerySmallSquare]  \"CKM\" - CKM matrix. \n"<>
"  \[FilledVerySmallSquare]  stateVec - this is a dim = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"  \[FilledVerySmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"  \[FilledVerySmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"  \[FilledVerySmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"  \[FilledVerySmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"  \[FilledVerySmallSquare]  YukMZValsO1 is an association containing the following keys:\n"<>
"  \[FilledVerySmallSquare]  \"ytMZVal\" - leading order diagonal Yukawa term in up.\n"<>
"  \[FilledVerySmallSquare]  \"ybMZVal\" - leading order diagonal Yukawa term in down.\n"<>
"  \[FilledVerySmallSquare]  \"PosUpLead\" - position along diagonal identifying the leading order eigenvalue in up.\n"<>
"  \[FilledVerySmallSquare]  \"PosDownLead\" - position along diagonal identifying the leading order eigenvalue in down.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadUpMat\" - 3x3 matrix giving the O(1) coefficients for the leading order terms in the up-yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"O1LeadDownMat\" - 3x3 matrix giving the O(1) coefficients for the leading order terms in the down-yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  This is nicely the output of the submodule FromAssocFormGetWOps[assocform].\n"<>
"\[FilledSmallSquare]  The module computes all the fitness functions relevant - Higgs,QMass,QMix,NonAb,GS and RG.\n"<>
"\[FilledSmallSquare]  The O(1) and texture contributions to the fitness is computed one-level up and not in this module.\n"<>
"\[FilledSmallSquare]  Output: An association with the following keys:\n"<>
"  \[FilledVerySmallSquare]  \"Fitness\" - the total fitness, a non-positive number.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHSec\" - total fitness from Higgs sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHiggs\" - fitness contributions from Higgs VEVs.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessTan\[Beta]\" - fitness contributions from tan\[Beta].\n"<>
"  \[FilledVerySmallSquare]  \"FitnessRGSec\" - total fitness from the non-abelian sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMassSec\" - total fitness from Quark Mass sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMixSec\" - total fitness from Quark Mix sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQmass\" - individual fitness contributions from quark masses, down families first.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessCKM\" - individual fitness contributions from CKM entries.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessGSSec\" - total fitness from the GS sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessNAFSec\" - total fitness from the non-abelian sector.";


(* ::Subsubsection::Closed:: *)
(*Auxiliary Module help - Miscellaneous*)


DetectBitlist::usage="DetectBitlist[object] detects whether the entered object is an appropriately-lengthed bitlist.\n"<>
"\[FilledSmallSquare]  The length of the bitlist is determined by the module GADetermDimStateVecwithVEVs.";

ChangeO1VecByAmount::usage="ChangeO1VecByAmount[O1Vec,n] changes the O1Vec by a fixed amount set by the Options at position n.\n"<>
"\[FilledSmallSquare]  Inputs:\n"<>
"  \[FilledVerySmallSquare]  O1Vec - this is an O1Vec, a typical example is the varible O1CLst.\n"<>
"  \[FilledVerySmallSquare]  n - a number smaller or equal to dim(O1Vec).\n"<>
"\[FilledSmallSquare]  Module adds \"O1CoeffVarDiff\" to the n-th position in O1Vec and returns the updated O1Vec.";

FitnessO1Coeff::usage="FitnessO1Coeff[varMassAssoc,massAssoc,O1val] computes the fitness contributions of O1 coefficients.\n"<>
"\[FilledSmallSquare]  Inputs: varMassAssoc, massAssoc are both mass associations and have the following keys:\n"<>
"\[FilledSmallSquare]  varMassAssoc is assumed to have one O(1) entry changed by an amount specified by \"O1CoeffVarDiff\" in the Options.\n"<>
"\[FilledSmallSquare]  O1val is the value of the O1 coefficient being changed.\n"<>
"\[FilledSmallSquare]  The module computes the numerical magnitude of the vector (\!\(\*SubscriptBox[\(dO\), \(i\)]\)/do) and returns this number as the key \"FitnessO1Sec\".\n"<>
"\[FilledSmallSquare]  The vector is returned as the key \"ListofdOdp\".";


(* ::Subsubsection:: *)
(*Auxiliary Module help - Data Analysis*)


LowestOrderTerms::usage="LowestOrderTerms[expression,vars] splits out the lowest order terms given variables and expressions.";

VecDuplSymm::usage="VecDuplSymm[vec] determines the duplicates and their positions in a given vector.\n"<>
"\[FilledSmallSquare]  Input: An integer vector.\n"<>
"\[FilledSmallSquare]  Returns association with \"RepEntry\" determining the entries that are repeated and \"RepPos\" their correponding positions in the input vector.";

RedunQSymm::usage="RedunQSymm[vec] determines the redundant symmetry in the state in statevector form given the settings in Options.\n"<>
"\[FilledSmallSquare]  Module returns association with key \"DirSymm\" indicating the directions of the remaining \!\(\*SubscriptBox[\(S\), \(2\)]\) symmetry.\n"<>
"\[FilledSmallSquare]  If {} is returned then it means there is no \!\(\*SubscriptBox[\(S\), \(2\)]\) symmetry.";

CalcQScoreVec::usage="CalcQScoreFix[fiveBarQ,oneQ,kahlQ] determines whether we need to exchange the charges by canonically defining a direction in the redundant \!\(\*SubscriptBox[\(S\), \(2\)]\) symmetry, if it exists.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  fiveBarQ - dim-6 vector specifying the charges in the five bar reps.\n"<>
"  \[FilledVerySmallSquare]  oneQ - dim-(2*\"nSinglet\") vector specifying the charges in the one reps.\n"<>
"  \[FilledVerySmallSquare]  kahlQ - dim-(\"nNonperts\"*\"nU1s\") vector specifying the line bundle integers.\n"<>
"\[FilledSmallSquare]  The module first uses RedunQSymm to find the directions where there is an \!\(\*SubscriptBox[\(S\), \(\(2\)\(\\\\\)\)]\)symmetry, and then determine a score for each of the directions if it is possible.\n"<>
"\[FilledSmallSquare]  It outputs an association with the following keys:\n"<>
"  \[FilledVerySmallSquare]  \"DirSymm\" - Direction where there is a symmetry. \n"<>
"  \[FilledVerySmallSquare]  \"DirSymmScores\" - Scores in each of the direction. This is the total charge of the Kahlers in that direction, plus the number of occurences it appears in the positive charges minus the number of negative charges in the 5 bar and one reps. \n"<>
"  \[FilledVerySmallSquare]  \"NeedReplace\" - A boolean determining whether a replacement rule is needed. \n"<>
"  \[FilledVerySmallSquare]  \"ReplRule\" - The replacement rule if \"NeedReplace\" is True to flip the directions to obtain the canonical ordering. \n"<>
"  \[FilledVerySmallSquare]  \"ResultDupli\" - If the scores are equal this is set to True, False otherwise. Saying this module does not work. ";

StateNormalForm::usage="StateNormalForm[stateVec] presents the normal form for the inputted statevector.\n"<>
"\[FilledSmallSquare]  stateVec must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"\[FilledSmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  The module finds whether there are extra symmetries using the module CalcQScoreVec and reorders the statevector into a Normal form to reduce all redundant symmetries.\n"<>
"\[FilledSmallSquare]  The output is the normal ordered form of the statevector.";

EvalStateVecSoloList::usage="EvalStateVecSoloList[inputList] evaluates the inputted list of bitlists and adjoins the \"StateVec\" association to the states.\n"<>
"\[FilledSmallSquare]  Input: The list must be a list of associations each having the key \"Bitlist\" which is the GA bitlist consistent with the Options. The key \"Terminal\" must also exists and determines whether the state is terminal.\n"<>
"\[FilledSmallSquare]  Output: The same list with the key \"StateVec\" adjoined to each of the correspondidng vectors. The states that are not terminal are eliminated.";


AppendSetNoToStates::usage="AppendSetNoToStates[inputList] tags the list of input states with the appropriate set number.\n"<>
"\[FilledSmallSquare]  Input: the inputList must be of the form {listOfAssoc,...} where each entry is a list of associations of states, {<||>,...}.\n"<>
"\[FilledSmallSquare]  The module then joins the key \"SetNo\" into each of the entries, where \"SetNo\" is the number of the entry in the list. If the entry is empty then it is skipped.\n"<>
"\[FilledSmallSquare]  The module outputs the amended list with the states attached with the numbers.";

WriteUniqueStates::usage="WriteUniqueStates[filename] reads the file and converts the content into a list of states that has been processed.\n"<>
"\[FilledSmallSquare]  Input: An accessible filename. Otherwise will return error.\n"<>
"\[FilledSmallSquare]  The states are read and evaluated from the file and adjoined with their positions of the run (so an accumulation graph can be generated.\n"<>
"\[FilledSmallSquare]  The list of states is then returned.";


SelectVecsSameCharge::usage="SelectVecsSameCharge[listOfStates] delete the duplicates of the list of states after being processed by the StateNormalForm function.\n"<>
"\[FilledSmallSquare]  Input: List of states in association form with the key \"StateVec\".\n"<>
"\[FilledSmallSquare]  Only unduplicated states will remain in the final list.\n"<>
"\[FilledSmallSquare]  This is used to produce the accumulation graph.";

SelectVecsSameModel::usage="SelectVecsSameModel[listOfStates] delete the duplicates of the list of states after being processed by the StateNormalForm function.\n"<>
"\[FilledSmallSquare]  Input: List of states in association form with the key \"StateVec\".\n"<>
"\[FilledSmallSquare]  Only unduplicated states will remain in the final list.\n"<>
"\[FilledSmallSquare]  This is used to produce the accumulation graph.";

ProduceReducedList::usage="ProduceReducedList[evalStates] runs through the inputted this and produces the list of states flattened and reduced using StateNormalForm.\n"<>
"\[FilledSmallSquare]  The input is a list of states that are being processed each containing the key \"StateVec\".\n"<>
"\[FilledSmallSquare]  Outputs the list with all redundant symmetries reduced, i.e. all states are now physical.";

ProduceReducedModelList::usage="ProduceReducedModelList[evalStates] runs through the inputted this and produces the list of states flattened and reduced using StateNormalForm.\n"<>
"\[FilledSmallSquare]  The input is a list of states that are being processed each containing the key \"StateVec\".\n"<>
"\[FilledSmallSquare]  Outputs the list with all redundant symmetries reduced, i.e. all states are now physical.";

ComputeNumModuli::usage="ComputeNumModuli[state] computes the number of moduli used in the leading order Yukawa textures.\n"<>
"\[FilledSmallSquare]  The input is a state in association form and must contain the key \"StateVec\" which encodes the state in statevector form.\n"<>
"\[FilledSmallSquare]  The state is evaluated using the module ComputeWOpsGivenStateVec and prints the number of fields used in the leading-order Yukawa texture via adjoining the keys \"numSing\" and \"numKahl\", the number of singlets and non-perturbative fields respectively.";

ProduceFullRedList::usage="ProduceFullRedList[inputFile,outputFile] takes in the input file of a GA output and prints out the list of reduced states in the outputFile.\n"<>
"\[FilledSmallSquare]  The inputs of the module must be file names accessible by the Mathematica kernel, the input file must exist.\n"<>
"\[FilledSmallSquare]  The module produces a file and prints the output list of states into the output file location.\n"<>
"\[FilledSmallSquare]  The number of runs is returned - this is the number of complete blocks of runs done by the GA.";

ProduceFullRedModelList::usage="ProduceFullRedModelList[inputFile,outputFile] takes in the input file of a GA output and prints out the list of reduced states in the outputFile.\n"<>
"\[FilledSmallSquare]  The inputs of the module must be file names accessible by the Mathematica kernel, the input file must exist.\n"<>
"\[FilledSmallSquare]  The module produces a file and prints the output list of states into the output file location.\n"<>
"\[FilledSmallSquare]  The number of runs is returned - this is the number of complete blocks of runs done by the GA.";

PlotAccumGraph::usage="PlotAccumGraph[accumVec,title] produces the GA accumulation file given the accumulation vector and the title.\n"<>
"\[FilledSmallSquare]  Inputs:\n"<>
"  \[FilledVerySmallSquare]  accumVec - this is the accumulated version of number of states, the y-axis intervals in the accumulation plot.\n"<>
"  \[FilledVerySmallSquare]  title - a string to be used as the plot title.\n"<>
"\[FilledSmallSquare]  Output is the accumulation plot.";

FromRedListGetPlot::usage="FromRedListGetPlot[filePath,num] produces the GA accumulation file given file path for the states.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  filePath - the path name for the file which conatins the process list of states, following the module ProduceFullRedList (as its second input).\n"<>
"  \[FilledVerySmallSquare]  num - The number of run cycles completed. This is the printed number in ProduceFullRedList.\n"<>
"\[FilledSmallSquare]  Module then produces an accumulation graph based on the data and prints it.";

GivenModNumFromRedListGetPlot::usage="GivenModNumFromRedListGetPlot[filePath,num,numSing,numKahl] produces the GA accumulation file given file path for the states specified to the number of moduli specified.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  filePath - the path name for the file which conatins the process list of states, following the module ProduceFullRedList (as its second input).\n"<>
"  \[FilledVerySmallSquare]  num - The number of run cycles completed. This is the printed number in ProduceFullRedList.\n"<>
"  \[FilledVerySmallSquare]  numSing - The number of target singlets in the Yukawa textures.\n"<>
"  \[FilledVerySmallSquare]  numKahl - The number of target non-perturbative fields in the Yukawa textures.\n"<>
"\[FilledSmallSquare]  Module then produces an accumulation graph based on the data and prints it, only picking the states that has a leading-order Yukawa texture with the number of moduli as specified.";


FromStateVecListGiveNum::usage="FromStateVecListGiveNum[filePath] reduces the list of states and gives the number of reduced states as output.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  filePath - the path name for the file which conatins the list of states with \"StateVec\" as one of their keys.\n"<>
"\[FilledSmallSquare]  Module then computes and finds the reduced number of states and prints that as output.";

ChargeStateNormalForm::usage="ChargeStateNormalForm[QstateVec] presents the normal form for the inputted charge part of statevector (without the VEVs).\n"<>
"\[FilledSmallSquare]  stateVec must be inputed as a vector, with dimension = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\".\n"<>
"\[FilledSmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"\[FilledSmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"\[FilledSmallSquare]  The module finds whether there are extra symmetries using the module CalcQScoreVec and reorders the statevector into a Normal form to reduce all redundant symmetries.\n"<>
"\[FilledSmallSquare]  The output is the normal ordered form of the Q-part of the statevector.";


(* ::Subsection::Closed:: *)
(*Main Module Help*)


CompleteStateGivenWOps::usage="CompleteStateGivenWOps[stateVec,WOpsOut,massHierAssoc,VEVsPowers,O1Vec] computes the state given the WOps.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  stateVec - this is a dim = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"  \[FilledVerySmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"  \[FilledVerySmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"  \[FilledVerySmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"  \[FilledVerySmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"  \[FilledVerySmallSquare]  WOpsOut - a dimension {18,1}-matrix with each entry either being a 0 (no allowed operators) or a list of allowed operators.\n"<>
"  \[FilledVerySmallSquare]  It is the output of the module ComputeWOpsGivenStateVec[stateVec].\n"<>
"  \[FilledVerySmallSquare]  massHierarchyAssoc - This is an association that contains the leading-order data of the Yukawa matrices.\n"<>
"  \[FilledVerySmallSquare]  This is naturally the output of the module CalcScaleAndHiggsO1Coeff[mHAssoc].\n"<>
"  \[FilledVerySmallSquare]  The required list of the keys for this association is as follows:\n"<>
"  \[FilledVerySmallSquare]  \"BestScale\" - This is best scale (geometric mean) by fixing the scales to the actual quark hierarchies.\n"<>
"  \[FilledVerySmallSquare]  \"BestO1Top\" - This is the best O(1) coefficient associated to the top mass (diagonal entry in Yukawa matrix that gives the top-mass).\n"<>
"  \[FilledVerySmallSquare]  \"PosOTop\" -> Diagonal entry position that contributes the top eigenvalue in the up-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare] VEVsPowers -> This is an association with the keys \"\[Phi]VEVPowers\" and \"\[CapitalPhi]VEVPowers\" that encodes the powers of the VEVs of the matrix.\n"<>
"  \[FilledVerySmallSquare] This is the output of the module ConvertStateVecToVEVPowers[statevec].\n"<>
"  \[FilledVerySmallSquare]  O1Vec - a list of O(1) coefficients which must have a larger dimension then the number of operator insertions.\n"<>
"\[FilledSmallSquare]  The module then computes the state and outputs the fullstate association with the following keys:\n"<>
"  \[FilledVerySmallSquare]  Fitness:\n"<>
"  \[FilledVerySmallSquare]  \"Fitness\" - the total fitness, a non-positive number.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHSec\" - total fitness from Higgs sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHiggs\" - fitness contributions from Higgs VEVs.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessTan\[Beta]\" - fitness contributions from tan\[Beta].\n"<>
"  \[FilledVerySmallSquare]  \"FitnessRGSec\" - total fitness from the non-abelian sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMassSec\" - total fitness from Quark Mass sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMixSec\" - total fitness from Quark Mix sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQmass\" - individual fitness contributions from quark masses, down families first.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessCKM\" - individual fitness contributions from CKM entries.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessGSSec\" - total fitness from the GS sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessNAFSec\" - total fitness from the non-abelian sector.\n"<>
"  \[FilledVerySmallSquare]  Mass Assoc:\n"<>
"  \[FilledVerySmallSquare]  \"Hd\" - VEV of down-Higgs (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"Hu\" - VEV of up-Higgs (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"Higgs\" - VEV of Higgs v (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"tan\[Beta]\" - the factor Hu/Hd.\n"<>
"  \[FilledVerySmallSquare]  \"CKMmat\" - the CKM matrix.\n"<>
"  \[FilledVerySmallSquare]  \"md\" - the down, strange and bottom quark masses in GeV.\n"<>
"  \[FilledVerySmallSquare]  \"mu\"- the up, charm and top quark masses in GeV.\n"<>
"  \[FilledVerySmallSquare]  Extras:\n"<>
"  \[FilledVerySmallSquare]  \"ValueLst\" - the values of H,QMass,QMix fitness contributions.\n"<>
"  \[FilledVerySmallSquare]  \"\[Phi]VEVs\" - the VEVs of singlet moduli.\n"<>
"  \[FilledVerySmallSquare]  \"\[CapitalPhi]VEVs\" - the VEVs of effective Kahler moduli.\n"<>
"  \[FilledVerySmallSquare]  \"YukdMatSt\" - leading order term of down-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"YukuMatSt\" - leading order term of up-Yukawa matrix.\n"<>
"\[FilledSmallSquare]  Module is used in one-level up in CompleteStateO1Coeffs.";

CompleteStateO1Coeffs::usage="CompleteStateO1Coeffs[stateVec,O1Vec] computes the full state given the O1 vector O1Vec.\n"<>
"\[FilledSmallSquare]  Inputs: \n"<>
"  \[FilledVerySmallSquare]  stateVec - this is a dim = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"  \[FilledVerySmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"  \[FilledVerySmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"  \[FilledVerySmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"  \[FilledVerySmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"  \[FilledVerySmallSquare]  O1Vec - a list of O(1) coefficients which must have a larger dimension then the number of operator insertions.\n"<>
"\[FilledSmallSquare]  The module then computes the state and outputs the fullstate association with the following keys (with O1 variation):\n"<>
"\[FilledSmallSquare]  The module then computes the state and outputs the fullstate association with the following keys:\n"<>
"  \[FilledVerySmallSquare]  \"Bits\" - the bitlist, if exists.\n"<>
"  \[FilledVerySmallSquare]  \"StateVec\" - the state-vector inputted. \n"<>
"  \[FilledVerySmallSquare]  \"Terminal\" - boolean whether the state is terminal according to bound set by \"FitRange\" of Options. \n"<>
"  \[FilledVerySmallSquare]  Fitness:\n"<>
"  \[FilledVerySmallSquare]  \"Fitness\" - the total fitness, a non-positive number.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHSec\" - total fitness from Higgs sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessHiggs\" - fitness contributions from Higgs VEVs.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessTan\[Beta]\" - fitness contributions from tan\[Beta].\n"<>
"  \[FilledVerySmallSquare]  \"FitnessRGSec\" - total fitness from the non-abelian sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMassSec\" - total fitness from Quark Mass sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQMixSec\" - total fitness from Quark Mix sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessQmass\" - individual fitness contributions from quark masses, down families first.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessCKM\" - individual fitness contributions from CKM entries.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessGSSec\" - total fitness from the GS sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessNAFSec\" - total fitness from the non-abelian sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessO1Sec\" - total fitness from the O(1) sector.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessO1List\" - individual contributions of 16 O(1)-number variations.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessTextSec\" -> This is the total fitness contribution from texture analysis.\n"<>
"  \[FilledVerySmallSquare]  \"FitnessPowers\" -> Fitness from negative powers in scales.\n"<>
"  \[FilledVerySmallSquare]  \"VarsScale\" -> Fitness from variance of calculated scales.\n"<>
"  \[FilledVerySmallSquare]  \"ScalePenalty\" -> Penalty from super-Planckian scales.\n"<>
"  \[FilledVerySmallSquare]  \"OTopPenalty\" -> Penalty from O(1).\n"<>
"  \[FilledVerySmallSquare]  Mass Assoc - this is the one without O(1) variation:\n"<>
"  \[FilledVerySmallSquare]  \"Hd\" - VEV of down-Higgs (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"Hu\" - VEV of up-Higgs (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"Higgs\" - VEV of Higgs v (in GeV).\n"<>
"  \[FilledVerySmallSquare]  \"tan\[Beta]\" - the factor Hu/Hd.\n"<>
"  \[FilledVerySmallSquare]  \"CKMmat\" - the CKM matrix.\n"<>
"  \[FilledVerySmallSquare]  \"md\" - the down, strange and bottom quark masses in GeV.\n"<>
"  \[FilledVerySmallSquare]  \"mu\"- the up, charm and top quark masses in GeV.\n"<>
"  \[FilledVerySmallSquare]  Extras:\n"<>
"  \[FilledVerySmallSquare]  \"ValueLst\" - the values of H,QMass,QMix fitness contributions.\n"<>
"  \[FilledVerySmallSquare]  \"\[Phi]VEVs\" - the VEVs of singlet moduli.\n"<>
"  \[FilledVerySmallSquare]  \"\[CapitalPhi]VEVs\" - the VEVs of effective Kahler moduli.\n"<>
"  \[FilledVerySmallSquare]  \"YukdMatSt\" - leading order term of down-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"YukuMatSt\" - leading order term of up-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  Charges:\n"<>
"  \[FilledVerySmallSquare]  \"Qq\" - a (3*\"nU1s\") matrix. This is the charge of the Q in the upstairs representation. \n"<>
"  \[FilledVerySmallSquare]  \"Qd\" - a (3*\"nU1s\") matrix, charge of the d in the upstairs representation. \n"<>
"  \[FilledVerySmallSquare]  \"Qu\" - a (3*\"nU1s\") matrix, charge of the u in the upstairs representation. \n"<>
"  \[FilledVerySmallSquare]  \"QL\" - a (3*\"nU1s\") matrix, charge of the L in the upstairs representation. \n"<>
"  \[FilledVerySmallSquare]  \"Qe\" - a (3*\"nU1s\") matrix, charge of the e in the upstairs representation. \n"<>
"  \[FilledVerySmallSquare]  \"QHd\" - a (1*\"nU1s\") matrix, charge of the Hd in the upstairs representation. \n"<>
"  \[FilledVerySmallSquare]  \"Q\[Phi]\" - a (\"nSinglet\"*\"nU1s\") matrix, charge of the bundle moduli in the upstairs representation. \n"<>
"  \[FilledVerySmallSquare]  \"Q\[CapitalPhi]\" - a (\"nNonperts\"*\"nU1s\") matrix, the charge of the Kahler moduli in the upstairs representation.\n"<>
"  \[FilledVerySmallSquare]  Texture Details:\n"<>
"  \[FilledVerySmallSquare]  \"mtPow\" - The leading-order eigenvalue from the up Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"mbPow\" - The leading-order eigenvalue from the down Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"mu/mcPow\" - The leading-order scale-difference between the up and charm quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"mc/mtPow\" - The leading-order scale-difference between the charm and top quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"md/msPow\" - The leading-order scale-difference between the down and strange quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"ms/mbPow\" - The leading-order scale-difference between the strange and bottom quarks (eigenvalues).\n"<>
"  \[FilledVerySmallSquare]  \"YukUpMatLead\" - This is a {3,3}-matrix listing the leading-order tmer of the up Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"YukDownMatLead\" - This is a {3,3}-matrix listing the leading-order tmer of the down Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"PosOTop\" - Diagonal entry position that contributes the top eigenvalue in the up-Yukawa matrix.\n"<>
"  \[FilledVerySmallSquare]  \"PosOBot\" - Diagonal entry position that contributes the bottom eigenvalue in the down-Yukawa matrix.\n";
CompleteStateQsector::usage="CompleteStateQsector[state] computes the full state.\n"<>
"\[FilledSmallSquare]  The input must be of state-vector or bitlist form.\n"<>
"  \[FilledVerySmallSquare]  bitlist form - in binary form, with dimension as determined by the module GADetermDimStateVecwithVEVs or GADetermDimStateVecFix.\n"<>
"  \[FilledVerySmallSquare]  state-vector form - this is a dim = 11+2*\"nSinglet\"+\"nNonperts\"*\"nU1s\"+2*(\"nSinglet\"+\"nNonperts\").\n"<>
"  \[FilledVerySmallSquare]  The first 11+2*\"nSinglet\" coordinates must be ranged from 1 to \"nU1s\", where \"nU1s\" is the number of U(1)s in the options.\n"<>
"  \[FilledVerySmallSquare]  The last \"nNonperts\"*\"nU1s\" coordinates of the first part must be ranged inside the range specified in the options in \"\[CapitalPhi]ChargeRange\".\n"<>
"  \[FilledVerySmallSquare]  The \[Phi]-VEVs coordinates are integer-valued and lie between {0,(2^(\"bitsfor\[Phi]\")-1)}.\n"<>
"  \[FilledVerySmallSquare]  The \[CapitalPhi]-VEVs coordinates are integer-valued and lie between {n,(2^(\"bitsfor\[CapitalPhi]\")-1+n)}, n=2.\n"<>
"\[FilledSmallSquare]  Module computes the entire full state, with an association with keys as sketched out in CompleteStateO1Coeffs.";
GAFitnessFuncQsector::usage="GAFitnessFuncQsector[state] computes the full state.\n"<>
"\[FilledSmallSquare]  The input must be of bitlist form.\n"<>
"\[FilledSmallSquare]  Module computes the entire full state, with an association with keys as sketched out in CompleteStateO1Coeffs.\n"<>
"\[FilledSmallSquare]  This module can be used as the FitnessFunction for GA.";
GAInitialiseStateQsector::usage="GAInitialiseStateQsector[] randomises a state and computes its full state.\n"<>
"\[FilledSmallSquare]  Module computes the entire full state, with an association with keys as sketched out in CompleteStateO1Coeffs.\n"<>
"\[FilledSmallSquare]  This module can be used as the InitFunction for GA.";



(* ::Subsubsection::Closed:: *)
(*Junk*)


(*
CompleteStateQsectorAssoc::usage="CompleteStateQsectorAssoc[stateassoc] computes the full state given an input state in association form.\n"<>
"\[FilledSmallSquare]  stateassoc -> must have the keys (\"Qq\",\"Qu\",\"Qe\") or (\"Q10\"); (\"Qd\",\"QL\") or (\"Q5b\"); (\"Q\[Phi]\") and/or (\"Q1\"); (\"Q\[CapitalPhi]\"); (\"\[Phi]VEVs\" and \"\[CapitalPhi]VEVs\").\n"<>
"\[FilledSmallSquare]  The description of the keys are as follows: \n"<>
"\[FilledSmallSquare]  \"Qq\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the Q in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qd\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the d in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qu\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the u in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"QL\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the L in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qe\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the e in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"QHd\" -> this must be a (1*\"nU1s\") matrix. This is the charge of the Hd in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Q\[Phi]\" -> this must be a (\"nSinglet\"*\"nU1s\") matrix. This is the charge of the bundle moduli in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Q\[CapitalPhi]\" this must be a (\"nNonperts\"*\"nU1s\") matrix. This is the charge of the Kahler moduli in the upstairs representation.\n"<>
"\[FilledSmallSquare]  \"\[Phi]VEVs\" -> this must be a dimension-(\"nSinglet\") vector. This is the VEVs for the bundle moduli. \n"<>
"\[FilledSmallSquare]  \"\[CapitalPhi]VEVs\" -> this must be a dimension-(\"nNonperts\") vector. This is the VEVs for the Kahler moduli. \n"<>
"\[FilledSmallSquare]  Note that \"Qq\", \"Qu\" and \"Qe\" should be equal since they are from the 10 rep in SU(5). Alternatively the package reads \"Q10\" as a 3-vector which just gives the charge-positions of 10s. \n"<>
"\[FilledSmallSquare]  Note that \"QL\" and \"Qd\" should be equal since they are from the 5 bar rep in SU(5). Alternatively the package reads \"Q5b\" as a (3*2)-matrix which gives the charge positions {e_a,e_b} for each 5 bars. \n"<>
"\[FilledSmallSquare]  Note that the 1 rep has charge pattern e_a-e_b. In this case there is a charge pattern penalty when the bundle moduli has zero charge in all U1s.\n"<>
"\[FilledSmallSquare]  In this case we will have to specify the charges using \"Q1\". This is a (\"nSinglet\"*2)-matrix which gives the charge positions {e_a,e_b} for each 1s. \n"<>
"\[FilledSmallSquare]  If the dimensions of the entered association is not matched the module will print an error message and exit.\n"<>
"\[FilledSmallSquare]  The O(1) coefficients are changed by manipulating SetO1coeff. See SetO1coeff and GetO1coeff for details.\n"<>
"\[FilledSmallSquare]  The output is a state association specifying the massrules and fitness of the state.";

FromAssocFormGetWOps::usage="FromAssocFormGetWOps[assocform] computes the WOps expansion from the inputted state in association form.\n"<>
"\[FilledSmallSquare]  stateassoc -> must have the keys (\"Qq\",\"Qu\",\"Qe\") or (\"Q10\"); (\"Qd\",\"QL\") or (\"Q5b\"); (\"Q\[Phi]\") and/or (\"Q1\"); (\"Q\[CapitalPhi]\"); (\"\[Phi]VEVs\" and \"\[CapitalPhi]VEVs\").\n"<>
"\[FilledSmallSquare]  The description of the keys are as follows: \n"<>
"\[FilledSmallSquare]  \"Qq\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the Q in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qd\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the d in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qu\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the u in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"QL\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the L in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Qe\" -> this must be a (3*\"nU1s\") matrix. This is the charge of the e in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"QHd\" -> this must be a (1*\"nU1s\") matrix. This is the charge of the Hd in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Q\[Phi]\" -> this must be a (\"nSinglet\"*\"nU1s\") matrix. This is the charge of the bundle moduli in the upstairs representation. \n"<>
"\[FilledSmallSquare]  \"Q\[CapitalPhi]\" this must be a (\"nNonperts\"*\"nU1s\") matrix. This is the charge of the Kahler moduli in the upstairs representation.\n"<>
"\[FilledSmallSquare]  \"\[Phi]VEVs\" -> this must be a dimension-(\"nSinglet\") vector. This is the VEVs for the bundle moduli. \n"<>
"\[FilledSmallSquare]  \"\[CapitalPhi]VEVs\" -> this must be a dimension-(\"nNonperts\") vector. This is the VEVs for the Kahler moduli. \n"<>
"\[FilledSmallSquare]  Note that \"Qq\", \"Qu\" and \"Qe\" should be equal since they are from the 10 rep in SU(5). Alternatively the package reads \"Q10\" as a 3-vector which just gives the charge-positions of 10s. \n"<>
"\[FilledSmallSquare]  Note that \"QL\" and \"Qd\" should be equal since they are from the 5 bar rep in SU(5). Alternatively the package reads \"Q5b\" as a (3*2)-matrix which gives the charge positions {e_a,e_b} for each 5 bars. \n"<>
"\[FilledSmallSquare]  Note that the 1 rep has charge pattern e_a-e_b. In this case there is a charge pattern penalty when the bundle moduli has zero charge in all U1s.\n"<>
"\[FilledSmallSquare]  In this case we will have to specify the charges using \"Q1\". This is a (\"nSinglet\"*2)-matrix which gives the charge positions {e_a,e_b} for each 1s. \n"<>
"\[FilledSmallSquare]  If the dimensions of the entered association is not matched the module will print an error message and exit.\n"<>
"\[FilledSmallSquare]  The module then ultilises the charge pattern to compute the allowed operators in the low-energy regime.\n"<>
"\[FilledSmallSquare]  The output is the expansion of the superpotential operators as a list (in each entry) in a {18,1} matrix, which gives the respective operator expansion with the corresponding coefficients for each expansion term.\n"<>
"\[FilledSmallSquare]  The ordering for the superpotential operators is the 9 down-Yukawa entries followed by the 9 up-Yukawa entries.\n"<>
"\[FilledSmallSquare]  If there are no operators then instead a 0 is outputted as the entry to the {18,1} matrix.";

CompleteStateFromAssocFormWOps::usage="Tests";
GAConvertStateVectoBitlistFix::usage="Need to add this in!";*)


(* ::Section:: *)
(*Initialise*)


(* ::Subsection::Closed:: *)
(*Start-up messages*)


Begin["`Private`"]
Print[Style["stringFNmodels: RL/GA Environment for Froggatt-Nielsen Mass Model Building with Non-pertubative Effects",Underlined,FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];
Print[Style["Constantin, Andrei and Fraser-Taliente, Cristofero and Harvey, Thomas R. and Leung, Lucas T. Y. and Lukas, Andre.",FontColor -> Black,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];
Print[Style["Last updated 19 Oct, 2024 by Lucas Leung.",FontColor -> Black,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];
Print[Style["Execute \!\(\*
StyleBox[\"?\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"stringFNmodels\",\nFontWeight->\"Bold\"]\) for help. Current stringFNmodels options:",FontColor -> Black,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];


(* ::Subsection::Closed:: *)
(*Options*)


Options[stringFNmodels]:={
(* Global Settings *)
(* global settings *)
"nU1s"->4,"nChargeVec"->{1,1,1,2},"nSinglet"->2,"nNonperts"->2,"\[CapitalPhi]ChargeRange"->{-7,8},
"DimOp"->5,"DimOpof\[CapitalPhi]"->2,"bitsfor\[Phi]"->3,"bitsfor\[CapitalPhi]"->3,
(* symbols *)
"BundleModuliFields"->Global`\[Phi],"KahlerEffFields"->Global`\[CapitalPhi],"VEVsScale"->Global`\[Epsilon],
(* vev settings *)
"\[Phi]VEVmin"->0.01,"\[Phi]VEVlogGen"->1.01,"\[CapitalPhi]VEVmin"->0.001,"\[CapitalPhi]VEVlogGen"->1.01,
(* fixed charges *)
"fixedCharges"->True,"tenQ"->{1,2,5},"fiveBarHiggsQ"->{4,5},
(* O1 settings *)
"O1CoeffRange"->{0.5,3},"numOfO1Coeff"->10000,"O1CoeffVarDiff"->0.001,
(* Method generation *)
"MethodGenOps"->"List",

(* Higgs settings: "vdmax"->174,"vumax"->174,"HiggsBound"->30,"HiggsFitScaleFactor"->10,*)
"vdmax"->174,"vumax"->174,

(* Fitness Settings *)
"FitRange"->-50,"factorFitnessH"->10,"factorFitnessRG"->1,"factorFitnessQ"->2,"factorFitnessCKM"->2,
"factorFitnessGS"->20,"factorFitnessNAF"->20,"factorFitnessTEX"->10,"factorFitnessO1"->5,
(* detail settings *)
"MaxLogPower"->10,
(* extra fitness penalties *)
"indPenal"->False,"diffPenal"->False,"indFitMax"->1,
(* individual fitness factors *)
"fitDownFac"->{1,1,1},"fitUpFac"->{1,1,1},"fitCKMFac"->{1,1,1,1,1,1,1,1,1},"fitTextFac"->{1,1,1,1},
(* higgs extra factors *)
(* "ExtraPenaltyforMaxHiggs"->1,
"Fitnesstan\[Beta]On"->True,"tan\[Beta]LowerBound"->0.6,"tan\[Beta]UpperBound"->10,"HiggsFitnessFactor"->1, *)
(* "OutputLong"->True, *)




(* RG Data *)
(*"\[Lambda]u33"->0.88,"\[Lambda]u3i"->0.88,"\[Lambda]uij"->1.49,"\[Lambda]d33"->2.11,"\[Lambda]d3i"->2.11,"\[Lambda]dij"->2.51,*)
(* Mass Data *)
(*"MZ"->80.379,*)
(*"Higgs"->174,*)
"v"->174,"vErr"->0.17,
"mu"->{0.00216,1.27,172.69},(*"muerror"->{0.00049,0.02,0.30},*)
"md"->{0.00467,0.0934,4.18},(*"mderror"->{0.00048,0.0086,0.03},*)
(*"sin12CKM"->0.2265,"sin13CKM"->0.0036,"sin23CKM"->0.0405,"deltaCKM"->1.196,*)
(*"me"->{0.00051099895000,0.1056583755,1.77686},"meerror"->{0.00000000000015,0.0000000023,0.00012},
"m\[Nu]"->Flatten[{225*10^(-9),190*10^(-9),18.2*10^(-6)}],"m\[Nu]error"->Flatten[{135*10^(-9),100*10^{-9},10*10^{-6}}],*)
"CKMmat"->Partition[Flatten[{{0.97373,0.2243,3.82*10^(-3)},{0.221,0.975,40.8*10^(-3)},{8.6*10^(-3),41.5*10^(-3),1.014}}],3]
(*"CKMmaterror"->Partition[Flatten[{{0.00031,0.0008,0.2*10^(-3)},{0.004,0.006,1.4*10^{-3}},{0.2*10^{-3},0.9*10^{-3},0.029}}],3],*)
(*"PMNSmatup3\[Sigma]bound"->{{0.845,0.579,0.155},{0.500,0.689,0.776},{0.525,0.694,0.756}},
"PMNSmatlow3\[Sigma]bound"->{{0.801,0.513,0.143},{0.234,0.471,0.637},{0.271,0.477,0.613}},
"tan\[Beta]"-> 0.03,"Jarlskog"->3.18*10^(-5)*)

};
(* These would be the bounds in which arbitrary choices are made *)
Print[Style[Options[stringFNmodels],"Output"]];


(* ::Section::Closed:: *)
(*Initialise O1 coefficients*)


(* ::Subsection::Closed:: *)
(*O1coefficients*)


stringFNmodels`O1CLst={};


(* ::Subsection::Closed:: *)
(*SetO1Coeff*)


SetO1Coeff[O1coeffvec_]:=Module[{O1list},
O1list=O1coeffvec;
O1CLst=O1list;
O1list
]


(* ::Subsection::Closed:: *)
(*GetO1Coeff*)


GetO1Coeff[]:=Module[{currentO1list},
currentO1list=O1CLst;
currentO1list
]


(* ::Subsection::Closed:: *)
(*NewO1Coeff*)


NewO1Coeff[]:=Module[{O1coeffvec},

(* This module generates a vector of O(1) coefficients. *)
O1coeffvec=RandomChoice[{-1,1},("numOfO1Coeff"/.Options[stringFNmodels])]*RandomReal[("O1CoeffRange"/.Options[stringFNmodels]),("numOfO1Coeff"/.Options[stringFNmodels])];
O1coeffvec

]
SetO1coeff[NewO1coeff[]];


(* ::Section:: *)
(*Auxiliary Modules - State Generation*)


(* ::Subsection::Closed:: *)
(*RandomStatewithVEVs*)


RandomStatewithVEVs[]:=Module[{dimchargelst,chrglstMSSM,singletchrglst,totalchrgvecT,num\[Phi]bit,\[Phi]VEVstate,randstate,
\[Phi]VEVsbits,\[Phi]VEVssigns,\[Phi]VEVsVec,\[CapitalPhi]VEVsbits,\[CapitalPhi]VEVssigns,\[CapitalPhi]VEVsVec,chrglstMSSMtemp},

(* First we find the dimension of the charge vector. *)
(* The number of charges is outlined in the notes*)
If[("fixedCharges"/.Options[stringFNmodels]),
dimchargelst=3*2+2*"nSinglet"/.Options[stringFNmodels];
,
dimchargelst=3*3+2+2*"nSinglet"/.Options[stringFNmodels];
];

(* Code positition of charges in the U(1)s of the MSSM particles and singlet fields *)
If[("fixedCharges"/.Options[stringFNmodels]),
chrglstMSSMtemp=RandomInteger[{1,"nU1s"}/.Options[stringFNmodels],dimchargelst];
chrglstMSSM=Join[("tenQ"/.Options[stringFNmodels]),Take[chrglstMSSMtemp,6],
("fiveBarHiggsQ"/.Options[stringFNmodels]),Take[chrglstMSSMtemp,-(2*("nSinglet"/.Options[stringFNmodels]))]]
,
chrglstMSSM=RandomInteger[{1,"nU1s"}/.Options[stringFNmodels],dimchargelst];
];

(* Code the charges taken by the non-perturbative fields separately. *)
singletchrglst=RandomInteger["\[CapitalPhi]ChargeRange"/.Options[stringFNmodels],{"nNonperts","nU1s"}/.Options[stringFNmodels]];
(* Now append the two charge vectors. *)
totalchrgvecT=Flatten[Append[chrglstMSSM,Flatten[singletchrglst]]];

(* Now add in the third part where the \[Phi]VEVs are determined. *)
num\[Phi]bit=("nSinglet"*"bitsfor\[Phi]"+"nNonperts"*"bitsfor\[CapitalPhi]")/.Options[stringFNmodels];
\[Phi]VEVstate=Table[RandomInteger[{0,1}],num\[Phi]bit];
(* Now determine the new chargevector for \[Phi]VEVs.*)
\[Phi]VEVsbits=Partition[Take[\[Phi]VEVstate,("nSinglet"*"bitsfor\[Phi]")/.Options[stringFNmodels]],"bitsfor\[Phi]"/.Options[stringFNmodels]];
\[Phi]VEVsVec=FromDigits[Transpose[\[Phi]VEVsbits],2];
If[\[Phi]VEVsVec==0,\[Phi]VEVsVec={}];
(* Repeat the same for \[CapitalPhi] part. *)
\[CapitalPhi]VEVsbits=Partition[Take[\[Phi]VEVstate,-("nNonperts"*"bitsfor\[CapitalPhi]")/.Options[stringFNmodels]],"bitsfor\[CapitalPhi]"/.Options[stringFNmodels]];
\[CapitalPhi]VEVsVec=FromDigits[Transpose[\[CapitalPhi]VEVsbits],2];
If[\[CapitalPhi]VEVsVec==0,\[CapitalPhi]VEVsVec={}];
\[Phi]VEVstate=Join[\[Phi]VEVsVec,\[CapitalPhi]VEVsVec];

randstate=Flatten[Append[totalchrgvecT,\[Phi]VEVstate]];
randstate

]



(* ::Subsection::Closed:: *)
(*ConvertStatewithVEVs*)


ConvertStatewithVEVs[state_]:=Module[{dimchargelst,chrglstMSSM,\[CapitalPhi]chrglst,Qcharges,dcharges,ucharges,Lcharges,echarges,Hdcharges,Hucharges,\[Phi]charges,convstateQ,
dim\[CapitalPhi]chrgbits,dim\[Phi]bits,\[Phi]VEVlst,\[Phi]VEVsbits,\[CapitalPhi]VEVsbits,\[Phi]VEVsVec,\[CapitalPhi]VEVsVec,convstateVEVs,convstate,\[Phi]VEVssigns,\[CapitalPhi]VEVssigns,\[Phi]VEVsVals,\[CapitalPhi]VEVsVals},

(* Chargevector to full association form. *)
dimchargelst=3*3+2+2*"nSinglet"/.Options[stringFNmodels];
dim\[Phi]bits=("nSinglet"+"nNonperts")/.Options[stringFNmodels];
dim\[CapitalPhi]chrgbits=Length[state]-dimchargelst-dim\[Phi]bits;

(* Now split the charge vector up. *)
chrglstMSSM=Take[state,dimchargelst];
chrglstMSSM=Table[UnitVector["nU1s",chrglstMSSM[[i]]],{i,dimchargelst}]/.Options[stringFNmodels];
\[CapitalPhi]chrglst=Take[state,{dimchargelst+1,dimchargelst+dim\[CapitalPhi]chrgbits}];
\[CapitalPhi]chrglst=Partition[\[CapitalPhi]chrglst,"nU1s"/.Options[stringFNmodels]];


(* Find all the charges with respect to the splitting suggested in the notes. *)
Qcharges=Take[chrglstMSSM,3];
dcharges=Take[chrglstMSSM,{4,6}]+Take[chrglstMSSM,{7,9}];
ucharges=Take[chrglstMSSM,3];
Lcharges=Take[chrglstMSSM,{4,6}]+Take[chrglstMSSM,{7,9}];
echarges=Take[chrglstMSSM,3];
Hdcharges=Take[chrglstMSSM,{10}]+Take[chrglstMSSM,{11}];
Hucharges=-Hdcharges;
\[Phi]charges=Take[chrglstMSSM,{12,11+"nSinglet"/.Options[stringFNmodels]}]-Take[chrglstMSSM,-"nSinglet"/.Options[stringFNmodels]];

(* Construct the association. *)
convstateQ=<|"Qq"-> Qcharges,"Qd"->dcharges,"Qu"->ucharges,"QL"->Lcharges,"Qe"->echarges,"QHd"->Hdcharges,"QHu"->Hucharges,"Q\[Phi]"->\[Phi]charges,"Q\[CapitalPhi]"->\[CapitalPhi]chrglst|>;

(* Now we will look at the moduli VEVs part. *)
(*\[Phi]VEVlst=Take[state,-dim\[Phi]bits];
\[Phi]VEVsbits=Take[\[Phi]VEVlst,("nSinglet")/.Options[stringFNmodels]];
\[Phi]VEVsVals=Abs[\[Phi]VEVsbits];
\[Phi]VEVsVec=Map[("\[Phi]VEVmin"/.Options[stringFNmodels])*("\[Phi]VEVlogGen"/.Options[stringFNmodels])^#&,\[Phi]VEVsVals];

\[CapitalPhi]VEVsbits=Take[\[Phi]VEVlst,-("nNonperts")/.Options[stringFNmodels]];
\[CapitalPhi]VEVsVals=Abs[\[CapitalPhi]VEVsbits];
\[CapitalPhi]VEVsVec=Map[("\[CapitalPhi]VEVmin"/.Options[stringFNmodels])*("\[CapitalPhi]VEVlogGen"/.Options[stringFNmodels])^#&,\[CapitalPhi]VEVsVals];

(*
\[Phi]VEVsVec=Map[#[[1]]+I #[[2]]&,Partition[\[Phi]VEVsVec,2]];
\[CapitalPhi]VEVsVec=Map[#[[1]]+I #[[2]]&,Partition[\[CapitalPhi]VEVsVec,2]];

*)
*)
(* Construct the association. *)
convstateQ
]


(* ::Subsection::Closed:: *)
(*StateVectoMatwithVEVs*)


StateVectoMatwithVEVs[state_]:=Module[{dimchargelst,chrglstMSSM,\[CapitalPhi]chrglst,Qcharges,dcharges,ucharges,Lcharges,echarges,Hdcharges,Hucharges,\[Phi]charges,chargemat,
dim\[Phi]bits,dim\[CapitalPhi]chrgbits},

(* Chargevector to full association form. *)
dimchargelst=3*3+2+2*"nSinglet"/.Options[stringFNmodels];
dim\[Phi]bits=("nSinglet"+"nNonperts")/.Options[stringFNmodels];
dim\[CapitalPhi]chrgbits=Length[state]-dimchargelst-dim\[Phi]bits;

(* Now split the charge vector up. *)
chrglstMSSM=Take[state,dimchargelst];
chrglstMSSM=Table[UnitVector["nU1s",chrglstMSSM[[i]]],{i,dimchargelst}]/.Options[stringFNmodels];
\[CapitalPhi]chrglst=Take[state,{dimchargelst+1,dimchargelst+dim\[CapitalPhi]chrgbits}];
\[CapitalPhi]chrglst=Partition[\[CapitalPhi]chrglst,"nU1s"/.Options[stringFNmodels]];


(* Find all the charges with respect to the splitting suggested in the notes. *)
Qcharges=Take[chrglstMSSM,3];
dcharges=Take[chrglstMSSM,{4,6}]+Take[chrglstMSSM,{7,9}];
ucharges=Take[chrglstMSSM,3];
Lcharges=Take[chrglstMSSM,{4,6}]+Take[chrglstMSSM,{7,9}];
echarges=Take[chrglstMSSM,3];
Hdcharges=Take[chrglstMSSM,{10}]+Take[chrglstMSSM,{11}];
Hucharges=-Hdcharges;
\[Phi]charges=Take[chrglstMSSM,{12,11+"nSinglet"/.Options[stringFNmodels]}]-Take[chrglstMSSM,-"nSinglet"/.Options[stringFNmodels]];

(* Construct the charge matrix. *)
chargemat=Partition[Flatten[Join[Qcharges,dcharges,ucharges,Lcharges,echarges,Hdcharges,Hucharges,\[Phi]charges,\[CapitalPhi]chrglst]],"nU1s"/.Options[stringFNmodels]]

]


(* ::Subsection::Closed:: *)
(*StateVectoVEVs*)


StateVectoVEVs[state_]:=Module[{dim\[Phi]bits,\[Phi]VEVlst,\[Phi]VEVsbits,\[Phi]VEVssigns,\[Phi]VEVsVals,\[Phi]VEVsVec,\[CapitalPhi]VEVsbits,\[CapitalPhi]VEVssigns,\[CapitalPhi]VEVsVals,\[CapitalPhi]VEVsVec,\[Phi]rules,\[CapitalPhi]rules,\[Phi]VEVrules},

(* First find the correponding vector to obtain the VEVs. *)
dim\[Phi]bits=("nSinglet"+"nNonperts")/.Options[stringFNmodels];
\[Phi]VEVlst=Take[state,-dim\[Phi]bits];
\[Phi]VEVsbits=Take[\[Phi]VEVlst,("nSinglet")/.Options[stringFNmodels]];
\[Phi]VEVsVals=Abs[\[Phi]VEVsbits];
\[Phi]VEVsVec=Map[("\[Phi]VEVmin"/.Options[stringFNmodels])*("\[Phi]VEVlogGen"/.Options[stringFNmodels])^#&,\[Phi]VEVsVals];

\[CapitalPhi]VEVsbits=Take[\[Phi]VEVlst,-("nNonperts")/.Options[stringFNmodels]];
\[CapitalPhi]VEVsVals=Abs[\[CapitalPhi]VEVsbits];
\[CapitalPhi]VEVsVec=Map[("\[CapitalPhi]VEVmin"/.Options[stringFNmodels])*("\[CapitalPhi]VEVlogGen"/.Options[stringFNmodels])^#&,\[CapitalPhi]VEVsVals];
(*
\[Phi]VEVsVec=Map[#[[1]]+I #[[2]]&,Partition[\[Phi]VEVsVec,2]];
\[CapitalPhi]VEVsVec=Map[#[[1]]+I #[[2]]&,Partition[\[CapitalPhi]VEVsVec,2]];
*)
(* Now obtain the \[Phi]VEV rules. *)
\[Phi]rules=Table[("BundleModuliFields"/.Options[stringFNmodels])[i]->\[Phi]VEVsVec[[i]],{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]rules=Table[("KahlerEffFields"/.Options[stringFNmodels])[i]->\[CapitalPhi]VEVsVec[[i]],{i,1,"nNonperts"/.Options[stringFNmodels]}];
\[Phi]VEVrules=Flatten[Join[\[Phi]rules,\[CapitalPhi]rules]];

\[Phi]VEVrules
]



(* ::Subsection::Closed:: *)
(*AssocFormToStateVec*)


AssocFormToStateVec[assoc_]:=Module[{Q10,Q5b,Q5bH,Qq,Qu,Qe,QL,Qd,QHd,QHu,sv10,sv5b,sv5bH,
numkahls,numperts,Q1,Q\[Phi],svbm,svkm,Q1vec,needQ1,QKah,numU1s,
bmlog0,kmlog0,bmloggen,kmloggen,bmvevs,kmvevs,svbmvev,svkmvev,sv
},
(* first read 10 charges. *)
If[KeyExistsQ[assoc,"Q10"],
Q10=assoc[["Q10"]];
sv10=Flatten[Q10];
If[Length[sv10]!=3,Print["Number of 10 charges incorrect; exiting.";Return[0]]];
,
If[KeyExistsQ[assoc,"Qq"],
Qq=assoc[["Qq"]],
Print["No Q10 associations, exiting."];
Return[0];
];
If[KeyExistsQ[assoc,"Qu"],Qu=assoc[["Qu"]];
If[Qq!=Qu,
Print["10 charges are not equal. Exiting Module."];
Return[0];
];
];
If[KeyExistsQ[assoc,"Qe"],Qe=assoc[["Qe"]];
If[Qq!=Qe,
Print["10 charges are not equal. Exiting Module."];
Return[0];
];
];
(* construct 10 part of sv *)
sv10=Flatten[Map[Position[#,1]&,Qq]];
];

(* read 5bar charge *)
If[KeyExistsQ[assoc,"Q5b"],
Q5b=assoc[["Q5b"]];
If[Dimensions[Q5b]!={3,2},Print["Q5b not in correct format, exiting."];Return[0]];
sv5b=Flatten[Transpose[Q5b]];
,
If[KeyExistsQ[assoc,"QL"],
QL=assoc[["QL"]],
Print["No Q5bar associations, exiting."];
Return[0];
];
If[KeyExistsQ[assoc,"Qd"],Qd=assoc[["Qd"]];
If[QL!=Qd,
Print["5 bar charges are not equal. Exiting Module."];
Return[0];
];
];
(* construct 5b part of sv *)
sv5b=Flatten[Transpose[Partition[Flatten[Map[If[Position[#,2]!={},
{Position[#,2],Position[#,2]},
Position[#,1]]&,QL]],2]]];

];

If[KeyExistsQ[assoc,"Q5bH"],
Q5bH=assoc[["Q5bH"]];
sv5bH=Flatten[Q5bH];
If[Length[sv5bH]!=2,Print["Number of 5bar charges for Higgs incorrect; exiting.";Return[0]]];
,
If[KeyExistsQ[assoc,"QHd"],
(* read 5barH charge *)
QHd=assoc[["QHd"]];
If[KeyExistsQ[assoc,"QHu"],
QHu=assoc[["QHu"]];
If[QHd!=-QHu,
Print["5 bar Higgs charges are not equal. Exiting Module."];
Return[0];
];
];
sv5bH=Flatten[Map[If[Position[#,2]!={},{Position[#,2],Position[#,2]},Position[#,1]]&,QHd]];
,
If[KeyExistsQ[assoc,"QHu"],
QHu=assoc[["QHu"]];
sv5bH=-Flatten[Map[If[Position[#,2]!={},{Position[#,2],Position[#,2]},Position[#,1]]&,QHu]];,
Print["No 5 bar Higgs charges entered; Exiting."];
Return[0];
]
];
];

numperts="nSinglet"/.Options[stringFNmodels];
numkahls="nNonperts"/.Options[stringFNmodels];
numU1s="nU1s"/.Options[stringFNmodels];

(* now do the bundle moduli charges *)
If[KeyExistsQ[assoc,"Q1"],
Q1=assoc[["Q1"]];
If[Dimensions[Q1]\[NotEqualTilde]{numperts,2},Print["Q1 charge vector dimensions are incorrect; exiting."];Return[0]];
svbm=Flatten[Transpose[Q1]];
,
(* ask if Q\[Phi] exists *)
If[KeyExistsQ[assoc,"Q\[Phi]"]&&(numperts!=0),
Q\[Phi]=assoc[["Q\[Phi]"]];
Q1vec=Map[If[Position[#,1]=={}&&Position[#,-1]=={},
{1,1},
(* otherwise we need to call Q1 *)
{Position[#,1],Position[#,-1]}
]
&,Q\[Phi]];
needQ1=Not[AllTrue[Q1vec,#!={1,1}&]];
If[needQ1,
Print["Q\[Phi] chargevector has ambiguity so therefore cannot be resolved. Please enter charge in assockey Q1 to clarify. Will assume e1."];
];
svbm=Flatten[Transpose[Q1vec]];
,
Print["No 1 charges entered; exiting."];
Return[0];
];
];

(* kahler charges *)
If[KeyExistsQ[assoc,"Q\[CapitalPhi]"],
QKah=assoc[["Q\[CapitalPhi]"]];
If[Dimensions[QKah]\[NotEqualTilde]{numkahls,numU1s},Print["Q\[CapitalPhi] charge vector dimensions are incorrect; exiting."];Return[0]];
svkm=Flatten[QKah];
,
If[numkahls==0,
svkm={};
,
Print["Q\[CapitalPhi] charge vector does not exist; exiting."];
Return[0];];
];

(* VEVs *)
(*bmlog0="\[Phi]VEVmin"/.Options[stringFNmodels];
kmlog0="\[CapitalPhi]VEVmin"/.Options[stringFNmodels];
bmloggen="\[Phi]VEVlogGen"/.Options[stringFNmodels];
kmloggen="\[CapitalPhi]VEVlogGen"/.Options[stringFNmodels];*)

(* check that there are entries *)
If[numperts!=0,
If[KeyExistsQ[assoc,"\[Phi]VEVs"],
bmvevs=assoc[["\[Phi]VEVs"]];
If[Length[bmvevs]!=numperts,Print["\[Phi]VEVs length do not match with number of bundle moduli; exiting."];Return[0]];
(*svbmvev=Flatten[Map[Round[Log[(#/bmlog0)]/Log[bmloggen]]&,bmvevs]];*)
svbmvev=bmvevs;
,
Print["\[Phi]VEVs does not exist; exiting."];Return[0];
];
,
svbmvev={};
];

If[numkahls!=0,
If[KeyExistsQ[assoc,"\[CapitalPhi]VEVs"],
kmvevs=assoc[["\[CapitalPhi]VEVs"]];
If[Length[kmvevs]!=numkahls,Print["\[CapitalPhi]VEVs length do not match with number of bundle moduli; exiting."];Return[0]];
(*svkmvev=Flatten[Map[Round[Log[(#/kmlog0)]/Log[kmloggen]]&,kmvevs]];*)
svkmvev=kmvevs;
,
Print["\[CapitalPhi]VEVs does not exist; exiting."];Return[0];
];
,
svkmvev={};
];
sv=Join[sv10,sv5b,sv5bH,svbm,svkm,svbmvev,svkmvev];
Return[sv]
]


(* ::Subsection:: *)
(*AssocFormToVEVs*)


AssocFormToVEVs[assoc_]:=Module[{\[Phi]VEVs,\[CapitalPhi]VEVs,numperts,numkahls,fields,len\[Phi]VEVs,len\[CapitalPhi]VEVs,\[Phi]rules,\[CapitalPhi]rules,\[Phi]VEVrules},
(* first check that the association keys exists. *)
If[Not[KeyExistsQ[assoc,"\[Phi]VEVs"]]||Not[KeyExistsQ[assoc,"\[CapitalPhi]VEVs"]]=={},
If[Not[KeyExistsQ[assoc,"\[Phi]VEVs"]],
Print["\[Phi]VEVs key does not exist. Exiting."];
];
If[Not[KeyExistsQ[assoc,"\[CapitalPhi]VEVs"]],
Print["\[CapitalPhi]VEVs key does not exist. Exiting."];
];
Return[0];
];

(* check the correct dimensions. *)
\[Phi]VEVs=assoc[["\[Phi]VEVs"]];
\[CapitalPhi]VEVs=assoc[["\[CapitalPhi]VEVs"]];
numperts="nSinglet"/.Options[stringFNmodels];
numkahls="nNonperts"/.Options[stringFNmodels];
len\[Phi]VEVs=Length[\[Phi]VEVs];
len\[CapitalPhi]VEVs=Length[\[CapitalPhi]VEVs];
If[len\[Phi]VEVs!=numperts,Print["\[Phi]VEVs dimensions incorrect, exiting."];Return[0]];
If[len\[CapitalPhi]VEVs!=numkahls,Print["\[CapitalPhi]VEVs dimensions incorrect, exiting."];Return[0]];

(* write the rules *)
\[Phi]rules=Table[("BundleModuliFields"/.Options[stringFNmodels])[i]->\[Phi]VEVs[[i]],{i,1,numperts}];
\[CapitalPhi]rules=Table[("KahlerEffFields"/.Options[stringFNmodels])[i]->\[CapitalPhi]VEVs[[i]],{i,1,numkahls}];

\[Phi]VEVrules=Flatten[Join[\[Phi]rules,\[CapitalPhi]rules]];

Return[\[Phi]VEVrules];
]


(* ::Subsection::Closed:: *)
(*ConvertStateVecToVEVPowers*)


ConvertStateVecToVEVPowers[stateVec_]:=Module[{dimVEVs,VEVlst,\[Phi]VEVslst,\[CapitalPhi]VEVslst,assocVEVList},
dimVEVs=("nSinglet"+"nNonperts")/.Options[stringFNmodels];
VEVlst=Take[stateVec,-dimVEVs];
\[Phi]VEVslst=Take[VEVlst,(("nSinglet")/.Options[stringFNmodels])];
\[CapitalPhi]VEVslst=Take[VEVlst,-(("nNonperts")/.Options[stringFNmodels])];
assocVEVList=<|"\[Phi]VEVPowers"->\[Phi]VEVslst,"\[CapitalPhi]VEVPowers"->\[CapitalPhi]VEVslst|>;
Return[assocVEVList]
]


(* ::Section:: *)
(*Auxiliary Modules - State Generation for GAbitlst*)


(* ::Subsection::Closed:: *)
(*GADetermDimStateVecwithVEVs*)


GADetermDimStateVecwithVEVs[]:=Module[{dimchargelst,num,dimension,\[CapitalPhi]Dim,chargerange,total\[CapitalPhi]dim,dimVEVbits},
(* This module determines the dimension of the binary vector to be generated in GA. *)

(* First we find the dimension of the charge vector. *)
dimchargelst=3*3+2+2*"nSinglet"/.Options[stringFNmodels];
num=("nU1s"/.Options[stringFNmodels])^dimchargelst;
dimension=Ceiling[Log[2,num]];

(* Now add in the dimension for \[CapitalPhi] sector. *)
chargerange="\[CapitalPhi]ChargeRange"/.Options[stringFNmodels];
\[CapitalPhi]Dim=(chargerange[[2]]-chargerange[[1]]+1);
\[CapitalPhi]Dim=Log[2,\[CapitalPhi]Dim];
total\[CapitalPhi]dim=\[CapitalPhi]Dim*("nNonperts"/.Options[stringFNmodels])*(("nU1s"/.Options[stringFNmodels])-1);
dimension=dimension+total\[CapitalPhi]dim;

(* Add in dimension for \[Phi]VEVs. *)
dimVEVbits=("nSinglet"*"bitsfor\[Phi]"+"nNonperts"*"bitsfor\[CapitalPhi]")/.Options[stringFNmodels];
dimension=dimVEVbits+dimension;

(* Returns dimension. *)
{dimension,\[CapitalPhi]Dim,dimVEVbits}
]


(* ::Subsection::Closed:: *)
(*GARandomStatewithVEVs*)


GARandomStatewithVEVs[]:=Module[{dim,bitlst},

(* Random 2 state for system. *)
dim=GADetermDimStateVecwithVEVs[];
bitlst=Table[RandomInteger[{0,1}],{dim[[1]]}];

bitlst
]


(* ::Subsection::Closed:: *)
(*GAConvertBitlstStateVecwithVEVs*)


GAConvertBitlstStateVecwithVEVs[bitlst_]:=Module[{dim,\[CapitalPhi]dim,tot\[CapitalPhi]dim,MSSMdim,\[CapitalPhi]chargevec,MSSMqvec,dimchargelst,MSSMqvecnew,fullQvec,
\[Phi]VEVdim,statevec\[Phi]VEV,fullstatevec,\[Phi]VEVstate,\[Phi]VEVsbits,\[Phi]VEVssigns,\[Phi]VEVsVec,\[CapitalPhi]VEVsbits,\[CapitalPhi]VEVssigns,\[CapitalPhi]VEVsVec},
(* This module converts a bitlst into the charge vector form to be used in the module. *)

(* Find dimension for \[CapitalPhi] sector. *)
dim=GADetermDimStateVecwithVEVs[];
\[CapitalPhi]dim=dim[[2]];
\[Phi]VEVdim=dim[[3]];
tot\[CapitalPhi]dim=\[CapitalPhi]dim*("nNonperts"/.Options[stringFNmodels])*("nU1s"/.Options[stringFNmodels]);
MSSMdim=dim[[1]]-tot\[CapitalPhi]dim-\[Phi]VEVdim;

(* \[CapitalPhi] sector. *)
If[("nNonperts"/.Options[stringFNmodels])==0,
\[CapitalPhi]chargevec={};
,
\[CapitalPhi]chargevec=Partition[Take[bitlst,{MSSMdim+1,MSSMdim+tot\[CapitalPhi]dim}],\[CapitalPhi]dim];
\[CapitalPhi]chargevec=FromDigits[Transpose[\[CapitalPhi]chargevec],2]-(2^(\[CapitalPhi]dim-1)-1);
];



(* MSSM sector. *)
MSSMqvec=Take[bitlst,MSSMdim];
dimchargelst=3*3+2+2*"nSinglet"/.Options[stringFNmodels];
MSSMqvec=IntegerDigits[FromDigits[MSSMqvec,2],"nU1s"/.Options[stringFNmodels]];
MSSMqvec=MSSMqvec+Table[1,Dimensions[MSSMqvec][[1]]];

If[Dimensions[MSSMqvec][[1]]<dimchargelst,
MSSMqvecnew=Join[Table[1,(dimchargelst-Dimensions[MSSMqvec][[1]])],MSSMqvec];,
	If[Dimensions[MSSMqvec][[1]]>dimchargelst,
	MSSMqvecnew=Take[MSSMqvec,-dimchargelst],
	MSSMqvecnew=MSSMqvec
	];
];


(* Constuct full charge vector. *)
fullQvec=Join[MSSMqvecnew,\[CapitalPhi]chargevec];

(* Construct the \[Phi]VEV sector. *)
\[Phi]VEVstate=Take[bitlst,-\[Phi]VEVdim];
(* Now determine the new chargevector for \[Phi]VEVs.*)
\[Phi]VEVsbits=Partition[Take[\[Phi]VEVstate,("nSinglet"*"bitsfor\[Phi]")/.Options[stringFNmodels]],"bitsfor\[Phi]"/.Options[stringFNmodels]];
\[Phi]VEVsVec=FromDigits[Transpose[\[Phi]VEVsbits],2];
If[\[Phi]VEVsVec==0,\[Phi]VEVsVec={}];
(* Repeat the same for \[CapitalPhi] part. *)
\[CapitalPhi]VEVsbits=Partition[Take[\[Phi]VEVstate,-("nNonperts"*"bitsfor\[CapitalPhi]")/.Options[stringFNmodels]],"bitsfor\[CapitalPhi]"/.Options[stringFNmodels]];
\[CapitalPhi]VEVsVec=FromDigits[Transpose[\[CapitalPhi]VEVsbits],2];
If[\[CapitalPhi]VEVsVec==0,\[CapitalPhi]VEVsVec={}];

statevec\[Phi]VEV=Join[\[Phi]VEVsVec,\[CapitalPhi]VEVsVec];

(* Constrcut full statevector. *)
fullstatevec=Join[fullQvec,statevec\[Phi]VEV]

]


(* ::Subsection::Closed:: *)
(*GADetermDimStateVecFix*)


GADetermDimStateVecFix[]:=Module[{dimchargelst,num,dimension,\[CapitalPhi]Dim,chargerange,total\[CapitalPhi]dim,dimVEVbits},
(* This module determines the dimension of the binary vector to be generated in GA. *)

(* First we find the dimension of the charge vector. *)
dimchargelst=3*2+2*"nSinglet"/.Options[stringFNmodels];
num=("nU1s"/.Options[stringFNmodels])^dimchargelst;
dimension=Ceiling[Log[2,num]];


(* Now add in the dimension for \[CapitalPhi] sector. *)
chargerange="\[CapitalPhi]ChargeRange"/.Options[stringFNmodels];
\[CapitalPhi]Dim=(chargerange[[2]]-chargerange[[1]]+1);
\[CapitalPhi]Dim=Log[2,\[CapitalPhi]Dim];
total\[CapitalPhi]dim=\[CapitalPhi]Dim*("nNonperts"/.Options[stringFNmodels])*(("nU1s"/.Options[stringFNmodels])-1);
dimension=dimension+total\[CapitalPhi]dim;

(* Add in dimension for \[Phi]VEVs. *)
dimVEVbits=("nSinglet"*"bitsfor\[Phi]"+"nNonperts"*"bitsfor\[CapitalPhi]")/.Options[stringFNmodels];
dimension=dimVEVbits+dimension;

(* Returns dimension. *)
{dimension,\[CapitalPhi]Dim,dimVEVbits}
]


(* ::Subsection::Closed:: *)
(*GARandomStateFix*)


GARandomStateFix[]:=Module[{dim,bitlst},

(* Random 2 state for system. *)
dim=GADetermDimStateVecFix[];
bitlst=Table[RandomInteger[{0,1}],{dim[[1]]}];

bitlst
]


(* ::Subsection::Closed:: *)
(*GAConvertBitlstStateVecFix*)


GAConvertBitlstStateVecFix[bitlst_]:=Module[{dim,\[CapitalPhi]dim,tot\[CapitalPhi]dim,MSSMdim,\[CapitalPhi]chargevec,MSSMqvec,dimchargelst,MSSMqvecnew,fullQvec,
\[Phi]VEVdim,statevec\[Phi]VEV,fullstatevec,\[Phi]VEVstate,\[Phi]VEVsbits,\[Phi]VEVssigns,\[Phi]VEVsVec,\[CapitalPhi]VEVsbits,\[CapitalPhi]VEVssigns,\[CapitalPhi]VEVsVec,
tenCharge,fiveBarHiggsCharge,MSSMqvecnewComb,statevec\[Phi]VEVNew,Max\[CapitalPhi]Charge,
lastDigit\[CapitalPhi]QVec},
(* This module converts a bitlst into the charge vector form to be used in the module. *)

(* Find dimension for \[CapitalPhi] sector. *)
If[("fixedCharges"/.Options[stringFNmodels]),
dim=GADetermDimStateVecFix[];
,
dim=GADetermDimStateVecwithVEVs[];
];

\[CapitalPhi]dim=dim[[2]];
\[Phi]VEVdim=dim[[3]];
tot\[CapitalPhi]dim=\[CapitalPhi]dim*("nNonperts"/.Options[stringFNmodels])*(("nU1s"/.Options[stringFNmodels])-1);
MSSMdim=dim[[1]]-tot\[CapitalPhi]dim-\[Phi]VEVdim;
Max\[CapitalPhi]Charge=-("\[CapitalPhi]ChargeRange"/.Options[stringFNmodels])[[1]];

(* \[CapitalPhi] sector. *)
If[("nNonperts"/.Options[stringFNmodels])===0,
\[CapitalPhi]chargevec={};
,
\[CapitalPhi]chargevec=Partition[Take[bitlst,{MSSMdim+1,MSSMdim+tot\[CapitalPhi]dim}],\[CapitalPhi]dim];
\[CapitalPhi]chargevec=FromDigits[Transpose[\[CapitalPhi]chargevec],2]-Max\[CapitalPhi]Charge;
(* determine last digit. *)
lastDigit\[CapitalPhi]QVec=Map[Total,Partition[\[CapitalPhi]chargevec,(("nU1s"/.Options[stringFNmodels])-1)]];
lastDigit\[CapitalPhi]QVec=Map[{-#}&,lastDigit\[CapitalPhi]QVec];
\[CapitalPhi]chargevec=Flatten[Flatten[Join[Partition[\[CapitalPhi]chargevec,(("nU1s"/.Options[stringFNmodels])-1)],lastDigit\[CapitalPhi]QVec,2]]];
];



(* MSSM sector. *)
MSSMqvec=Take[bitlst,MSSMdim];
dimchargelst=3*2+2*"nSinglet"/.Options[stringFNmodels];
MSSMqvec=IntegerDigits[FromDigits[MSSMqvec,2],"nU1s"/.Options[stringFNmodels]];
MSSMqvec=MSSMqvec+Table[1,Dimensions[MSSMqvec][[1]]];

If[Dimensions[MSSMqvec][[1]]<dimchargelst,
MSSMqvecnew=Join[Table[1,(dimchargelst-Dimensions[MSSMqvec][[1]])],MSSMqvec];,
	If[Dimensions[MSSMqvec][[1]]>dimchargelst,
	MSSMqvecnew=Take[MSSMqvec,-dimchargelst],
	MSSMqvecnew=MSSMqvec
	];
];

(* now add in the 10 and 5 bar charges *)
tenCharge=("tenQ"/.Options[stringFNmodels]);
fiveBarHiggsCharge=("fiveBarHiggsQ"/.Options[stringFNmodels]);
MSSMqvecnewComb=Join[tenCharge,Take[MSSMqvecnew,6],fiveBarHiggsCharge,Take[MSSMqvecnew,-(2*"nSinglet"/.Options[stringFNmodels])]];

(* Constuct full charge vector. *)
fullQvec=Join[MSSMqvecnewComb,\[CapitalPhi]chargevec];

(* Construct the \[Phi]VEV sector. *)
\[Phi]VEVstate=Take[bitlst,-\[Phi]VEVdim];
(* Now determine the new chargevector for \[Phi]VEVs.*)
\[Phi]VEVsbits=Partition[Take[\[Phi]VEVstate,("nSinglet"*"bitsfor\[Phi]")/.Options[stringFNmodels]],"bitsfor\[Phi]"/.Options[stringFNmodels]];
\[Phi]VEVsVec=FromDigits[Transpose[\[Phi]VEVsbits],2];
If[\[Phi]VEVsVec==0,\[Phi]VEVsVec={}];
(* Repeat the same for \[CapitalPhi] part. *)
\[CapitalPhi]VEVsbits=Partition[Take[\[Phi]VEVstate,-("nNonperts"*"bitsfor\[CapitalPhi]")/.Options[stringFNmodels]],"bitsfor\[CapitalPhi]"/.Options[stringFNmodels]];
\[CapitalPhi]VEVsVec=FromDigits[Transpose[\[CapitalPhi]VEVsbits],2];
If[\[CapitalPhi]VEVsVec==0,\[CapitalPhi]VEVsVec={}];
\[CapitalPhi]VEVsVec=\[CapitalPhi]VEVsVec+Table[2,Length[\[CapitalPhi]VEVsVec]];

statevec\[Phi]VEV=Join[\[Phi]VEVsVec,\[CapitalPhi]VEVsVec];

(* add 1s to statevector *)
statevec\[Phi]VEVNew=statevec\[Phi]VEV+Table[1,Length[statevec\[Phi]VEV]];

(* Constrcut full statevector. *)
fullstatevec=Join[fullQvec,statevec\[Phi]VEVNew]

]


(* ::Subsection::Closed:: *)
(*GAConvertStateVectoBitlst*)


GAConvertStateVectoBitlst[stateVec_]:=Module[{smQPart,VEVPart,numPerts,numKahls,numVEVs,numSmQ,
\[Phi]VEVpart,\[CapitalPhi]VEVpart,bits\[Phi]VEVs,bits\[CapitalPhi]VEVs,bitsVEVs,
numSmQPart,numU1,smQPartToConv,smVecBits,lenBitsVEVs,reqLength,bitsSMVec,allBits,SmQPartNew,
\[CapitalPhi]ChargePart,num\[CapitalPhi]Qs,bits\[CapitalPhi]Qs,num\[CapitalPhi]QBitsPer,\[CapitalPhi]ChargePartToConv,Max\[CapitalPhi]Charge
},
numPerts="nSinglet"/.Options[stringFNmodels];
numKahls="nNonperts"/.Options[stringFNmodels];
numVEVs=numPerts+numKahls;
numSmQ=11+2*numPerts;
smQPart=Take[stateVec,numSmQ];
VEVPart=Take[stateVec,-numVEVs]-Table[1,numVEVs];
\[Phi]VEVpart=Take[VEVPart,numPerts];
\[CapitalPhi]VEVpart=Take[VEVPart,-numKahls];
\[CapitalPhi]VEVpart=\[CapitalPhi]VEVpart-Table[2,Length[\[CapitalPhi]VEVpart]];

Max\[CapitalPhi]Charge=-("\[CapitalPhi]ChargeRange"/.Options[stringFNmodels])[[1]];

(* convert the VEV parts *)
bits\[Phi]VEVs=Flatten[(Map[PadLeft[IntegerDigits[#,2],("bitsfor\[Phi]"/.Options[stringFNmodels])]&,\[Phi]VEVpart])];
bits\[CapitalPhi]VEVs=Flatten[(Map[PadLeft[IntegerDigits[#,2],("bitsfor\[CapitalPhi]"/.Options[stringFNmodels])]&,\[CapitalPhi]VEVpart])];

bitsVEVs=Join[bits\[Phi]VEVs,bits\[CapitalPhi]VEVs];
lenBitsVEVs=Length[bitsVEVs];

(* now convert the sm charge part *)
If[("fixedCharges"/.Options[stringFNmodels]),
SmQPartNew=Join[Take[smQPart,{4,9}],Take[smQPart,-(2*numPerts)]];
numSmQPart=Length[SmQPartNew];
smQPartToConv=SmQPartNew-Table[1,numSmQPart];
,
numSmQPart=Length[smQPart];
smQPartToConv=smQPart-Table[1,numSmQPart];
];

numU1=("nU1s"/.Options[stringFNmodels]);
smVecBits=IntegerDigits[(FromDigits[smQPartToConv,numU1]),2];

(* \[CapitalPhi]ChargePart conversion *)
If[numKahls===0,
bits\[CapitalPhi]Qs={};
,
num\[CapitalPhi]Qs=numU1*numKahls;
\[CapitalPhi]ChargePart=Take[stateVec,{numSmQ+1,numSmQ+num\[CapitalPhi]Qs}];
\[CapitalPhi]ChargePart=\[CapitalPhi]ChargePart+Table[Max\[CapitalPhi]Charge,Length[\[CapitalPhi]ChargePart]];
\[CapitalPhi]ChargePartToConv=Flatten[Flatten[Map[Delete[#,numU1]&,Partition[\[CapitalPhi]ChargePart,numU1]]]];
num\[CapitalPhi]QBitsPer=Ceiling[(Log[(("\[CapitalPhi]ChargeRange"/.Options[stringFNmodels])[[2]]-("\[CapitalPhi]ChargeRange"/.Options[stringFNmodels])[[1]])]/Log[2])];
bits\[CapitalPhi]Qs=Flatten[(Map[PadLeft[IntegerDigits[#,2],num\[CapitalPhi]QBitsPer]&,\[CapitalPhi]ChargePartToConv])];
];
If[("fixedCharges"/.Options[stringFNmodels]),
reqLength=GADetermDimStateVecFix[][[1]]-lenBitsVEVs-Length[bits\[CapitalPhi]Qs],
reqLength=GADetermDimStateVecwithVEVs[][[1]]-lenBitsVEVs-Length[bits\[CapitalPhi]Qs]
];
bitsSMVec=PadLeft[smVecBits,reqLength];

(* construct full bitlist *)
allBits=Join[bitsSMVec,bits\[CapitalPhi]Qs,bitsVEVs];

Return[allBits]

]


(* ::Section::Closed:: *)
(*Auxiliary Modules - Finding Operators (Quark Sector)*)


(* ::Subsection::Closed:: *)
(*QMatRepChange*)


QMatRepChange[Qchargevec_,nchargevec_]:=Module[{Fnominus1,Qreducedvec,n0,q0,nredvec,newvec},
Fnominus1="nU1s"-1/.Options[stringFNmodels];
Qreducedvec=Take[Qchargevec,-Fnominus1];
nredvec=Take[nchargevec,-Fnominus1];
q0=Qchargevec[[1]];
n0=nchargevec[[1]];
newvec=n0*Qreducedvec-q0*nredvec;
newvec
]


(* ::Subsection::Closed:: *)
(*FromPowertoOps*)


FromPowertoOps[pvec_]:=Module[{fieldvec,output},
(* This is a small function to change from the power of the fields to the required operator form. *)
fieldvec=Flatten[Join[Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}],Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}]]];
output=Inner[Power,fieldvec,pvec,Times]
]


(* ::Subsection::Closed:: *)
(*FindOpsforWIJaListwith\[Phi]VEVsQsector*)


FindOpsforWIJaListwith\[Phi]VEVsQsector[state_]:=Module[{nFields,\[Phi]fields,\[CapitalPhi]fields,fields,solvqmat,Neqns,QtableQHd,QtableQHu,QtableLHe,QtableLH,Qtablesolve,eqntosolve,powersolsforW,opsforW,
j,dimOplist,O1coeffneed,newO1coeffneed,Oplist,newO1count,raggedOplist,i,Convertfieldslst\[CapitalPhi],Wexpansion,Adjoin\[Kappa]Rule,Adjoin\[Sigma]Rule,nchargevec,
powersphi,powersPhi,powerslst,chargelst,ruleslst,viablerules,goodpowers},

(* First calcualte the total number of fields. *)
nFields= "nSinglet"+"nNonperts"/.Options[stringFNmodels];

(* Run StateVectoMat to get the charge matrix. *)
solvqmat=StateVectoMatwithVEVs[state];

(* We now want to do this as described in the notes. *)
(* Find the negative of the charge of the respective SM fields in Yukawa terms in the superpotential. *)
QtableQHd=Join[Table[solvqmat[[1]]+solvqmat[[i+3]]+solvqmat[[16]],{i,1,3}],Table[solvqmat[[2]]+solvqmat[[i+3]]+solvqmat[[16]],{i,1,3}],Table[solvqmat[[3]]+solvqmat[[i+3]]+solvqmat[[16]],{i,1,3}]];
QtableQHu=Join[Table[solvqmat[[1]]+solvqmat[[i+6]]+solvqmat[[17]],{i,1,3}],Table[solvqmat[[2]]+solvqmat[[i+6]]+solvqmat[[17]],{i,1,3}],Table[solvqmat[[3]]+solvqmat[[i+6]]+solvqmat[[17]],{i,1,3}]];
(* Append the matrix and solve for the system. *)
Qtablesolve=Join[QtableQHd,QtableQHu];

(* The fields to solve for the charges is the \[Phi] and \[CapitalPhi] expansions. *)
solvqmat=Take[solvqmat,-nFields];
(* Produce a list of fields. *)
fields=Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,nFields,1}];
\[Phi]fields=Take[fields,"nSinglet"/.Options[stringFNmodels]];
\[CapitalPhi]fields=Take[fields,-"nNonperts"/.Options[stringFNmodels]];

(* Charge matrix reduction to representations of J. *)
nchargevec="nChargeVec"/.Options[stringFNmodels];
Qtablesolve=Map[QMatRepChange[#,nchargevec]&,Qtablesolve,1];
solvqmat=Map[QMatRepChange[#,nchargevec]&,solvqmat,1];

(* Calculate the number of eqns. *)

(* In this module we will list out all the operators and compare them to the minus the operator charge. *)
(* First compute all powers of fields. *)
(*powersphi=Tuples[Table[i,{i,0,("DimOp"/.Options[stringFNmodels])}],("nSinglet"/.Options[stringFNmodels])];
powersphi=Select[powersphi,Total[#]<=("DimOp"/.Options[stringFNmodels])&];
powersPhi=Tuples[Table[i,{i,0,("DimOpof\[CapitalPhi]"/.Options[stringFNmodels])}],("nNonperts"/.Options[stringFNmodels])];
powersPhi=Select[powersPhi,Total[#]<=("DimOpof\[CapitalPhi]"/.Options[stringFNmodels])&];
powerslst=Map[Flatten,Tuples[{powersphi,powersPhi}]];
*)
powersphi=Tuples[Table[i,{i,0,("DimOp"/.Options[stringFNmodels])}],nFields];
powersphi=Map[Reverse,powersphi];
powersphi=Select[powersphi,Total[Take[#,("nSinglet"/.Options[stringFNmodels])]]<=("DimOp"/.Options[stringFNmodels])&];
powerslst=Select[powersphi,Total[Take[#,-("nNonperts"/.Options[stringFNmodels])]]<=("DimOpof\[CapitalPhi]"/.Options[stringFNmodels])&];

(* generate the list of all charges. *)
chargelst=powerslst . solvqmat;
ruleslst=Thread[Rule[powerslst,chargelst]];

(* Now select the cases which match the required operator charge. *)
viablerules=Function[x,Select[ruleslst,#[[2]]==-x&]]/@Qtablesolve;
goodpowers=If[#=={},{},Map[First,#]]&/@viablerules;
Oplist=Map[If[#=={},0,Map[FromPowertoOps,#]]&,goodpowers];

Oplist
]



(* ::Subsection::Closed:: *)
(*UpdateVectorWithMinimum*)


UpdateVectorWithMinimum[vec_,posList_]:=Module[{minValue},
Fold[(minValue=Min[vec[[#2]]];ReplacePart[#,Thread[#2->minValue]])&,vec,posList]];


(* ::Subsection::Closed:: *)
(*MinReplaceListRule*)


MinReplaceListRule[vec_,posList_]:=Module[{replaceZeroList,replaceOneList},
replaceZeroList=Flatten[Map[Complement[#,{#[[Position[vec[[#]],Min[vec[[#]]]][[1]][[1]]]]}]&,posList]];
replaceOneList=Flatten[Complement[Range[Length[vec]],replaceZeroList]];
Total[IdentityMatrix[Length[vec]][[replaceOneList]]]
(*Range[Length[vec]]/.Join[Thread[replaceZeroList->0],Thread[replaceOneList->1]]*)
]


(* ::Subsection::Closed:: *)
(*FindOpsFromGraphMethod*)


FindOpsFromGraphMethod[smoduliCharges_,smoduliVEVs_,chargesToSolve_]:=Module[{nSings,nKahls,numU1s,singQList,
duplicates,dupliElements,dupliPos,smoduliVEVsPos,smoduliVEVsEdge,
graf,edgesGraphRules,flowMat,sparseFlowMat,maxDim,powersS,
zeroVecQs},
nSings="nSinglet"/.Options[stringFNmodels];
nKahls="nNonperts"/.Options[stringFNmodels];
numU1s="nU1s"/.Options[stringFNmodels];
singQList=Transpose[Partition[smoduliCharges,nSings]];
maxDim="DimOp"/.Options[stringFNmodels];

(* We need to check the chargeToSolve - if it already is leading order power then we need to return the powers. *)
zeroVecQs=Table[0,numU1s];
If[chargesToSolve==zeroVecQs,
powersS=Table[0,nSings];
Return[powersS];
];

edgesGraphRules=Map[#[[2]]->#[[1]]&,singQList];

(* locate duplicates of the singlets and their positions *)
duplicates=Select[Tally[singQList],Last[#]>1&];
dupliElements=duplicates[[All,1]];
dupliPos=Map[#->Position[singQList,#]&,dupliElements];
If[dupliPos!={},
smoduliVEVsPos=Map[Flatten,(dupliElements/.dupliPos)];
smoduliVEVsEdge=UpdateVectorWithMinimum[smoduliVEVs,smoduliVEVsPos];
,
smoduliVEVsEdge=smoduliVEVs;
];

(* now replace the VEVs *)
graf=Graph[Table[i,{i,1,numU1s}],edgesGraphRules,EdgeCost->smoduliVEVsEdge,VertexLabels->Automatic,EdgeLabels->"EdgeTag"];
flowMat=FindMinimumCostFlow[graf,chargesToSolve,"FlowMatrix",EdgeCapacity->ConstantArray[10000(*this is just a large number*),Length[EdgeList[graf]]]];
sparseFlowMat=If[Total[Transpose[flowMat]-flowMat]===chargesToSolve//Normal,flowMat,Null];
(* return the operator as required. *)
If[sparseFlowMat===Null,
Return[Null]
];
If[
dupliPos!={},
powersS=(Map[{#[[2]],#[[1]]}&,singQList]/.ArrayRules[sparseFlowMat])*MinReplaceListRule[smoduliVEVs,smoduliVEVsPos];
If[Total[powersS]>maxDim,powersS=Null];
Return[powersS];
,
powersS=Map[{#[[2]],#[[1]]}&,singQList]/.ArrayRules[sparseFlowMat];
If[Total[powersS]>maxDim,powersS=Null];
Return[powersS]
]
]


(* ::Subsection::Closed:: *)
(*GenerateOrderedTuples*)


GenerateOrderedTuples[sum_,length_]:=Module[{combs,sortCombs},
combs=Select[Tuples[Range[0,sum],length],Total[#]<=sum&];
sortCombs=SortBy[combs,{Total[#],-#[[#]]}&];
sortCombs
]


(* ::Subsection::Closed:: *)
(*SumChargesByPowers*)


SumChargesByPowers[kPows_,kCharges_]:=Module[{kChargeRes},
kChargeRes=Transpose[kCharges] . kPows;
Return[kChargeRes]
]


(* ::Subsection::Closed:: *)
(*ProcessModPowers*)


ProcessModPowers[vec1_,vec2_]:=Module[{posNull,newVec1,newVec2,result},
posNull=Position[vec1,Null];
newVec1=Delete[vec1,posNull];
newVec2=Delete[vec2,posNull];
result=MapThread[Join,{newVec1,newVec2}];
result]


(* ::Subsection::Closed:: *)
(*FindOpsFromGraphAdd\[CapitalPhi]*)


FindOpsFromGraphAdd\[CapitalPhi][smoduliCharges_,smoduliVEVs_,chargesToSolve_,kCharges_]:=Module[{nSings,nKahls,numU1s,singQList,
dimKahl,listKPows,newQsToSolve,listSPows,listModPows,nVec,newkCharges},
nSings="nSinglet"/.Options[stringFNmodels];
nKahls="nNonperts"/.Options[stringFNmodels];
numU1s="nU1s"/.Options[stringFNmodels];
singQList=Transpose[Partition[smoduliCharges,nSings]];
dimKahl="DimOpof\[CapitalPhi]"/.Options[stringFNmodels];

(* Find correct representation of the KahlCharges *)
nVec="nChargeVec"/.Options[stringFNmodels];
newkCharges=Map[#-(Total[#]/5)*nVec&,kCharges];
listKPows=GenerateOrderedTuples[dimKahl,nKahls];
newQsToSolve=Map[SumChargesByPowers[#,newkCharges]+chargesToSolve&,listKPows];

listSPows=Map[FindOpsFromGraphMethod[smoduliCharges,smoduliVEVs,#]&,newQsToSolve];

listModPows=ProcessModPowers[listSPows,listKPows];
Return[listModPows];

]


(* ::Subsection::Closed:: *)
(*FindOpsforWIJaGraphwith\[Phi]VEVsQsector*)


FindOpsforWIJaGraphwith\[Phi]VEVsQsector[statevec_]:=Module[{nSings,nKahls,nFields,solvqmat,QtableQHd,QtableQHu,QTableSolve,nVec,NewQTableSolve,
svSingCharges,edgesRules,nU1s,graphToCompute,flowMatrix,svVEVs,KQTable,
listOfOps,Oplist},
(* number of total moduli fields *)
nSings="nSinglet"/.Options[stringFNmodels];
nKahls="nNonperts"/.Options[stringFNmodels];
nU1s="nU1s"/.Options[stringFNmodels];
nFields= nSings+nKahls;
(* charge matrix *)
solvqmat=StateVectoMatwithVEVs[statevec];

(* We now want to do this as described in the notes. *)
(* Find the negative of the charge of the respective SM fields in Yukawa terms in the superpotential. *)QtableQHd=Join[Table[solvqmat[[1]]+solvqmat[[i+3]]+solvqmat[[16]],{i,1,3}],Table[solvqmat[[2]]+solvqmat[[i+3]]+solvqmat[[16]],{i,1,3}],Table[solvqmat[[3]]+solvqmat[[i+3]]+solvqmat[[16]],{i,1,3}]];QtableQHu=Join[Table[solvqmat[[1]]+solvqmat[[i+6]]+solvqmat[[17]],{i,1,3}],Table[solvqmat[[2]]+solvqmat[[i+6]]+solvqmat[[17]],{i,1,3}],Table[solvqmat[[3]]+solvqmat[[i+6]]+solvqmat[[17]],{i,1,3}]];(* Append the matrix and solve for the system. *)
QTableSolve=Join[QtableQHd,QtableQHu];

(* transform each row of the charge matrix to be solved into the representative such that it adds up to zero. *)
nVec="nChargeVec"/.Options[stringFNmodels];
NewQTableSolve=Map[#-(Total[#]/5)*nVec&,QTableSolve];

(* extract singlet moduli charges *)
svSingCharges=Take[statevec,{12,11+(2*nSings)}];
svVEVs=Take[Take[statevec,-nFields],nSings];
edgesRules=Flatten[Map[#[[2]]->#[[1]]&,Transpose[Partition[svSingCharges,nSings]]]];
(* assume no kahler *)
If[nKahls==0
,
listOfOps=Map[FindOpsFromGraphMethod[svSingCharges,svVEVs,#]&,NewQTableSolve];
Oplist=Map[If[#===Null,0,{FromPowertoOps[#]}]&,listOfOps];
Return[Oplist]
,
KQTable=2*Take[solvqmat,-nKahls];
listOfOps=Map[FindOpsFromGraphAdd\[CapitalPhi][svSingCharges,svVEVs,#,KQTable]&,NewQTableSolve];
(*Print[listOfOps];*)
Oplist=Map[If[#==={},0,Map[FromPowertoOps,#]]&,listOfOps];
Return[Oplist]
]

]


(* ::Subsection::Closed:: *)
(*ComputeWOpsGivenStateVec*)


ComputeWOpsGivenStateVec[statevec_]:=Module[{WIJaOps},

(* Step 1: Find the relevant operators. *)
(* Note O1CLst is the O1vec required. *)
(* In this step we have now used the newer modules using the Smith decomposition. *)
(* Setting to turn on what type of method we are generating the operators. *)
If[(("MethodGenOps"/.Options[stringFNmodels])=="List")||(("MethodGenOps"/.Options[stringFNmodels])=="Graph")||(("MethodGenOps"/.Options[stringFNmodels])=="Optimal"),

If[("MethodGenOps"/.Options[stringFNmodels])=="List",
WIJaOps=FindOpsforWIJaListwith\[Phi]VEVsQsector[statevec];
];
If[("MethodGenOps"/.Options[stringFNmodels])=="Graph"||("MethodGenOps"/.Options[stringFNmodels])=="Optimal",
WIJaOps=FindOpsforWIJaGraphwith\[Phi]VEVsQsector[statevec];
];
,
Return["Error"];
];

WIJaOps

]


(* ::Section::Closed:: *)
(*Auxiliary Modules - Texture Analysis*)


(* ::Subsection::Closed:: *)
(*FromWOpsGetLeadingYukawaMatrix*)


FromWOpsGetLeadingYukawaMatrix[WOps_,VEVPowersListAssoc_]:=Module[{\[Phi]Fields,\[CapitalPhi]Fields,Scale,\[Phi]VEVsList,\[CapitalPhi]VEVsList,\[Phi]VEVsReplaceRule,\[CapitalPhi]VEVsReplaceRule,WOpswithVEVs,WOpsLeading,YukDownLeadPowers,YukUpLeadPowers,assocYukLeadPowers,WOpsLeadingRemZero},
Scale="VEVsScale"/.Options[stringFNmodels];
\[Phi]Fields=Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]Fields=Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}];
\[Phi]VEVsList=Map[Scale^#&,VEVPowersListAssoc[["\[Phi]VEVPowers"]]];
\[CapitalPhi]VEVsList=Map[Scale^#&,VEVPowersListAssoc[["\[CapitalPhi]VEVPowers"]]];
\[Phi]VEVsReplaceRule=Table[(\[Phi]Fields[[i]]->\[Phi]VEVsList[[i]]),{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]VEVsReplaceRule=Table[(\[CapitalPhi]Fields[[i]]->\[CapitalPhi]VEVsList[[i]]),{i,1,"nNonperts"/.Options[stringFNmodels]}];
WOpswithVEVs=WOps/.\[Phi]VEVsReplaceRule/.\[CapitalPhi]VEVsReplaceRule;
(* Now extract the lowest order term and return the matrices. *)
WOpsLeading=
Map[
If[#===0,
#,
LowestOrderOfMonomialListOneVar[#,Scale]
]&,WOpswithVEVs];

WOpsLeadingRemZero=Map[If[Exponent[#,Scale]>40||#===0,Scale^40,#]&,WOpsLeading];

YukDownLeadPowers=Partition[Take[WOpsLeadingRemZero,9],3];
YukUpLeadPowers=Partition[Take[WOpsLeadingRemZero,-9],3];
assocYukLeadPowers=<|"YukdLead"->YukDownLeadPowers,"YukuLead"->YukUpLeadPowers|>;
Return[assocYukLeadPowers]
]


(* ::Subsection::Closed:: *)
(*LowestOrderOfMonomialListOneVar*)


LowestOrderOfMonomialListOneVar[monomialList_,field_]:=Module[{lowestOrder,term},
lowestOrder=Min[Exponent[#,field]&/@monomialList];
term=Select[monomialList,Exponent[#,field]==lowestOrder&];
Return[term[[1]]];
]


(* ::Subsection::Closed:: *)
(*RemoveCoeffFromMonomial*)


RemoveCoeffFromMonomial[monomial_]:=Module[{term,termvar},
termvar=Exponent[monomial,Variables[monomial]];
term=monomial/(Coefficient[monomial,Variables[monomial]^termvar]);
Return[term]
]


(* ::Subsection::Closed:: *)
(*MassHierarchyFromLeadingOrders*)


MassHierarchyFromLeadingOrders[YukLeadPowAssoc_]:=Module[{YukdLead,YukuLead,\[Lambda],Yukd\[Chi]Pol,Yuku\[Chi]Pol,YukDown\[Chi]Diag,YukDown\[Chi]Lin,YukDown\[Chi]Const,YukDown\[Chi]DiagLeadTerm,YukDown\[Chi]LinLeadTerm,YukDown\[Chi]ConstLeadTerm,YukDownEigens,
YukUp\[Chi]Diag,YukUp\[Chi]Lin,YukUp\[Chi]Const,YukUp\[Chi]DiagLeadTerm,YukUp\[Chi]LinLeadTerm,YukUp\[Chi]ConstLeadTerm,YukUpEigens,
mUpCharmRatio,mCharmTopRatio,mDownStrangeRatio,mStrangeBottomRatio,mBottomTopRatio,
AssocToReturn,PosTopEntry,YukUp\[Chi]DiagLeadTermPowers,YukUp\[Chi]DiagPowers,YukdLeadYYt,YukdLeadList,YukdLeadMon,
YukDownEigensPowers,YukUpEigensPowers,YukDown\[Chi]DiagLeadTermPowers,YukDown\[Chi]DiagPowers,PosBottomEntry
},
(* first calculate the characteristic polynomial *)
YukdLeadYYt=YukLeadPowAssoc[["YukdLead"]] . Transpose[YukLeadPowAssoc[["YukdLead"]]];
YukdLeadList=Flatten[Map[LowestOrderTerms[#,{("VEVsScale"/.Options[stringFNmodels])}]&,Flatten[YukdLeadYYt]]];
YukdLeadMon=Map[RemoveCoeffFromMonomial[#]&,YukdLeadList]/.{{}->{1}};
YukdLead=Partition[Flatten[YukdLeadMon],3];
YukDown\[Chi]Diag={YukdLead[[1]][[1]],YukdLead[[2]][[2]],YukdLead[[3]][[3]]};
YukDown\[Chi]Lin={(YukdLead[[1]][[3]]*YukdLead[[3]][[1]]),(YukdLead[[1]][[2]]*YukdLead[[2]][[1]]),(YukdLead[[2]][[3]]*YukdLead[[3]][[2]])};
YukDown\[Chi]Const={(YukdLead[[1]][[1]]*YukdLead[[2]][[2]]*YukdLead[[3]][[3]]),(YukdLead[[1]][[2]]*YukdLead[[2]][[3]]*YukdLead[[3]][[1]]),(YukdLead[[1]][[3]]*YukdLead[[3]][[2]]*YukdLead[[2]][[1]]),(YukdLead[[1]][[3]]*YukdLead[[3]][[1]]*YukdLead[[2]][[2]]),(YukdLead[[1]][[2]]*YukdLead[[2]][[1]]*YukdLead[[3]][[3]]),(YukdLead[[2]][[3]]*YukdLead[[3]][[2]]*YukdLead[[1]][[1]])};
YukDown\[Chi]DiagLeadTerm=LowestOrderOfMonomialListOneVar[DeleteElements[YukDown\[Chi]Diag,{0}],("VEVsScale"/.Options[stringFNmodels])];
YukDown\[Chi]LinLeadTerm=LowestOrderOfMonomialListOneVar[DeleteElements[YukDown\[Chi]Lin,{0}],("VEVsScale"/.Options[stringFNmodels])];
YukDown\[Chi]ConstLeadTerm=LowestOrderOfMonomialListOneVar[DeleteElements[YukDown\[Chi]Const,{0}],("VEVsScale"/.Options[stringFNmodels])];
YukDownEigens={Sqrt[YukDown\[Chi]DiagLeadTerm],Sqrt[(YukDown\[Chi]LinLeadTerm/YukDown\[Chi]DiagLeadTerm)],Sqrt[(YukDown\[Chi]ConstLeadTerm/YukDown\[Chi]LinLeadTerm)]};
YukDownEigensPowers=Map[Exponent[#,("VEVsScale"/.Options[stringFNmodels])]&,YukDownEigens];
(* repeat from Up *)
YukuLead=YukLeadPowAssoc[["YukuLead"]];
YukUp\[Chi]Diag={YukuLead[[1]][[1]],YukuLead[[2]][[2]],YukuLead[[3]][[3]]};
YukUp\[Chi]Lin={(YukuLead[[1]][[3]]*YukuLead[[3]][[1]]),(YukuLead[[1]][[2]]*YukuLead[[2]][[1]]),(YukuLead[[2]][[3]]*YukuLead[[3]][[2]])};
YukUp\[Chi]Const={(YukuLead[[1]][[1]]*YukuLead[[2]][[2]]*YukuLead[[3]][[3]]),(YukuLead[[1]][[2]]*YukuLead[[2]][[3]]*YukuLead[[3]][[1]]),(YukuLead[[1]][[3]]*YukuLead[[3]][[2]]*YukuLead[[2]][[1]]),(YukuLead[[1]][[3]]*YukuLead[[3]][[1]]*YukuLead[[2]][[2]]),(YukuLead[[1]][[2]]*YukuLead[[2]][[1]]*YukuLead[[3]][[3]]),(YukuLead[[2]][[3]]*YukuLead[[3]][[2]]*YukuLead[[1]][[1]])};
YukUp\[Chi]DiagLeadTerm=LowestOrderOfMonomialListOneVar[DeleteElements[YukUp\[Chi]Diag,{0}],("VEVsScale"/.Options[stringFNmodels])];
YukUp\[Chi]LinLeadTerm=LowestOrderOfMonomialListOneVar[DeleteElements[YukUp\[Chi]Lin,{0}],("VEVsScale"/.Options[stringFNmodels])];
YukUp\[Chi]ConstLeadTerm=LowestOrderOfMonomialListOneVar[DeleteElements[YukUp\[Chi]Const,{0}],("VEVsScale"/.Options[stringFNmodels])];
YukUpEigens={YukUp\[Chi]DiagLeadTerm,(YukUp\[Chi]LinLeadTerm/YukUp\[Chi]DiagLeadTerm),(YukUp\[Chi]ConstLeadTerm/YukUp\[Chi]LinLeadTerm)};
YukUp\[Chi]DiagLeadTermPowers=Exponent[YukUp\[Chi]DiagLeadTerm,("VEVsScale"/.Options[stringFNmodels])];
YukUp\[Chi]DiagPowers=Map[Exponent[#,("VEVsScale"/.Options[stringFNmodels])]&,YukUp\[Chi]Diag];
YukUpEigensPowers=Map[Exponent[#,("VEVsScale"/.Options[stringFNmodels])]&,YukUpEigens];
PosTopEntry=Position[YukUp\[Chi]DiagPowers,YukUp\[Chi]DiagLeadTermPowers][[1]][[1]];
(* Find Pos Down Entry too *)
YukDown\[Chi]DiagLeadTermPowers=Exponent[YukDown\[Chi]DiagLeadTerm,("VEVsScale"/.Options[stringFNmodels])];
YukDown\[Chi]DiagPowers=Map[Exponent[#,("VEVsScale"/.Options[stringFNmodels])]&,YukDown\[Chi]Diag];
PosBottomEntry=Position[YukDown\[Chi]DiagPowers,YukDown\[Chi]DiagLeadTermPowers][[1]][[1]];
(* calculate the scales *)
mUpCharmRatio=YukUpEigensPowers[[3]]-YukUpEigensPowers[[2]];
mCharmTopRatio=YukUpEigensPowers[[2]]-YukUpEigensPowers[[1]];
mDownStrangeRatio=YukDownEigensPowers[[3]]-YukDownEigensPowers[[2]];
mStrangeBottomRatio=YukDownEigensPowers[[2]]-YukDownEigensPowers[[1]];
AssocToReturn=<|"YukUpEigens"->YukUpEigensPowers,"YukDownEigens"->YukDownEigensPowers,"mtPow"->YukUpEigensPowers[[1]],"mbPow"->YukDownEigensPowers[[1]],"mu/mcPow"->mUpCharmRatio,"mc/mtPow"->mCharmTopRatio,"md/msPow"->mDownStrangeRatio,"ms/mbPow"->mStrangeBottomRatio,"YukUpMatLead"->YukLeadPowAssoc[["YukuLead"]],"YukDownMatLead"->YukLeadPowAssoc[["YukdLead"]],"PosOTop"->PosTopEntry,"PosOBot"->PosBottomEntry|>;
Return[AssocToReturn];
]


(* ::Subsection::Closed:: *)
(*IdentifyLowestPowerPos*)


IdentifyLowestPowerPos[ListExp_]:=Module[{usedScale,expList,posNum},
If[ListExp===0,Return[1]];
usedScale="VEVsScale"/.Options[stringFNmodels];
expList=Exponent[ListExp,usedScale];
posNum=Position[expList,Min[expList]][[1]][[1]];
Return[posNum]
]


(* ::Subsection::Closed:: *)
(*LocateExtractO1Coeffs*)


LocateExtractO1Coeffs[WOpsOut_,VEVPowersAssoc_,MHAssoc_,O1coefflist_]:=Module[{diagBnum,diagTnum,usedScale,\[Phi]Fields,\[CapitalPhi]Fields,\[Phi]VEVsList,\[CapitalPhi]VEVsList,
\[Phi]VEVsReplaceRule,\[CapitalPhi]VEVsReplaceRule,scalesList,
dimsList,startPosListC,startPosList,leadPosListC,O1PosList,O1LeadingList,
O1UpMat,O1DownMat,assocReturn,
ReplaceBotPosRule,ReplaceTopPosRule,OTopPosInVec,OBotPosInVec,
numSing,numKahl},

(* create global variables *)
numSing="nSinglet"/.Options[stringFNmodels];
numKahl="nNonperts"/.Options[stringFNmodels];

(* create the ScalesList *)
usedScale="VEVsScale"/.Options[stringFNmodels];
\[Phi]Fields=Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,numSing}];
\[CapitalPhi]Fields=Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,numKahl}];
\[Phi]VEVsList=Map[usedScale^#&,VEVPowersAssoc[["\[Phi]VEVPowers"]]];
\[CapitalPhi]VEVsList=Map[usedScale^#&,VEVPowersAssoc[["\[CapitalPhi]VEVPowers"]]];
\[Phi]VEVsReplaceRule=Table[(\[Phi]Fields[[i]]->\[Phi]VEVsList[[i]]),{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]VEVsReplaceRule=Table[(\[CapitalPhi]Fields[[i]]->\[CapitalPhi]VEVsList[[i]]),{i,1,"nNonperts"/.Options[stringFNmodels]}];
scalesList=WOpsOut/.\[Phi]VEVsReplaceRule/.\[CapitalPhi]VEVsReplaceRule;

(* First we extract the dimensions. *)
dimsList=Map[Length[#]&,(scalesList/.{0->{0}})];

(* Extract starting positions. *)
startPosListC=Delete[Join[{0},Accumulate[dimsList]],-1];
startPosList=startPosListC+Table[1,18];

(* Extract leading power positions. *)
leadPosListC=Map[IdentifyLowestPowerPos[#]&,scalesList]-Table[1,18];
O1PosList=leadPosListC+startPosList;

(* Extract Leading O1 Coeff Lists *)
O1LeadingList=Map[O1coefflist[[#]]&,O1PosList];

O1UpMat=Partition[Take[O1LeadingList,-9],3];
O1DownMat=Partition[Take[O1LeadingList,9],3];

(* Now find the relevant O1 Bottom positions and O1 Top positions. *)
ReplaceBotPosRule={1->1,2->5,3->9};
ReplaceTopPosRule={1->10,2->14,3->18};

(* Extract the relevant positions. *)
diagTnum=MHAssoc[["PosOTop"]]/.ReplaceTopPosRule;
diagBnum=MHAssoc[["PosOBot"]]/.ReplaceBotPosRule;

OTopPosInVec=O1PosList[[diagTnum]];
OBotPosInVec=O1PosList[[diagBnum]];


assocReturn=<|"O1StartPosList"->startPosList,"O1LeadPosList"->O1PosList,"O1LeadList"->O1LeadingList,
"O1LeadUpMat"->O1UpMat,"O1LeadDownMat"->O1DownMat,"OTopVecPos"->OTopPosInVec,"OBotVecPos"->OBotPosInVec|>;
Return[assocReturn]

]


(* ::Subsection::Closed:: *)
(*CalcScaleAndHiggsO1Coeff*)


CalcScaleAndHiggsO1Coeff[massHierarchyAssoc_]:=Module[{ScaleSym,\[Epsilon]1,\[Epsilon]2,\[Epsilon]3,\[Epsilon]4,\[Epsilon]5,ListofScales,Best\[Epsilon],BestO1Top,mTopScale,mBottomScale,
mTScal,mBScal,HiggsTopSq,HiggsBottomSq,AllAssoc,ListofScalesToTest,
NumDiag,NumDiagNumYuk,WOpsListBot,Oplist,dimOplist,O1coeffneed,newO1count,raggedOplist,a,j,newO1coeffneed,
\[Phi]Fields,\[CapitalPhi]Fields,\[Phi]VEVsList,\[CapitalPhi]VEVsList,\[Phi]VEVsReplaceRule,\[CapitalPhi]VEVsReplaceRule,ScalesList,NumDiagNumEntry,ScalesPowerList,O1CoeffBottom,BestO1TopSq
},
ScaleSym="VEVsScale"/.Options[stringFNmodels];
\[Epsilon]1=If[massHierarchyAssoc[["mu/mcPow"]]===0,3,
(("mu"/.Options[stringFNmodels])[[1]]/("mu"/.Options[stringFNmodels])[[2]])^(1/massHierarchyAssoc[["mu/mcPow"]])];
\[Epsilon]2=If[massHierarchyAssoc[["mc/mtPow"]]===0,3,
(("mu"/.Options[stringFNmodels])[[2]]/("mu"/.Options[stringFNmodels])[[3]])^(1/massHierarchyAssoc[["mc/mtPow"]])];
\[Epsilon]3=If[massHierarchyAssoc[["md/msPow"]]===0,3,
(("md"/.Options[stringFNmodels])[[1]]/("md"/.Options[stringFNmodels])[[2]])^(1/massHierarchyAssoc[["md/msPow"]])];
\[Epsilon]4=If[massHierarchyAssoc[["ms/mbPow"]]===0,3,
(("md"/.Options[stringFNmodels])[[2]]/("md"/.Options[stringFNmodels])[[3]])^(1/massHierarchyAssoc[["ms/mbPow"]])];
ListofScalesToTest={\[Epsilon]1,\[Epsilon]2,\[Epsilon]3,\[Epsilon]4};
ListofScales=Map[If[#>3,3,#]&,ListofScalesToTest];
Best\[Epsilon]=GeometricMean[ListofScales];
mTopScale=Best\[Epsilon]^(massHierarchyAssoc[["mtPow"]]);
mBottomScale=Best\[Epsilon]^(massHierarchyAssoc[["mbPow"]]);

O1CoeffBottom=1;

HiggsTopSq=(("mu"/.Options[stringFNmodels])[[3]])^2/((mTopScale)^2);
HiggsBottomSq=(("md"/.Options[stringFNmodels])[[3]])^2/((mBottomScale*O1CoeffBottom)^2);
BestO1TopSq=HiggsTopSq/(("v"/.Options[stringFNmodels])^2-HiggsBottomSq);
BestO1Top=Sqrt[BestO1TopSq];
AllAssoc=Append[massHierarchyAssoc,<|"ListOfScales"->ListofScales,"BestScale"->Best\[Epsilon],"HiggsSquared"->{HiggsTopSq,HiggsBottomSq},"BestO1Top"->BestO1Top,"O1t^2"->BestO1TopSq|>];
Return[AllAssoc]
]


(* ::Subsection::Closed:: *)
(*GeometricStandardDeviation*)


GeometricStandardDeviation[listOfValues_]:=Module[{num,geometricMean,GSD},
num=Length[listOfValues];
geometricMean=GeometricMean[listOfValues];
GSD=Exp[Sqrt[Total[Map[(Log[#/geometricMean])^2&,listOfValues]]/num]];
Return[N[GSD]]
]


(* ::Subsection::Closed:: *)
(*FitnessTexture*)


FitnessTexture[MassHierarchyAssoc_,VEVsPowerListAssoc_]:=Module[{VarianceOfScale,AssocToReturn,FitnessTextureTotal,ListMassRatios,ListPowers,FitnessPowers,extraFactors,O1TopRange,O1CoeffTop,
BestCalcScale,OrderOneCoeffPenal,ScalePenal,\[Phi]VEVPows,\[CapitalPhi]VEVPows,listOf\[Phi]VEVs,listOf\[CapitalPhi]VEVs,
fixedMaxSing,fixedMinSing,fixedMaxKahl,fixedMinKahl,max\[Phi]VEVCont,min\[Phi]VEVCont,max\[CapitalPhi]VEVCont,min\[CapitalPhi]VEVCont,
HiggsSMVEV,otSq,HiggsUpSq,HiggsDownSq,
numSing,numKahl},

(* external setting options for later use. *)
numSing="nSinglet"/.Options[stringFNmodels];
numKahl="nNonperts"/.Options[stringFNmodels];

extraFactors="fitTextFac"/.Options[stringFNmodels];

(* Variance contributions. *)
VarianceOfScale=GeometricStandardDeviation[MassHierarchyAssoc[["ListOfScales"]]];

(* Power contributions. *)
ListPowers={MassHierarchyAssoc[["mu/mcPow"]],MassHierarchyAssoc[["mc/mtPow"]],MassHierarchyAssoc[["md/msPow"]],MassHierarchyAssoc[["ms/mbPow"]]};
FitnessPowers=Total[Map[If[#>=0,0,Abs[#]]&,ListPowers]];

(* O1Top contributions. *)
HiggsSMVEV="v"/.Options[stringFNmodels];
{HiggsUpSq,HiggsDownSq}=MassHierarchyAssoc[["HiggsSquared"]];
otSq=HiggsUpSq/(HiggsSMVEV^2-HiggsDownSq);

O1TopRange="O1CoeffRange"/.Options[stringFNmodels];
O1CoeffTop=Abs[MassHierarchyAssoc[["BestO1Top"]]];
BestCalcScale=MassHierarchyAssoc[["BestScale"]];
If[otSq>0,
	(* Higgs can be solved. *) 
	If[O1CoeffTop<O1TopRange[[1]]||O1CoeffTop>O1TopRange[[2]],
	If[O1CoeffTop<O1TopRange[[1]],OrderOneCoeffPenal=Abs[Log10[O1CoeffTop/O1TopRange[[1]]]]];
	If[O1CoeffTop>O1TopRange[[2]],OrderOneCoeffPenal=Abs[Log10[O1CoeffTop/O1TopRange[[2]]]]];
	,
	OrderOneCoeffPenal=0;
	];
	,
	(* Higgs can be solved. *) 
	OrderOneCoeffPenal=Abs[otSq];
];

(* scale penalities *)
(* added in stuff so it works for different types of moduli being set to zero *)
If[numSing>0,
\[Phi]VEVPows=VEVsPowerListAssoc[["\[Phi]VEVPowers"]];
listOf\[Phi]VEVs=Map[BestCalcScale^#&,\[Phi]VEVPows];
fixedMaxSing=0.6;
fixedMinSing=0.005;
max\[Phi]VEVCont=If[Max[listOf\[Phi]VEVs]>=fixedMaxSing,Abs[Log[(Max[listOf\[Phi]VEVs]/fixedMaxSing)]],0];
min\[Phi]VEVCont=If[Min[listOf\[Phi]VEVs]<=fixedMinSing,Abs[Log[(Min[listOf\[Phi]VEVs]/fixedMinSing)]],0];
,
max\[Phi]VEVCont=0;
min\[Phi]VEVCont=0;
];
If[numKahl>0,
\[CapitalPhi]VEVPows=VEVsPowerListAssoc[["\[CapitalPhi]VEVPowers"]];
listOf\[CapitalPhi]VEVs=Map[BestCalcScale^#&,\[CapitalPhi]VEVPows];
fixedMaxKahl=0.1;
fixedMinKahl=0.0001;
max\[CapitalPhi]VEVCont=If[Max[listOf\[CapitalPhi]VEVs]>=fixedMaxKahl,Abs[Log[(Max[listOf\[CapitalPhi]VEVs]/fixedMaxKahl)]],0];
min\[CapitalPhi]VEVCont=If[Min[listOf\[CapitalPhi]VEVs]<=fixedMinKahl,Abs[Log[(Min[listOf\[CapitalPhi]VEVs]/fixedMinKahl)]],0];
,
max\[CapitalPhi]VEVCont=0;
min\[CapitalPhi]VEVCont=0;
];

ScalePenal=max\[Phi]VEVCont+min\[Phi]VEVCont+max\[CapitalPhi]VEVCont+min\[CapitalPhi]VEVCont;

FitnessTextureTotal=FitnessPowers*extraFactors[[1]]+VarianceOfScale*extraFactors[[2]]+ScalePenal*extraFactors[[3]]+OrderOneCoeffPenal*extraFactors[[4]];
AssocToReturn=<|"FitnessTextSec"->FitnessTextureTotal,"FitnessPowers"->FitnessPowers,"VarsScale"->VarianceOfScale,"ScalePenalty"->ScalePenal,"OTopPenalty"->OrderOneCoeffPenal|>;
Return[AssocToReturn]
]


(* ::Section::Closed:: *)
(*Auxiliary Modules - Operators to Matrices*)


(* ::Subsection::Closed:: *)
(*AdjoinO1coeffWOpsTextQSec*)


AdjoinO1coeffWOpsTextQSec[WOplistoutput_,MHAssoc_,VEVPowersListAssoc_,O1Assoc_,O1vec_]:=Module[
{Oplist,dimOplist,O1coeffneed,raggedOplist,a,j,newO1coeffneed,newO1Vec,Wexpansion,
usedScale,\[Phi]Fields,\[CapitalPhi]Fields,\[Phi]VEVsList,\[CapitalPhi]VEVsList,\[Phi]VEVsReplaceRule,\[CapitalPhi]VEVsReplaceRule,ReplaceScaleRule,
\[Phi]VEVsRepRule,\[CapitalPhi]VEVsRepRules,\[CapitalPhi]VEVsRepRule,
WIJaOps,WExpVals,YdMatVals,YuMatVals,outAssoc,YdMatOps,YuMatOps,oldO1LeadList,
posToReplaceO1T,posToReplaceO1B,newO1LeadList,newO1LeadUp,newO1LeadDown,ytMZVal,ybMZVal},

usedScale="VEVsScale"/.Options[stringFNmodels];
\[Phi]Fields=Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]Fields=Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}];
\[Phi]VEVsList=Map[usedScale^#&,VEVPowersListAssoc[["\[Phi]VEVPowers"]]];
\[CapitalPhi]VEVsList=Map[usedScale^#&,VEVPowersListAssoc[["\[CapitalPhi]VEVPowers"]]];
\[Phi]VEVsReplaceRule=Table[(\[Phi]Fields[[i]]->\[Phi]VEVsList[[i]]),{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]VEVsReplaceRule=Table[(\[CapitalPhi]Fields[[i]]->\[CapitalPhi]VEVsList[[i]]),{i,1,"nNonperts"/.Options[stringFNmodels]}];

(* line to update the newO1coeff *)
newO1Vec=O1vec;
newO1Vec[[O1Assoc[["OTopVecPos"]]]]=Abs[MHAssoc[["BestO1Top"]]];
newO1Vec[[O1Assoc[["OBotVecPos"]]]]=1.0;

(* Now append O(1) coefficients. *)
Oplist=Map[Flatten[{#}]&,WOplistoutput];
dimOplist=Total[Map[Dimensions,Oplist,1]][[1]];
O1coeffneed=Take[newO1Vec,{1,dimOplist}];

(* Write O1coefflist as the required ragged array. *)
raggedOplist=Oplist/.Table[("BundleModuliFields"/.Options[stringFNmodels])[i]->a,{i,1,"nSinglet"/.Options[stringFNmodels],1}]/.Table[("KahlerEffFields"/.Options[stringFNmodels])[i]->a,{i,1,"nNonperts"/.Options[stringFNmodels],1}];
raggedOplist=raggedOplist/.a->1;
j=1;
newO1coeffneed=Map[O1coeffneed[[j++]]&,raggedOplist,{-1}];

Wexpansion=newO1coeffneed*Oplist;
Wexpansion=Map[Flatten[{Total[#]}]&,Wexpansion,1];

(* Now put in the VEVs of the fields. *)
ReplaceScaleRule={usedScale->MHAssoc[["BestScale"]]};
\[Phi]VEVsRepRule=\[Phi]VEVsReplaceRule/.ReplaceScaleRule;
\[CapitalPhi]VEVsRepRule=\[CapitalPhi]VEVsReplaceRule/.ReplaceScaleRule;
WIJaOps=Map[Normal[#]&,Wexpansion];
WExpVals=Map[Normal[#]&,Wexpansion]/.\[Phi]VEVsRepRule/.\[CapitalPhi]VEVsRepRule;
YdMatVals=Partition[Take[(Flatten[WExpVals]),9],3];
YuMatVals=Partition[Take[(Flatten[WExpVals]),-9],3];

YdMatOps=Partition[Take[(Flatten[WIJaOps]),9],3];
YuMatOps=Partition[Take[(Flatten[WIJaOps]),-9],3];

outAssoc=<|"YdMatVal"->YdMatVals,"YuMatVal"->YuMatVals,"YdMatOps"->YdMatOps,"YuMatOps"->YuMatOps|>;

(* Update O1 Assoc with the new O1 coefficients. *)
oldO1LeadList=O1Assoc[["O1LeadList"]];
posToReplaceO1T=MHAssoc[["PosOTop"]]/.{1->10,2->14,3->18};
posToReplaceO1B=MHAssoc[["PosOBot"]]/.{1->1,2->5,3->9};
newO1LeadList=oldO1LeadList;
newO1LeadList[[posToReplaceO1T]]=Abs[MHAssoc[["BestO1Top"]]];
newO1LeadList[[posToReplaceO1B]]=1.0;

newO1LeadUp=Partition[Take[newO1LeadList,-9],3];
newO1LeadDown=Partition[Take[newO1LeadList,9],3];

ytMZVal=YuMatVals[[MHAssoc[["PosOTop"]]]][[MHAssoc[["PosOTop"]]]];
ybMZVal=YdMatVals[[MHAssoc[["PosOBot"]]]][[MHAssoc[["PosOBot"]]]];

outAssoc=<|"YdMatVal"->YdMatVals,"YuMatVal"->YuMatVals,"YdMatOps"->YdMatOps,"YuMatOps"->YuMatOps,"ytMZVal"->ytMZVal,"ybMZVal"->ybMZVal|>;
outAssoc=Join[outAssoc,O1Assoc];
AssociateTo[outAssoc,{"O1LeadList"->newO1LeadList,"O1LeadUpMat"->newO1LeadUp,"O1LeadDownMat"->newO1LeadDown,"PosUpLead"->MHAssoc[["PosOTop"]],"PosDownLead"->MHAssoc[["PosOBot"]]}];

Return[outAssoc]
]


(* ::Section::Closed:: *)
(*Auxiliary Modules - RG Evolution*)


(* ::Subsection::Closed:: *)
(*RGGivenYukMZFindYukMGUT*)


RGGivenYukMZFindYukMGUT[ytMZ_,ybMZ_]:=Module[{g3OneLoop,g3,g30,b3,t,
ytEQN,yt,yb,ybEQN,SaveSoln,ytMGUT,ybMGUT
},
g3OneLoop=g3->(Sqrt[g30^2/(1-(2*g30^2*b3*t))]/.{b3->-3,g30->Sqrt[4*\[Pi]*0.1155]});
ytEQN=yt'[t]==yt[t]*(6*yt[t]^2+yb[t]^2-(16/3)*g3^2)/.g3OneLoop;
ybEQN=yb'[t]==yb[t]*(6*yb[t]^2+yt[t]^2-(16/3)*g3^2)/.g3OneLoop;
SaveSoln=NDSolve[{ytEQN,ybEQN,yt[0.0]==ytMZ,yb[0.0]==ybMZ},{yt[t],yb[t]},{t,0.0,0.201}];
ytMGUT=(yt[t]/.SaveSoln)/.t->0.201;
ybMGUT=(yb[t]/.SaveSoln)/.t->0.201;
Return[Flatten[{ytMGUT,ybMGUT}]]
]


(* ::Subsection::Closed:: *)
(*RGGivenYukMGUTFindYukMZ*)


RGGivenYukMGUTFindYukMZ[ytMGUT_,ybMGUT_]:=Module[{g3OneLoopSq,g3Sq,g30Sq,b3,t,
ytEQN,yt,yb,ybEQN,SaveSoln,ytMZ,ybMZ},
g3OneLoopSq=g3Sq->(g30Sq/(1-(2*g30Sq*b3*t))/.{b3->-3,g30Sq->4*\[Pi]*0.1155});
ytEQN=yt'[t]==yt[t]*(6*yt[t]^2+yb[t]^2-(16/3)*g3Sq)/.g3OneLoopSq;
ybEQN=yb'[t]==yb[t]*(6*yb[t]^2+yt[t]^2-(16/3)*g3Sq)/.g3OneLoopSq;
SaveSoln=NDSolve[{ytEQN,ybEQN,yt[0.201]==ytMGUT,yb[0.201]==ybMGUT},{yt[t],yb[t]},{t,0.0,0.201}];
ytMZ=(yt[t]/.SaveSoln)/.t->0.0;
ybMZ=(yb[t]/.SaveSoln)/.t->0.0;
Return[Flatten[{ytMZ,ybMZ}]]
]


(* ::Subsection::Closed:: *)
(*RGGivenYukMZFindYukRatiosAll*)


RGGivenYukMZFindYukRatiosAll[ytMZ_,ybMZ_]:=Module[{SaveSoln,ytMGUT,ybMGUT,
yu3iRatio,yuijRatio,yd3iRatio,ydijRatio,g3OneLoop,g3,g30,b3,t,
ytEQN,yt,yb,ybEQN,yu3iEQN,yuijEQN,yd3iEQN,ydijEQN,returnAssoc},
g3OneLoop=g3->(Sqrt[g30^2/(1-(2*g30^2*b3*t))]/.{b3->-3,g30->Sqrt[4*\[Pi]*0.1155]});
ytEQN=yt'[t]==yt[t]*(6*yt[t]^2+yb[t]^2-(16/3)*g3^2)/.g3OneLoop;
ybEQN=yb'[t]==yb[t]*(6*yb[t]^2+yt[t]^2-(16/3)*g3^2)/.g3OneLoop;

yu3iEQN=6yt[t]^2+yb[t]^2-(16/3)*(g3/.g3OneLoop)^2;
yuijEQN=3yt[t]^2-(16/3)(g3/.g3OneLoop)^2;
yd3iEQN=6yb[t]^2+yt[t]^2-(16/3)*(g3/.g3OneLoop)^2;
ydijEQN=3yb[t]^2-(16/3)(g3/.g3OneLoop)^2;

SaveSoln=NDSolve[{ytEQN,ybEQN,yt[0.0]==ytMZ,yb[0.0]==ybMZ},{yt[t],yb[t]},{t,0.0,0.201}];
ytMGUT=(yt[t]/.SaveSoln[[1]])/.t->0.201;
ybMGUT=(yb[t]/.SaveSoln[[1]])/.t->0.201;
yu3iRatio=Exp[NIntegrate[yu3iEQN/.SaveSoln[[1]],{t,0,0.201}]];
yuijRatio=Exp[NIntegrate[yuijEQN/.SaveSoln[[1]],{t,0,0.201}]];
yd3iRatio=Exp[NIntegrate[yd3iEQN/.SaveSoln[[1]],{t,0,0.201}]];
ydijRatio=Exp[NIntegrate[ydijEQN/.SaveSoln[[1]],{t,0,0.201}]];
returnAssoc=<|"ytMZ"->ytMZ,"ybMZ"->ybMZ,"yu33R"->ytMGUT/ytMZ,"yu3iR"->yu3iRatio,"yuijR"->yuijRatio,"yd33R"->ybMGUT/ybMZ,"yd3iR"->yd3iRatio,"ydijR"->ydijRatio|>;
Return[returnAssoc]
]


(* ::Subsection:: *)
(*RGGivenYukGUTFindYukRatiosAll*)


RGGivenYukGUTFindYukRatiosAll[ytMGUT_,ybMGUT_]:=Module[{SaveSoln,ytMZ,ybMZ,
yu3iRatio,yuijRatio,yd3iRatio,ydijRatio,g3OneLoop,g3,g30,b3,t,
ytEQN,yt,yb,ybEQN,yu3iEQN,yuijEQN,yd3iEQN,ydijEQN,returnAssoc},
g3OneLoop=g3->(Sqrt[g30^2/(1-(2*g30^2*b3*t))]/.{b3->-3,g30->Sqrt[4*\[Pi]*0.1155]});
ytEQN=yt'[t]==yt[t]*(6*yt[t]^2+yb[t]^2-(16/3)*g3^2)/.g3OneLoop;
ybEQN=yb'[t]==yb[t]*(7*yb[t]^2+yt[t]^2-(16/3)*g3^2)/.g3OneLoop;

yu3iEQN=6yt[t]^2+yb[t]^2-(16/3)*(g3/.g3OneLoop)^2;
yuijEQN=3yt[t]^2-(16/3)(g3/.g3OneLoop)^2;
yd3iEQN=7yb[t]^2+yt[t]^2-(16/3)*(g3/.g3OneLoop)^2;
ydijEQN=4yb[t]^2-(16/3)(g3/.g3OneLoop)^2;

SaveSoln=NDSolve[{ytEQN,ybEQN,yt[0.201]==ytMGUT,yb[0.201]==ybMGUT},{yt[t],yb[t]},{t,0.0,0.201}];
ytMZ=(yt[t]/.SaveSoln[[1]])/.t->0.0;
ybMZ=(yb[t]/.SaveSoln[[1]])/.t->0.0;
yu3iRatio=Exp[NIntegrate[-yu3iEQN/.SaveSoln[[1]],{t,0,0.201}]];
yuijRatio=Exp[NIntegrate[-yuijEQN/.SaveSoln[[1]],{t,0,0.201}]];
yd3iRatio=Exp[NIntegrate[-yd3iEQN/.SaveSoln[[1]],{t,0,0.201}]];
ydijRatio=Exp[NIntegrate[-ydijEQN/.SaveSoln[[1]],{t,0,0.201}]];
returnAssoc=<|"ytMGUT"->ytMGUT,"ybMGUT"->ybMGUT,"ytMZ"->ytMZ,"ybMZ"->ybMZ,
"yu33R"->ytMZ/ytMGUT,"yu3iR"->yu3iRatio,"yuijR"->yuijRatio,
"yd33R"->ybMZ/ybMGUT,"yd3iR"->yd3iRatio,"ydijR"->ydijRatio|>;
Return[returnAssoc]
]


(* ::Subsection:: *)
(*RGGivenYukMZFindYukRatiosAllFixedPlane*)


RGGivenYukMZFindYukRatiosAllFixedPlane[ytMZ_,ybMZ_]:=Module[{yu33Factor,yu3iFactor,yuijFactor,
yd33Factor,yd3iFactor,ydijFactor,returnAssoc},
(* fixed plane *)
yu33Factor=-94.4478+96.4698*ytMZ+2.99524*ybMZ;
yu3iFactor=-94.4476+96.4697*ytMZ+2.99523*ybMZ;
yuijFactor=-17.4073+18.3248*ytMZ+0.496986*ybMZ;
yd33Factor=-2.88346+3.40746*ytMZ+0.242808*ybMZ;
yd3iFactor=-2.88346+3.40746*ytMZ+0.242808*ybMZ;
ydijFactor=0.396947+0.00647511*ytMZ+0.0550254*ybMZ;
returnAssoc=<|"ytMZ"->ytMZ,"ybMZ"->ybMZ,"yu33R"->yu33Factor,"yu3iR"->yu3iFactor,"yuijR"->yuijFactor,"yd33R"->yd33Factor,"yd3iR"->yd3iFactor,"ydijR"->ydijFactor|>;
Return[returnAssoc]
]


(* ::Subsection::Closed:: *)
(*ConstructRGFactors*)


ConstructRGFactors[y33Factor_,y3iFactor_,yijFactor_,leadPos_]:=Module[{outputMat},
If[leadPos==1||leadPos==2||leadPos==3,
	If[leadPos==1,
		outputMat={y33Factor,y3iFactor,y3iFactor,y3iFactor,yijFactor,yijFactor,y3iFactor,yijFactor,yijFactor};];
	If[leadPos==2,
		outputMat={yijFactor,y3iFactor,yijFactor,y3iFactor,y33Factor,y3iFactor,yijFactor,y3iFactor,yijFactor};];
	If[leadPos==3,
		outputMat={yijFactor,yijFactor,y3iFactor,yijFactor,yijFactor,y3iFactor,y3iFactor,y3iFactor,y33Factor};];
	,
	Print["leadPos is out of range. Check."];
	Return[0];	
];
Return[outputMat]		
]



(* ::Subsection::Closed:: *)
(*GivenYukMatsFindRGMatricesPlane*)


GivenYukMatsFindRGMatricesPlane[yukMZAssoc_]:=Module[{leadPosOfUp,leadPosOfDown,
yu33Ratio,yu3iRatio,yuijRatio,yd33Ratio,yd3iRatio,ydijRatio,
RGUpMatPlane,RGDownMatPlane,outAssoc,
posToExtractUp,posToExtractDown,yukTopMZ,yukBotMZ,yukRatiosAssoc},
leadPosOfUp=yukMZAssoc[["PosUpLead"]];
leadPosOfDown=yukMZAssoc[["PosDownLead"]];
yukTopMZ=yukMZAssoc[["ytMZVal"]];
yukBotMZ=yukMZAssoc[["ybMZVal"]];
(*Print[yukTopMZ,yukBotMZ];*)
(* Compute the ratios *)
yukRatiosAssoc=RGGivenYukMZFindYukRatiosAllFixedPlane[yukTopMZ,yukBotMZ];
yu33Ratio=yukRatiosAssoc[["yu33R"]];
yu3iRatio=yukRatiosAssoc[["yu3iR"]];
yuijRatio=yukRatiosAssoc[["yuijR"]];
yd33Ratio=yukRatiosAssoc[["yd33R"]];
yd3iRatio=yukRatiosAssoc[["yd3iR"]];
ydijRatio=yukRatiosAssoc[["ydijR"]];
RGUpMatPlane=ConstructRGFactors[yu33Ratio,yu3iRatio,yuijRatio,leadPosOfUp];
RGDownMatPlane=ConstructRGFactors[yd33Ratio,yd3iRatio,ydijRatio,leadPosOfDown];
outAssoc=<|"RGUpMat"->Partition[RGUpMatPlane,3],"RGDownMat"->Partition[RGDownMatPlane,3]|>;
(*Print[outAssoc];*)
Return[outAssoc]
]


(* ::Subsection::Closed:: *)
(*O1LogDevPenal*)


O1LogDevPenal[O1numIn_]:=Module[{O1RangeMin,O1RangeMax,O1Pen,tol},
{O1RangeMin,O1RangeMax}="O1CoeffRange"/.Options[stringFNmodels];
tol=10^(-10);
If[Abs[O1numIn]>=O1RangeMax||Abs[O1numIn]<=O1RangeMin,
	If[Abs[O1numIn]>=O1RangeMax,O1Pen=Abs[Log10[O1RangeMax/Abs[O1numIn]]]];
	If[Abs[O1numIn]<=O1RangeMin&&Abs[O1numIn]>=tol,O1Pen=Abs[Log10[O1RangeMin/Abs[O1numIn]]]];
	If[Abs[O1numIn]<tol,O1Pen=Abs[Log10[O1RangeMin/tol]]];
,
	O1Pen=0;
];
(*Print["Input number : "<>ToString[O1numIn]<>" -> Penalty: "<>ToString[Abs[O1Pen]]];*)
Return[Abs[O1Pen]]
]


(* ::Subsection::Closed:: *)
(*FitnessRG*)


FitnessRG[YukO1Assoc_]:=Module[{leadO1Assoc,RGMatrixAssoc,
O1UpMat,O1DownMat,O1UpRatioToDiv,O1DownRatioToDiv,O1UpGUTMat,O1DownGUTMat,O1GUTList,RGFitCont},

RGMatrixAssoc=GivenYukMatsFindRGMatricesPlane[YukO1Assoc];
O1UpMat=YukO1Assoc[["O1LeadUpMat"]];
O1DownMat=YukO1Assoc[["O1LeadDownMat"]];
O1UpRatioToDiv=RGMatrixAssoc[["RGUpMat"]];
O1DownRatioToDiv=RGMatrixAssoc[["RGDownMat"]];
O1UpGUTMat=O1UpMat*O1UpRatioToDiv;
O1DownGUTMat=O1DownMat*O1DownRatioToDiv;
(*Print[O1DownGUTMat];
Print[O1UpGUTMat];*)
O1GUTList=Join[Flatten[O1UpGUTMat],Flatten[O1DownGUTMat]];

RGFitCont=Total[Map[O1LogDevPenal,O1GUTList]];
Return[<|"FitnessRGSec"->RGFitCont|>]
]


(* ::Section::Closed:: *)
(*Auxiliary Modules - Computing Mass*)


(* ::Subsection::Closed:: *)
(*ComputeSMQuantities*)


(* This submodule computes all the SM quantities. *)
ComputeSMQuantities[Yukd_,Yuku_]:=Module[{Yddec,Yudec,vuval,vdval,vuvalmax,vdvalmax,HuRule,HdRule,Tan\[Beta]Rule,
transmat,dmass,umass,Uu,Ud,CKMmatval,mdrule,murule,CKMrule,HiggsVEVRule,
maxTan\[Beta],tol,tan\[Beta]Temp},

(*Print["Matrices before decomposition."];
Print["Yukd:"];
Print[Yukd//MatrixForm];
Print["Yuku:"];
Print[Yuku//MatrixForm];*)

(* First SVD both Yukawa matrices. *)
Yddec=SingularValueDecomposition[Yukd];
Yudec=SingularValueDecomposition[Yuku];

(* Find if vu is solvable. *)
vuvalmax=("vumax"/.Options[stringFNmodels]);
If[Yudec[[2]][[1]][[1]]==0,
vuval=vuvalmax,
vuval=(("mu"/.Options[stringFNmodels])[[3]])/(Yudec[[2]][[1]][[1]]);
If[vuval<vuvalmax,vuval,vuval=vuvalmax];
];
HuRule=<|"Hu"->vuval|>;

(* We might as well just find vd using the bottom quark. *)
vdvalmax=("vdmax"/.Options[stringFNmodels]);
If[Yddec[[2]][[1]][[1]]==0,
vdval=vdvalmax,
vdval=(("md"/.Options[stringFNmodels])[[3]])/(Yddec[[2]][[1]][[1]]);
If[vdval<vdvalmax,vdval,vdval=vdvalmax];
];

HdRule=<|"Hd"->vdval|>;

(* Now compute tan\[Beta]. *)
maxTan\[Beta]=1000;
tan\[Beta]Temp=vuval/vdval;
If[tan\[Beta]Temp>maxTan\[Beta],
Tan\[Beta]Rule=<|"tan\[Beta]"->maxTan\[Beta]|>,
Tan\[Beta]Rule=<|"tan\[Beta]"->tan\[Beta]Temp|>
];

HiggsVEVRule=<|"Higgs"->Sqrt[vuval^2+vdval^2]|>;

(* Output of Higgs*)
Join[HdRule,HuRule,HiggsVEVRule,Tan\[Beta]Rule];

(* Now want to canonically order the matrices. *)
(* Note SingularValueDecomposition in Mathematica gives inverse the order we want. *)
transmat={{0,0,1},{0,1,0},{1,0,0}};
dmass=transmat . Yddec[[2]] . Transpose[transmat]*vdval;
umass=transmat . Yudec[[2]] . Transpose[transmat]*vuval;

Ud=Yddec[[1]];
Uu=Yudec[[1]];

(* Compute CKM Matrix. *)
CKMmatval=ConjugateTranspose[transmat] . ConjugateTranspose[Uu] . Ud . transmat;

(* Note that transmat is orthogonal so the Hermitian conjugate *)

(* Take real value *)
(*CKMmatval=Re[CKMmatval];*)

(* Convert to substitution rules. *)
mdrule=<|"md"->Diagonal[dmass]|>;
murule=<|"mu"->Diagonal[umass]|>;

CKMrule=<|"CKMmat"->CKMmatval|>;

(* Return: *)
Join[HdRule,HuRule,HiggsVEVRule,Tan\[Beta]Rule,CKMrule,mdrule,murule]

]



(* ::Section::Closed:: *)
(*Auxiliary Modules - Computing Fitness*)


(* ::Subsection::Closed:: *)
(*NormalLog*)


NormalLog[numerator_,denominator_]:=Module[{result,maxLog},
maxLog="MaxLogPower"/.Options[stringFNmodels];
If[numerator==0,
result=maxLog;
,
result=Abs[Log10[numerator/denominator]];
If[result>=maxLog,result=maxLog];
];
result
]


(* ::Subsection::Closed:: *)
(*FitnessHiggs*)


FitnessHiggs[mass_]:=Module[{HiggsSM,HiggsMeasured,HiggsFactor,HiggsBound,Fitvev,FitH,Higgsmax,fitHextra,fittan\[Beta]},

(* This modules concerns the calculation of the fitness of the Higgs sector. *)
(*
(* 1. Contribution of deviation away from <H>. *)
HiggsSM=("HVEV"/.mass)*("MplinGeV"/.Options[stringFNmodels]);
HiggsBound=("HiggsBound"/.Options[stringFNmodels]);
HiggsFactor=("HiggsFitScaleFactor"/.Options[stringFNmodels]);
HiggsMeasured=("Higgs"/.Options[stringFNmodels]);
(*Higgsmax=Sqrt[("vumax"/.Options[stringFNmodels])^2+("vdmax"/.Options[stringFNmodels])^2]/(("MplinGeV"/.Options[stringFNmodels]));
(* New Change : to a linear scaling as the EW-breaking scale is well-fixed. *)
Fitvev=If[Abs[HiggsSM-HiggsMeasured]<=(("vErr"/.Options[stringFNmodels])/("MplinGeV"/.Options[stringFNmodels])),
0
,
Log10[1+(Abs[HiggsSM-HiggsMeasured]/(("vErr"/.Options[stringFNmodels])/("MplinGeV"/.Options[stringFNmodels])))]
];*)

(*Abs[HiggsSM-HiggsMeasured]*("HiggsFitnessFactor"/.Options[stringFNmodels]);*)
(* Update: Apr 2024 - linear fitness for penalising. *)
Fitvev=If[Abs[HiggsSM-HiggsMeasured]<=HiggsBound,
Abs[(HiggsSM-HiggsMeasured)/HiggsBound],
If[HiggsSM>=HiggsMeasured+HiggsBound,
Abs[HiggsFactor*(HiggsSM-HiggsMeasured-HiggsBound)/HiggsBound]+1,
Abs[HiggsFactor*(HiggsSM-HiggsMeasured+HiggsBound)/HiggsBound]+1
]
];


(* 2. Extra penalty for hitting "vdmax" or "vumax". *)
fitHextra=0;
If[(("vdmax"/.Options[stringFNmodels])==("MplinGeV"/.Options[stringFNmodels])*("Hd"/.mass)),fitHextra=fitHextra+("ExtraPenaltyforMaxHiggs"/.Options[stringFNmodels]),fitHextra=fitHextra+0];
If[(("vumax"/.Options[stringFNmodels])==("MplinGeV"/.Options[stringFNmodels])*("Hu"/.mass)),fitHextra=fitHextra+("ExtraPenaltyforMaxHiggs"/.Options[stringFNmodels]),fitHextra=fitHextra+0];


(* 3. Tangent \[Beta] contributions. *)
If[("Fitnesstan\[Beta]On"/.Options[stringFNmodels])&&(("tan\[Beta]"/.mass)>("tan\[Beta]UpperBound")/.Options[stringFNmodels]),
fittan\[Beta]=NormalLog[("tan\[Beta]"/.mass),(("tan\[Beta]UpperBound")/.Options[stringFNmodels]),1000];
,
fittan\[Beta]=0];
(* Total Fitness: *)
FitH=Fitvev+fitHextra+fittan\[Beta];

*)
FitH=0;
Fitvev=0;
fitHextra=0;
fittan\[Beta]=0;
<|"FitnessHSec"->FitH,"FitnessHiggs"->Fitvev+fitHextra,"FitnessTan\[Beta]"->fittan\[Beta]|>

]


(* ::Subsection::Closed:: *)
(*FitnessQuarksMassMix*)


FitnessQuarksMassMix[mass_]:=Module[{massdexp,massuexp,massdmodel,massumodel,complicatedmassd,complicatedmassu,maxQMassCont,output,output1,output2,
CKMmatcomplicated,output3,CKMcmplx,Jarlskogmeasured,Jarlscont,extrafit,extrafitCKM,outputQmass,outputQmix,fitQtotal,fitCKMtotal},
(* This module concerns the calculation of the fitness of the quark sector. *)

(* 1. Mass Differences *)
massdexp=("md"/.Options[stringFNmodels]);
massuexp=("mu"/.Options[stringFNmodels]);
massdmodel=("md"/.mass);
massumodel=("mu"/.mass);
complicatedmassd=Transpose[Partition[Join[massdmodel,massdexp],3]];
complicatedmassu=Transpose[Partition[Join[massumodel,massuexp],3]];
output1=Map[NormalLog[#[[1]],#[[2]]]&,complicatedmassd]*("fitDownFac"/.Options[stringFNmodels]);
output2=Map[NormalLog[#[[1]],#[[2]]]&,complicatedmassu]*("fitUpFac"/.Options[stringFNmodels]);

(* 2. Also output the indiivial CKM contributions. *)
CKMmatcomplicated=Transpose[Partition[Flatten[Join[Abs[("CKMmat"/.mass)],("CKMmat"/.Options[stringFNmodels])]],9]];
output3=Map[NormalLog[#[[1]],#[[2]]]&,CKMmatcomplicated]*("fitCKMFac"/.Options[stringFNmodels]);
CKMcmplx=("CKMmat"/.mass);


(* Total Fitness: *)
outputQmass=Join[output1,output2];
outputQmix=output3;
fitQtotal=Total[outputQmass];
fitCKMtotal=Total[outputQmix];

<|"FitnessQMassSec"->fitQtotal,"FitnessQMixSec"->fitCKMtotal,"FitnessQmass"->outputQmass,"FitnessCKM"->outputQmix|>
]



(* ::Subsection::Closed:: *)
(*FitnessGSAnomaly*)


FitnessGSAnomaly[Qmatrix_]:=Module[{A5anomaly,A3anomaly,A2anomaly,A1anomaly,fitanomaly,qHsum,qHsumdown,\[Beta]vec,\[CapitalPhi]Chargemat,MatrixEqn,sol,\[CapitalPhi]Chargematlow,
A5anomalylow,\[Beta],realsol,roundingerror,
\[Beta]vars,x,image\[CapitalPhi]mat,distance,regioneqs},
(* This module calculates the fitness contributing from whether the U(1) factors can be GS cancelled. *)

(* Compute JSU(5)^2 anomaly. *)
If[("nNonperts"/.Options[stringFNmodels])!=0,
A5anomaly=3*(Qmatrix[[1]]+Qmatrix[[2]]+Qmatrix[[3]])+Qmatrix[[4]]+Qmatrix[[5]]+Qmatrix[[6]];
\[Beta]vec=Table[\[Beta][i],{i,1,("nNonperts"/.Options[stringFNmodels])}];
\[CapitalPhi]Chargemat=Take[Qmatrix,-("nNonperts"/.Options[stringFNmodels])];
(* Change to downstairs rep. *)
\[CapitalPhi]Chargematlow=Map[QMatRepChange[#,("nChargeVec"/.Options[stringFNmodels])]&,\[CapitalPhi]Chargemat];
A5anomalylow=QMatRepChange[A5anomaly,("nChargeVec"/.Options[stringFNmodels])];
MatrixEqn=Transpose[\[CapitalPhi]Chargematlow] . \[Beta]vec==A5anomalylow;
sol=FindInstance[MatrixEqn,\[Beta]vec,Integers];
(*Print[sol];*)
If[sol=={},
(* At the moment this means that there are no solutions that would satisfy the constraint. *)
(* Compute the solution in the field of Reals. *)
realsol=FindInstance[MatrixEqn,\[Beta]vec,Reals];
	If[realsol!={},
	(* Case 1: There is a solution in the field of Reals. *)
	roundingerror=Total[Map[Min[{Abs[(#-Floor[#])],Abs[(Ceiling[#]-#)]}]&,\[Beta]vec/.realsol[[1]]]];
	fitanomaly=roundingerror;
	,
	(* Case 2 : There is no solution in the field of reals. We calculate the image of A, and then find the minimum distance between \[Beta] and imA. *)
	\[Beta]vars=Table[x[i],{i,1,(Dimensions[\[CapitalPhi]Chargematlow][[1]])}];
	regioneqs=(Transpose[\[CapitalPhi]Chargematlow] . (\[Beta]vars));
	image\[CapitalPhi]mat=Map[ParametricRegion[regioneqs,#]&,{\[Beta]vars}][[1]];
	distance=N[RegionDistance[image\[CapitalPhi]mat,A5anomalylow]];
	fitanomaly=Max[Length[A5anomalylow],distance];
	];
,
fitanomaly=0];
,
fitanomaly=0;
];

<|"FitnessGSSec"->fitanomaly|>
]


(* ::Subsection::Closed:: *)
(*CheckCharge*)


CheckCharge[vector_]:=Module[{element1,element2,num},
element1=vector[[1]];
element2=vector[[2]];
If[element1==element2,
num=element1,
num=0
];
num
]


(* ::Subsection::Closed:: *)
(*FitnessNonAbelianFactors*)


FitnessNonAbelianFactors[state_]:=Module[{Qchargepos,Hchargepos,\[Phi]chargepos,chargepos,idenchargevec,replaceposrule,fitnessNAF,numSing},
(* The module does a consistency check on whether the charge pattern can be generated by the given charge vector n inputed as Options in the code. *)
(* The idea is sketched out in the note, essentially this checks whether the system can generate non-abelian charge patterns. *)
Qchargepos=Transpose[Partition[Flatten[Take[state,{4,9}]],3]];
Hchargepos=Partition[Flatten[Take[state,{10,11}]],2];
numSing="nSinglet"/.Options[stringFNmodels];
If[numSing>0,
\[Phi]chargepos=Transpose[Partition[Flatten[Take[state,{12,11+2*("nSinglet"/.Options[stringFNmodels])}]],("nSinglet"/.Options[stringFNmodels])]];
chargepos=Join[Qchargepos,Hchargepos,\[Phi]chargepos],
chargepos=Join[Qchargepos,Hchargepos]
];
idenchargevec=Map[CheckCharge[#]&,chargepos];
replaceposrule=Join[Map[#->1&,Flatten[Position["nChargeVec"/.Options[stringFNmodels],1]],1],
Map[#->0&,Complement[Table[i,{i,1,"nU1s"/.Options[stringFNmodels]}],Flatten[Position["nChargeVec"/.Options[stringFNmodels],1]]],1]];
fitnessNAF=Abs[Total[idenchargevec/.replaceposrule]];
<|"FitnessNAFSec"->fitnessNAF|>
]


(* ::Subsection::Closed:: *)
(*FitnessCalcQsectorOnly*)


FitnessCalcQsectorOnly[massrules_,stateVec_,yukO1Assoc_]:=Module[{TotalFit,FitHsector,FitQsector,FitnessCKMsector,FitGSsector,FitNAFsector,FitDiffExtra,ThreeFitSecs,
FitRGSector,FitQCKMsector,FitIndivExtra,Qmat,
QmassIndPen,fitnessQmass,fitnessCKM,QmixIndPen,totalfitassoc},

Qmat=StateVectoMatwithVEVs[stateVec];

FitHsector=FitnessHiggs[massrules];
FitGSsector=FitnessGSAnomaly[Qmat];
FitNAFsector=FitnessNonAbelianFactors[stateVec];
FitQCKMsector=FitnessQuarksMassMix[massrules];
FitRGSector=FitnessRG[yukO1Assoc];

(* 1. Sector Differences. *)
ThreeFitSecs={("FitnessHSec"/.FitHsector),("FitnessQMassSec"/.FitQCKMsector),("FitnessQMixSec"/.FitQCKMsector)};
If[("diffPenal"/.Options[stringFNmodels]),
FitDiffExtra=Abs[Max[ThreeFitSecs]-Min[ThreeFitSecs]],
FitDiffExtra=0];

(* 2. Individual Differences. *)
If["indPenal"/.Options[stringFNmodels],
fitnessQmass="FitnessQmass"/.FitQCKMsector;
fitnessCKM="FitnessCKM"/.FitQCKMsector;
QmassIndPen=Map[If[#<=0,0,#]&,(fitnessQmass-ConstantArray[("indFitMax"/.Options[stringFNmodels]),Length[fitnessQmass]])];
QmixIndPen=Map[If[#<=0,0,#]&,(fitnessCKM-ConstantArray[("indFitMax"/.Options[stringFNmodels]),Length[fitnessCKM]])];
FitIndivExtra=Total[QmassIndPen]+Total[QmixIndPen];
,
FitIndivExtra=0
];


TotalFit=("FitnessHSec"/.FitHsector)*("factorFitnessH"/.Options[stringFNmodels])+("FitnessQMassSec"/.FitQCKMsector)*("factorFitnessQ"/.Options[stringFNmodels])
+("FitnessQMixSec"/.FitQCKMsector)*("factorFitnessCKM"/.Options[stringFNmodels])
+("FitnessGSSec"/.FitGSsector)*("factorFitnessGS"/.Options[stringFNmodels])
+("FitnessNAFSec"/.FitNAFsector)*("factorFitnessNAF"/.Options[stringFNmodels])
+("FitnessRGSec"/.FitRGSector)*("factorFitnessRG"/.Options[stringFNmodels])
+FitDiffExtra+FitIndivExtra;

totalfitassoc=<|"Fitness"->-TotalFit|>;
Join[totalfitassoc,FitHsector,FitRGSector,FitQCKMsector,FitGSsector,FitNAFsector]


]


(* ::Section::Closed:: *)
(*Auxiliary Modules - Complete State Modules*)


(* ::Subsection::Closed:: *)
(*DetectBitlist*)


DetectBitlist[object_]:=Module[{lengthObj,requiredlength,lengthboolean,onlybitsboolean},
(* Detect whether the object is a bitlist. *)

lengthObj=Length[object];
If[("fixedCharges"/.Options[stringFNmodels]),
requiredlength=GADetermDimStateVecFix[][[1]];
,
requiredlength=GADetermDimStateVecwithVEVs[][[1]];
];
lengthboolean=(lengthObj==requiredlength);
onlybitsboolean=MemberQ[object,{2,3,4,5,6,7,8,9}];
lengthboolean&&Not[onlybitsboolean]

]


(* ::Subsection::Closed:: *)
(*ChangeO1VecByAmount*)


ChangeO1VecByAmount[O1VecToChange_,Pos_]:=Module[{NewO1Vec},
NewO1Vec=O1VecToChange;
If[Length[NewO1Vec]<Pos||Pos<=0,Print["Error. Position of O1 is out-of-range."];Return[0]];
NewO1Vec[[Pos]]=O1VecToChange[[Pos]]+("O1CoeffVarDiff"/.Options[stringFNmodels]);
Return[NewO1Vec]
]


(* ::Subsection::Closed:: *)
(*FitnessO1Coeff*)


FitnessO1Coeff[variedmassassoc_,massassoc_,O1val_]:=Module[{HiggsVEV,HiggsVEVvaried,tan\[Beta],tan\[Beta]varied,massDown,massDownVaried,massUp,massUpVaried,CKMentries,CKMentriesVaried,
physQuantsList,physQuantsListVaried,physQuantsDiff,O1DiffScale,dOdpList,assocToReturn,FitnessO1SecNum,physQuantsListInv},
HiggsVEV=massassoc[["Higgs"]];
HiggsVEVvaried=variedmassassoc[["Higgs"]];
tan\[Beta]=massassoc[["tan\[Beta]"]];
tan\[Beta]varied=variedmassassoc[["tan\[Beta]"]];
massDown=massassoc[["md"]];
massDownVaried=variedmassassoc[["md"]];
massUp=massassoc[["mu"]];
massUpVaried=variedmassassoc[["mu"]];
CKMentries=Flatten[massassoc[["CKMmat"]]];
CKMentriesVaried=Flatten[variedmassassoc[["CKMmat"]]];
physQuantsList=Join[{HiggsVEV},{tan\[Beta]},massUp,massDown,CKMentries];
physQuantsListVaried=Flatten[Join[{HiggsVEVvaried},{tan\[Beta]varied},massUpVaried,massDownVaried,CKMentriesVaried]];
(*Print[physQuantsListVaried];*)
physQuantsDiff=(physQuantsListVaried-physQuantsList);
O1DiffScale=("O1CoeffVarDiff"/.Options[stringFNmodels]);
(*Print[physQuantsList];*)
physQuantsListInv=Map[If[Abs[#]<10^-8,10^-8,Abs[#]]&,physQuantsList];
(*Print[physQuantsListInv];*)
dOdpList=Abs[O1val*physQuantsDiff/O1DiffScale*(1/physQuantsListInv)];
FitnessO1SecNum =Total[Map[Abs[#]&,dOdpList]];
assocToReturn=<|"FitnessO1Sec"->FitnessO1SecNum,"ListofdOdp"->dOdpList|>;
Return[assocToReturn]
]


(* ::Section::Closed:: *)
(*Auxiliary Modules - Finding Operators for Quark Sector*)


(* ::Subsection::Closed:: *)
(*AdjoinO1coeffWOps*)


AdjoinO1coeffWOps[WOplistoutput_,O1vec_,O1count_]:=Module[{Oplist,dimOplist,O1coeffneed,newO1count,raggedOplist,newO1coeffneed,j,Wexpansion,a,Adjoin\[Kappa]Rule,Adjoin\[Sigma]Rule},
(* Now append O(1) coefficients. *)
Oplist=Map[Flatten[{#}]&,WOplistoutput];
dimOplist=Total[Map[Dimensions,Oplist,1]][[1]];
O1coeffneed=Take[O1vec,{O1count,O1count+dimOplist-1}];
newO1count=O1count+dimOplist;
(* Write O1coefflist as the required ragged array. *)
raggedOplist=Oplist/.Table[("BundleModuliFields"/.Options[stringFNmodels])[i]->a,{i,1,"nSinglet"/.Options[stringFNmodels],1}]/.Table[("KahlerEffFields"/.Options[stringFNmodels])[i]->a,{i,1,"nNonperts"/.Options[stringFNmodels],1}];
raggedOplist=raggedOplist/.a->1;
j=1;
newO1coeffneed=Map[O1coeffneed[[j++]]&,raggedOplist,{-1}];
Wexpansion=newO1coeffneed*Oplist;
Wexpansion=Map[Flatten[{Total[#]}]&,Wexpansion,1];

(* Adjoining \[Kappa] factors. Only add to \[Phi] fields. *)
Adjoin\[Kappa]Rule=Table[("BundleModuliFields"/.Options[stringFNmodels])[i]->\[Kappa]*("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels],1}];
Wexpansion=Wexpansion/.Adjoin\[Kappa]Rule;
Wexpansion=Wexpansion+O[\[Kappa]]^(1+"DimOp"/.Options[stringFNmodels]);

(* Adjoining \[Sigma] factors for \[CapitalPhi] fields. *)
Adjoin\[Sigma]Rule=Table[("KahlerEffFields"/.Options[stringFNmodels])[i]->\[Sigma]*("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels],1}];
Wexpansion=Wexpansion/.Adjoin\[Sigma]Rule;
Wexpansion=Wexpansion+O[\[Sigma]]^(1+"DimOpof\[CapitalPhi]"/.Options[stringFNmodels]);

(* Returns the expansion and True. *)
{Wexpansion,newO1count}

]


(* ::Subsection::Closed:: *)
(*AdjoinO1coeffWOpsText*)


AdjoinO1coeffWOpsText[WOplistoutput_,MassRatioAssoc_,VEVPowersListAssoc_,O1vec_,O1count_]:=Module[
{Oplist,dimOplist,O1coeffneed,newO1count,raggedOplist,newO1coeffneed,j,Wexpansion,a,Adjoin\[Kappa]Rule,Adjoin\[Sigma]Rule,
NumDiag,WOpsListTop,\[Phi]Fields,\[CapitalPhi]Fields,\[Phi]VEVsList,\[CapitalPhi]VEVsList,\[Phi]VEVsReplaceRule,\[CapitalPhi]VEVsReplaceRule,Scale,ScalesList,
NumDiagNumYuk,NumDiagNumEntry,ScalesPowerList},
(* determine which O(1) coefficient to replace. *)
NumDiag=MassRatioAssoc[["PosOTop"]];
If[NumDiag==1||NumDiag==2||NumDiag==3,
If[NumDiag==1,NumDiagNumYuk=10;WOpsListTop=WOplistoutput[[10]]];
If[NumDiag==2,NumDiagNumYuk=14;WOpsListTop=WOplistoutput[[14]]];
If[NumDiag==3,NumDiagNumYuk=18;WOpsListTop=WOplistoutput[[18]]];
,
Print["PosOTop is Wrong. Exiting."];
Return[0];
];

Scale="VEVsScale"/.Options[stringFNmodels];
\[Phi]Fields=Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]Fields=Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}];
\[Phi]VEVsList=Map[Scale^#&,VEVPowersListAssoc[["\[Phi]VEVPowers"]]];
\[CapitalPhi]VEVsList=Map[Scale^#&,VEVPowersListAssoc[["\[CapitalPhi]VEVPowers"]]];
\[Phi]VEVsReplaceRule=Table[(\[Phi]Fields[[i]]->\[Phi]VEVsList[[i]]),{i,1,"nSinglet"/.Options[stringFNmodels]}];
\[CapitalPhi]VEVsReplaceRule=Table[(\[CapitalPhi]Fields[[i]]->\[CapitalPhi]VEVsList[[i]]),{i,1,"nNonperts"/.Options[stringFNmodels]}];
ScalesList=WOpsListTop/.\[Phi]VEVsReplaceRule/.\[CapitalPhi]VEVsReplaceRule;
If[ScalesList===0,
NumDiagNumEntry=1;
,
ScalesPowerList=Map[Exponent[#,Scale]&,ScalesList];
NumDiagNumEntry=Flatten[Position[ScalesPowerList,Min[ScalesPowerList]]][[1]];
];

(* Now append O(1) coefficients. *)
Oplist=Map[Flatten[{#}]&,WOplistoutput];
dimOplist=Total[Map[Dimensions,Oplist,1]][[1]];
O1coeffneed=Take[O1vec,{O1count,O1count+dimOplist-1}];
newO1count=O1count+dimOplist;
(* Write O1coefflist as the required ragged array. *)
raggedOplist=Oplist/.Table[("BundleModuliFields"/.Options[stringFNmodels])[i]->a,{i,1,"nSinglet"/.Options[stringFNmodels],1}]/.Table[("KahlerEffFields"/.Options[stringFNmodels])[i]->a,{i,1,"nNonperts"/.Options[stringFNmodels],1}];
raggedOplist=raggedOplist/.a->1;
j=1;
newO1coeffneed=Map[O1coeffneed[[j++]]&,raggedOplist,{-1}];
(* line to update the newO1coeff *)
newO1coeffneed[[NumDiagNumYuk]][[NumDiagNumEntry]]=Abs[MassRatioAssoc[["BestO1Top"]]];
Wexpansion=newO1coeffneed*Oplist;
Wexpansion=Map[Flatten[{Total[#]}]&,Wexpansion,1];

(* Adjoining \[Kappa] factors. Only add to \[Phi] fields. *)
Adjoin\[Kappa]Rule=Table[("BundleModuliFields"/.Options[stringFNmodels])[i]->\[Kappa]*("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels],1}];
Wexpansion=Wexpansion/.Adjoin\[Kappa]Rule;
Wexpansion=Wexpansion+O[\[Kappa]]^(1+"DimOp"/.Options[stringFNmodels]);

(* Adjoining \[Sigma] factors for \[CapitalPhi] fields. *)
Adjoin\[Sigma]Rule=Table[("KahlerEffFields"/.Options[stringFNmodels])[i]->\[Sigma]*("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels],1}];
Wexpansion=Wexpansion/.Adjoin\[Sigma]Rule;
Wexpansion=Wexpansion+O[\[Sigma]]^(1+"DimOpof\[CapitalPhi]"/.Options[stringFNmodels]);

(* Returns the expansion and True. *)
{Wexpansion,newO1count}

]


(* ::Subsection::Closed:: *)
(*FindYukExp*)


(* This submodule computes the Yukawa term expanions. *)
(* We are ignoring all the allowed operators in the exponentials - we are basically saying that they go into the scaling of the VEVs. *)
FindYukExp[WOps_]:=Module[{WOpsn,KOps,nFields,fields,Cvtlst\[CapitalPhi],add\[Kappa]rule,SUGRAterms,YukExp},
(* KOps=FindOpsforKmod[state,seed]; *)
WOpsn=WOps;
WOpsn=Map[Total[#]&,WOps,1];
YukExp=Normal[WOpsn];

YukExp
]

(* Here we have ignored the exponential scaling - we treat that as contributing to the VEVs. *)


(* ::Subsection::Closed:: *)
(*FindYukMatExpwith\[Phi]VEVsQsector*)


(* This module basically separates the Yukawa Expansion into Yukawa matrices. *)
FindYukMatExpwith\[Phi]VEVsQsector[Yexp_]:=Module[{Yd,Yu,Ydn,Yun},

(* Form the Yukawa matrices. *)
Yd=Partition[Take[Yexp,9],3];
Yu=Partition[Take[Yexp,{10,18}],3];

(* Convert list. *)
Ydn=Normal[Yd]/.{\[Kappa]->"Mpl/Mcy"/.Options[stringFNmodels],\[Sigma]->1};
Yun=Normal[Yu]/.{\[Kappa]->"Mpl/Mcy"/.Options[stringFNmodels],\[Sigma]->1};

{Ydn,Yun}
]


(* ::Subsection::Closed:: *)
(*YukMatsToAssoc*)


YukMatsToAssoc[YdMat_,YuMat_]:=Module[{yukMatAssoc},
yukMatAssoc=<|"YdMatNum"->YdMat,"YuMatNum"->YuMat|>;
Return[yukMatAssoc]
]


(* ::Section:: *)
(*Auxiliary Modules - Data Analysis*)


(* ::Subsection::Closed:: *)
(*LowestOrderTerms*)


LowestOrderTerms[expr_,vars_]:=Module[{polylist,orderlist,mindegree,express,return},
If[expr===0,
return=0,
polylist=Flatten[{MonomialList[expr]}];
orderlist=Map[Total[Exponent[#,vars]]&,polylist];
mindegree=Min[orderlist];
express=Select[polylist,Total[Exponent[#,vars]]==mindegree&];
return=Total[express]];
Return[return]
]


(* ::Subsection::Closed:: *)
(*VecDuplSymm*)


VecDuplSymm[vec_]:=Module[{dups,pos,outAssoc},
dups=Select[Tally[vec],#[[2]]>1&][[All,1]];
pos=Flatten[Position[vec,#]]&/@dups;
outAssoc=<|"RepEntry"->dups,"RepPos"->pos|>;
Return[outAssoc];
]


(* ::Subsection::Closed:: *)
(*RedunQSymm*)


RedunQSymm[]:=Module[{nVec,fixedQBoo,tenFQ,fiveBarHFQ,RepeatQPos,assocSymm},
nVec="nChargeVec"/.Options[stringFNmodels];
fixedQBoo="fixedCharges"/.Options[stringFNmodels];
tenFQ="tenQ"/.Options[stringFNmodels];
fiveBarHFQ="fiveBarHiggsQ"/.Options[stringFNmodels];

(* suppose the charge is fixed. *)
If[fixedQBoo,
	(* first check the available sectors. *)
	(* completely split case *)
	If[nVec=={1,1,1,1,1},
		RepeatQPos=Complement[tenFQ,Intersection[tenFQ,fiveBarHFQ]];
		Return[<|"SymmDir"->RepeatQPos|>]
	  ];
	(* n = (1,1,1,2) case *)
	If[nVec=={1,1,1,2},
		RepeatQPos=Complement[tenFQ,Intersection[tenFQ,fiveBarHFQ]];
		If[MemberQ[RepeatQPos,4],RepeatQPos={}];
		Return[<|"SymmDir"->RepeatQPos|>]
	];
	(* n = (1,1,3) case *)
	If[nVec=={1,1,3},
		RepeatQPos=Complement[tenFQ,Intersection[tenFQ,fiveBarHFQ]];
		Return[<|"SymmDir"->RepeatQPos|>]
	  ];
	(* n = (1,2,2) case *)
	If[nVec=={1,2,2},
		RepeatQPos={};
		Return[<|"SymmDir"->RepeatQPos|>]
	  ];
	Print["nVec or fixedQs entered incorrectly; please check."];
	Return[0];
,

(* now charges are not fixed. *)
assocSymm=VecDuplSymm[nVec];
Return[AssociateTo[assocSymm,"SymmDir"->assocSymm[["RepPos"]]]]
];

]


(* ::Subsection::Closed:: *)
(*CalcQScoreVec*)


CalcQScoreVec[fiveBarQ_,oneQ_,kahlQ_]:=Module[{numSing,numKahl,numU1s,nVec, symmAssoc,symmDir,onePosQs,oneNegQs,
fixedQBoo,symmDirOne,symmDirTwo,countOne,countTwo},
(* Module is used to calculate the score between the directions where the charges are equal and has an Subscript[S, 2] symmetry. *)

numSing="nSinglet"/.Options[stringFNmodels];
numKahl="nNonperts"/.Options[stringFNmodels];
numU1s="nU1s"/.Options[stringFNmodels];
nVec="nChargeVec"/.Options[stringFNmodels];
fixedQBoo="fixedCharges"/.Options[stringFNmodels];

(* Extract the relevant positions. *)
If[fixedQBoo,
	symmAssoc=RedunQSymm[];
	symmDir=symmAssoc[["SymmDir"]];
	If[symmDir=={},
	(* this means we dont have any replacement rules. *)
	AssociateTo[symmAssoc,"NeedReplace"->False];
	AssociateTo[symmAssoc,"ReplRule"->{}];
	AssociateTo[symmAssoc,"ResultDupli"->False];
	Return[symmAssoc]
	,
	symmDirOne=symmDir[[1]];
	symmDirTwo=symmDir[[2]];
	];
,
	Print["Not implemented yet. Will come back to this."];
	symmAssoc=RedunQSymm[];
	symmDir=symmAssoc[["SymmDir"]][[1]];
	If[symmDir=={},
	(* this means we dont have any replacement rules. *)
	AssociateTo[symmAssoc,"NeedReplace"->False];
	AssociateTo[symmAssoc,"ReplRule"->{}];
	Return[symmAssoc]
	,
	symmDirOne=symmDir[[1]];
	symmDirTwo=symmDir[[2]];
	];
	
];

(*Print[{symmDirOne,symmDirTwo}];*)

If[numKahl===0,
	onePosQs=Take[oneQ,numSing];
	oneNegQs=Take[oneQ,-numSing];
	countOne=Count[fiveBarQ,symmDirOne]+Count[onePosQs,symmDirOne]-Count[oneNegQs,symmDirOne];
	countTwo=Count[fiveBarQ,symmDirTwo]+Count[onePosQs,symmDirTwo]-Count[oneNegQs,symmDirTwo];
,
	onePosQs=Take[oneQ,numSing];
	oneNegQs=Take[oneQ,-numSing];
	countOne=Count[fiveBarQ,symmDirOne]+Count[onePosQs,symmDirOne]-Count[oneNegQs,symmDirOne]+Total[Partition[kahlQ,numU1s][[All,symmDirOne]]];
	countTwo=Count[fiveBarQ,symmDirTwo]+Count[onePosQs,symmDirTwo]-Count[oneNegQs,symmDirTwo]+Total[Partition[kahlQ,numU1s][[All,symmDirTwo]]];
];

If[countOne>countTwo,
	AssociateTo[symmAssoc,"NeedReplace"->False];
	AssociateTo[symmAssoc,"SymmDirScore"->{countOne,countTwo}];
	AssociateTo[symmAssoc,"ReplRule"->{}];
	AssociateTo[symmAssoc,"ResultDupli"->False];
Return[symmAssoc]
,
	If[countOne<countTwo,
		AssociateTo[symmAssoc,"NeedReplace"->True];
		AssociateTo[symmAssoc,"SymmDirScore"->{countOne,countTwo}];
		AssociateTo[symmAssoc,"ReplRule"->{symmDirTwo->symmDirOne,symmDirOne->symmDirTwo}];
		AssociateTo[symmAssoc,"ResultDupli"->False];
		Return[symmAssoc];
	,
		AssociateTo[symmAssoc,"NeedReplace"->False];
		AssociateTo[symmAssoc,"SymmDirScore"->{countOne,countTwo}];
		AssociateTo[symmAssoc,"ReplRule"->{}];
		AssociateTo[symmAssoc,"ResultDupli"->True];
		Return[symmAssoc];
	]
];

]


(* ::Subsection::Closed:: *)
(*StateNormalForm*)


StateNormalForm[stateVec_]:=Module[{numSing,numKahl,numU1s,tenQs,fiveBarQs,HdQ,\[Phi]Qs,\[CapitalPhi]Qs,VEVpart,\[Phi]VEVs,\[CapitalPhi]VEVs,
fixedQBoo,symmAssoc,
newFiveBarHiggsQs,newFiveBarQs,newTenQs,new\[Phi]Qs,new\[CapitalPhi]Qs,
new\[CapitalPhi]Qmat,old\[CapitalPhi]Qmat,
sorted\[Phi]data,new\[Phi]VEVs,old\[CapitalPhi]QmatTr,sorted\[CapitalPhi]data,new\[CapitalPhi]VEVs,newSV
},
(* This module provides the normal form of a statevector *)
numSing="nSinglet"/.Options[stringFNmodels];
numKahl="nNonperts"/.Options[stringFNmodels];
numU1s="nU1s"/.Options[stringFNmodels];
fixedQBoo="fixedCharges"/.Options[stringFNmodels];

(* extract the components from the stateVec *)
tenQs=Take[stateVec,3];
fiveBarQs=Take[stateVec,{4,9}];
HdQ=Take[stateVec,{10,11}];
\[Phi]Qs=Take[stateVec,{12,11+2*numSing}];
VEVpart=Take[stateVec,-(numSing+numKahl)];
\[Phi]VEVs=Take[VEVpart,numSing];
\[CapitalPhi]Qs={};
If[numKahl>0,
\[CapitalPhi]Qs=Take[stateVec,{12+2*numSing,11+2*numSing+numKahl*numU1s}];
\[CapitalPhi]VEVs=Take[VEVpart,-numKahl];
];

(* If the fix charge option is on, then we relabel the state in the following manner *)

If[Not[fixedQBoo],
Print["Not implemented. Exiting temporarily."];
Return[0];
,
symmAssoc=CalcQScoreVec[fiveBarQs,\[Phi]Qs,\[CapitalPhi]Qs];
];

(* check if we need replacing. *)
If[symmAssoc[["NeedReplace"]],
	newTenQs=tenQs;
	newFiveBarQs=fiveBarQs/.symmAssoc[["ReplRule"]];
	newFiveBarHiggsQs=HdQ/.symmAssoc[["ReplRule"]];
	new\[Phi]Qs=\[Phi]Qs/.symmAssoc[["ReplRule"]];
	
	If[numKahl>0,
		old\[CapitalPhi]Qmat=Partition[\[CapitalPhi]Qs,numU1s];	
		old\[CapitalPhi]QmatTr=Transpose[old\[CapitalPhi]Qmat];
		new\[CapitalPhi]Qmat=Transpose[ReplacePart[old\[CapitalPhi]QmatTr,{(symmAssoc[["SymmDir"]][[1]])->old\[CapitalPhi]QmatTr[[(symmAssoc[["SymmDir"]][[2]])]],(symmAssoc[["SymmDir"]][[2]])->old\[CapitalPhi]QmatTr[[(symmAssoc[["SymmDir"]][[1]])]]}]];
	];
,
	newTenQs=tenQs;
	newFiveBarQs=fiveBarQs;
	newFiveBarHiggsQs=HdQ;
	new\[Phi]Qs=\[Phi]Qs;
	If[numKahl>0,
		new\[CapitalPhi]Qmat=Partition[\[CapitalPhi]Qs,numU1s];
	];
];

(* now rearrange so we get a nice rearrangement *)
(* sort 5bar + Higgs charges: *)
newFiveBarQs=Flatten[Transpose[SortBy[Map[Sort,Transpose[Partition[newFiveBarQs,3]]],First]]];
newFiveBarHiggsQs=Sort[newFiveBarHiggsQs];

(* sort \[Phi]s *)
If[numSing>0,
sorted\[Phi]data=SortBy[Transpose[Join[Partition[new\[Phi]Qs,numSing],Transpose[Partition[\[Phi]VEVs,1]]]],First];
new\[Phi]Qs=Flatten[Take[Transpose[sorted\[Phi]data],2]];
new\[Phi]VEVs=Flatten[Transpose[Take[Transpose[sorted\[Phi]data],-1]]];
,
new\[Phi]Qs={};
new\[Phi]VEVs={};
];


(* sort the Kahler effective charges *)
	If[numKahl>0,
	sorted\[CapitalPhi]data=SortBy[Transpose[Join[Transpose[new\[CapitalPhi]Qmat],Transpose[Partition[\[CapitalPhi]VEVs,1]]]],First];
	(*Print[sorted\[CapitalPhi]data];*)
	new\[CapitalPhi]Qs=Flatten[Take[Transpose[sorted\[CapitalPhi]data],numU1s]];
	new\[CapitalPhi]VEVs=Flatten[Transpose[Take[Transpose[sorted\[CapitalPhi]data],-1]]];
	,
	new\[CapitalPhi]Qs={};
	new\[CapitalPhi]VEVs={};
	];

(* join the snippets to give the overall matrix. *)
newSV=Flatten[Join[newTenQs,newFiveBarQs,newFiveBarHiggsQs,Flatten[new\[Phi]Qs],Flatten[new\[CapitalPhi]Qs],Flatten[new\[Phi]VEVs],Flatten[new\[CapitalPhi]VEVs]]];
Return[newSV]
]


(* ::Subsection::Closed:: *)
(*EvalStateVecSoloList*)


EvalStateVecSoloList[inputlist_]:=Module[{assoclist},
assoclist=Map[
Append[#,<|"StateVec"->(GAConvertBitlstStateVecFix[#[["Bitlist"]]])|>]
&,inputlist
];
assoclist=DeleteCases[assoclist,#[["Terminal"]]==False&];
Return[assoclist]
]


(* ::Subsection::Closed:: *)
(*AppendSetNoToStates*)


AppendSetNoToStates[states_]:=Module[{indexedStates},
indexedStates=MapIndexed[Function[{sublist,index},If[Length[sublist]==0,sublist,Map[Append[#,"SetNo"->First[index]]&,sublist]]],states];
indexedStates]


(* ::Subsection::Closed:: *)
(*WriteUniqueStates*)


WriteUniqueStates[filename_]:=Module[{text,splittext,runstates,evalstates,setNumRule},
text=ToExpression[ReadList[filename]];
splittext=text/.{stringFNmodels`BREAK->{}};
runstates=Map[ToExpression,splittext];
evalstates=Map[EvalStateVecSoloList,runstates];
evalstates=AppendSetNoToStates[evalstates];
Return[evalstates]
]


(* ::Subsection::Closed:: *)
(*SelectVecsSameCharge*)


SelectVecsSameCharge[liststates_]:=Module[{dimQ,outputlist},
dimQ=11+("nSinglet"/.Options[stringFNmodels])*2+("nNonperts"/.Options[stringFNmodels])*("nU1s"/.Options[stringFNmodels]);
outputlist=DeleteDuplicatesBy[liststates,Take[StateNormalForm[#[["StateVec"]]],dimQ]&];
Return[outputlist];
]


(* ::Subsection::Closed:: *)
(*SelectVecsSameModel*)


SelectVecsSameModel[liststates_]:=Module[{dimQ,outputlist},
dimQ=11+("nSinglet"/.Options[stringFNmodels])*2+("nNonperts"/.Options[stringFNmodels])*("nU1s"/.Options[stringFNmodels]);
outputlist=DeleteDuplicatesBy[liststates,StateNormalForm[#[["StateVec"]]]&];
Return[outputlist];
]


(* ::Subsection:: *)
(*ProduceReducedList*)


ProduceReducedList[evalstates_]:=Module[{evalStatesList,redEvalStsLst},
evalStatesList=Flatten[evalstates];
redEvalStsLst=SelectVecsSameCharge[evalStatesList];
Return[redEvalStsLst]
]


(* ::Subsection:: *)
(*ProduceReducedModelList*)


ProduceReducedModelList[evalstates_]:=Module[{evalStatesList,redEvalStsLst},
evalStatesList=Flatten[evalstates];
redEvalStsLst=SelectVecsSameModel[evalStatesList];
Return[redEvalStsLst]
]


(* ::Subsection::Closed:: *)
(*ComputeNumModuli*)


ComputeNumModuli[state_]:=Module[{stateVec,superOps,vars,countSing,countKahl,assocNum,newStateVec},
stateVec=state[["StateVec"]];
superOps=ComputeWOpsGivenStateVec[stateVec];
vars=Variables[superOps];

countSing=Count[vars,("BundleModuliFields"/.Options[stringFNmodels])[_Integer?Positive]];
countKahl=Count[vars,("KahlerEffFields"/.Options[stringFNmodels])[_Integer?Positive]];
assocNum=<|"numSing"->countSing,"numKahl"->countKahl|>;
newStateVec=Append[state,assocNum];
Return[newStateVec]
]


(* ::Subsection:: *)
(*ProduceFullRedList*)


ProduceFullRedList[filename_,outputFile_]:=Module[{evalstates,redevalstates,statesWithNumMod},
evalstates=WriteUniqueStates[filename];
redevalstates=Flatten[ProduceReducedList[evalstates]];
statesWithNumMod=Map[ComputeNumModuli,redevalstates];
Put[statesWithNumMod,outputFile];
Print[Length[evalstates]];
]


(* ::Subsection:: *)
(*ProduceFullRedModelList*)


ProduceFullRedModelList[filename_,outputFile_]:=Module[{evalstates,redevalstates,statesWithNumMod},
evalstates=WriteUniqueStates[filename];
redevalstates=Flatten[ProduceReducedModelList[evalstates]];
statesWithNumMod=Map[ComputeNumModuli,redevalstates];
Put[statesWithNumMod,outputFile];
Print[Length[evalstates]];
]


(* ::Subsection::Closed:: *)
(*PlotAccumGraph*)


PlotAccumGraph[line_,title_]:=Module[{length,tableofnum,DataToPlot},
length=Length[line];
tableofnum=Table[i,{i,0,length-1}];
tableofnum=tableofnum*20*40000;
DataToPlot=Transpose[{tableofnum,line}];
ListPlot[DataToPlot,PlotRange->All,Joined->True,PlotLabel->Style[title,20,Bold],FrameLabel->{Style["Number of states visited by GA",12],Style["Number of valid models",12]},GridLines->Automatic,Frame->True,Axes->True,PlotRangePadding->{{.75,Automatic},{5,Automatic}}]
]


(* ::Subsection::Closed:: *)
(*FromRedListGetPlot*)


FromRedListGetPlot[filePath_,numTotRuns_]:=Module[{lengthRuns,setNoList,listOfNums,countsVector,redevalstates,graph2},
redevalstates=Get[filePath];
lengthRuns=numTotRuns;
setNoList=Range[lengthRuns];
listOfNums=Counts[redevalstates[[All,"SetNo"]]];
countsVector=ConstantArray[0,lengthRuns];
countsVector=ReplacePart[countsVector,KeyValueMap[#1->#2&,listOfNums]];
graph2=PlotAccumGraph[Accumulate[countsVector],"Number of terminal states reduced by symmetry against GA runs"];
Print[graph2];
]


(* ::Subsection::Closed:: *)
(*GivenModNumFromRedListGetPlot*)


GivenModNumFromRedListGetPlot[filePath_,numTotRuns_,numSing_,numKahl_]:=Module[{lengthRuns,setNoList,listOfNums,countsVector,redevalstates,graph2,redevalstatesmod},
redevalstates=Get[filePath];
redevalstatesmod=Select[redevalstates,(#[["numSing"]]==numSing&&#[["numKahl"]]==numKahl)&];
(*Print[redevalstatesmod];*)
lengthRuns=numTotRuns;
setNoList=Range[lengthRuns];
listOfNums=Counts[redevalstatesmod[[All,"SetNo"]]];
countsVector=ConstantArray[0,lengthRuns];
countsVector=ReplacePart[countsVector,KeyValueMap[#1->#2&,listOfNums]];
graph2=PlotAccumGraph[Accumulate[countsVector],"Number of terminal states reduced by symmetry against GA runs"];
Print[graph2];
]


(* ::Subsection:: *)
(*FromStateVecListGiveNum*)


FromStateVecListGiveNum[filename_]:=Module[{listStates,listRedStates,lenList,dimQ},
(* module extracts the statevector, and reduces it and give the number *)
listStates=Get[filename];
dimQ=11+("nSinglet"/.Options[stringFNmodels])*2+("nNonperts"/.Options[stringFNmodels])*("nU1s"/.Options[stringFNmodels]);
listRedStates=DeleteDuplicatesBy[listStates,ChargeStateNormalForm[(Take[#[["StateVec"]],dimQ])]&];
(*Print[listRedStates];*)
lenList=Length[listRedStates];
Return[lenList]
]


(* ::Subsection::Closed:: *)
(*ChargeStateNormalForm*)


ChargeStateNormalForm[QstateVec_]:=Module[{numSing,numKahl,numU1s,fixedQBoo,
tenQs,fiveBarQs,HdQ,\[Phi]Qs,\[CapitalPhi]Qs,symmAssoc,newTenQs,newFiveBarQs,newFiveBarHiggsQs,new\[Phi]Qs,
old\[CapitalPhi]Qmat,old\[CapitalPhi]QmatTr,new\[CapitalPhi]Qmat,sorted\[Phi]data,sorted\[CapitalPhi]data,new\[CapitalPhi]Qs,newSV

},
(* This module provides the normal form of a statevector *)
numSing="nSinglet"/.Options[stringFNmodels];
numKahl="nNonperts"/.Options[stringFNmodels];
numU1s="nU1s"/.Options[stringFNmodels];
fixedQBoo="fixedCharges"/.Options[stringFNmodels];

(* extract the components from the stateVec *)
tenQs=Take[QstateVec,3];
fiveBarQs=Take[QstateVec,{4,9}];
HdQ=Take[QstateVec,{10,11}];
\[Phi]Qs=Take[QstateVec,{12,11+2*numSing}];
\[CapitalPhi]Qs={};
If[numKahl>0,
\[CapitalPhi]Qs=Take[QstateVec,-numKahl*numU1s];
];

(* If the fix charge option is on, then we relabel the state in the following manner *)

If[Not[fixedQBoo],
Print["Not implemented. Exiting temporarily."];
Return[0];
,
symmAssoc=CalcQScoreVec[fiveBarQs,\[Phi]Qs,\[CapitalPhi]Qs];
];

(* check if we need replacing. *)
If[symmAssoc[["NeedReplace"]],
	newTenQs=tenQs;
	newFiveBarQs=fiveBarQs/.symmAssoc[["ReplRule"]];
	newFiveBarHiggsQs=HdQ/.symmAssoc[["ReplRule"]];
	new\[Phi]Qs=\[Phi]Qs/.symmAssoc[["ReplRule"]];
	
	If[numKahl>0,
		old\[CapitalPhi]Qmat=Partition[\[CapitalPhi]Qs,numU1s];	
		old\[CapitalPhi]QmatTr=Transpose[old\[CapitalPhi]Qmat];
		new\[CapitalPhi]Qmat=Transpose[ReplacePart[old\[CapitalPhi]QmatTr,{(symmAssoc[["SymmDir"]][[1]])->old\[CapitalPhi]QmatTr[[(symmAssoc[["SymmDir"]][[2]])]],(symmAssoc[["SymmDir"]][[2]])->old\[CapitalPhi]QmatTr[[(symmAssoc[["SymmDir"]][[1]])]]}]];
	];
,
	newTenQs=tenQs;
	newFiveBarQs=fiveBarQs;
	newFiveBarHiggsQs=HdQ;
	new\[Phi]Qs=\[Phi]Qs;
	If[numKahl>0,
		new\[CapitalPhi]Qmat=Partition[\[CapitalPhi]Qs,numU1s];
	];
];

(* now rearrange so we get a nice rearrangement *)
(* sort 5bar + Higgs charges: *)
newFiveBarQs=Flatten[Transpose[SortBy[Map[Sort,Transpose[Partition[newFiveBarQs,3]]],First]]];
newFiveBarHiggsQs=Sort[newFiveBarHiggsQs];

(* sort \[Phi]s *)
If[numSing>0,
sorted\[Phi]data=SortBy[Transpose[Partition[new\[Phi]Qs,numSing]],First];
new\[Phi]Qs=Flatten[Transpose[sorted\[Phi]data]];
,
new\[Phi]Qs={};
];


(* sort the Kahler effective charges *)
	If[numKahl>0,
	sorted\[CapitalPhi]data=SortBy[Transpose[Transpose[new\[CapitalPhi]Qmat]],First];
	(*Print[sorted\[CapitalPhi]data];*)
	new\[CapitalPhi]Qs=Flatten[Transpose[sorted\[CapitalPhi]data]];
	,
	new\[CapitalPhi]Qs={};
	];

(* join the snippets to give the overall matrix. *)
newSV=Flatten[Join[newTenQs,newFiveBarQs,newFiveBarHiggsQs,Flatten[new\[Phi]Qs],Flatten[new\[CapitalPhi]Qs]]];
Return[newSV]
]


(* ::Section:: *)
(*Main Modules*)


(* ::Subsection::Closed:: *)
(*CompleteStateGivenWOps*)


CompleteStateGivenWOps[stateVec_,WOps_,assocMassHierScale_,powerList_,O1CoeffVec_]:=Module[{O1Assoc,YukMZValsO1Assoc,QMassAssoc,
fitnessAssoc,\[Phi]Vals,\[CapitalPhi]Vals,fitness3Secs,YukdOps,YukuOps,fields,YukdSt,YukuSt,stateAssoc
},

(* we have now changed this function to introduce the new fitness Texture part. *)
(* Step 1: Determine the O(1) coefficient details . *)

O1Assoc=LocateExtractO1Coeffs[WOps,powerList,assocMassHierScale,O1CoeffVec];

(* Step 2: Determine the \[Phi]VEVs, and substitute into the corresponding operators to obtain the valued matrices. *)
YukMZValsO1Assoc=AdjoinO1coeffWOpsTextQSec[WOps,assocMassHierScale,powerList,O1Assoc,O1CoeffVec];
(*If[("ComputeKahlerOn"/.Options[stringFNmodels]),
{KmodVal,KSMVal,KHiggsVal}={KmodOps,KSMOps,KHiggsOps}/.\[Phi]VEVsVal;];*)


(* Step 3: Compute the SM Quantities *)
QMassAssoc=ComputeSMQuantities[YukMZValsO1Assoc[["YdMatVal"]],YukMZValsO1Assoc[["YuMatVal"]]];

(* Step 4: Detemine the fitness of the system. *)
fitnessAssoc=FitnessCalcQsectorOnly[QMassAssoc,stateVec,YukMZValsO1Assoc];

(* Step 5: Finally convert everthing into association form. *)
\[Phi]Vals=Map[(assocMassHierScale[["BestScale"]])^#&,powerList[["\[Phi]VEVPowers"]]];
\[CapitalPhi]Vals=Map[(assocMassHierScale[["BestScale"]])^#&,powerList[["\[CapitalPhi]VEVPowers"]]];

YukdOps=YukMZValsO1Assoc[["YdMatOps"]];
YukuOps=YukMZValsO1Assoc[["YuMatOps"]];
fitness3Secs=Join[{("FitnessHSec"/.fitnessAssoc)},{("FitnessQMassSec"/.fitnessAssoc)},{("FitnessQMixSec"/.fitnessAssoc)}];

fields=Join[Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}],Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}]];
YukdSt=Partition[Flatten[Map[LowestOrderTerms[#,fields]&,Flatten[YukdOps]]],3];
YukuSt=Partition[Flatten[Map[LowestOrderTerms[#,fields]&,Flatten[YukuOps]]],3];


stateAssoc=Join[fitnessAssoc,QMassAssoc,<|"ValueLst"->fitness3Secs,"\[Phi]VEVs"->\[Phi]Vals,"\[CapitalPhi]VEVs"->\[CapitalPhi]Vals,"YukdMatSt"->(YukdSt//MatrixForm),"YukuMatSt"->(YukuSt//MatrixForm)|>];
stateAssoc
]



(* ::Subsection::Closed:: *)
(*CompleteStateO1Coeffs*)


CompleteStateO1Coeffs[state_,O1CoeffVec_,opt___]:=Module[{stateBitList,bitExists,stateVec,LengthSV,
WijOp,powerList,assocLeadPowers,assocMassHier,assocMassHierScale,fitnessTextureAssoc,UnchangedO1VecNum,
UnchangedO1VecNumB,UnchangedO1VecNumT,O1CoeffPosChange,O1Assoc,O1coeffVecList,oldO1Assoc,
newO1AssocList,O1FitnessAssocList,O1FitList,O1FitTot,newFitness,stateAssoc,terminal,
partState,Qassoc,fullState,O1CoeffChangeLst,O1FitListMat},

If[AssociationQ[state],
If[KeyExistsQ[state,"Bits"],
stateBitList=state["Bits"];
bitExists=True,
bitExists=False
];
stateVec=If[KeyExistsQ[state,"StateVec"],state["StateVec"],
(If[("fixedCharges"/.Options[stringFNmodels]),GAConvertBitlstStateVecFix[stateBitList],GAConvertBitlstStateVecwithVEVs[stateBitList]])];
,
	(* not an association *)
	If[DetectBitlist[state],
	stateBitList=state;
	stateVec=(If[("fixedCharges"/.Options[stringFNmodels]),GAConvertBitlstStateVecFix[stateBitList],GAConvertBitlstStateVecwithVEVs[stateBitList]]);
	bitExists=True;
	,
		LengthSV=3*3+2+(2*"nSinglet"/.Options[stringFNmodels])+("nNonperts"/.Options[stringFNmodels])*("nU1s"/.Options[stringFNmodels])+("nSinglet"+"nNonperts")/.Options[stringFNmodels];
		If[(Length[state]==LengthSV),
		stateVec=state;
		bitExists=False;
		,
		Return["Error, invalid input."];
		];
	];
];


WijOp=ComputeWOpsGivenStateVec[stateVec];

(* Computation of the texture quantities. *)
powerList=ConvertStateVecToVEVPowers[stateVec];
assocLeadPowers=FromWOpsGetLeadingYukawaMatrix[WijOp,powerList];
assocMassHier=MassHierarchyFromLeadingOrders[assocLeadPowers];
assocMassHierScale=CalcScaleAndHiggsO1Coeff[assocMassHier];
O1Assoc=LocateExtractO1Coeffs[WijOp,powerList,assocMassHierScale,O1CoeffVec];
fitnessTextureAssoc=FitnessTexture[assocMassHierScale,powerList];

(* Computation of the O(1) Coefficients Variations. *)
If[assocMassHierScale[["PosOTop"]]===1||assocMassHierScale[["PosOTop"]]===2||assocMassHierScale[["PosOTop"]]===3,
If[assocMassHierScale[["PosOTop"]]===1,UnchangedO1VecNumT=10];
If[assocMassHierScale[["PosOTop"]]===2,UnchangedO1VecNumT=14];
If[assocMassHierScale[["PosOTop"]]===3,UnchangedO1VecNumT=18];
,
Print["Error, PosOTop is incorrect."];Return[0];
];

If[assocMassHierScale[["PosOBot"]]===1||assocMassHierScale[["PosOBot"]]===2||assocMassHierScale[["PosOBot"]]===3,
If[assocMassHierScale[["PosOBot"]]===1,UnchangedO1VecNumB=1];
If[assocMassHierScale[["PosOBot"]]===2,UnchangedO1VecNumB=5];
If[assocMassHierScale[["PosOBot"]]===3,UnchangedO1VecNumB=9];
,
Print["Error, PosOBot is incorrect."];Return[0];
];
(* We have fixed both the top and the bottom. *)
UnchangedO1VecNum={{UnchangedO1VecNumB},{UnchangedO1VecNumT}};

O1CoeffPosChange=O1Assoc[["O1LeadPosList"]];
O1CoeffChangeLst=O1Assoc[["O1LeadList"]];
O1CoeffPosChange=Delete[O1CoeffPosChange,UnchangedO1VecNum];
O1CoeffChangeLst=Delete[O1CoeffChangeLst,UnchangedO1VecNum];
O1coeffVecList=Map[ChangeO1VecByAmount[O1CoeffVec,#]&,O1CoeffPosChange];



oldO1Assoc=CompleteStateGivenWOps[stateVec,WijOp,assocMassHierScale,powerList,O1CoeffVec];
newO1AssocList=Map[CompleteStateGivenWOps[stateVec,WijOp,assocMassHierScale,powerList,#]&,O1coeffVecList];
(*Print[newO1AssocList];*)

(* Now Compute the O1 variation procedure. *)
O1FitnessAssocList=Map[FitnessO1Coeff[#[[1]],oldO1Assoc,#[[2]]]&,Thread[{newO1AssocList,O1CoeffChangeLst}]];
O1FitListMat=Map[#[["ListofdOdp"]]&,O1FitnessAssocList];
(*Print[O1FitListMat];*)
O1FitList=Map[Max,Transpose[O1FitListMat]];
O1FitTot=Total[O1FitList];

(*
fitnesslst=Map[#[["Fitness"]]&,assoclist];
bestfitness=Max[fitnesslst];
meanfitness=Mean[fitnesslst];
stddevfitness=StandardDeviation[fitnesslst];
O1vecMaximise=Flatten[Position[fitnesslst,Max[fitnesslst]]][[1]];
bestassoc=assoclist[[O1vecMaximise]];
newfitness=-(Abs[meanfitness]+Abs[stddevfitness]*("factorFitnessO1"/.Options[stringFNmodels])+fitnessTextureAssoc[["FitnessTextSec"]]*("factorFitnessTEX"/.Options[stringFNmodels]));
*)
newFitness=-(Abs[oldO1Assoc[["Fitness"]]]
+Abs[O1FitTot]*("factorFitnessO1"/.Options[stringFNmodels])
+Abs[fitnessTextureAssoc[["FitnessTextSec"]]]*("factorFitnessTEX"/.Options[stringFNmodels]));

(* Update fitness *)
stateAssoc=oldO1Assoc;
stateAssoc=Insert[stateAssoc,"Fitness"->newFitness,Key["Fitness"]];
stateAssoc=Insert[stateAssoc,"FitnessO1Sec"->O1FitTot,Key["ValueLst"]];
stateAssoc=Insert[stateAssoc,"FitnessO1List"->O1FitList,Key["FitnessO1Sec"]];
(*stateassoc=Insert[stateassoc,"BestO1vec"->O1vecMaximise,Key["JarlskogInv"]];*)

terminal=(newFitness>=("FitRange"/.Options[stringFNmodels]));
If[bitExists,
partState=Association[{"Bits"->stateBitList,"StateVec"->stateVec,"Terminal"->terminal}],
partState=Association[{"StateVec"->stateVec,"Terminal"->terminal}]
];
Qassoc=ConvertStatewithVEVs[stateVec];

fullState=Join[stateAssoc,partState,fitnessTextureAssoc,Qassoc,assocMassHierScale];

Return[fullState]
]


(* ::Subsection::Closed:: *)
(*CompleteStateQsector*)


CompleteStateQsector[state_,opt___]:=Module[{fullstate},
fullstate=CompleteStateO1Coeffs[state,O1CLst];
Return[fullstate]
]


(* ::Subsection::Closed:: *)
(*GAFitnessFuncQsector*)


GAFitnessFuncQsector[bitlst_,opt___]:=Module[{fullstate,stateassoc,state,fitness,terminal,partstate,Qassoc},

fullstate=CompleteStateO1Coeffs[bitlst,O1CLst];

If[fullstate[["Fitness"]]<=-1000,
AssociateTo[fullstate,"Fitness"->-1000];
];

Return[fullstate]
]


(* ::Subsection::Closed:: *)
(*GAInitialiseStateQsector*)


GAInitialiseStateQsector[opt___]:=Module[{randombitlst,stateassoc},

If[("fixedCharges"/.Options[stringFNmodels]),
randombitlst=GARandomStateFix[];,
randombitlst=GARandomStatewithVEVs[];
];
stateassoc=CompleteStateQsector[randombitlst];

stateassoc
]


(* ::Section::Closed:: *)
(*Recycling Bin...*)


(* ::Subsection::Closed:: *)
(*CompleteStateQsectorAssoc*)


CompleteStateQsectorAssoc[stateassocin_,opt___]:=Module[
{state,\[Phi]VEVsVal,O1vec,O1coeffcounter,WIJaOps,KmodOps,KSMOps,KHiggsOps,WIJaVal,KmodVal,KSMVal,KHiggsVal,YukVal,YdVal,YuVal,YeVal,
Qbarmat,dmat,umat,Lbarmat,emat,Hdscale,Huscale,Ydnew,Yunew,HiggsRules,CKMRules,Qmassrules,chargematfromstate,fitnessval,
HdVal,HuVal,tan\[Beta]Val,CKMmatVal,mdVal,muVal,\[Phi]Vals,\[CapitalPhi]Vals,stateasssoc,conversionvec,YukdOps,YukuOps,
fitnessHval,fitnessQval,fitnessGSval,fitnessCKMval,fitnessNAFval,fitnesstestval,fields,YukdSt,YukuSt,
Jarlskogmeasured,CKMcmplx},

(* first check this is a stateassoc *)
If[Not[AssociationQ[stateassocin]],Print["Input not an association."];Return[0]];

state=AssocFormToStateVec[stateassocin];

(* Step 1: Find the relevant operators. *)
(* Note O1coefflist is the O1vec required. *)
(* In this step we have now used the newer modules using the Smith decomposition. *)
(* Setting to turn on what type of method we are generating the operators. *)
If[(("MethodGenOps"/.Options[stringFNmodels])=="Solve")||(("MethodGenOps"/.Options[stringFNmodels])=="Smith")||
(("MethodGenOps"/.Options[stringFNmodels])=="List")||(("MethodGenOps"/.Options[stringFNmodels])=="SolveA")||(("MethodGenOps"/.Options[stringFNmodels])=="Optimal"),

If[("MethodGenOps"/.Options[stringFNmodels])=="Solve",
{WIJaOps,O1coeffcounter}=AdjoinO1coeffWOps[(FindOpsforWIJawith\[Phi]VEVsQsector[state]),O1coefflist,1];
];
If[("MethodGenOps"/.Options[stringFNmodels])=="Smith",
{WIJaOps,O1coeffcounter}=AdjoinO1coeffWOps[(FindOpsforWIJaSmithwith\[Phi]VEVsQsector[state]),O1coefflist,1];
];
If[("MethodGenOps"/.Options[stringFNmodels])=="SolveA",
{WIJaOps,O1coeffcounter}=AdjoinO1coeffWOps[(FindOpsforWIJaSolvewith\[Phi]VEVsQsector[state]),O1coefflist,1];
];
If[("MethodGenOps"/.Options[stringFNmodels])=="List",
{WIJaOps,O1coeffcounter}=AdjoinO1coeffWOps[(FindOpsforWIJaListwith\[Phi]VEVsQsector[state]),O1coefflist,1];
];

If[("MethodGenOps"/.Options[stringFNmodels])=="Optimal",
If[("nSinglet"+"nNonperts")/.Options[stringFNmodels]<=3,
{WIJaOps,O1coeffcounter}=AdjoinO1coeffWOps[FindOpsforWIJaListwith\[Phi]VEVsQsector[state],O1coefflist,1],
{WIJaOps,O1coeffcounter}=AdjoinO1coeffWOps[FindOpsforWIJaSmithwith\[Phi]VEVsQsector[state],O1coefflist,1]]
];
,
Return["Error"];
];


If[("ComputeKahlerOn"/.Options[stringFNmodels]),
{KmodOps,O1coeffcounter}=FindOpsforKmodSmithwith\[Phi]VEVs[state,O1coefflist,O1coeffcounter];
{KSMOps,O1coeffcounter}=FindOpsforKSMSmithwith\[Phi]VEVs[state,O1coefflist,O1coeffcounter];
KHiggsOps=FindOpsforKHiggsSmithwith\[Phi]VEVs[state,O1coefflist,O1coeffcounter];];

(* Step 2: Determine the \[Phi]VEVs, and substitute into the corresponding operators to obtain the valued matrices. *)
\[Phi]VEVsVal=AssocFormToVEVs[stateassocin];
WIJaVal=WIJaOps/.\[Phi]VEVsVal;
If[("ComputeKahlerOn"/.Options[stringFNmodels]),
{KmodVal,KSMVal,KHiggsVal}={KmodOps,KSMOps,KHiggsOps}/.\[Phi]VEVsVal;];
YukVal=FindYukExp[WIJaVal];
{YdVal,YuVal}=FindYukMatExpwith\[Phi]VEVsQsector[YukVal];


(* Step 3: Find the Canonical basis for the system and compute the transformation needed to find the Yukawa terms. *)
If[("ComputeKahlerOn"/.Options[stringFNmodels]),
{Qbarmat,dmat,umat,Lbarmat,emat}=CanonicalBasis[\[Phi]VEVsVal,KSMVal];
{Hdscale,Huscale}=RenormHiggs[\[Phi]VEVsVal,KHiggsVal];
Ydnew=Qbarmat . YdVal . dmat/.\[Phi]VEVsVal;
Yunew=Qbarmat . YuVal . umat/.\[Phi]VEVsVal;
YdVal=Ydnew*Hdscale;
YuVal=Yunew*Huscale;
];

(* Step 4: Run RG. First find out the largest operator along the diagonal in Yu and then run RG according to that scale. *)
{YdVal,YuVal}=RenormYuk[YdVal,YuVal];
If[("ComputeKahlerOn"/.Options[stringFNmodels]),
If[Chop[Abs[KmodVal]]==0,
{YdVal,YuVal},
YdVal=MultiplyExpKmod[\[Phi]VEVsVal,YdVal,KmodOps];
YuVal=MultiplyExpKmod[\[Phi]VEVsVal,YuVal,KmodOps];
];
];

(* Step 5: Find Higgs. First we find the Hu expected value. Depending on how the things are constructed, we want to perhaps find Hd using some rules defined in the module. *)
Qmassrules=ComputeSMQuantities[YdVal,YuVal];

(* Step 6: Detemine the fitness of the system. *)
chargematfromstate=StateVectoMatwithVEVs[state];
fitnessval=FitnessCalcQsectorOnly[Qmassrules,chargematfromstate,state];

(* Step 7: Finally convert everthing into association form. *)
HdVal=("Hd"/.Qmassrules)*(("MplinGeV"/.Options[stringFNmodels]));
HuVal=("Hu"/.Qmassrules)*(("MplinGeV"/.Options[stringFNmodels]));
tan\[Beta]Val="tan\[Beta]"/.Qmassrules;
CKMmatVal=Abs[("CKMmat"/.Qmassrules)];
conversionvec={1,1,1}*(("MplinGeV"/.Options[stringFNmodels]));
mdVal=("md"/.Qmassrules)*conversionvec;
muVal=("mu"/.Qmassrules)*conversionvec;
\[Phi]Vals=Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}]/.\[Phi]VEVsVal;
\[CapitalPhi]Vals=Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}]/.\[Phi]VEVsVal;
CKMcmplx="CKMmat"/.Qmassrules;
Jarlskogmeasured=Im[CKMcmplx[[1]][[2]]*CKMcmplx[[2]][[3]]*Conjugate[CKMcmplx[[1]][[3]]]*Conjugate[CKMcmplx[[2]][[2]]]];

YukdOps=Partition[Normal[Take[WIJaOps,9]],3]/.{\[Kappa]->"Mpl/Mcy"/.Options[stringFNmodels],\[Sigma]->1};
YukuOps=Partition[Normal[Take[WIJaOps,-9]],3]/.{\[Kappa]->"Mpl/Mcy"/.Options[stringFNmodels],\[Sigma]->1};

fitnesstestval=Join[{("FitnessHSec"/.fitnessval),("FitnessQMassSec"/.fitnessval),("FitnessQMixSec"/.fitnessval)}];

fields=Join[Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}],Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}]];
YukdSt=Partition[Flatten[Map[LowestOrderTerms[#,fields]&,Flatten[YukdOps]]],3];
YukuSt=Partition[Flatten[Map[LowestOrderTerms[#,fields]&,Flatten[YukuOps]]],3];

stateasssoc=Join[fitnessval,<|"ValueLst"->fitnesstestval,
"Hd"->HdVal,"Hu"->HuVal,"tan\[Beta]"->tan\[Beta]Val,"md"->mdVal,"mu"->muVal,"CKMmat"->CKMmatVal,"\[Phi]VEVs"->\[Phi]Vals,"\[CapitalPhi]VEVs"->\[CapitalPhi]Vals,"YukdMat"->(YukdOps//MatrixForm),"YukuMat"->(YukuOps//MatrixForm),"YukdMatSt"->(YukdSt//MatrixForm),"YukuMatSt"->(YukuSt//MatrixForm),"JarlskogInv"->Jarlskogmeasured|>];

Return[stateasssoc]
]


(* ::Subsection::Closed:: *)
(*FromAssocFormGetWOps*)


FromAssocFormGetWOps[stateassocin_]:=Module[{statevec,WIJaOps},

(* first check this is a stateassoc *)
If[Not[AssociationQ[stateassocin]],Print["Input not an association."];Return[0]];

statevec=AssocFormToStateVec[stateassocin];

(* Step 1: Find the relevant operators. *)
(* Note O1coefflist is the O1vec required. *)
(* In this step we have now used the newer modules using the Smith decomposition. *)
(* Setting to turn on what type of method we are generating the operators. *)
If[(("MethodGenOps"/.Options[stringFNmodels])=="Solve")||(("MethodGenOps"/.Options[stringFNmodels])=="Smith")||
(("MethodGenOps"/.Options[stringFNmodels])=="List")||(("MethodGenOps"/.Options[stringFNmodels])=="SolveA")||(("MethodGenOps"/.Options[stringFNmodels])=="Optimal"),

If[("MethodGenOps"/.Options[stringFNmodels])=="Solve",
WIJaOps=FindOpsforWIJawith\[Phi]VEVsQsector[statevec];
];
If[("MethodGenOps"/.Options[stringFNmodels])=="Smith",
WIJaOps=FindOpsforWIJaSmithwith\[Phi]VEVsQsector[statevec];
];
If[("MethodGenOps"/.Options[stringFNmodels])=="SolveA",
WIJaOps=FindOpsforWIJaSolvewith\[Phi]VEVsQsector[statevec];
];
If[("MethodGenOps"/.Options[stringFNmodels])=="List",
WIJaOps=FindOpsforWIJaListwith\[Phi]VEVsQsector[statevec];
];

If[("MethodGenOps"/.Options[stringFNmodels])=="Optimal",
If[("nSinglet"+"nNonperts")/.Options[stringFNmodels]<=3,
WIJaOps=FindOpsforWIJaListwith\[Phi]VEVsQsector[statevec],
WIJaOps=FindOpsforWIJaSmithwith\[Phi]VEVsQsector[statevec]]
];
,
Return["Error"];
];

WIJaOps
]


(* ::Subsection:: *)
(*CompleteStateFromAssocFormWOps*)


CompleteStateFromAssocFormWOps[stateassocin_,WOps_,O1coeffvec_]:=Module[{statevec,\[Phi]VEVsVal,O1vec,O1coeffcounter,WIJaOps,KmodOps,KSMOps,KHiggsOps,WIJaVal,KmodVal,KSMVal,KHiggsVal,YukVal,YdVal,YuVal,YeVal,
Qbarmat,dmat,umat,Lbarmat,emat,Hdscale,Huscale,Ydnew,Yunew,HiggsRules,CKMRules,Qmassrules,chargematfromstate,fitnessval,
HdVal,HuVal,tan\[Beta]Val,CKMmatVal,mdVal,muVal,\[Phi]Vals,\[CapitalPhi]Vals,stateasssoc,conversionvec,YukdOps,YukuOps,
fitnessHval,fitnessQval,fitnessGSval,fitnessCKMval,fitnessNAFval,fitnesstestval,
Jarlskogmeasured,CKMcmplx,
fields,YukdSt,YukuSt
},

(* first check this is a stateassoc *)
If[Not[AssociationQ[stateassocin]],Print["Input not an association."];Return[0]];

statevec=AssocFormToStateVec[stateassocin];

{WIJaOps,O1coeffcounter}=AdjoinO1coeffWOps[WOps,O1coeffvec,1];


If[("ComputeKahlerOn"/.Options[stringFNmodels]),
{KmodOps,O1coeffcounter}=FindOpsforKmodSmithwith\[Phi]VEVs[statevec,O1coefflist,O1coeffcounter];
{KSMOps,O1coeffcounter}=FindOpsforKSMSmithwith\[Phi]VEVs[statevec,O1coefflist,O1coeffcounter];
KHiggsOps=FindOpsforKHiggsSmithwith\[Phi]VEVs[statevec,O1coefflist,O1coeffcounter];];


(* Step 2: Determine the \[Phi]VEVs, and substitute into the corresponding operators to obtain the valued matrices. *)
\[Phi]VEVsVal=StateVectoVEVs[statevec];
WIJaVal=WIJaOps/.\[Phi]VEVsVal;
If[("ComputeKahlerOn"/.Options[stringFNmodels]),
{KmodVal,KSMVal,KHiggsVal}={KmodOps,KSMOps,KHiggsOps}/.\[Phi]VEVsVal;];
YukVal=FindYukExp[WIJaVal];
{YdVal,YuVal}=FindYukMatExpwith\[Phi]VEVsQsector[YukVal];


(* Step 3: Find the Canonical basis for the system and compute the transformation needed to find the Yukawa terms. *)
If[("ComputeKahlerOn"/.Options[stringFNmodels]),
{Qbarmat,dmat,umat,Lbarmat,emat}=CanonicalBasis[\[Phi]VEVsVal,KSMVal];
{Hdscale,Huscale}=RenormHiggs[\[Phi]VEVsVal,KHiggsVal];
Ydnew=Qbarmat . YdVal . dmat/.\[Phi]VEVsVal;
Yunew=Qbarmat . YuVal . umat/.\[Phi]VEVsVal;
YdVal=Ydnew*Hdscale;
YuVal=Yunew*Huscale;
];

(* Step 4: Run RG. First find out the largest operator along the diagonal in Yu and then run RG according to that scale. *)
{YdVal,YuVal}=RenormYuk[YdVal,YuVal];
If[("ComputeKahlerOn"/.Options[stringFNmodels]),
If[Chop[Abs[KmodVal]]==0,
{YdVal,YuVal},
YdVal=MultiplyExpKmod[\[Phi]VEVsVal,YdVal,KmodOps];
YuVal=MultiplyExpKmod[\[Phi]VEVsVal,YuVal,KmodOps];
];
];

(* Step 5: Find Higgs. First we find the Hu expected value. Depending on how the things are constructed, we want to perhaps find Hd using some rules defined in the module. *)
Qmassrules=ComputeSMQuantities[YdVal,YuVal];

(* Step 6: Detemine the fitness of the system. *)
chargematfromstate=StateVectoMatwithVEVs[statevec];
fitnessval=FitnessCalcQsectorOnly[Qmassrules,chargematfromstate,statevec];



(* Step 7: Finally convert everthing into association form. *)
HdVal=("Hd"/.Qmassrules)*(("MplinGeV"/.Options[stringFNmodels]));
HuVal=("Hu"/.Qmassrules)*(("MplinGeV"/.Options[stringFNmodels]));
tan\[Beta]Val="tan\[Beta]"/.Qmassrules;
CKMmatVal=Abs[("CKMmat"/.Qmassrules)];
conversionvec={1,1,1}*(("MplinGeV"/.Options[stringFNmodels]));
mdVal=("md"/.Qmassrules)*conversionvec;
muVal=("mu"/.Qmassrules)*conversionvec;
\[Phi]Vals=Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}]/.\[Phi]VEVsVal;
\[CapitalPhi]Vals=Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}]/.\[Phi]VEVsVal;
CKMcmplx="CKMmat"/.Qmassrules;
Jarlskogmeasured=Im[CKMcmplx[[1]][[2]]*CKMcmplx[[2]][[3]]*Conjugate[CKMcmplx[[1]][[3]]]*Conjugate[CKMcmplx[[2]][[2]]]];

YukdOps=Partition[Normal[Take[WIJaOps,9]],3]/.{\[Kappa]->"Mpl/Mcy"/.Options[stringFNmodels],\[Sigma]->1};
YukuOps=Partition[Normal[Take[WIJaOps,-9]],3]/.{\[Kappa]->"Mpl/Mcy"/.Options[stringFNmodels],\[Sigma]->1};

fitnesstestval=Join[{("FitnessHSec"/.fitnessval)},{("FitnessQMassSec"/.fitnessval)},{("FitnessQMixSec"/.fitnessval)}];

fields=Join[Table[("BundleModuliFields"/.Options[stringFNmodels])[i],{i,1,"nSinglet"/.Options[stringFNmodels]}],Table[("KahlerEffFields"/.Options[stringFNmodels])[i],{i,1,"nNonperts"/.Options[stringFNmodels]}]];
YukdSt=Partition[Flatten[Map[LowestOrderTerms[#,fields]&,Flatten[YukdOps]]],3];
YukuSt=Partition[Flatten[Map[LowestOrderTerms[#,fields]&,Flatten[YukuOps]]],3];

stateasssoc=Join[fitnessval,<|"ValueLst"->fitnesstestval,
"Hd"->HdVal,"Hu"->HuVal,"tan\[Beta]"->tan\[Beta]Val,"md"->mdVal,"mu"->muVal,"CKMmat"->CKMmatVal,"\[Phi]VEVs"->\[Phi]Vals,"\[CapitalPhi]VEVs"->\[CapitalPhi]Vals,"YukdMat"->(YukdOps//MatrixForm),
"YukuMat"->(YukuOps//MatrixForm),"YukdMatSt"->(YukdSt//MatrixForm),"YukuMatSt"->(YukuSt//MatrixForm),"JarlskogInv"->Jarlskogmeasured|>];

stateasssoc


]




(* ::Section:: *)
(*End Package*)


End[]
EndPackage[]
