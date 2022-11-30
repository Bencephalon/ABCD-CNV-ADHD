import pandas as pd
import numpy as np
import re
from tqdm import tqdm
from glob import glob
import os
import random

wd = "/data/NCR_SBRB/jungbt/CNV/ABCD"
file_data = f"{wd}/data/assessments/orig/big_data_abridged.csv"

def open_abcd(file):
    """
    Open ABCD file and extract descriptive info
    """
    df = pd.read_csv(file,sep="\t",low_memory = False)
    description = pd.DataFrame({"description":df.iloc[0,:]},index=df.columns)
    df = df.iloc[1:,:]
    return df, description

#Load KSAD
data_df = pd.read_csv(file_data)
data_df.ID = [x.replace("NDAR","NDAR_") for x in data_df.ID]
#Drop irrelevant columns
data_df = data_df.loc[:,~data_df.isna().all()]
data_df = data_df.drop([x for x in data_df.columns if "run" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "raw" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "Unnamed" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "Scan" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "Past" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "Present" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "Diagnosis" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "Easily" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if "Difficulty" in x],axis=1)
data_df = data_df.drop([x for x in data_df.columns if x.startswith("N_")],axis=1)
data_df = data_df.drop(["Unique ID provided by KSADS for dataset","Session.2","ID.1","sx_source","Sex of the subject","ever_ADHD.1","Session","Session.1","alt_ID","longitudinal_subject","matched_group","Age_round","siblings_twins","The event name for which the data was collected.1","The NDAR Global Unique Identifier (GUID) for research subject","The NDAR Global Unique Identifier (GUID) for research subject.1"],axis=1,errors="ignore")


#Add neurocognitive measures
measures = dict()
measures["abcd_tbss01.txt"] = ["nihtbx_picvocab_fc","nihtbx_flanker_fc","nihtbx_list_fc","nihtbx_cardsort_fc","nihtbx_pattern_fc","nihtbx_picture_fc","nihtbx_reading_fc","nihtbx_fluidcomp_fc","nihtbx_cryst_fc","nihtbx_totalcomp_fc",
                               "nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_list_uncorrected","nihtbx_cardsort_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected","nihtbx_reading_uncorrected","nihtbx_fluidcomp_uncorrected","nihtbx_cryst_uncorrected","nihtbx_totalcomp_uncorrected"]
measures["abcd_ps01.txt"] = ["pea_wiscv_trs","pea_wiscv_tss","RAVLT_items_learned","RAVLT_immediate_recall","RAVLT_delayed_recall"]
measures["abcd_gdss01.txt"] = ["gdt_scr_values_risky"]



measures_flat = [item for sublist in list(measures.values()) for item in sublist]


#Add neurocog measures
for file in measures.keys():
    #Load neurocog file
    neurocog_df, _ = open_abcd(f"data/assessments/orig/{file}")
    cols = ['src_subject_id',"eventname"]
    if file == "abcd_ps01.txt":
        neurocog_df["RAVLT_items_learned"] = neurocog_df.pea_ravlt_sd_trial_v_tc
        neurocog_df["RAVLT_immediate_recall"] = neurocog_df.pea_ravlt_sd_trial_vi_tc
        neurocog_df["RAVLT_delayed_recall"] = neurocog_df.pea_ravlt_ld_trial_vii_tc
    elif file == "abcd_gdss01.txt":
        #Only available year 2, so remove the event criteria
        cols = ['src_subject_id']
    neurocog_df = neurocog_df.loc[:,cols + measures[file]]
    neurocog_df[measures[file]] = neurocog_df[measures[file]].apply(pd.to_numeric, errors='coerce')
    #Merge neurocog into data
    data_df = data_df.merge(neurocog_df, left_on=['ID',"The event name for which the data was collected"], right_on=['src_subject_id',"eventname"],how="left",suffixes = ["_x",None])

data_df.eventname = data_df["The event name for which the data was collected"]

#Make a dataframe containing all subjects with CNV probe data
sub_df = pd.DataFrame({"File":glob(f"{wd}/data/penncnv/*.txt")})
#Add Genomic information
sub_df["ID"] = sub_df.File.apply(lambda x: "_".join(x.split("_")[-2:]).replace(".txt",""))
sub_df["ID2"] = sub_df.File.apply(lambda x: "_".join(x.split("_")[-3:-2]))
sub_df["BATCH"] = sub_df.File.apply(lambda x: "_".join(os.path.basename(x).split("_")[0:2]))
sub_df["METHOD"] = sub_df.File.apply(lambda x: x.split("_")[3])
#Recode index
sub_df.index = sub_df["ID"]

#Load PLINK fam file
plink_df = pd.read_csv(f"{wd}/data/PLINK/ABCD_release_3.0_QCed.filtered.fam",sep=" ")
plink_df.columns = ["FID","IID","PID","MID","SEX","PHENO"]
plink_df.index = plink_df.IID

#Convert plink sex to M/F
sex_dict = {1:"M",2:"F"}
plink_df["SEX"] = plink_df.SEX.apply(lambda x: sex_dict[x])

#Define values to store subject-level sex and family ID
sub_df["sex"] = np.nan
sub_df["rel_family_id"] = np.nan

# Assign sex for each subject
for ind in sub_df.index:
    try:
        #Assign sex to each subject
        sub_df.loc[ind,"sex"]           = plink_df.loc[sub_df.loc[ind,"ID"],"SEX"]
        #Assign family ID to each subject
        sub_df.loc[ind,"rel_family_id"] = plink_df.loc[sub_df.loc[ind,"ID"],"FID"]
    except KeyError:
        print(f"{sub_df.loc[ind,'ID']} not found in PLINK file")

#Add maximum scores for ADHD and matrix reasoning scores across all time points
max_adhd_df = data_df[["ID","sx_inatt_present","sx_hi_present","sx_total"] + measures_flat].groupby(by="ID").max()
max_adhd_df.columns = [f"{x}_max" for x in max_adhd_df.columns]
sub_df = pd.concat([sub_df,max_adhd_df],axis=1)

#Add baseline scores
baseline_df = data_df.loc[data_df.eventname == "baseline_year_1_arm_1",:]
base_cols = [f"{x}_base" for x in ["sx_inatt_present","sx_hi_present","sx_total"] + measures_flat]
for sub in tqdm(sub_df.index):
    try:
        sub_df.loc[sub,base_cols] = list(baseline_df.loc[baseline_df.ID == sub,["sx_inatt_present","sx_hi_present","sx_total"] + measures_flat].iloc[0])
    except IndexError:
        pass


# Load KSADS-COMP
ksad_df, _ = open_abcd(f"{wd}/data/assessments/orig/abcd_ksad01.txt")

# Load ADHD raw questions
adhd_raw, _ = open_abcd(f"{wd}/data/assessments/orig/attn_deficit_hyperactiv_p01.txt")

#Merge KSADS data
ksad_df = ksad_df.merge(adhd_raw,on = ["subjectkey","eventname"], how = "left")

# Load KSADS-COMP child
ksad_child_df, _ = open_abcd(f"{wd}/data/assessments/orig/abcd_ksad501.txt")


# Relevant codes:
    # Current ADHD = ksads_14_853_p
    # ADHD-NOS    = ksads_14_856_p
    # ADHD in partial remissions = ksads_14_855_p
    # past ADHD = ksads_14_854_p (at baseline)
    # Difficulties = ksads_adhd_raw_990_p , ksads_adhd_raw_991_p, ksads_adhd_raw_992_p, ksads_adhd_raw_994_p
    # Difficulties = ksads_adhd_raw_1012_p , ksads_adhd_raw_1013_p, ksads_adhd_raw_1014_p, ksads_adhd_raw_1016_p

################################ ADHD DIAGNOSIS ################################
#### STEP 1: Make a simplified defition of ADHD
#Define subjects with ADHD as any of the following:
#   current ADHD, ADHD-NOS, ADHD in partial remission
ksad_df["ADHD_base"] = (ksad_df.ksads_14_853_p.isin([1,"1"])) | (ksad_df.ksads_14_856_p.isin([1,"1"])) | (ksad_df.ksads_14_855_p.isin([1,"1"]))

#### STEP 2: Generate list of ADHD cases with impairment in two or more settings
#Calculate if difficulties are present in two or more settings - present
difficulty_codes_present = ["ksads_adhd_raw_990_p","ksads_adhd_raw_991_p","ksads_adhd_raw_992_p","ksads_adhd_raw_993_p","ksads_adhd_raw_994_p"]
ksad_df["Difficulties_Present"] = ksad_df[difficulty_codes_present].astype(float).sum(axis=1) >= 2
#Define ever_ADHD subjects: Has ADHD and present difficulties
ksad_df["ADHD"] = False
ksad_df.loc[ksad_df["ADHD_base"] & ksad_df["Difficulties_Present"],"ADHD"] = True
adhd_cases = list(ksad_df.loc[ksad_df["ADHD_base"] & ksad_df["Difficulties_Present"],"subjectkey"].unique())

#### STEP 3: Define past only ADHD
#Calculate if difficulties are present in two or more settings - past
difficulty_codes_past = ["ksads_adhd_raw_1012_p","ksads_adhd_raw_1013_p","ksads_adhd_raw_1014_p","ksads_adhd_raw_1015_p","ksads_adhd_raw_1016_p"]
ksad_df["Difficulties_Past"] = ksad_df[difficulty_codes_past].astype(float).sum(axis=1) >= 2
#Define past-only ADHD: Positive past ADHD dx with no current ADHD dx
ksad_df["ADHD_past_only"] = ksad_df.ksads_14_854_p.isin([1,"1"]) & ~ksad_df.subjectkey.isin(adhd_cases)

#### STEP 4: Generate list of past only ADHD not taking medication and ADHD w/o impairment
#Find subjects taking ADHD meds during baseline
adhd_meds = data_df.loc[data_df[["stimulants_during_study","other_ADHD_meds_during_study"]].any(axis=1) & (ksad_df.eventname == "baseline_year_1_arm_1"),"ID"].unique()
#Exclude past-only ADHD with no meds (don't take into account interference)
adhd_exclude = list(ksad_df.loc[~ksad_df.subjectkey.isin(adhd_meds) & ksad_df.ADHD_past_only & (ksad_df.eventname == "baseline_year_1_arm_1"),"subjectkey"].unique())
#Exclude past-only ADHD with meds but no interference
adhd_exclude += list(ksad_df.loc[ksad_df["ADHD_past_only"] & ksad_df.subjectkey.isin(adhd_meds) & ~ksad_df["Difficulties_Past"],"subjectkey"].unique())
#Exclude ADHD cases without impairment
adhd_exclude += [x for x in ksad_df.loc[ksad_df.ADHD_base == True,"subjectkey"].unique() if not x in adhd_cases]
adhd_exclude = set(adhd_exclude)

#### STEP 5: Generate list of past only ADHD taking medication
#Define past-only ADHD w/ meds and impairment
ksad_df.loc[ksad_df["ADHD_past_only"] & ksad_df.subjectkey.isin(adhd_meds) & ksad_df["Difficulties_Past"],"ADHD"] = True
past_only_meds = list(ksad_df.loc[ksad_df["ADHD_past_only"] & ksad_df.subjectkey.isin(adhd_meds) & ksad_df["Difficulties_Past"],"subjectkey"].unique())
adhd_cases += past_only_meds

#### STEP 6: Define ADHD Cases and Controls
sub_df["ever_ADHD"] = ""
#Define ADHD from list of cases
sub_df.loc[sub_df.ID.isin(adhd_cases),"ever_ADHD"] = "ADHD"
#Exclude those lacking impairment and past-only with no meds
sub_df.loc[sub_df.ID.isin(adhd_exclude),"ever_ADHD"] = "EXCLUDE_Ambiguous_ADHD"
#Define control subjects
sub_df.loc[sub_df.ever_ADHD == "","ever_ADHD"] = "non-ADHD"

#### STEP 7: Incorporate CBCL T-Scores to exclude controls
#Exclude control subjects with subthreshold scores
high_cbcl_tscores = data_df.loc[data_df["ADHD CBCL DSM5 Scale (t-score)"] >= 65,"ID"].unique()
high_sx_counts    = data_df.loc[data_df.sx_total >= 5,"ID"].unique()
sub_df.loc[(sub_df.ever_ADHD == "non-ADHD") & (sub_df.ID.isin(high_cbcl_tscores) | sub_df.ID.isin(high_sx_counts)),"ever_ADHD"] = "EXCLUDE_Subthreshold"

#### STEP 8: Exclude psychosis and bipolar
#Psychosis
ever_psychosis = ksad_df.loc[ksad_df.ksads_4_851_p.isin([1,"1"]),"subjectkey"].unique()
sub_df.loc[sub_df.ID.isin(ever_psychosis),"ever_ADHD"] = "EXCLUDE_Psychosis"
#Create psychosis phenotype
sub_df["ever_psychosis"] = 0
sub_df.loc[sub_df.ID.isin(ever_psychosis),"ever_psychosis"] = 1
#Bipolar
#Calculate depression dx - Adult
depression_cols = ["ksads_1_840_p","ksads_1_841_p","ksads_1_842_p","ksads_1_843_p","ksads_1_844_p" ,"ksads_1_845_p","ksads_1_846_p","ksads_1_847_p"]
ksad_df["depression"] = (ksad_df.loc[:,depression_cols] == 1).any(axis=1)
#Calculate bipolar-1
ever_bipolar = list(ksad_df.loc[ksad_df.ksads_2_830_p.isin([1,"1"]) | (ksad_df.ksads_2_833_p.isin([1,"1"]) & ksad_df["depression"]),"subjectkey"].unique())
#Calculate depression dx - Child
depression_cols = ["ksads_1_840_t","ksads_1_841_t","ksads_1_842_t","ksads_1_843_t","ksads_1_844_t" ,"ksads_1_845_t","ksads_1_846_t","ksads_1_847_t"]
ksad_child_df["depression"] = (ksad_child_df.loc[:,depression_cols] == 1).any(axis=1)
#Add child assessments of bipolar-1 to exclusion list
ever_bipolar += list(ksad_child_df.loc[ksad_child_df.ksads_2_830_t.isin([1,"1"]) | (ksad_child_df.ksads_2_833_t.isin([1,"1"]) & ksad_child_df["depression"]),"subjectkey"].unique())
#Exclude subjects with bipolar-1
sub_df.loc[sub_df.ID.isin(ever_bipolar),"ever_ADHD"] = "EXCLUDE_Bipolar"
#Create bipolar phenotype
sub_df["ever_bipolar"] = 0
sub_df.loc[sub_df.ID.isin(ever_bipolar),"ever_bipolar"] = 1



# Exclude all subjects from plate 461
exclude_df = pd.read_csv(f"{wd}/data/assessments/exclude_plate_461.csv",names=["ID2","ID","plate","BATCH"])
sub_df.loc[sub_df.ID.isin(exclude_df.ID),"ever_ADHD"]  = "EXCLUDE_461"

#Extract current age (i.e., age of last assessment)
ksad_df.interview_age_x = ksad_df.interview_age_x.astype(int)
max_ages = ksad_df.loc[ksad_df.sort_values(['interview_age_x']).drop_duplicates(['subjectkey'],keep='last').index,["subjectkey","interview_age_x"]]
max_ages.index = max_ages.subjectkey
max_ages.columns = ["subjectkey","interview_age_max"]
sub_df = sub_df.merge(max_ages.loc[:,"interview_age_max"],left_index = True, right_index = True,how="left")

#Extract age of first DX (i.e., study timepoint with first ADHD dx)
ksad_df_adhd = ksad_df.loc[ksad_df.ADHD,:]
min_ADHD_ages = ksad_df_adhd.loc[ksad_df_adhd.sort_values(['interview_age_x']).drop_duplicates(['subjectkey'],keep='last').index,["subjectkey","interview_age_x"]]
min_ADHD_ages.index = min_ADHD_ages.subjectkey
min_ADHD_ages.columns = ["subjectkey","interview_age_ADHD_min"]
sub_df = sub_df.merge(min_ADHD_ages.loc[:,"interview_age_ADHD_min"],left_index = True, right_index = True,how="left")

#Load Family and race/ethnicity data
fam_df, fam_description = open_abcd(f"{wd}/data/assessments/orig/acspsw03.txt")
#Only interested in baseline
fam_df = fam_df.loc[fam_df.eventname == "baseline_year_1_arm_1",:]

#Convert to numeric
for col in fam_df.columns:
    fam_df[col] = pd.to_numeric(fam_df[col],errors="ignore")


fam_df.index = fam_df.subjectkey
#Add race
sub_df["race"] = fam_df.race_ethnicity



#Add CBCL Scores (Max; T-Scores only)
max_cbcl_df = data_df[["ID"] + [x for x in data_df.columns if "CBCL" in x]].groupby(by="ID").max()
max_cbcl_df.columns = [x.replace("(t-score)","t_max").replace(" ","_") for x in max_cbcl_df.columns]
sub_df = pd.concat([sub_df,max_cbcl_df],axis=1)

# Add Ever DX of Comorbidities
ksad_df, _ = open_abcd(f"{wd}/data/assessments/orig/abcd_ksad01.txt")

comorbid_dict = {"ODD":["ksads_15_901_p"],
                 "CD":["ksads_16_897_p","ksads_16_898_p"],
                 "MDD":["ksads_1_840_p","ksads_1_841_p"],
                 "GAD":["ksads_10_869_p"],
                 "OCD":["ksads_11_917_p"]}

flat_list = [item for sublist in list(comorbid_dict.values()) for item in sublist]
ksad_df[flat_list] = ksad_df[flat_list].astype(int).replace({555:0, 888:0})
#Loop through each diagnosis and define status for each subject
for dx in tqdm(list(comorbid_dict.keys())):
    #Find ever dx for a given timepoint
    dx_max = ksad_df[comorbid_dict[dx]].astype(int).max(axis=1)
    # Assume no DX for codes 555 and 888 in ABCD
    dx_max[dx_max == 555] = 0
    dx_max[dx_max == 888] = 0
    #Create dataframe
    dx_df = pd.DataFrame({"subjectkey": ksad_df["subjectkey"], "dx": dx_max})
    #Find max dx score
    dx_df_max = (dx_df.groupby(['subjectkey'],as_index=False)['dx'].max()
                    .groupby('subjectkey')['dx'].agg(['sum']))
    sub_df[dx] = 0
    #Assign each subject the relevant dx
    for n in sub_df.index:
        try:
            sub_df.loc[n,dx] = dx_df_max.loc[sub_df.loc[n,"ID"],"sum"] > 0
        except KeyError:
            pass

# Define broad comorbidities
sub_df.loc[:,"externalizing"] = (sub_df.ODD == 1) | (sub_df.CD == 1)
sub_df.loc[:,"anxiety"]       = sub_df.GAD == 1
sub_df.loc[:,"mood"]          = sub_df.MDD == 1
sub_df.loc[:,"OCD"]          = sub_df.OCD == 1

#Save subject data
sub_df.to_csv(f"{wd}/data/assessments/subject_data.csv")

########################## PFB File Generation Prep ##########################
#Remove subjects lacking SNP files
sub_df = sub_df.loc[~sub_df.File.isnull(),:]
#Save a random family member from each family for PFB generation
rand_fam = []
for fam in tqdm(sub_df.rel_family_id.unique()):
    #Subset the df to extract the family
    tmp_fam_df = sub_df.loc[sub_df.rel_family_id == fam,:].copy()
    #Only subjects with PennCNV files
    tmp_fam_df = tmp_fam_df.loc[~(tmp_fam_df.File == ""),:]
    #Add a random family member to the list (if one exists)
    if tmp_fam_df.shape[0] > 0:
        rand_fam.append(tmp_fam_df.sample().File.iloc[0])

#Randomize entries in case order matters
random.shuffle(rand_fam)

# Save the random family members for PFB generation
with open(f"{wd}/data/model/pfb_subjects.txt", 'w') as f:
    for item in rand_fam:
        f.write("%s\n" % item)
