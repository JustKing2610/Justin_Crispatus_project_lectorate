import pandas as pd

tsv_file = "/mnt/StudentFiles/2023/Justin/mlst_crispatus/tsv_MLST_scheme/cgMLST_profiles.tsv"

df = pd.read_csv(tsv_file, delimiter='\t')

dict_GCA = {}

#for each row, take the GCA ID from the column named "file" and join the values of the row together with "-"

for index, row in df.iterrows():
    GCA_ID = row["FILE"]

    row_values = list(row[1:])

    row_string = '-'.join(map(str, row_values))

    dict_GCA[GCA_ID] = row_string

results_df = pd.DataFrame(list(dict_GCA.items()), columns=["GCA", 'Allele_values'])


results_df.to_excel("/mnt/StudentFiles/2023/Justin/mlst_crispatus/output_python_mlst_scheme/output_dict_Alleleschema.xlsx", index=False)



values = list(dict_GCA.values())

labels = pd.factorize(values)[0] + 1

ST_types_label = [f'ST{label}' for label in labels]

final_dict = dict(zip(dict_GCA.keys(), ST_types_label))


print("Check final dictonary for GCA numbers and ST types:", final_dict)

ST_dataframe = pd.DataFrame(list(final_dict.items()), columns=["Identifier", "Sequence_types"])

ST_dataframe.to_excel("/mnt/StudentFiles/2023/Justin/mlst_crispatus/output_python_mlst_scheme/ST_types_Crispatus_chewbacca.xlsx", index=False)
