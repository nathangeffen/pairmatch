require(xtable)
df000 = read.csv("output5k_ra_000.csv", comment.char = "#")
df025 = read.csv("output5k_ra_025.csv", comment.char = "#")
df050 = read.csv("output5k_ra_050.csv", comment.char = "#")
df075 = read.csv("output5k_ra_075.csv", comment.char = "#")
df100 = read.csv("output5k_ra_100.csv", comment.char = "#")

# Distance

results = df000
result_rows = results[grepl("mean distance", results$Measure),]
df = data.frame(Algorithm=result_rows$Algorithm) 
df["Mean 0.00"] = result_rows$Value
df["Ratio to best 0.00"] = result_rows$Value / min(result_rows$Value)

results = df025
result_rows = results[grepl("mean distance", results$Measure),]
df["Mean 0.25"] = result_rows$Value
df["Ratio to best 0.25"] = result_rows$Value / min(result_rows$Value)

results = df050
result_rows = results[grepl("mean distance", results$Measure),]
df["Mean 0.50"] = result_rows$Value
df["Ratio to best 0.50"] = result_rows$Value / min(result_rows$Value)

results = df075
result_rows = results[grepl("mean distance", results$Measure),]
df["Mean 0.75"] = result_rows$Value
df["Ratio to best 0.75"] = result_rows$Value / min(result_rows$Value)

results = df100
result_rows = results[grepl("mean distance", results$Measure),]
df["Mean 1.00"] = result_rows$Value
df["Ratio to best 1.00"] = result_rows$Value / min(result_rows$Value)

df = df[order(df["Ratio to best 0.50"]),]
print(df)
xtable(df)

# Rankings

results = df000
result_rows = results[grepl("mean rank", results$Measure),]
df = data.frame(Algorithm=result_rows$Algorithm) 
df["Mean 0.00"] = result_rows$Value
df["Ratio to best 0.00"] = result_rows$Value / min(result_rows$Value)

results = df025
result_rows = results[grepl("mean rank", results$Measure),]
df["Mean 0.25"] = result_rows$Value
df["Ratio to best 0.25"] = result_rows$Value / min(result_rows$Value)

results = df050
result_rows = results[grepl("mean rank", results$Measure),]
df["Mean 0.50"] = result_rows$Value
df["Ratio to best 0.50"] = result_rows$Value / min(result_rows$Value)

results = df075
result_rows = results[grepl("mean rank", results$Measure),]
df["Mean 0.75"] = result_rows$Value
df["Ratio to best 0.75"] = result_rows$Value / min(result_rows$Value)

results = df100
result_rows = results[grepl("mean rank", results$Measure),]
df["Mean 1.00"] = result_rows$Value
df["Ratio to best 1.00"] = result_rows$Value / min(result_rows$Value)

df = df[order(df["Ratio to best 0.50"]),]
print(df)
xtable(df)