import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_vcf(file_path):
    """Reads a VCF file into a Pandas DataFrame.

    Args:
        file_path (str): The path to the VCF file.

    Returns:
        pd.DataFrame: A Pandas DataFrame containing the VCF data.
    """
    with open(file_path, 'r') as f:
        header_lines = 0
        for line in f:
            if line.startswith('##'):
                header_lines += 1
            elif line.startswith('#CHROM'):
                column_names = line.strip().split('\t')
                break

    df = pd.read_csv(file_path, sep='\t', skiprows=header_lines, names=column_names)
    return df

# create pandas df of each filtered vcf
vcf_LM_df = read_vcf('LolMulSortedFilter.8.025.vcf')
vcf_LP_df = read_vcf('LolPeSortedFilter.8.025.vcf')
vcf_LR_df = read_vcf('LolRiSortedFilter.8.025.vcf')

# manually check columns
print('columns', vcf_LR_df.columns)

# extract just the sample names
vcf_LM_sample_names = vcf_LM_df.columns.to_list()[9:]
vcf_LP_sample_names = vcf_LP_df.columns.to_list()[9:]
vcf_LR_sample_names = vcf_LR_df.columns.to_list()[9:]

# print how many reads we have in each alignment
print('num readings LM', len(vcf_LM_df))
print('num readings LP', len(vcf_LP_df))
print('num readings LR', len(vcf_LR_df))

# find overlap in what readings aligned to each ref - uses ID - not super sure if this is best way to think abt this?
print('overlap ID of readings LM and LP', len(set(vcf_LM_df['ID'].to_list()[1:]).intersection(set(vcf_LP_df['ID'].to_list()[1:]))))
print('overlap ID of readings LM and LR', len(set(vcf_LM_df['ID'].to_list()[1:]).intersection(set(vcf_LR_df['ID'].to_list()[1:]))))
print('overlap ID of readings LP and LR', len(set(vcf_LP_df['ID'].to_list()[1:]).intersection(set(vcf_LR_df['ID'].to_list()[1:]))))

#make tuples of (chrom, pos)
chrom_pos_tuples_LP = vcf_LP_df.apply(lambda row: (row['#CHROM'], row['POS']), axis=1).tolist()[1:]
chrom_pos_tuples_LR = vcf_LR_df.apply(lambda row: (row['#CHROM'], row['POS']), axis=1).tolist()[1:]
chrom_pos_tuples_LM = vcf_LM_df.apply(lambda row: (row['#CHROM'], row['POS']), axis=1).tolist()[1:]
# find overlap in (chrom, pos) between diff alignments
print('overlap in (chrom, pos) pairs between LM and LP', len(set(chrom_pos_tuples_LM).intersection(set(chrom_pos_tuples_LP))))
print('overlap in (chrom, pos) pairs between LM and LR', len(set(chrom_pos_tuples_LM).intersection(set(chrom_pos_tuples_LR))))
print('overlap in (chrom, pos) pairs between LP and LR', len(set(chrom_pos_tuples_LP).intersection(set(chrom_pos_tuples_LR))))

# print how many plant samples (cloumns of plants) there are in each alignment
print('num plant samples LM', len(vcf_LM_sample_names))
print('num plant samples LP', len(vcf_LP_sample_names))
print('num plant samples LR', len(vcf_LR_sample_names))

# find overlap in what plant samples alligned to each ref
print('overlap plant samples LM and LP', len(set(vcf_LM_sample_names).intersection(set(vcf_LP_sample_names))))
print('overlap plant samples LM and LR', len(set(vcf_LM_sample_names).intersection(set(vcf_LR_sample_names))))
print('overlap plant samples LP and LR', len(set(vcf_LP_sample_names).intersection(set(vcf_LR_sample_names))))

print('samples found in LP and LR but not LM', set(vcf_LP_sample_names).difference(set(vcf_LM_sample_names)))

# make pandas df out of qc metrics
qc_df = pd.read_csv("GBS summary(Sheet1).csv")

# figuring out format / how to extract information
# these lines of code extract column values of a given qc metric if the sample aligned (aka if the sample name can be found in the list of vcf sample names)
fastqc_mean_quality_read_aligned_to_LP = qc_df.loc[qc_df['general-rawsamplename'].isin(vcf_LP_sample_names), 'fastqc-meanreadqualityR1']
print(fastqc_mean_quality_read_aligned_to_LP)
print(sum(fastqc_mean_quality_read_aligned_to_LP)/len(fastqc_mean_quality_read_aligned_to_LP))

# the addition of the ~ makes this do the opposite - extracts information on what does not align
fastqc_mean_quality_read_not_aligned_to_LP = qc_df.loc[~qc_df['general-rawsamplename'].isin(vcf_LP_sample_names), 'fastqc-meanreadqualityR1']
print(fastqc_mean_quality_read_not_aligned_to_LP)
print(sum(fastqc_mean_quality_read_not_aligned_to_LP)/len(fastqc_mean_quality_read_not_aligned_to_LP))


'''
version 1 of making graphs
this made 1 graph per qc metric - with 6 columns (2 per vcf)
however the graphs weren't super descriptive since many of the bars ended up around same height
so i switched how i wanted to visualize this info (seen further down in code)
'''

# all_alignment_sample_names = [vcf_LM_sample_names, vcf_LP_sample_names, vcf_LR_sample_names]

# for column in columns_to_compare:
#     y_values_graph = []
#     for sample_names in all_alignment_sample_names:
#         fastqc_quality_aligned = qc_df.loc[qc_df['general-rawsamplename'].isin(sample_names), column]
#         y_values_graph.append(sum(fastqc_quality_aligned)/len(fastqc_quality_aligned))
#         fastqc_quality_not_aligned = qc_df.loc[~qc_df['general-rawsamplename'].isin(sample_names), column]
#         y_values_graph.append(sum(fastqc_quality_not_aligned)/len(fastqc_quality_not_aligned))
    
#     x_values_graph = ['LM aligned', 'LM didn\'t align', 'LP aligned', 'LP didn\'t align', 'LR aligned', 'LR didn\'t align']
#     plt.bar(x_values_graph, y_values_graph)
#     plt.title('average ' + column + ' across different alignments')
#     plt.show()

# x_values_graph = ['LM aligned', 'LM didn\'t align', 'LP aligned', 'LP didn\'t align', 'LR aligned', 'LR didn\'t align']
# plt.bar(columns_to_compare, y_values_graph)
# plt.xticks(rotation='vertical')
# plt.tight_layout()
# plt.title('quality control difference between samples that aligned and those that didn\'t ')
# plt.show()


'''
version 2
this makes one bar chart per vcf (i manually changed lines 125, 135, 151, & 168 to make different graphs)
each plots the difference in several qc metrics of what did and didn't align to
'''
# list of columns (found in GBS summary) to compare
columns_to_compare = ['fastqc-meanreadqualityR1','fastqc-pct-adapter','fastqc-pct-deduplicated','fastqc-pct-dimer','fastqc-pct-gc','fastqc-sequencelength','dedupmetric-PCRDuplication5M','dedupmetric-PCRduplication','dedupmetric-TileDuplication5M','dedupmetric-Tileduplication','fastqspeciesscreen-pct-Adapters','fastqspeciesscreen-pct-Chloroplast','fastqspeciesscreen-pct-Mitochondria','fastqspeciesscreen-pct-NoHit','fastqspeciesscreen-pct-PhiX','fastqspeciesscreen-pct-Vectors','fastqspeciesscreen-pct-ecoli','fastqspeciesscreen-pct-human','fastqspeciesscreen-pct-mouse','fastqspeciesscreen-pct-rRNA']
y_values_graph = []
for column in columns_to_compare:
    # find and extract values of those that aligned
    fastqc_quality_aligned = qc_df.loc[qc_df['general-rawsamplename'].isin(vcf_LR_sample_names), column]
    # sum - manually ignore NaNs since they're in unorthodox format
    sum_aligned = 0
    for entry in fastqc_quality_aligned:
        if entry != 'undef':
            sum_aligned += float(entry)
    # take the average
    avg_align = sum_aligned/len(fastqc_quality_aligned)

    # repeat for those that did not align
    fastqc_quality_not_aligned = qc_df.loc[~qc_df['general-rawsamplename'].isin(vcf_LR_sample_names), column]
    sum_not_aligned = 0
    for entry in fastqc_quality_not_aligned:
        if entry != 'undef':
            sum_not_aligned += float(entry)
    avg_didnt_align = sum_not_aligned/len(fastqc_quality_not_aligned)

    # calculate and store the average
    y_values_graph.append(round(avg_align - avg_didnt_align, 3))

# create bar chart out of above information 
# https://stackoverflow.com/questions/28931224/how-to-add-value-labels-on-a-bar-chart
plt.figure(figsize=(12, 9))
plt.tight_layout()
y_series = pd.Series(y_values_graph)
ax = y_series.plot(kind="bar")
ax.set_title("difference in quality control metrics between samples that did and did not align to L. rigidum")
ax.set_xlabel("fastqc metrics")
ax.set_ylabel("difference (aligned - not aligned) in average")
ax.set_xticklabels(columns_to_compare)
rects = ax.patches
labels = y_values_graph

# add labels to bar chart so easier to understand - code from the above stack overflow post ! 
for rect, label in zip(rects, labels):
    height = rect.get_height()
    adjust = .1 if abs(height) == height else -1.3
    ax.text(
        rect.get_x() + rect.get_width() / 2, height + adjust, label, ha="center", va="bottom"
    )

plt.subplots_adjust(bottom=0.4) # made the Long x axis labels readable
plt.savefig('LR_bar_chart.png', dpi=300)
plt.show()



# messed around with these but popmap.txt not really helpful for what I was trying to think about
# and SampleInfo would've been helpful but the naming convention is different 
# and my attempts to match up the names were unsuccessful (not shown)

sampleinfo_df = pd.read_csv('SampleInfo.txt', sep='\s+', header=None)
sampleinfo_df.columns = ['sample name', 'info']
print(sampleinfo_df)
sampleinfo_aligned_to_LP = sampleinfo_df.loc[sampleinfo_df['sample name'].isin(vcf_LP_sample_names), 'info']
print(sampleinfo_aligned_to_LP)
print(len(sampleinfo_aligned_to_LP))

popmap_df = pd.read_csv('popmap.txt', sep='\t', header=None)
popmap_df.columns = ['sample name', 'pop']
popmap_aligned_to_LP = popmap_df.loc[popmap_df['sample name'].isin(vcf_LP_sample_names), 'pop']
print(popmap_aligned_to_LP)
print(len(popmap_aligned_to_LP))




'''
same logic to above graph but now compare the 4 that did not align to LM vs everything else
'''
# list of columns (found in GBS summary) to compare
columns_to_compare = ['fastqc-meanreadqualityR1','fastqc-pct-adapter','fastqc-pct-deduplicated','fastqc-pct-dimer','fastqc-pct-gc','fastqc-sequencelength','dedupmetric-PCRDuplication5M','dedupmetric-PCRduplication','dedupmetric-TileDuplication5M','dedupmetric-Tileduplication','fastqspeciesscreen-pct-Adapters','fastqspeciesscreen-pct-Chloroplast','fastqspeciesscreen-pct-Mitochondria','fastqspeciesscreen-pct-NoHit','fastqspeciesscreen-pct-PhiX','fastqspeciesscreen-pct-Vectors','fastqspeciesscreen-pct-ecoli','fastqspeciesscreen-pct-human','fastqspeciesscreen-pct-mouse','fastqspeciesscreen-pct-rRNA']
y_values_graph = []
the_outliers = ['VLR691', 'VLR692', 'VLR695', 'VLR696']
for column in columns_to_compare:
    # find and extract values of those that aligned
    fastqc_quality_aligned = qc_df.loc[qc_df['general-rawsamplename'].isin(the_outliers), column]
    # sum - manually ignore NaNs since they're in unorthodox format
    sum_aligned = 0
    for entry in fastqc_quality_aligned:
        if entry != 'undef':
            sum_aligned += float(entry)
    # take the average
    avg_align = sum_aligned/len(fastqc_quality_aligned)

    # repeat for those that did not align
    fastqc_quality_not_aligned = qc_df.loc[~qc_df['general-rawsamplename'].isin(the_outliers), column]
    sum_not_aligned = 0
    for entry in fastqc_quality_not_aligned:
        if entry != 'undef':
            sum_not_aligned += float(entry)
    avg_didnt_align = sum_not_aligned/len(fastqc_quality_not_aligned)

    # calculate and store the average
    y_values_graph.append(round(avg_align - avg_didnt_align, 3))

# create bar chart out of above information 
# https://stackoverflow.com/questions/28931224/how-to-add-value-labels-on-a-bar-chart
plt.figure(figsize=(12, 9))
plt.tight_layout()
y_series = pd.Series(y_values_graph)
ax = y_series.plot(kind="bar")
ax.set_title("difference in quality control metrics between VLR69 samples and everything else")
ax.set_xlabel("fastqc metrics")
ax.set_ylabel("difference (VLR69 - everything else) in average")
ax.set_xticklabels(columns_to_compare)
rects = ax.patches
labels = y_values_graph

# add labels to bar chart so easier to understand - code from the above stack overflow post ! 
for rect, label in zip(rects, labels):
    height = rect.get_height()
    adjust = .1 if abs(height) == height else -1.6
    ax.text(
        rect.get_x() + rect.get_width() / 2, height + adjust, label, ha="center", va="bottom"
    )

plt.subplots_adjust(bottom=0.4) # made the Long x axis labels readable
plt.savefig('outlier_bar_chart.png', dpi=500)
plt.show()