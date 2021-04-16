
import statistics
import sys
import matplotlib
import matplotlib.pyplot as plt

infile_samples = open(sys.argv[1])
infile_coverage = open(sys.argv[2])
infile_normal_samples = open(sys.argv[3])
infile_normal_coverage = open(sys.argv[4])
infile_normal_cnv1 = open(sys.argv[5])
outfile_cnv = open(sys.argv[6], "w")
outplotsdir = sys.argv[7]

#infile_normal_samples = open("/home/jonas/Investigations/CNV_WP3/samples_TE6_37.txt")
#infile_normal_coverage = open("/home/jonas/Investigations/CNV_WP3/TE6_37_all_panels.cov")
#infile_normal_cnv1 = open("/home/jonas/Investigations/CNV_WP3/Known_CNV_data/WW_25m_CNV_ChAS3.0.aed")
#infile_normal_samples = open("/projects/wp3/nobackup/TWIST/OUTBOX/TE41_210122/Jonas_CNV/samples_TE6_37.txt")
#infile_normal_coverage = open("/projects/wp3/nobackup/TWIST/OUTBOX/TE41_210122/Jonas_CNV/TE6_37_all_panels.cov")
#infile_normal_cnv1 = open("/projects/wp3/nobackup/TWIST/OUTBOX/TE41_210122/Jonas_CNV/WW_25m_CNV_ChAS3.0.aed")
#outfile_cnv = open("/home/jonas/Investigations/CNV_WP3/TE6_37_all_panels_cnv.csv", "w")

outfile_cnv.write("CNV_type\trun\tsample\tchrom\tstart\tend\texons\tnr_of_exons\tcoverage\tgene_corrected_coverage\tmedian_coverage\tnr_std\tgene_corrected_nr_std\tCN\tgene_corrected_CN\tgene_median\tgene_stdev\tgene_coverage\tstd_from_gene_median\tQC_total_std\thet_variants\tCNV_in_Normals\n")

bad_samples = ["D20-01807", "D20-05905", "D20-00636", "D20-02964", "47202"]
#bad_runs = ["TE19", "TE999"]
bad_runs = ["TE6","TE8","TE9","TE10","TE11","TE12","TE13","TE14","TE15","TE19","TE999"]


def Get_sample_and_run_names(sample_infile, sample_type) :
    sample_dict = {}
    i = 0
    for line in sample_infile :
        sample = line.strip().split("/")[-1].split("-ready")[0]
        if sample_type == "normal" :
            run = line.strip().split("/")[-3].split("_")[0]
            complete_run = line.strip().split("/")[-3]
        else :
            run = line.strip().split("/")[-5].split("_")[0]
            complete_run = line.strip().split("/")[-5]
        sample_dict[i] = [run, sample, complete_run]
        i += 1
    return sample_dict

def Get_normal_cnv1() :
    Normal_cnv_dict = {}
    header = True
    for line in infile_normal_cnv1 :
        if header :
            if line[:3] == "chr" :
                header = False
            else :
                continue
        lline = line.strip().split("\t")
        chrom = lline[0]
        start = int(lline[1])
        stop = int(lline[2])
        CN = lline[6]
        nr_markers = lline[7]
        confidence = lline[9]
        if chrom not in Normal_cnv_dict :
            Normal_cnv_dict[chrom] = []
        Normal_cnv_dict[chrom].append([chrom, start, stop, CN, nr_markers, confidence])
    return Normal_cnv_dict

def CNV_in_Normal(chrom, start, stop, Normal_cnv_dict) :
    start = int(start)
    stop = int(stop)
    Normal_CNV = ""
    for CNV in Normal_cnv_dict[chrom] :
        if (start >= CNV[1] and start <= CNV[2]) or (stop >= CNV[1] and stop <= CNV[2]) :
            if Normal_CNV == "" :
                    Normal_CNV += CNV[0] + ":" + str(CNV[1]) + "-" + str(CNV[2]) + "," + CNV[3] + "," + CNV[4] + "," + CNV[5]
            else :
                Normal_CNV += ";" + CNV[0] + ":" + str(CNV[1]) + "-" + str(CNV[2]) + "," + CNV[3] + "," + CNV[4] + "," + CNV[5]

    return Normal_CNV



def Read_coverage_file(exon_coverage, gene_coverage, sample_included_list, gene_chr, sample_dict, coverage_file) :
    i = 0
    for line in coverage_file :
        lline = line.strip().split("\t")
        gene = lline[3].split("_")[0]
        transcript = lline[3].split("_[")[0]
        if gene not in gene_coverage :
            gene_coverage[gene] = []
            gene_chr[gene] = lline[0]
            k = 0
            for coverage in lline[4:] :
                run, sample, run_complete = sample_dict[k]
                if run in bad_runs or sample in bad_samples :# or transcript in bad_transcript :
                    k += 1
                    continue
                gene_coverage[gene].append(0)
                k += 1
        exon_coverage.append([[lline[0], lline[1], lline[2], lline[3]],[]])
        k = 0
        j = 0
        for coverage in lline[4:] :
            run, sample, run_complete = sample_dict[k]
            if run in bad_runs or sample in bad_samples :# or transcript in bad_transcript :
                k += 1
                continue
            if i == 0 :
                sample_included_list.append([run, sample, run_complete])
            exon_coverage[i][1].append(int(coverage))
            gene_coverage[gene][j] += int(coverage)
            k += 1
            j += 1
        i += 1

def Calculate_avg_coverage_per_sample(avg_cov, avg_cov_X, exon_coverage, tc, tc_len, tcx, tcx_len):
    total_cov = 0
    total_cov_X = 0
    for exon in exon_coverage :
        if exon[0][0] == "chrX" :
            if avg_cov_X == [] :
                i = 0
                for c in exon[1] :
                    avg_cov_X.append(c)
                    total_cov_X += c
                    i += 1
            else :
                i = 0
                for c in exon[1] :
                    avg_cov_X[i] += c
                    total_cov_X += c
                    i += 1
        else :
            if avg_cov == [] :
                i = 0
                for c in exon[1] :
                    avg_cov.append(c)
                    total_cov += c
                    i += 1
            else :
                i = 0
                for c in exon[1] :
                    avg_cov[i] += c
                    total_cov += c
                    i += 1
    i = 0
    for ac in avg_cov :
        if tc > 0 :
            avg_cov[i] = (ac / float(tc)) * tc_len
        else :
            avg_cov[i] = (ac / float(total_cov)) * len(avg_cov)
        i += 1
    i = 0
    for ac in avg_cov_X :
        if tcx > 0 :
            avg_cov_X[i] = (ac / float(tcx)) * tcx_len
        else :
            avg_cov_X[i] = (ac / float(total_cov_X)) * len(avg_cov_X)
        i += 1
    return total_cov, len(avg_cov), total_cov_X, len(avg_cov_X)


def Normalize_exons(exon_coverage, avg_cov, avg_cov_X) :
    i = 0
    for exon in exon_coverage :
        if exon[0][0] == "chrX" :
            j = 0
            for c in exon[1] :
                exon_coverage[i][1][j] = c / avg_cov_X[j]
                j += 1
        else :
            j = 0
            for c in exon[1] :
                exon_coverage[i][1][j] = c / avg_cov[j]
                j += 1
        i += 1

def Normalize_genes(gene_coverage, gene_chr, avg_cov, avg_cov_X) :
    i = 0
    for gene in gene_coverage :
        if gene_chr[gene] == "chrX" :
            j = 0
            for c in gene_coverage[gene] :
                gene_coverage[gene][j] = c / avg_cov_X[j]
                j += 1
        else :
            j = 0
            for c in gene_coverage[gene]  :
                gene_coverage[gene][j] = c / avg_cov[j]
                j += 1
        i += 1

def Calculate_sample_QC(sample_exon_coverage, exon_median, exon_stdev) :
    j = 0
    sample_std = {}
    sample_std_sum = []
    Sample_QC = []
    for exon in sample_exon_coverage :
        i = 0
        for coverage in exon[1] :
            Std_from_exon_median = 0
            if exon_stdev[j] > 0 :
                Std_from_exon_median = (coverage - exon_median[j]) / exon_stdev[j]
            if j == 0 :
                sample_std[i] = [Std_from_exon_median]
                sample_std_sum.append(abs(Std_from_exon_median))
            else :
                sample_std[i].append(Std_from_exon_median)
                sample_std_sum[i] += abs(Std_from_exon_median)
            i += 1
        j += 1
    for sample in sample_std :
        i = 0
        max_i = len(sample_std[sample])
        std_sum = 0.0
        while i < max_i - 10 :
            std_sum += abs(statistics.median(sample_std[sample][i:i+10]))
            i += 5
        Sample_QC.append([std_sum, sample_std_sum[sample], statistics.median(sample_std[sample]), statistics.stdev(sample_std[sample])])
        if std_sum > 750 :
            outfile_cnv.write("Warning: " + sample_dict[sample][1] + " has " + str(int(std_sum)) + " in QC measure: total standard deviations. This is higher than the threshold set to 750\n")
    return Sample_QC

def Call_single_exon_CNV(sample_exon_coverage, exon_median, exon_stdev, sample_gene_coverage, gene_median, gene_stdev, sample_included_list, Sample_QC) :
    j = 0
    results = []
    for exon in sample_exon_coverage :
        gene = exon[0][3].split("_")[0]
        i = 0
        for coverage in exon[1] :
            if exon_stdev[j] == 0 or exon_median[j] == 0 or sample_gene_coverage[gene][i] == 0 :
                continue
            Std_from_gene_median = (sample_gene_coverage[gene][i] - gene_median[gene]) / gene_stdev[gene]
            Std_from_exon_median = (coverage - exon_median[j]) / exon_stdev[j]
            cn = coverage / exon_median[j] * 2
            #Correct exon coverage based on the coverage of the gene in relation to the median coverage of the gene
            gene_corrected_coverage = coverage * gene_median[gene] / sample_gene_coverage[gene][i]
            gene_corrected_cn = gene_corrected_coverage  / exon_median[j] * 2
            gene_corrected_std = (gene_corrected_coverage - exon_median[j]) / exon_stdev[j]
            #Deletion
            if (gene_corrected_cn < 1.2 and gene_corrected_std < -4.5 and Std_from_exon_median < -4.5) :
                #Variant annotation
                variants_in_region = ""
                #variants_in_region = Annotate_with_vcf(sample_dict[i][2], sample_dict[i][1], exon[0][0],  exon[0][1],  exon[0][2])
                #CNV Normal annotation
                CNV_Normal = CNV_in_Normal(exon[0][0], exon[0][1], exon[0][2], Normal_cnv1_dict)
                #outfile_cnv.write("Deletion\t" + sample_included_list[i][0] + "\t" + sample_included_list[i][1] + "\t" + exon[0][0] + "\t" + exon[0][1] + "\t" + exon[0][2] + "\t" + exon[0][3] + "\t1\t" + str(coverage) + "\t" + str(gene_corrected_coverage) + "\t" + str(exon_median[j]) + "\t" + str(Std_from_exon_median) + "\t" + str(gene_corrected_std) + "\t" + str(cn) + "\t" + str(gene_corrected_cn) + "\t" + str(gene_median[gene]) + "\t" + str(gene_stdev[gene]) + "\t" + str(sample_gene_coverage[gene][i]) + "\t" + str(Std_from_gene_median) + "\t" + str(Sample_QC[i][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                results.append(["Deletion", sample_included_list[i][0], sample_included_list[i][1], exon[0][0], exon[0][1], exon[0][2], exon[0][3], 1, coverage, gene_corrected_coverage, exon_median[j], Std_from_exon_median, gene_corrected_std, cn, gene_corrected_cn, gene_median[gene], gene_stdev[gene], sample_gene_coverage[gene][i], Std_from_gene_median, Sample_QC[i][0], variants_in_region, CNV_Normal])
            #Duplication
            elif gene_corrected_cn > 2.8 and gene_corrected_std > 4.5 and Std_from_exon_median > 4.5 :
                #Variant annotation
                variants_in_region = ""
                #variants_in_region = Annotate_with_vcf(sample_dict[i][2], sample_dict[i][1], exon[0][0],  exon[0][1],  exon[0][2])
                #CNV Normal annotation
                CNV_Normal = CNV_in_Normal(exon[0][0], exon[0][1], exon[0][2], Normal_cnv1_dict)
                #outfile_cnv.write("Duplication\t" + sample_included_list[i][0] + "\t" + sample_included_list[i][1] + "\t" + exon[0][0] + "\t" + exon[0][1] + "\t" + exon[0][2] + "\t" + exon[0][3] + "\t1\t" + str(coverage) + "\t" + str(gene_corrected_coverage) + "\t" + str(exon_median[j]) + "\t" + str(Std_from_exon_median) + "\t" + str(gene_corrected_std) + "\t" + str(cn) + "\t" + str(gene_corrected_cn) + "\t" + str(gene_median[gene]) + "\t" + str(gene_stdev[gene]) + "\t" + str(sample_gene_coverage[gene][i]) + "\t" + str(Std_from_gene_median) + "\t" + str(Sample_QC[i][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                results.append(["Duplication", sample_included_list[i][0], sample_included_list[i][1], exon[0][0], exon[0][1], exon[0][2], exon[0][3], 1, coverage, gene_corrected_coverage, exon_median[j], Std_from_exon_median, gene_corrected_std, cn, gene_corrected_cn, gene_median[gene], gene_stdev[gene], sample_gene_coverage[gene][i], Std_from_gene_median, Sample_QC[i][0], variants_in_region, CNV_Normal])
            i += 1
        j += 1
    return results

def Get_non_median_exons(sample_exon_coverage, exon_candidates_del, exon_candidates_dup, exon_median, exon_stdev, std) :
    j = 0
    for exon in sample_exon_coverage :
        gene = exon[0][3].split("_")[0]
        i = 0
        for coverage in exon[1] :
            Std_from_exon_median = 0
            if exon_stdev[j] > 0 :
                Std_from_exon_median = (coverage - exon_median[j]) / exon_stdev[j]
            if Std_from_exon_median < -1 * std :
                if gene not in exon_candidates_del :
                    exon_candidates_del[gene] = {}
                if i in exon_candidates_del[gene] :
                    exon_candidates_del[gene][i].append([coverage,j,Std_from_exon_median])
                else :
                    exon_candidates_del[gene][i] = [[coverage,j,Std_from_exon_median]]
            if Std_from_exon_median > 1 * std :
                if gene not in exon_candidates_dup :
                    exon_candidates_dup[gene] = {}
                if i in exon_candidates_dup[gene] :
                    exon_candidates_dup[gene][i].append([coverage,j,Std_from_exon_median])
                else :
                    exon_candidates_dup[gene][i] = [[coverage,j,Std_from_exon_median]]
            i += 1
        j += 1

def Merge_adjecent_exons(exon_candidates, size) :
    exon_candidates2 = {}
    i = 0
    for gene in exon_candidates :
        for sample in exon_candidates[gene] :
            if len(exon_candidates[gene][sample]) <= size :
                continue
            start_j = 0
            stop_j = 0
            cov_list = []
            std_list = []
            for exon in exon_candidates[gene][sample] :
                #if size > 1 :
                    #print(start_j, stop_j, cov_list, std_list)
                if start_j == 0 :
                    start_j = exon[1]
                    stop_j = exon[1]
                    cov_list.append(exon[0])
                    std_list.append(exon[2])
                elif exon[1] == stop_j + 1 :
                    stop_j += 1
                    cov_list.append(exon[0])
                    std_list.append(exon[2])
                else :
                    #Trim low variation end exons
                    nr_exons = stop_j - start_j + 1
                    while nr_exons > size :
                        if abs(std_list[0]) < 1.5 :
                            start_j += 1
                            std_list = std_list[1:]
                            cov_list = cov_list[1:]
                            nr_exons -= 1
                        elif abs(std_list[-1]) < 1.5 :
                            stop_j -= 1
                            std_list = std_list[:-1]
                            cov_list = cov_list[:-1]
                            nr_exons -= 1
                        else :
                            break
                    #Keep exon regions > size after trimming
                    #print("!!!", start_j, stop_j, cov_list, std_list)
                    if nr_exons > size :
                        if gene not in exon_candidates2 :
                            exon_candidates2[gene] = {}
                        if sample not in exon_candidates2[gene] :
                            exon_candidates2[gene][sample] = []
                        exon_candidates2[gene][sample].append([sum(cov_list), start_j, stop_j])
                    start_j = exon[1]
                    stop_j = exon[1]
                    cov_list = [exon[0]]
            #Trim low variation end exons
            nr_exons = stop_j - start_j + 1
            while nr_exons > size :
                if abs(std_list[0]) < 1.5 :
                    start_j += 1
                    std_list = std_list[1:]
                    cov_list = cov_list[1:]
                    nr_exons -= 1
                elif abs(std_list[-1]) < 1.5 :
                    stop_j -= 1
                    std_list = std_list[:-1]
                    cov_list = cov_list[:-1]
                    nr_exons -= 1
                else :
                    break
            #Keep exon regions > size after trimming
            if nr_exons > size :
                if gene not in exon_candidates2 :
                    exon_candidates2[gene] = {}
                if sample not in exon_candidates2[gene] :
                    exon_candidates2[gene][sample] = []
                exon_candidates2[gene][sample].append([sum(cov_list), start_j, stop_j])
    return exon_candidates2


def Get_non_median_exons_large_CNV(sample_exon_coverage, exon_candidates_del, exon_candidates_dup, exon_median, exon_stdev, std) :
    j = 0
    for exon in sample_exon_coverage :
        gene = exon[0][3].split("_")[0]
        chrom = gene_chr[gene]
        i = 0
        for coverage in exon[1] :
            Std_from_exon_median = 0
            if exon_stdev[j] > 0:
                Std_from_exon_median = (coverage - exon_median[j]) / exon_stdev[j]
            if Std_from_exon_median < -1 * std :
                if chrom not in exon_candidates_del :
                    exon_candidates_del[chrom] = {}
                if i in exon_candidates_del[chrom] :
                    exon_candidates_del[chrom][i].append([coverage,j,Std_from_exon_median])
                else :
                    exon_candidates_del[chrom][i] = [[coverage,j,Std_from_exon_median]]
            if Std_from_exon_median > 1 * std :
                if chrom not in exon_candidates_dup :
                    exon_candidates_dup[chrom] = {}
                if i in exon_candidates_dup[chrom] :
                    exon_candidates_dup[chrom][i].append([coverage,j,Std_from_exon_median])
                else :
                    exon_candidates_dup[chrom][i] = [[coverage,j,Std_from_exon_median]]
            i += 1
        j += 1


def Call_multi_exon_CNV(CNV_type, exon_candidates2, normal_exon_coverage, sample_exon_coverage, gene_median, sample_gene_coverage, gene_stdev, sample_included_list, Sample_QC, CNV_size_type) :
    results = []
    for gene in exon_candidates2 :
        for sample in exon_candidates2[gene] :
            for exon_region in exon_candidates2[gene][sample] :
                exon_region_coverage = []
                start_exon = exon_region[1]
                stop_exon = exon_region[2]
                i = 0
                while start_exon <= stop_exon :
                    j = 0
                    for coverage in normal_exon_coverage[start_exon][1]:
                        if i == 0 :
                            exon_region_coverage.append(coverage)
                        else :
                            exon_region_coverage[j] += coverage
                        j += 1
                    i = 1
                    start_exon += 1
                exon_region_median = statistics.median(exon_region_coverage)
                exon_region_std = statistics.stdev(exon_region_coverage)
                if exon_region_std == 0 or exon_region_median == 0 :
                    continue
                std_diff = (exon_region[0] - exon_region_median) / exon_region_std
                cn = exon_region[0] / exon_region_median * 2
                nr_exons = exon_region[2] - exon_region[1] + 1
                if CNV_size_type == "gene" :
                    if sample_gene_coverage[gene][sample] == 0:
                        continue
                    Std_from_gene_median = (sample_gene_coverage[gene][sample] - gene_median[gene]) / gene_stdev[gene]
                    if exon_region[0] / sample_gene_coverage[gene][sample] > 0.75 :
                        gene_corrected_exon_region_coverage = exon_region[0]
                    else :
                        gene_corrected_exon_region_coverage = exon_region[0] * ((gene_median[gene] - exon_region_median) / (sample_gene_coverage[gene][sample] - exon_region[0]))
                    gene_corrected_std_diff = (gene_corrected_exon_region_coverage - exon_region_median) / exon_region_std
                    gene_corrected_cn = gene_corrected_exon_region_coverage / exon_region_median * 2
                    if nr_exons <= 3 :
                        if CNV_type == "Deletion" :
                            if (std_diff < -4.0 and gene_corrected_std_diff < -4.0 and gene_corrected_cn < 1.2) :
                                #Variant annotation
                                variants_in_region = ""
                                #variants_in_region = Annotate_with_vcf(sample_included_list[sample][2], sample_included_list[sample][1], sample_exon_coverage[exon_region[1]][0][0],  exon_coverage[exon_region[1]][0][1],  exon_coverage[exon_region[2]][0][2])
                                #CNV Normal annotation
                                CNV_Normal = CNV_in_Normal(sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[1]][0][2], Normal_cnv1_dict)
                                outfile_cnv.write("Deletion\t" + sample_included_list[sample][0] + "\t" + sample_included_list[sample][1] + "\t" + sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[2]][0][2] + "\t" + sample_exon_coverage[exon_region[1]][0][3] + " - " + sample_exon_coverage[exon_region[2]][0][3] + "\t" + str(nr_exons) + "\t" + str(exon_region[0]) + "\t" + str(gene_corrected_exon_region_coverage) + "\t" + str(exon_region_median) + "\t" + str(std_diff) + "\t" + str(gene_corrected_std_diff) + "\t" + str(cn) + "\t" + str(gene_corrected_cn) + "\t" + str(gene_median[gene]) + "\t" + str(gene_stdev[gene]) + "\t" + str(sample_gene_coverage[gene][sample]) + "\t" + str(Std_from_gene_median) + "\t" + str(Sample_QC[sample][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                                results.append([sample_exon_coverage[exon_region[1]][0][0], sample_exon_coverage[exon_region[1]][0][1], sample_exon_coverage[exon_region[2]][0][2], cn, gene_corrected_std_diff, sample_included_list[sample][1]])
                                outfile_sample_cnv = open("CNV/" + sample_included_list[sample][1] + "_calls.aed", "a")
                                outfile_sample_cnv.write(sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[1]][0][2] + "\t" + str(gene_corrected_cn) + "\trgb(166,0,0)" + "\n")
                                outfile_sample_cnv.close()
                        elif CNV_type == "Duplication" :
                            if std_diff > 4.0 and gene_corrected_std_diff > 4.0 and gene_corrected_cn > 2.8 :
                                #Variant annotation
                                variants_in_region = ""
                                #variants_in_region = Annotate_with_vcf(sample_included_list[sample][2], sample_included_list[sample][1], sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[2]][0][2])
                                #CNV Normal annotation
                                CNV_Normal = CNV_in_Normal(sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[1]][0][2], Normal_cnv1_dict)
                                outfile_cnv.write("Duplication\t" + sample_included_list[sample][0] + "\t" + sample_included_list[sample][1] + "\t" + sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[2]][0][2] + "\t" + sample_exon_coverage[exon_region[1]][0][3] + " - " + sample_exon_coverage[exon_region[2]][0][3] + "\t" + str(nr_exons) + "\t" + str(exon_region[0]) + "\t" + str(gene_corrected_exon_region_coverage) + "\t" + str(exon_region_median) + "\t" + str(std_diff) + "\t" + str(gene_corrected_std_diff) + "\t" + str(cn) + "\t" + str(gene_corrected_cn) + "\t" + str(gene_median[gene]) + "\t" + str(gene_stdev[gene]) + "\t" + str(sample_gene_coverage[gene][sample]) + "\t" + str(Std_from_gene_median) + "\t" + str(Sample_QC[sample][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                                results.append([sample_exon_coverage[exon_region[1]][0][0], sample_exon_coverage[exon_region[1]][0][1], sample_exon_coverage[exon_region[2]][0][2], cn, gene_corrected_std_diff, sample_included_list[sample][1]])
                                outfile_sample_cnv = open("CNV/" + sample_included_list[sample][1] + "_calls.aed", "a")
                                outfile_sample_cnv.write(sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[1]][0][2] + "\t" + str(gene_corrected_cn) + "\trgb(0,0,166)" + "\n")
                                outfile_sample_cnv.close()
                    else :
                        if CNV_type == "Deletion" :
                            if (std_diff < -4.0 and gene_corrected_std_diff < -4.0 and gene_corrected_cn < 1.3) :
                                #Variant annotation
                                variants_in_region = ""
                                #variants_in_region = Annotate_with_vcf(sample_included_list[sample][2], sample_included_list[sample][1], sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[2]][0][2])
                                #CNV Normal annotation
                                CNV_Normal = CNV_in_Normal(sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[1]][0][2], Normal_cnv1_dict)
                                outfile_cnv.write("Deletion\t" + sample_included_list[sample][0] + "\t" + sample_included_list[sample][1] + "\t" + sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[2]][0][2] + "\t" + sample_exon_coverage[exon_region[1]][0][3] + " - " + sample_exon_coverage[exon_region[2]][0][3] + "\t" + str(nr_exons) + "\t" + str(exon_region[0]) + "\t" + str(gene_corrected_exon_region_coverage) + "\t" + str(exon_region_median) + "\t" + str(std_diff) + "\t" + str(gene_corrected_std_diff)+ "\t" + str(cn) + "\t" + str(gene_corrected_cn) + "\t" + str(gene_median[gene]) + "\t" + str(gene_stdev[gene]) + "\t" + str(sample_gene_coverage[gene][sample]) + "\t" + str(Std_from_gene_median) + "\t" + str(Sample_QC[sample][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                                results.append([sample_exon_coverage[exon_region[1]][0][0], sample_exon_coverage[exon_region[1]][0][1], sample_exon_coverage[exon_region[2]][0][2], cn, gene_corrected_std_diff, sample_included_list[sample][1]])
                                outfile_sample_cnv = open("CNV/" + sample_included_list[sample][1] + "_calls.aed", "a")
                                outfile_sample_cnv.write(sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[1]][0][2] + "\t" + str(cn) + "\trgb(166,0,0)" + "\n")
                                outfile_sample_cnv.close()
                        elif CNV_type == "Duplication" :
                            if std_diff > 4.0 and gene_corrected_std_diff > 4.0 and gene_corrected_cn > 2.7 :
                                #Variant annotation
                                variants_in_region = ""
                                #variants_in_region = Annotate_with_vcf(sample_included_list[sample][2], sample_included_list[sample][1], sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[2]][0][2])
                                #CNV Normal annotation
                                CNV_Normal = CNV_in_Normal(sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[1]][0][2], Normal_cnv1_dict)
                                outfile_cnv.write("Duplication\t" + sample_included_list[sample][0] + "\t" + sample_included_list[sample][1] + "\t" + sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[2]][0][2] + "\t" + sample_exon_coverage[exon_region[1]][0][3] + " - " + sample_exon_coverage[exon_region[2]][0][3] + "\t" + str(nr_exons) + "\t" + str(exon_region[0]) + "\t" + str(gene_corrected_exon_region_coverage) + "\t" + str(exon_region_median) + "\t" + str(std_diff) + "\t" + str(gene_corrected_std_diff) + "\t" + str(cn) + "\t" + str(gene_corrected_cn) + "\t" + str(gene_median[gene]) + "\t" + str(gene_stdev[gene]) + "\t" + str(sample_gene_coverage[gene][sample]) + "\t" + str(Std_from_gene_median) + "\t" + str(Sample_QC[sample][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                                results.append([sample_exon_coverage[exon_region[1]][0][0], sample_exon_coverage[exon_region[1]][0][1], sample_exon_coverage[exon_region[2]][0][2], cn, gene_corrected_std_diff, sample_included_list[sample][1]])
                                outfile_sample_cnv = open("CNV/" + sample_included_list[sample][1] + "_calls.aed", "a")
                                outfile_sample_cnv.write(sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[1]][0][2] + "\t" + str(cn) + "\trgb(0,0,166)" + "\n")
                                outfile_sample_cnv.close()
                elif CNV_size_type == "large" :
                    if CNV_type == "Deletion" :
                        if std_diff < -4.5 and cn < 1.7 :
                            #Variant annotation
                            variants_in_region = ""
                            #variants_in_region = Annotate_with_vcf(sample_included_list[sample][2], sample_included_list[sample][1], sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[2]][0][2])
                            #CNV Normal annotation
                            CNV_Normal = CNV_in_Normal(sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[1]][0][2], Normal_cnv1_dict)
                            outfile_cnv.write("Deletion\t" + sample_included_list[sample][0] + "\t" + sample_included_list[sample][1] + "\t" + sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[2]][0][2] + "\t" + sample_exon_coverage[exon_region[1]][0][3] + " - " + sample_exon_coverage[exon_region[2]][0][3] + "\t" + str(nr_exons) + "\t" + str(exon_region[0]) + "\t\t" + str(exon_region_median) + "\t" + str(std_diff) + "\t\t" + str(cn) + "\t\t\t\t\t\t" + str(Sample_QC[sample][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                            results.append([sample_exon_coverage[exon_region[1]][0][0], sample_exon_coverage[exon_region[1]][0][1], sample_exon_coverage[exon_region[2]][0][2], cn, std_diff, sample_included_list[sample][1]])
                            outfile_sample_cnv = open("CNV/" + sample_included_list[sample][1] + "_calls.aed", "a")
                            outfile_sample_cnv.write(sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[1]][0][2] + "\t" + str(cn) + "\trgb(166,0,0)" + "\n")
                            outfile_sample_cnv.close()
                    elif CNV_type == "Duplication" :
                        if std_diff > 4.5 and cn > 2.35 and Sample_QC[sample][0] < 7500 :#750 :
                            #Variant annotation
                            variants_in_region = ""
                            #variants_in_region = Annotate_with_vcf(sample_included_list[sample][2], sample_included_list[sample][1], sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[2]][0][2])
                            #CNV Normal annotation
                            CNV_Normal = CNV_in_Normal(sample_exon_coverage[exon_region[1]][0][0],  sample_exon_coverage[exon_region[1]][0][1],  sample_exon_coverage[exon_region[1]][0][2], Normal_cnv1_dict)
                            outfile_cnv.write("Duplication\t" + sample_included_list[sample][0] + "\t" + sample_included_list[sample][1] + "\t" + sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[2]][0][2] + "\t" + sample_exon_coverage[exon_region[1]][0][3] + " - " + sample_exon_coverage[exon_region[2]][0][3] + "\t" + str(nr_exons) + "\t" + str(exon_region[0]) + "\t\t" + str(exon_region_median) + "\t" + str(std_diff) + "\t\t" + str(cn) + "\t\t\t\t\t\t" + str(Sample_QC[sample][0]) + "\t" + variants_in_region + "\t" + CNV_Normal + "\n")
                            results.append([sample_exon_coverage[exon_region[1]][0][0], sample_exon_coverage[exon_region[1]][0][1], sample_exon_coverage[exon_region[2]][0][2], cn, std_diff, sample_included_list[sample][1]])
                            outfile_sample_cnv = open("CNV/" + sample_included_list[sample][1] + "_calls.aed", "a")
                            outfile_sample_cnv.write(sample_exon_coverage[exon_region[1]][0][0] + "\t" + sample_exon_coverage[exon_region[1]][0][1] + "\t" + sample_exon_coverage[exon_region[1]][0][2] + "\t" + str(cn) + "\trgb(0,0,166)" + "\n")
                            outfile_sample_cnv.close()
    return results

def Annotate_with_vcf(run, sample, chrom, start, end) :
    if run == "TE6_200213" :
        return "NA"
    #print(run,sample,chrom,start,end)
    start = int(start)
    end = int(end)
    #in_vcf = open("/home/jonas/Investigations/CNV_WP3/VCF/" + run + "." + sample + ".All_panels_het_SNP_DP20.vcf")
    in_vcf = open("/projects/wp3/nobackup/TWIST/OUTBOX/" + run  + "/VCF/" + sample + "-gatk-haplotype_cg_1.vcf")
    header = True
    variants = ""
    for line in in_vcf :
        if header :
            header = False
            continue
        lline = line.strip().split("\t")
        vcf_chrom = lline[0]
        if chrom != vcf_chrom :
            continue
        pos = int(lline[1])
        if pos >= start and pos <= end :
            AF = lline[7].split("AF=")[1].split(";")[0]
            DP = lline[7].split("DP=")[1].split(";")[0]
            if variants != "" :
                variants += ";"
            variants += AF + ":" + DP
            #print(variants)
    in_vcf.close()
    return variants



def Plot_CNV_figure(sample_included_list, sample_exon_coverage, exon_median, exon_stdev, Single_exon_calls, Multi_exon_calls) :
    sample_nr = 0
    for info in sample_included_list :
        run = info[0]
        sample = info[1]

        j = 0
        plot_data = []
        for exon in sample_exon_coverage :
            coverage = exon[1][sample_nr]
            Std_from_exon_median = (coverage - exon_median[j]) / exon_stdev[j]
            cn = coverage / exon_median[j] * 2
            Std_from_exon_median2 = 0
            cn2 = 0
            found = False
            for single_call in Single_exon_calls :
                if exon[0][0] == single_call[3] and int(exon[0][1]) <= int(single_call[4]) and int(exon[0][2]) >= int(single_call[4]) and sample == single_call[3] :
                    Std_from_exon_median2 = single_call[12]
                    cn2 = single_call[14]
                    found = True
            for multi_call in Multi_exon_calls :
                if exon[0][0] == multi_call[0] and int(exon[0][1]) >= int(multi_call[1]) and int(exon[0][1]) <= int(multi_call[2]) and sample == multi_call[5] :
                    Std_from_exon_median2 = multi_call[4]
                    cn2 = multi_call[3]
                    found = True
            int_chr = exon[0][0][3:]
            if int_chr == "X" :
                int_chr = 23
            elif int_chr == "Y" :
                int_chr = 24
            elif int_chr == "MT" :
                int_chr = 25
            elif len(int_chr) > 2 :
                int_chr = 26
            else :
                int_chr = int(int_chr)
            pos = int(exon[0][1])
            if found :
                plot_data.append([int_chr, pos, Std_from_exon_median2, cn2, 1])
            plot_data.append([int_chr, pos, Std_from_exon_median, cn, 0])
            j += 1
        plot_data.sort()

        j = 0
        std_data = []
        cn_data = []
        std_data2 = []
        cn_data2 = []
        x = []
        x2 = []
        chr_boundaries = []
        chrom = 1
        for data in plot_data :
            std_data.append(data[2])
            cn_data.append(data[3])
            x.append(j)
            if data[4] == 1 :
                std_data2.append(data[2])
                cn_data2.append(data[3])
                x2.append(j)
            if data[0] != chrom :
                chr_boundaries.append(j)
                chrom = data[0]
            j += 1

        fig, ax = plt.subplots(2)
        ax[0].plot([0,len(x)], [4,4], "r")
        ax[0].plot([0,len(x)], [-4,-4], "r")
        for cb in chr_boundaries:
            ax[0].axvline(x=cb, color='y', linestyle='dashed', linewidth=1)
        ax[0].scatter(x, std_data, s=1, color='b')
        ax[0].scatter(x2, std_data2, s=1, color='r')
        ax[0].set(ylabel='Stdev from median',
               title='CNV plot of ' + sample + " in run " + run)
        ax[0].set_xticks([250,950,1550,1900,2200,2500,2820,3030,3250,3550,3860,4250,4600,4880,5150,5550,5930,6200,6580,6800])
        ax[0].set_xticklabels(["1","2","3","4","5","6","7","8","9","10","11","12","14","15","16","17","18","19","22","X"])
        #ax[0].grid()

        ax[1].plot([0,len(x)], [2.8,2.8], "r")
        ax[1].plot([0,len(x)], [1.2,1.2], "r")
        for cb in chr_boundaries:
            ax[1].axvline(x=cb, color='y', linestyle='dashed', linewidth=1)
        ax[1].scatter(x, cn_data, s=1, color='b')
        ax[1].scatter(x2, cn_data2, s=1, color='r')
        ax[1].set(xlabel='Chromosomes', ylabel='CN')
        ax[1].set_xticks([250,950,1550,1900,2200,2500,2820,3030,3250,3550,3860,4250,4600,4880,5150,5550,5930,6200,6580,6800])
        ax[1].set_xticklabels(["1","2","3","4","5","6","7","8","9","10","11","12","14","15","16","17","18","19","22","X"])
        #ax[1].grid()

        fig.savefig(outplotsdir + run + "_" + sample + ".png")
        #plt.show()
        plt.close(fig)

        sample_nr += 1



"""Start of main script"""
#Get normals and run names
normal_sample_dict = Get_sample_and_run_names(infile_normal_samples, "normal")

#Get samples and run names
sample_dict = Get_sample_and_run_names(infile_samples, "sample")

#Read Normal CNV file for annotation
Normal_cnv1_dict = Get_normal_cnv1()


#Read normal coverage file
normal_exon_coverage = []
normal_gene_coverage = {}
normal_sample_included_list = []
normal_gene_chr = {}
Read_coverage_file(normal_exon_coverage, normal_gene_coverage, normal_sample_included_list, normal_gene_chr, normal_sample_dict, infile_normal_coverage)

#Read run coverage file
exon_coverage = []
gene_coverage = {}
sample_included_list = []
gene_chr = {}
Read_coverage_file(exon_coverage, gene_coverage, sample_included_list, gene_chr, sample_dict, infile_coverage)


#Calculate avg coverage per normal
normal_avg_cov = []
normal_avg_cov_X = []
tc, tc_len, tcx, tcx_len = Calculate_avg_coverage_per_sample(normal_avg_cov, normal_avg_cov_X, normal_exon_coverage, 0, 0, 0, 0)

#Calculate avg coverage per sample
avg_cov = []
avg_cov_X = []
Calculate_avg_coverage_per_sample(avg_cov, avg_cov_X, exon_coverage, tc, tc_len, tcx, tcx_len)


#Normalize exons against total coverage per sample
Normalize_exons(normal_exon_coverage, normal_avg_cov, normal_avg_cov_X)
Normalize_exons(exon_coverage, avg_cov, avg_cov_X)

#Normalize genes against total coverage per sample
Normalize_genes(normal_gene_coverage, normal_gene_chr, normal_avg_cov, normal_avg_cov_X)
Normalize_genes(gene_coverage, gene_chr, avg_cov, avg_cov_X)


#Calculate exon statistics
normal_exon_median = []
normal_exon_stdev = []
for exon in normal_exon_coverage :
    normal_exon_median.append(statistics.median(exon[1]))
    normal_exon_stdev.append(statistics.stdev(exon[1]))


#Calculate gene statistics
normal_gene_median = {}
normal_gene_stdev = {}
for gene in normal_gene_coverage :
    normal_gene_median[gene] = statistics.median(normal_gene_coverage[gene])
    normal_gene_stdev[gene] = statistics.stdev(normal_gene_coverage[gene])


#Calculate sample waveiness and noise
Sample_QC = Calculate_sample_QC(exon_coverage, normal_exon_median, normal_exon_stdev)


#Create .aed output files with header
for sample_i in sample_dict :
    sample = sample_dict[sample_i][1]
    outfile_sample_cnv = open("CNV/" + sample + "_calls.aed", "w")
    outfile_sample_cnv.write("bio:sequence(aed:String)\tbio:start(aed:Integer)\tbio:end(aed:Integer)\taed:name(aed:String)\taed:value(aed:String)\tstyle:color(aed:Color)\n")
    outfile_sample_cnv.close()


#Call single exon CNV:s
#3 = chrom, 4 = pos
Single_exon_calls = Call_single_exon_CNV(exon_coverage, normal_exon_median, normal_exon_stdev, normal_gene_coverage, normal_gene_median, normal_gene_stdev, sample_included_list, Sample_QC)


#Multi-exon CNV step1: Collect exons with some variation and split them into deletions and duplications
exon_candidates_del = {}
exon_candidates_dup = {}
Get_non_median_exons(exon_coverage, exon_candidates_del, exon_candidates_dup, normal_exon_median, normal_exon_stdev, 1.0)

#Multi-exon CNV step2: Merge adjecent exons, keep only consecutive exons (> 1 exons)
exon_candidates_del2 = Merge_adjecent_exons(exon_candidates_del, 1)
exon_candidates_dup2 = Merge_adjecent_exons(exon_candidates_dup, 1)
#print(exon_candidates_dup2)

#Multi-exon CNV step3: Call multi-exon CNV within one gene
Multi_exon_calls = Call_multi_exon_CNV("Deletion", exon_candidates_del2, normal_exon_coverage, exon_coverage, normal_gene_median, gene_coverage, normal_gene_stdev, sample_included_list, Sample_QC, "gene")
Multi_exon_calls += (Call_multi_exon_CNV("Duplication", exon_candidates_dup2, normal_exon_coverage, exon_coverage, normal_gene_median, gene_coverage, normal_gene_stdev, sample_included_list, Sample_QC, "gene"))

#Large CNV step1: Collect exons with some variation and split them into deletions and duplications
large_CNV_candidates_del = {}
large_CNV_candidates_dup = {}
Get_non_median_exons_large_CNV(exon_coverage, large_CNV_candidates_del, large_CNV_candidates_dup, normal_exon_median, normal_exon_stdev, 0.0)

#Multi-exon CNV step2: Merge adjecent exons, keep only large CNVs (> 15 exons)
large_CNV_candidates_del2 = Merge_adjecent_exons(large_CNV_candidates_del, 15)
large_CNV_candidates_dup2 = Merge_adjecent_exons(large_CNV_candidates_dup, 15)

#Multi-exon CNV step3: Call large CNVs
Multi_exon_calls += (Call_multi_exon_CNV("Deletion", large_CNV_candidates_del2, normal_exon_coverage, exon_coverage, normal_gene_median, gene_coverage, normal_gene_stdev, sample_included_list, Sample_QC, "large"))
Multi_exon_calls += (Call_multi_exon_CNV("Duplication", large_CNV_candidates_dup2, normal_exon_coverage, exon_coverage, normal_gene_median, gene_coverage, normal_gene_stdev, sample_included_list, Sample_QC, "large"))

#Output single exon CNVs not overlapping with larger regions
#Output .aed files
for single_call in Single_exon_calls :
    found = False
    single_chr = single_call[3]
    single_pos = single_call[4]
    for multi_call in Multi_exon_calls :
        if single_chr == multi_call[0] and int(single_pos) >= int(multi_call[1]) and int(single_pos) <= int(multi_call[2]) and single_call[2] == multi_call[5] :
            found = True
            break
    if not found :
        outfile_cnv.write(single_call[0])
        for output in single_call[1:] :
            outfile_cnv.write("\t" + str(output))
        outfile_cnv.write("\n")
        sample = single_call[2]
        outfile_sample_cnv = open("CNV/" + sample + "_calls.aed", "a")
        if single_call[12] > 2:
            #Blue
            outfile_sample_cnv.write(str(single_call[3]) + "\t" + str(single_call[4]) + "\t" + str(single_call[5]) + "\t" + str(single_call[12]) + "\trgb(0,0,166)" + "\n")
        else :
            #Red
            outfile_sample_cnv.write(str(single_call[3]) + "\t" + str(single_call[4]) + "\t" + str(single_call[5]) + "\t" + str(single_call[12]) + "\trgb(166,0,0)" + "\n")
        outfile_sample_cnv.close()
outfile_cnv.close()

#Plot CNV figures
Plot_CNV_figure(sample_included_list, exon_coverage, normal_exon_median, normal_exon_stdev, Single_exon_calls, Multi_exon_calls)
