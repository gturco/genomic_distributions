from intersection import Intersecter, Feature
from flatfeature import Bed
from collections import defaultdict

def loadintointersect(bed_file):
        query_list = {}
        feature_list = Bed(bed_file)
        for feature in feature_list:
            if float(feature['accn']) < .4: continue
            if feature['seqid'] not in list(query_list):
                query_list[feature['seqid']] = Intersecter()
            query_list[feature['seqid']].add_interval(Feature(int(feature['start']),
            int(feature['end']),name=feature['strand']))
        return query_list

def main(features_bed, results_bed, interval,padding):
    header = "start\tfreq\tchr\tType"
    print header
    features = Bed(features_bed)
    query_list = loadintointersect(results_bed)
    for feature in features:
        region = range(int(feature["start"]),int(feature["end"]))
        gaps = region[::interval-1]
        for start in gaps:
            end = start + interval
            freq = query_list[feature['seqid']].find(start, end)
            #print feature['seqid'], start, len(freq)
            print "{0}\t{1}\t{2}\t{3}".format(start,len(freq),feature['seqid'],"type of methylation group")

def maketable(dic,filename,mtype,TSS_type):
    #outfile = open(filename, "wb")
    header = list(dic.keys())
    values = []
    #outfile.write("pos\tfreq\n")
    for key in dic.keys():
        value = sum(dic[key])
        values.append(value)
        line = "{0}\t{1}\t{2}\t{3}".format(key,value,mtype,TSS_type)
        #outfile.write(line)
    	print line
    #outfile.close()


def genespace(results_bed,genelistbed,window_size,mtype):
    genelist = Bed(genelistbed)
    searchrange = range(-2000,2000)
    TSS_sites = defaultdict(list)
    TE_sites = defaultdict(list)
    query_list = loadintointersect(results_bed)
    for gene in genelist:
        TSS = gene["start"]
        TE = gene["end"]
        for region in searchrange[::window_size]:
            TSS_start = region + TSS
            TSS_end = region + window_size + TSS
            TE_start = region + TE
            TE_end = region + window_size + TE
            tss_matches = query_list[gene['seqid']].find(TSS_start, TSS_end)
            TSS_freq = len(list(tssm.name for tssm in tss_matches if tssm.name == gene['strand']))
            TE_matches = query_list[gene['seqid']].find(TE_start, TE_end)
            TE_freq = len(list(tem.name for tem in TE_matches if tem.name == gene['strand']))
            TSS_sites[str(region)].append(TSS_freq)
            TE_sites[str(region)].append(TE_freq)
    maketable(TSS_sites,"{0}_TSS_sites.txt".format(results_bed),mtype,"TSS")
    maketable(TE_sites,"{0}_TE_sites.txt".format(results_bed),mtype,"TE")
#going to need to add in methylation type in later
#numb  numb methlyation_group(H3K4Me3) str

#main("lib7_rmdup-W200-G600-FDR.001-islandfiltered.bed", "tair.bed", 100,1000)
#genespace("H3K27_rmdup-W200-G600-FDR.005-islandfiltered.bed", "tair.bed", 20, "H3K27")
#genespace("H3K4_rmdup-W200-G600-FDR.005-islandfiltered.bed", "tair.bed", 20, "H3K4")
#
#qplot(pos,freq, data=m, geom="line")

### download from http://www.phytozome.net/biomart/martview/bb80d3ad2d21beb8f457c27a1eeb3c87
## change ^C/c

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-f", dest="file", help="bedfile containing sig hits islandfiltered.bed")
    parser.add_option("--ref", dest="ref", help="ref genome bed (sorg_v2)")
    parser.add_option("--window", dest="window", help="window for grouping values")
    parser.add_option("--mtype", dest="mtype", help="methylation type")
    (options, _) = parser.parse_args()

    genespace(options.file, options.ref, int(options.window), options.mtype)
