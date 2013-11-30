from intersection import Intersecter, Feature
from flatfeature import Bed
from collections import defaultdict

def loadintointersect(bed_file):
        """loads features of bedfile into a tree for finding their freq within genomic regions"""
        query_list = {}
        feature_list = Bed(bed_file)
        for feature in feature_list:
            if feature['seqid'] not in list(query_list):
                query_list[feature['seqid']] = Intersecter()
            query_list[feature['seqid']].add_interval(Feature(int(feature['start']),
            int(feature['end']),name=['name']))
        return query_list

def maketable(dic,filename,mtype):
    outfile = open(filename, "wb")
    header = list(dic.keys())
    values = []
    outfile.write("pos\tfreq\tmtype\n")
    for key in dic.keys():
        value = sum(dic[key])
        values.append(value)
        line = "{0}\t{1}\t{2}\n".format(key,value,mtype)
        outfile.write(line)
    outfile.close()


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
            TSS_freq = len(query_list[gene['seqid']].find(TSS_start, TSS_end))
            TE_freq = len(query_list[gene['seqid']].find(TE_start, TE_end))
            TSS_sites[str(region)].append(TSS_freq)
            TE_sites[str(region)].append(TE_freq)
    maketable(TSS_sites,"{0}_TSS_sites.txt".format(results_bed),mtype)
    maketable(TE_sites,"{0}_TE_sites.txt".format(results_bed),mtype)
#going to need to add in methylation type in later
#numb  numb methlyation_group(H3K4Me3) str

#genespace("H3K27_rmdup-W200-G600-FDR.005-islandfiltered.bed", "tair.bed", 20, "H3K27")
#genespace("H3K4_rmdup-W200-G600-FDR.005-islandfiltered.bed", "tair.bed", 20, "H3K4")
#qplot(pos,freq, data=m, geom="line")

### download from http://www.phytozome.net/biomart/martview/bb80d3ad2d21beb8f457c27a1eeb3c87
## change ^C/c
#
if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--genomic", dest="genomic_bed", help="bedfile containing gene information")
    parser.add_option("--data", dest="data_bed", help="ref genome bed (sorg_v2)")
    parser.add_option("--window", dest="window_size",type="int", help="window for geting freq of sites along gene")
    parser.add_option("--type", dest="mtype", help="methylation type")
    (options, _) = parser.parse_args() 
    
    genespace(options.data_bed, options.genomic_bed, options.window_size, options.mtype)
