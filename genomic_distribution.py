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



def get_freq(features,data_points,interval,freq_file):
    for feature in features:
        region = range(int(feature["start"]),int(feature["end"]))
        gaps = region[::interval-1]
        for start in gaps:
            end = start + interval
            matches = data_points[feature['seqid']].find(start, end)
            freq = [(min(m.stop,end)-max(start,m.start))for m in matches]
            circos_line = "at{0}\t{1}\t{2}\t{3}\n".format(feature['seqid'].strip("chr"),start,end,sum(freq))
            freq_file.write(circos_line)

def main(genomic_features_bed, data_bed,interval,outfile):
    header = "chr\tstart\tstop\tfreq\ttype\n"
    freq_file = open(outfile,"wb")
    freq_file.write(header)
    features = Bed(genomic_features_bed)
    data_points = loadintointersect(data_bed)
    get_freq(features,data_points,interval,freq_file)
    freq_file.close()

#main("karyotype.tair.bed","md/H3K27.islandfiltered.bed", 50000,"H327_tracks.txt")
#main("karyotype.tair.bed","md/H3K4.islandfiltered.bed", 50000,"H34_tracks.txt")

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--genomic", dest="genomic_bed", help="bedfile containing genomic info ie: chr start stops")
    parser.add_option("--data", dest="data_bed", help="bedfile containing data information ie histone mark mapping")
    parser.add_option("--window", dest="window",type='int', help="window for calc freqs")
    parser.add_option("--out", dest="outfile", help="name of outfile")
    (options, _) = parser.parse_args()

    main(options.genomic_bed, options.data_bed, options.window, options.outfile)
