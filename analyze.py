from ibidas import *
from ibidas.utils import util
import numpy
import sys
def load_data(subreads,ccs):
    d = Read(subreads,delimiter='\t')
    e = Read(ccs,delimiter='\t')
    ex = e.Get(_.f0.SplitOnPattern('/')[1].Cast(int)/'name',
                _.f2/'ccs_chrom',
                _.f3.Cast(int)/'ccs_pos',
                _.f1.Cast(int)/'ccs_mapflags',
                _.f4.Cast(float)/'ccs_mapping_quality',
                _.f5/'ccs_cigar',
                _.f6.Cast(int)/'ccs_length',
                    _.f7.SplitOnPattern(':')[2].Cast(float)/'coverage', 
                    _.f8.SplitOnPattern(':')[2].Cast(int)/'complete_coverage', 
                    _.f9.SplitOnPattern(':')[2].Cast(float)/'quality', 
                    _.f10.SplitOnPattern(':')[2].SplitOnPattern(',')[1:].Cast(float)/'ccs_snr',
                    _.f12.SplitOnPattern(':').Array().Each(lambda x: float(x[2]) if len(x) == 3 else Missing,dtype='real64?')/'ccs_mapping_concensus').Copy()

    dx = d.Get(_.f0.SplitOnPattern('/').Get(_[1].Cast(int)/'name', 
                _[2].SplitOnPattern('_').Get(_[0].Cast(int)/'start',_[1].Cast(int)/'stop')), 
                    _.f2/'subread_chrom',
                    _.f3.Cast(int)/'subread_pos',
                    _.f1.Cast(int)/'subread_mapflags',
                    _.f4.Cast(float)/'subread_mapping_quality',
                    _.f5/'subread_cigar',
                    _.f6.Cast(int)/'subread_length',
                _.f12.SplitOnPattern(':')[2].SplitOnPattern(',')[1:].Cast(float)/'subread_snr',
                _.f7.SplitOnPattern(':')[2].Cast(int)/'subread_flags',
                _.f15.SplitOnPattern(':').Array().Each(lambda x: float(x[2]) if len(x) == 3 else Missing,dtype='real64?')/'subread_mapping_concensus').Copy()

    gex = ex.GroupBy(_.name).Get(_.name, ((_.ccs_mapflags & 2048)>0).Sum()/'ccs_supalign')
    ex = (ex|Match(_.name, merge_same=True)| gex)

    dx = dx.Get(_.Get((_.name, _.start, _.stop))/'identifier','*')
    gdx = dx.GroupBy(_.identifier).Get(_.identifier, ((_.subread_mapflags & 2048)>0).Sum()/'subread_supalign', _[(_.subread_mapflags & 2048) == 0].Get(_.subread_length.Sum()/'zmw_length', _.subread_length.Count()/'zmw_cycles'))
    dx = dx|Match(_.identifier, merge_same=True)| gdx

    dx = dx % ('subreads','snr_sub')
    ex = ex % ('ccs','snr_ccs')

    d = (dx |Match(_.name, jointype='left', merge_same=True)| ex).Sort(_.name, _.start).Copy()
    return d



def run_stats(d,filename):
    #statistics
    f = open(filename,'w')
    dx = d[(_.subread_mapflags & 2048) == 0][IsMissing(_.ccs_mapflags)|((_.ccs_mapflags & 2048) == 0)].Copy() #drop supplementary alignments
    ccs = dx[~IsMissing(_.coverage)].Copy()
    nonccs = dx[IsMissing(_.coverage)].Copy()

    ccs_reads = ccs.name.Unique().Count()()
    nonccs_reads = nonccs.name.Unique().Count()()
    f.write('Number of ZWM-reads (percentage)\n')
    f.write('CCS:     %d (%.1f%%)\n' % (ccs_reads, (float(ccs_reads) / (ccs_reads + nonccs_reads)) * 100.0))
    f.write('not-CCS: %d (%.1f%%)\n\n' % (nonccs_reads, (float(nonccs_reads) / (ccs_reads + nonccs_reads)) * 100.0))

    ccs_10k_reads = ccs[_.zmw_length > 10000].name.Unique().Count()()
    nonccs_10k_reads = nonccs[_.zmw_length > 10000].name.Unique().Count()()
    f.write('Number of ZMW-reads (ZMW length > 10000) (percentage)\n')
    f.write('CCS:     %d (%.1f%%)\n' % (ccs_10k_reads, (ccs_10k_reads / (float(ccs_10k_reads) + nonccs_10k_reads)) * 100.0))
    f.write('not-CCS: %d (%.1f%%)\n\n' % (nonccs_10k_reads, (nonccs_10k_reads / (float(ccs_10k_reads) + nonccs_10k_reads)) * 100.0))

    ccs_input_gb = ccs.subread_length.Sum()()
    nonccs_input_gb = nonccs.subread_length.Sum()()

    ccs_output_gb = ccs.Get(_.name, _.ccs_length).Unique().ccs_length.Sum()()
    ccs_average_passes = float(ccs_input_gb) / float(ccs_output_gb)
    f.write('Input read Gb\n')
    f.write('CCS:     %.1f Gb (%.1f%%)\n' % (ccs_input_gb / 1000000000.0, (float(ccs_input_gb) / (ccs_input_gb + nonccs_input_gb)) * 100.0))
    f.write('not-CCS: %.1f Gb (%.1f%%)\n\n' % (nonccs_input_gb / 1000000000.0, (float(nonccs_input_gb) / (ccs_input_gb + nonccs_input_gb)) * 100.0))
    
    f.write('Output read Gb\n')
    ccs_output_median_subread = ccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length').subread_length.Sum()()
    ccs_output_ccs_read = ccs.GroupBy(_.name,flat=_.ccs_length).Get(_.name, _.ccs_length).ccs_length.Sum()()
    nonccs_output_median_subread = nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length').subread_length.Sum()()
    nonccs_output_median_subread_qual = nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length', _.subread_snr.Mean(), _.subread_length.Sum()/'fr_length')[_.subread_snr > 7.5][_.fr_length > 10000].subread_length.Sum()()
    f.write('CCS (median subread):      %.1f Mb (%.1f%%)\n' % (ccs_output_median_subread / 1000000.0, (float(ccs_output_median_subread) / (ccs_output_median_subread + nonccs_output_median_subread)) * 100.0))
    f.write('CCS (ccs read):            %.1f Mb (for reference)\n' % (ccs_output_ccs_read / 1000000.0))
    f.write('not-CCS (median subread):  %.1f Mb (%.1f%%)\n' % (nonccs_output_median_subread / 1000000.0, (float(nonccs_output_median_subread) / (ccs_output_median_subread + nonccs_output_median_subread)) * 100.0))
    f.write('not-CCS (median subread):  %.1f Mb (quality filter: snr>7.5, zmw-length>10000)\n\n' % (nonccs_output_median_subread_qual / 1000000.0,))

   
    f.write('Output read Gb (insert size length > 10000)\n')
    ccs_output_median_subread_10k = ccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length')[_.subread_length > 10000].subread_length.Sum()()
    ccs_output_ccs_read_10k = ccs.GroupBy(_.name,flat=_.ccs_length).Get(_.name, _.ccs_length)[_.ccs_length > 10000].ccs_length.Sum()()
    nonccs_output_median_subread_10k = nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length')[_.subread_length > 10000].subread_length.Sum()()
    nonccs_output_median_subread_qual = nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length', _.subread_snr.Mean(), _.subread_length.Sum()/'fr_length')[_.subread_snr > 7.5][_.subread_length > 10000].subread_length.Sum()()

    f.write('CCS (median subread):      %.1f Mb (%.1f%%)\n' % (ccs_output_median_subread_10k / 1000000.0, (float(ccs_output_median_subread_10k) / (ccs_output_median_subread_10k + nonccs_output_median_subread_10k)) * 100.0))
    f.write('CCS (ccs read):            %.1f Mb (for reference)\n' % (ccs_output_ccs_read_10k / 1000000.0))
    f.write('not-CCS (median subread):  %.1f Mb (%.1f%%)\n' % (nonccs_output_median_subread_10k / 1000000.0, (float(nonccs_output_median_subread_10k) / (ccs_output_median_subread_10k + nonccs_output_median_subread_10k)) * 100.0))
    f.write('not-CCS (median subread):  %.1f Mb (quality filter: snr>7.5)\n\n' % (nonccs_output_median_subread_qual / 1000000.0,))

    ccs_avg_fullread_lengths= numpy.percentile(ccs.GroupBy(_.name).Get(_.name, _.subread_length.Sum()).subread_length(), [0,5, 25,50, 75,95,100])
    nonccs_avg_fullread_lengths = numpy.percentile(nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Sum()).subread_length(), [0,5,25,50,75,95,100])
    f.write('ZMW-read length percentiles (0,5,25,(50),75,95,100)\n')
    f.write('CCS:     %d - %d - %d - (%d) - %d - %d - %d\n' % tuple(ccs_avg_fullread_lengths))
    f.write('not-CCS: %d - %d - %d - (%d) - %d - %d - %d\n\n' % tuple(nonccs_avg_fullread_lengths))

   

    f.write('Median subread length percentiles (0,5,25,(50),75,95,100)\n')
    ccs_avg_median_subread_lengths= numpy.percentile(ccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()).subread_length(), [0,5,25,50,75,95,100])
    nonccs_avg_median_subread_lengths= numpy.percentile(nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()).subread_length(), [0,5,25,50,75,95,100])
    f.write('CCS:     %d - %d - %d - (%d) - %d - %d - %d\n' % tuple(ccs_avg_median_subread_lengths))
    f.write('not-CCS: %d - %d - %d - (%d) - %d - %d - %d\n\n' % tuple(nonccs_avg_median_subread_lengths))

    f.write('CCS read length percentiles (0,5,25,(50),75,95,100)\n')
    ccs_avg_ccs_length = numpy.percentile(ccs.GroupBy(_.name).Get(_.name, _.ccs_length.Median()).ccs_length(), [0,5,25,50,75,95,100])
    f.write('CCS:     %d - %d - %d - (%d) - %d - %d - %d\n\n' % tuple(ccs_avg_ccs_length))

    ccs_median_coverage = numpy.percentile(ccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length', _.subread_length.Sum()/'full_length').Get(_.full_length / _.subread_length)(), [0,5,25,50,75,95,100])
    ccs_max_coverage = numpy.percentile(ccs.GroupBy(_.name).Get(_.name, _.subread_length.Max()/'subread_length', _.subread_length.Sum()/'full_length').Get(_.full_length / _.subread_length)(), [0,5,25,50,75,95,100])
    ccs_actual_coverage = numpy.percentile(ccs.GroupBy(_.name,flat=_.coverage).Get(_.name, _.coverage).coverage(), [0,5,25,50,75,95,100])
    nonccs_median_coverage = numpy.percentile(nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Median()/'subread_length', _.subread_length.Sum()/'full_length').Get(_.full_length / _.subread_length)(), [0,5,25,50,75,95,100])
    nonccs_max_coverage = numpy.percentile(nonccs.GroupBy(_.name).Get(_.name, _.subread_length.Max()/'subread_length', _.subread_length.Sum()/'full_length').Get(_.full_length / _.subread_length)(), [0,5,25,50,75,95,100])
    f.write('Coverage percentiles (0,5,25,(50),75,95,100)\n')
    f.write('CCS (wrt. median subread):     %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(ccs_median_coverage))
    f.write('CCS (wrt. max subread):     %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(ccs_max_coverage))
    f.write('CCS (wrt. ccs-read):     %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(ccs_actual_coverage))
    f.write('not-CCS (wrt. median subread):     %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(nonccs_median_coverage))
    f.write('not-CCS (wrt. max subread):     %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n\n' % tuple(nonccs_max_coverage))
    

    ccs_std = ccs.GroupBy(_.name)[_.subread_length.Count() >= 5].Get(_.name, _.subread_length.Median()/'median_subread_length', _.subread_length[1:-1]).\
                                    Get( _.subread_length.Count()/'read_count', _[(_.subread_length > 0.9 * _.median_subread_length) & (_.subread_length < 1.1 * _.median_subread_length)].subread_length.Count()/'read_fil_count').\
                                    Get((_.read_count - 1).Sum(), (_.read_fil_count - 1).Sum())()
    nonccs_std = nonccs.GroupBy(_.name)[_.subread_length.Count() >= 5].Get(_.name, _.subread_length.Median()/'median_subread_length', _.subread_length[1:-1]).\
                                    Get(_.subread_length.Count()/'read_count', _[(_.subread_length > 0.9 * _.median_subread_length) & (_.subread_length < 1.1 * _.median_subread_length)].subread_length.Count()/'read_fil_count').\
                                    Get((_.read_count - 1).Sum(), (_.read_fil_count - 1).Sum())()
    f.write('Length variation (fraction of subreads within 10% of length of median subread)\n')
    f.write(' (>=5 subreads, first/last subread excluded, 1 read subtracted to account for median subread)\n')
    f.write('CCS:     %.2f (%d/%d enclosed subreads)\n' % (ccs_std[1] / float(ccs_std[0]), ccs_std[1], ccs_std[0]))
    f.write('not-CCS: %.2f (%d/%d enclosed subreads)\n\n' % (nonccs_std[1] / float(nonccs_std[0]), nonccs_std[1], nonccs_std[0]))



    ccs_snr = numpy.percentile(ccs.GroupBy(_.name, flat=_.subread_snr).Get(_.name, _.subread_snr.Mean()).subread_snr(),[0,5,25,50,75,95,100])
    nonccs_snr = numpy.percentile(nonccs.GroupBy(_.name, flat=_.subread_snr).Get(_.name, _.subread_snr.Mean()).subread_snr(),[0,5,25,50,75,95,100])
    nonccs_snr_10k = numpy.percentile(nonccs.GroupBy(_.name, flat=_.subread_snr).Get(_.name, _.subread_snr.Mean(), _.subread_length.Sum())[_.subread_length > 10000].subread_snr(), [0,5,25,50,75,95,100])

    f.write('Subread Signal/Noise ratio percentiles (average ACGT)\n')
    f.write('CCS:                    %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(ccs_snr))
    f.write('not-CCS:                %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(nonccs_snr))
    f.write('not-CCS (fr > 10000):   %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n\n' % tuple(nonccs_snr_10k))
    
    
    
    ccs_aligned = ccs.Get(((_.subread_mapflags & 4) > 0).Sum()/'unaligned', _.subread_mapflags.Count()/'all')()
    nonccs_aligned = nonccs.Get(((_.subread_mapflags & 4) > 0).Sum()/'unaligned', _.subread_mapflags.Count()/'all')()
    nonccs_aligned_qual = nonccs[_.zmw_length > 10000][_.subread_snr.Mean() > 7.5].Get(((_.subread_mapflags & 4) > 0).Sum()/'unaligned', _.subread_mapflags.Count()/'all')()
    f.write('Subreads unaligned\n')
    f.write('CCS:                          %.1f%%\n' % ((ccs_aligned[0] / float(ccs_aligned[1])) * 100.0))
    f.write('not-CCS:                      %.1f%%\n' % ((nonccs_aligned[0] / float(nonccs_aligned[1])) * 100.0))
    f.write('not-CCS (snr>7.5,zl>10000):   %.1f%%\n\n' % ((nonccs_aligned_qual[0] / float(nonccs_aligned_qual[1])) * 100.0))
    
    ccs_mapconc = numpy.percentile(ccs[(_.subread_mapflags & 4) == 0].subread_mapping_concensus(),[0,5,25,50,75,95,100])
    ccs_mapconc_ref = numpy.percentile(ccs[(_.ccs_mapflags & 4) == 0].ccs_mapping_concensus(),[0,5,25,50,75,95,100])
    nonccs_mapconc = numpy.percentile(nonccs[(_.subread_mapflags & 4) == 0].subread_mapping_concensus(),[0,5,25,50,75,95,100])
    nonccs_mapconc_qual = numpy.percentile(nonccs[_.zmw_length > 10000][_.subread_snr.Mean() > 7.5][(_.subread_mapflags & 4) == 0].subread_mapping_concensus(),[0,5,25,50,75,95,100])
    f.write('Subreads mapping concensus percentiles (0,5,25,(50),75,95,100)\n')
    f.write('CCS (concensus reference):   %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(ccs_mapconc_ref))
    f.write('CCS:                         %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(ccs_mapconc))
    f.write('not-CCS:                     %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n' % tuple(nonccs_mapconc))
    f.write('not-CCS (snr>7.5,zl>10000):  %.1f - %.1f - %.1f - (%.1f) - %.1f - %.1f - %.1f\n\n' % tuple(nonccs_mapconc_qual))

    f.close()



d = load_data(sys.argv[1],sys.argv[2])
run_stats(d, sys.argv[3])
