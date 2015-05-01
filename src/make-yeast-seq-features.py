#! /usr/local/bin/python

import sys, os, math, string, random, argparse
import translate, cai, stats, biofile, util, protprop, na

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Calculate basic features of coding sequences")
	parser.add_argument(dest="cds_in_fname", type=str, help="FASTA file containing coding sequences")
	parser.add_argument(dest="prot_in_fname", type=str, help="FASTA file containing protein sequences")
	parser.add_argument(dest="feature_fname", type=str, help="SGD file containing sequence features")
	parser.add_argument(dest="paralog_fname", type=str, help="Yeast Gene Order Browser formatted file of paralog identifications")
	parser.add_argument("--aa", dest="do_aa", default=False, action="store_true", help="compute amino-acid frequencies?")
	parser.add_argument("--gc", dest="do_gc", default=False, action="store_true", help="compute GC frequencies?")
	parser.add_argument("--mw", dest="do_mw", default=False, action="store_true", help="compute molecular weights?")
	parser.add_argument("--target-aas", dest="target_aas", type=str, default=translate.AAs(), help="amino acids (e.g. ACDEF) for frequency analysis")
	parser.add_argument("-p", "--pseudo", dest="pseudocount", type=float, default=0.0, help="pseudocount to add to all frequencies")
	parser.add_argument("-o", "--out", dest="out_fname", type=str, default=None, help="output filename")
	options = parser.parse_args()

	cdna_dict = biofile.readFASTADict(os.path.expanduser(options.cds_in_fname))
	prot_dict = biofile.readFASTADict(os.path.expanduser(options.prot_in_fname))

	# Read paralog data from Yeast Gene Order Browser file
	ygob_data = util.readTable(file(os.path.expanduser(options.paralog_fname),'r'))
	paralog_dict = {}
	for flds in ygob_data.dictrows:
		scer1 = flds['scer1'].strip()
		scer2 = flds['scer2'].strip()
		if not (na.isNA(scer1) or na.isNA(scer2)):
			paralog_dict[scer1] = scer2
			paralog_dict[scer2] = scer1

	# Read SGD data
	sgd_features = util.readTable(file(os.path.expanduser(options.feature_fname),'r'), header=False)
	'''
	http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.README
	1.   Primary SGDID (mandatory)
	2.   Feature type (mandatory)
	3.   Feature qualifier (optional)
	4.   Feature name (optional)
	5.   Standard gene name (optional)
	6.   Alias (optional, multiples separated by |)
	7.   Parent feature name (optional)
	8.   Secondary SGDID (optional, multiples separated by |)
	9.   Chromosome (optional)
	10.  Start_coordinate (optional)
	11.  Stop_coordinate (optional)
	12.  Strand (optional)
	13.  Genetic position (optional)
	14.  Coordinate version (optional)
	15.  Sequence version (optional)
	16.  Description (optional)
	'''
	fld_names = ['sgdid','feature.type','feature.qualifier','orf','gene','alias','parent.feature','secondary.sgdid','chromosome','start','stop','strand','genetic.position','coordinate.version','sequence.version','description']
	fld_inds = dict([(fld_names[i],i) for i in range(len(fld_names))])
	sgd_dict = {}
	for flds in sgd_features.rows:
		if flds[fld_inds['feature.type']] == 'ORF':
			# Save this row
			orf = flds[fld_inds['orf']]
			sgd_dict[orf] = flds

	# Start up output
	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Assemble header information
	aa_freq_flds = []
	if options.do_aa:
		aa_freq_dict = dict([(aa, "f.%s"%aa) for aa in options.target_aas])
		aa_freq_flds = aa_freq_dict.values()
		aa_freq_flds.sort()
	gc_freq_flds = []
	if options.do_gc:
		gc_freq_flds = ["gc1","gc2","gc3","gc"]
	mw_flds = []
	if options.do_mw:
		mw_flds = ["mw"]
	#freq_flds = aa_freq_flds + gc_freq_flds + fop_flds + cai_flds + mw_flds

	aas = translate.AAs()
	pprop = protprop.ProteinProperties()
	#aa_line = '\t'.join(["f.%s"%aa for aa in aas])
	#outf.write("orf\tlen\tmw\tcai\tfop\tgc\tgc1\tgc2\tgc3\t%s\n" % aa_line)


	dout = util.DelimitedOutput()
	dout.addHeader('orf','S.cerevisiae systematic open reading frame identifier','s')
	dout.addHeader('gene','Common gene name','s')
	dout.addHeader('pair.orf','Paralog according to YGOB (http://ygob.ucd.ie)','s')
	dout.addHeader('feature.qualifier','Feature qualifier (Verified, Uncharacterized, Dubious)','s')
	dout.addHeader('chromosome','Chromosome','s')
	dout.addHeader('strand','Strand (W = Watson, C = Crick)','s')
	dout.addHeader('length.nt','Length of coding sequence in nucleotides','d')
	dout.addHeader('length.aa','Length of encoded protein sequence in amino acids','d')
	dout.addHeader('mw','Molecular weight of protein (Da)','f')
	dout.addHeader('pi','Isoelectric point','f')
	dout.addHeader('charge','Net charge at pH 7.0','f')
	for aa in aas:
		dout.addHeader('n.{}'.format(aa), 'Number of {} in protein sequence'.format(aa), 'd')
	dout.addHeader('alias','Aliases','s')
	dout.addHeader('description','Description','s')

	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	n_written = 0
	for orf in sorted(sgd_dict.keys()):
		feat = sgd_dict[orf]
		seq = cdna_dict[orf]
		prot = prot_dict[orf]
		#p1 = translate.translate(seq)
		#if p1 is None:
		#	print orf
		#	print translate.translateRaw(seq)
		fdict = {}

		for fld in ['orf','gene','alias','chromosome','strand','feature.qualifier','description']:
			fdict[fld.replace(".",'_')] = feat[fld_inds[fld]]
		fdict['pair_orf'] = paralog_dict.get(orf,'')
		fdict['mw'] = pprop.getMolecularWeight(prot)
		fdict['pi'] = pprop.getIsoelectricPoint(prot)
		fdict['charge'] = pprop.getCharge(prot, pH=7.0)
		fdict['length_nt'] = len(seq)
		fdict['length_aa'] = len(prot)
		fdict['chromosome'] = str(fdict['chromosome'])
		for aa in aas:
			fdict['n_{}'.format(aa)] = prot.count(aa)

		line = dout.getFormat().format(**fdict)
		data_outs.write(line)
		#print line,
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()
