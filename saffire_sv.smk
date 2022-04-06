import pandas as pd
from pybedtools import BedTool
import numpy as np
import networkx as nx
import more_itertools as mit


configfile: 'config/saffire_sv.yaml'

PARTS=config.get('PARTS', 15)

MANIFEST = config.get('MANIFEST', 'config/manifest.tab')
SV_SIZE = config.get('SV_SIZE', '30000')
REF = config.get('REFERENCE')
FAI = REF+'.fai'

manifest_df = pd.read_csv(MANIFEST, sep='\t', index_col=['SAMPLE'])

def merge_cat(neighborhood, iloc):
	neighborhood[iloc].append(iloc+1)
	if iloc+1 in neighborhood:
		neighborhood[iloc+1].append(iloc)
	else:
		neighborhood[iloc+1] = [iloc]
	return neighborhood

def search_index(search_df, inv_index):
	search_list = []
	if inv_index > 0:
		search_list.append(search_df.at[inv_index-1, 'q_chr'])
	if inv_index < len(search_df)-1:
		search_list.append(search_df.at[inv_index+1, 'q_chr'])
	if search_df.at[inv_index, 'q_chr'] in search_list:
		return True
	else:
		return False

def find_qpos(row_check, broken_check):
	row_df = broken_check.loc[(broken_check['r_chr'] == row_check['r_chr']) & (broken_check['q_chr'] == row_check['q_chr'])]
	if len(row_df.loc[row_df['r_end'] == row_check['r_pos']]) > 0:
		return row_df.loc[row_df['r_end'] == row_check['r_pos']]['q_end'].values[0]
	elif len(row_df.loc[row_df['r_pos'] == row_check['r_end']]) > 0:
		return row_df.loc[row_df['r_pos'] == row_check['r_end']]['q_pos'].values[0]
	elif len(row_df.loc[row_df['r_pos'] == row_check['r_pos']]) > 0:
		return row_df.loc[row_df['r_pos'] == row_check['r_pos']]['q_pos'].values[0]
	elif len(row_df.loc[row_df['r_end'] == row_check['r_end']]) > 0:
		return row_df.loc[row_df['r_end'] == row_check['r_end']]['q_end'].values[0]
	else:
		return broken_check.loc[(broken_check['r_chr'] == row_check['r_chr']) & (broken_check['q_chr'] == broken_check['q_chr']) & ((broken_check['r_end'] >= broken_check['r_pos']) | (broken_check['r_pos'] >= row_check['r_end']))]['q_pos'].values[0]

def check_trans(calls):
	if 'SYNTENIC' in calls:
		return 'SYNTENIC'
	else:
		return calls [0]

def find_paf(wildcards):
	return manifest_df.at[wildcards.sample, 'PAF']

def find_contigs(wildcards):
	return gather.split('tmp/{sample}.{scatteritem}-{status}.paf', sample=wildcards.sample, status=wildcards.status)

scattergather:
    split=PARTS

wildcard_constraints:
	status='broken|orient',
	sample='|'.join(manifest_df.index)

localrules: all


rule all:
	input:
		expand('results/raw/{sample}/{sample}_var_all.bed', sample=manifest_df.index), expand('results/merged/{sample}/{sample}_var_all.bed', sample=manifest_df.index)


rule split_paf:
	input:
		paf = find_paf
	output:
		flag = temp(scatter.split("tmp/{{sample}}.{scatteritem}.paf"))
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.paf, sep='\t', low_memory=False, header=None)
		col_out = df.columns.values
		for i, contig in enumerate(df[0].unique()):
			out_num = (i % PARTS) + 1 
			df.loc[df[0] == contig][col_out].to_csv(f'tmp/{wildcards.sample}.{out_num}-of-{PARTS}.paf', sep='\t', index=False, header=False, mode='a+')

rule trim_break_orient_paf:
	input:
		paf = 'tmp/{sample}.{scatteritem}.paf'
	output:
		contig = temp('tmp/{sample}.{scatteritem}-orient.paf'),
		broken = temp('tmp/{sample}.{scatteritem}-broken.paf')
	threads: 1
	resources: 
		mem = 24,
		hrs = 24
	shell:
		'''
		rustybam orient {input.paf} | rustybam trim-paf | rb filter --paired-len 1000000 > {output.contig}
		rustybam break-paf --max-size {SV_SIZE} {output.contig} > {output.broken}
		'''

rule combine_paf:
	input:
		paf = find_contigs,
		flag = rules.split_paf.output.flag
	output:
		paf = 'tmp/{sample}-{status}.paf'
	threads: 1
	resources: 
		mem = 4,
		hrs = 24
	shell:
		'''
		cat {input.paf} > {output.paf} 
		'''


rule call_gaps:
	input:
		contig = find_paf,
		fai = FAI
	output:
		gap = 'tmp/{sample}_gap.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	shell:
		'''
		awk -vOFS="\\t" '{{print $6,$8,$9}}' {input.contig} | tail -n +2 | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools sort -i - -g {input.fai} | bedtools complement -i - -g {input.fai} > {output.gap}
		'''

rule call_simple_insdel:
	input:
		contig = 'tmp/{sample}-orient.paf',
		broken = 'tmp/{sample}-broken.paf'
	output:
		tsv = 'tmp/{sample}_sv_insdel.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		contig_df = pd.read_csv(input.contig, sep='\t', header=None, usecols=[0,1,2,3,4,5,7,8], names=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_pos', 'r_end'])
		broken_df = pd.read_csv(input.broken, sep='\t', header=None, usecols=[0,1,2,3,4,5,7,8], names=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_pos', 'r_end'])
		broken_df = broken_df.loc[broken_df['q_len'] > 1000000]
		contig_df = contig_df.loc[contig_df['q_len'] > 1000000]
		# Where alignments are normally contiguous, remove small alignments
		if np.median(broken_df['q_end'] - broken_df['q_pos']) > 100000:
			broken_df_short = broken_df.loc[broken_df['q_end'] - broken_df['q_pos'] < 20000].copy()
			drop_short = []
			for index in range(len(broken_df_short)):
				if index == 0 or index == len(broken_df)-1:
					continue
				else:
					check_up = index+1
					check_down = index-1
					contig = broken_df.iloc[index]['q_chr']
					if contig == broken_df.iloc[check_up]['q_chr'] and contig == broken_df.iloc[check_down]['q_chr']:
						drop_short.append(broken_df_short.iloc[index].name)
			broken_df = broken_df.drop(labels=drop_short).reset_index(drop=True).reset_index()
		else:
			broken_df = broken_df.reset_index(drop=True).reset_index()

		contig_df = contig_df.reset_index(drop=True).reset_index()
		broken_ins = BedTool.from_dataframe(broken_df[['q_chr', 'q_pos', 'q_end', 'index']])
		contig_ins = BedTool.from_dataframe(contig_df[['q_chr', 'q_pos', 'q_end', 'index']])
		broken_del = BedTool.from_dataframe(broken_df[['r_chr', 'r_pos', 'r_end', 'index']])
		contig_del = BedTool.from_dataframe(contig_df[['r_chr', 'r_pos', 'r_end', 'index', 'q_chr']])
		del_df = contig_del.subtract(broken_del).to_dataframe(names=['r_chr', 'r_pos', 'r_end', 'index', 'q_chr'])
		ins_df = contig_ins.subtract(broken_ins).to_dataframe(names=['q_chr', 'q_end', 'q_pos', 'index'])
		ins_df = ins_df.rename(columns={'index' : 'ins_index'})

		if len(ins_df) > 0:
			ins_df['SVLEN'] = ins_df['q_pos'] - ins_df['q_end']
			ins_df = ins_df.loc[ins_df['SVLEN'] > 30000]
			ins_uniq_start = broken_df[['r_chr', 'r_end', 'q_chr', 'q_end', 'strand']].merge(ins_df[['q_chr', 'q_end', 'ins_index', 'SVLEN']]).groupby(['q_chr', 'q_end', 'ins_index', 'SVLEN', 'r_chr', 'r_end']).count().reset_index()[['q_chr', 'q_end', 'ins_index', 'SVLEN', 'r_end', 'r_chr', 'strand']]
			ins_uniq_start = ins_uniq_start.loc[ins_uniq_start['strand'] == 1]
			ins_uniq_end = broken_df[['r_chr', 'r_pos', 'q_chr', 'q_pos', 'strand']].merge(ins_df[['q_chr', 'q_pos', 'ins_index', 'SVLEN']]).groupby(['q_chr', 'q_pos', 'ins_index', 'SVLEN', 'r_chr', 'r_pos']).count().reset_index()[['q_chr', 'q_pos', 'ins_index', 'SVLEN', 'r_pos', 'r_chr', 'strand']]
			ins_uniq_end = ins_uniq_end.loc[ins_uniq_end['strand'] == 1]
			ins_ref_locs = ins_uniq_end.merge(ins_uniq_start)
			ins_ref_locs['SVTYPE'] = 'DUP'
			ins_ref_locs['r_end'] = ins_ref_locs['r_end']+1
		else:
			ins_ref_locs = pd.DataFrame(columns=['SVLEN', 'SVTYPE', 'r_chr', 'r_pos', 'r_end', 'q_chr', 'q_pos', 'q_end'])


		if len(del_df) > 0:
			del_df['SVTYPE'] = 'DEL'
			del_df['SVLEN'] = del_df['r_end'] - del_df['r_pos']
			del_df = del_df.loc[del_df['SVLEN'] > 30000]
			del_df['q_pos'] = del_df.apply(lambda row: find_qpos(row, broken_df), axis=1)
			del_df['q_end'] = del_df['q_pos']+1
		else:
			del_df = pd.DataFrame(columns=['SVLEN', 'SVTYPE', 'r_chr', 'r_pos', 'r_end', 'q_chr', 'q_pos', 'q_end'])


		sv_out = ins_ref_locs.append(del_df)
		sv_out['ID'] = sv_out['r_chr']+'-'+sv_out['r_pos'].astype(str)+'-'+sv_out['SVTYPE']+'-'+sv_out['SVLEN'].astype(str)
		sv_out['CONTIG_SEQ'] = sv_out['q_chr']+":"+sv_out['q_pos'].astype(str)+"-"+sv_out['q_end'].astype(str)
		sv_out = sv_out.rename(columns={'r_chr' : '#CHROM', 'r_pos' : 'POS', 'r_end' : 'END'})

		sv_out = sv_out.loc[sv_out['POS'] <= sv_out['END']]

		sv_check = sv_out.copy()

		sv_check['POS'] = sv_check.apply(lambda row: max(0, row['POS']-50000), axis=1)
		sv_check['END'] = sv_check.apply(lambda row: row['POS']+50000, axis=1)

		contig_bed = BedTool.from_dataframe(sv_check[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID']].sort_values(['#CHROM', 'POS']))

		contig_check = contig_bed.intersect(BedTool.from_dataframe(contig_df[['r_chr', 'r_pos', 'r_end', 'q_chr']]), wa=True, wb=True).to_dataframe().groupby(['strand', 'blockCount']).count().reset_index().groupby('strand').count().reset_index()

		contig_check = contig_check.loc[contig_check['chrom'] == 1][['strand']].rename(columns={'strand' : 'ID'})

		sv_out = sv_out.merge(contig_check)

		sv_out[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID', 'CONTIG_SEQ']].sort_values(['#CHROM', 'POS']).to_csv(output.tsv, sep='\t', index=False)



rule call_inv:
	input:
		broken = 'tmp/{sample}-broken.paf'
	output:
		inv = 'tmp/{sample}_sv_inv.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.broken, sep='\t', header=None, usecols=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_pos', 'r_end', 'nmatch'], names=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_len', 'r_pos', 'r_end', 'nmatch', 'alen', 'mapq', 'id', 'cg'])

		df = df.loc[df['q_len'] > 1000000].sort_values(['q_chr', 'q_pos', 'r_pos']).reset_index(drop=True)

		df_inv = df.loc[(df['strand'] == '-') & (df['r_end'] - df['r_pos'] > 30000)].copy()

		for inv_index in range(len(df_inv)):
			orig_index = df_inv.iloc[inv_index].name
			inv_contig = df_inv.iloc[inv_index]['q_chr']
			df_inv.at[orig_index, 'INV_STATUS'] = search_index(df, orig_index)


		df_inv = df_inv.loc[df_inv['INV_STATUS'] == True]

		df_inv_out = pd.DataFrame()

		for grouping in [list(group) for group in mit.consecutive_groups(df_inv.sort_values(['q_chr', 'r_pos', 'q_pos'], ascending=False).index)]:
			df_inv_subset = df_inv.loc[grouping].copy()
			if len(df_inv_subset['q_chr'].unique()) == 1 and len(df_inv_subset['r_chr'].unique()) == 1:
				df_inv_qchr = df_inv_subset['q_chr'].unique()[0]
				df_inv_rchr = df_inv_subset['r_chr'].unique()[0]
				df_inv_rpos = min(df_inv_subset['r_pos'])
				df_inv_rend = max(df_inv_subset['r_end'])
				df_inv_qpos = min(df_inv_subset['q_pos'])
				df_inv_qend = max(df_inv_subset['q_end'])
				df_inv_out = df_inv_out.append(pd.DataFrame.from_dict({'q_chr' : [df_inv_qchr], 'q_pos' : [df_inv_qpos], 'q_end' : [df_inv_qend], 'r_chr' : [df_inv_rchr], 'r_pos' : [df_inv_rpos], 'r_end' : [df_inv_rend]}))
			else:
				for chrom in df_inv_subset['r_chr'].unique():
					for contig in df_inv_subset['q_chr'].unique():
						df_inv_subset_split = df_inv_subset.loc[(df_inv_subset['r_chr'] == chrom) & (df_inv_subset['q_chr'] == contig)]
						if len(df_inv_subset_split) > 0:
							df_inv_qchr = df_inv_subset_split['q_chr'].unique()[0]
							df_inv_rchr = df_inv_subset_split['r_chr'].unique()[0]
							df_inv_rpos = min(df_inv_subset_split['r_pos'])
							df_inv_rend = max(df_inv_subset_split['r_end'])				
							df_inv_qpos = min(df_inv_subset_split['q_pos'])
							df_inv_qend = max(df_inv_subset_split['q_end'])
							df_inv_out = df_inv_out.append(pd.DataFrame.from_dict({'q_chr' : [df_inv_qchr], 'q_pos' : [df_inv_qpos], 'q_end' : [df_inv_qend], 'r_chr' : [df_inv_rchr], 'r_pos' : [df_inv_rpos], 'r_end' : [df_inv_rend]}))


		# df_inv['INV_ID'] = df_inv['q_chr']+"_"+df_inv['q_pos'].astype(str)+'_'+df_inv['q_end'].astype(str)

		df_inv_out['SVLEN'] = df_inv_out['r_end']-df_inv_out['r_pos']
		df_inv_out['ID'] = df_inv_out['r_chr']+'-'+df_inv_out['r_pos'].astype(str)+'-INV-'+df_inv_out['SVLEN'].astype(str)
		df_inv_out['SVTYPE'] = 'INV'
		df_inv_out['CONTIG_SEQ'] = df_inv_out['q_chr']+":"+df_inv_out['q_pos'].astype(str)+"-"+df_inv_out['q_end'].astype(str)

		df_inv_out[['r_chr', 'r_pos', 'r_end', 'ID', 'SVLEN', 'SVTYPE', 'CONTIG_SEQ']].rename(columns={'r_chr' : '#CHROM', 'r_pos' : 'POS', 'r_end' : 'END'}).to_csv(output.inv, sep='\t',index=False)


rule call_transposition:
	input:
		broken = 'tmp/{sample}-broken.paf',
		sv = rules.call_simple_insdel.output.tsv,
		gap = rules.call_gaps.output.gap,
		inv = rules.call_inv.output.inv
	output:
		trans = 'tmp/{sample}_transp-{direct}.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		dir_dict = {'rev' : {'q_pos' : 'q_end', 'q_end' : 'q_pos', 'r_pos' : 'r_end', 'r_end' : 'r_pos'}, 'for' : {'q_pos' : 'q_pos', 'q_end' : 'q_end', 'r_pos' : 'r_pos', 'r_end' : 'r_end'}}
		sort_dict = {'for' : True, 'rev' : False} 

		df = pd.read_csv(input.broken, sep='\t', header=None, names=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_len', 'r_pos', 'r_end', 'nmatch', 'alen', 'mapq', 'id', 'cg'])

		df = df.sort_values(['r_chr', 'r_pos']).reset_index(drop=True)

		# df['browse'] = df['r_chr']+":"+df['r_pos'].astype(str)+"-"+df['r_end'].astype(str)
		# sv_df['browse'] = sv_df['#CHROM']+":"+sv_df['POS'].astype(str)+"-"+sv_df['END'].astype(str)

		sv_df = pd.concat( [ pd.read_csv(x, sep='\t') for x in [input.sv, input.inv] ] )

		gap_df = pd.read_csv(input.gap, sep='\t', header=None, names=['#CHROM', 'POS', 'END'])

		sv_bed = BedTool.from_dataframe(sv_df[['#CHROM', 'POS', 'END']].append(gap_df))

		for contig in df['q_chr'].unique():
			last_trans = False
			trans_len = 0
			last_syn_id = 0
			prev_inter = False
			contig_df = df.loc[df['q_chr'] == contig].sort_values(['q_pos'], ascending=sort_dict[wildcards.direct]).copy()
			contig_chr = contig_df.groupby('r_chr').sum().reset_index().sort_values(['alen'], ascending=False).iloc[0]['r_chr']
			for i, index in enumerate(contig_df.index):
				if contig_df.iloc[i]['q_len'] < 1000000:
					df.at[index, f'TRANS_{wildcards.direct}'] = 'SMALL_CTG'
					prev_inter = False
				elif (contig_df.iloc[i]['r_chr'] != contig_chr):
					df.at[index, f'TRANS_{wildcards.direct}'] = 'INTER'
					prev_inter = True
				elif df.at[index, 'strand'] == '-':
					df.at[index, f'TRANS_{wildcards.direct}'] = 'INV'
					prev_inter = False
				else:
					if prev_inter == True and i != len(contig_df)-1 and i != 0:
						if len(contig_df) == 1:
							df.at[index, f'TRANS_{wildcards.direct}'] = 'SYNTENIC'
							prev_inter = False
						else:
							check_up = i+1
							end_pred = contig_df.iloc[check_up][dir_dict[wildcards.direct]['r_pos']]
							start_pred = max(0, end_pred-contig_df.iloc[i]['nmatch'])
							if np.abs(start_pred-contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']]) < 30000:
								df.at[index, f'TRANS_{wildcards.direct}'] = 'SYNTENIC'
								prev_inter = False
								last_trans = False
								last_syn_id = i
							else:
								gap_bed = BedTool.from_dataframe(pd.DataFrame.from_dict({'#chr' : [contig_chr], 'start' : [min(contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']], start_pred)], 'end' : [max(contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']], start_pred)]}))
								cov_df = gap_bed.coverage(sv_bed).to_dataframe(names=['#chr', 'start', 'end', 'events', 'bases', 'len', 'perc'])
								if cov_df.iloc[0]['perc'] > .90:
									df.at[index, f'TRANS_{wildcards.direct}'] = 'SYNTENIC'
									prev_inter = False
									last_trans = False
									last_syn_id = i
								else:
									df.at[index, f'TRANS_{wildcards.direct}'] = 'TRANSPOSE'
									prev_inter = False
									last_trans = True
									trans_len = df.at[index, 'alen']
					elif i == len(contig_df)-1 or i == 0:
						df.at[index, f'TRANS_{wildcards.direct}'] = 'EDGE'
					elif i != 0 and i != len(contig_df)-1:
						if last_trans == False:
							check_down = i-1
							start_pred = contig_df.iloc[check_down][dir_dict[wildcards.direct]['r_end']]
						else:
							check_down = last_syn_id
							start_pred = contig_df.iloc[check_down][dir_dict[wildcards.direct]['r_end']]+trans_len
						check_up = i+1
						end_pred = contig_df.iloc[check_up][dir_dict[wildcards.direct]['r_pos']]
						if np.abs(start_pred-contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']]) < 30000 or contig_df.iloc[i]['alen'] > 1000000:
							df.at[index, f'TRANS_{wildcards.direct}'] = 'SYNTENIC'
							prev_inter = False
							last_trans = False
							last_syn_id = i
						elif (contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']] < contig_df.iloc[check_down][dir_dict[wildcards.direct]['r_end']] and contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']] > contig_df.iloc[check_down][dir_dict[wildcards.direct]['r_end']]) or (contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']] < contig_df.iloc[check_up][dir_dict[wildcards.direct]['r_end']] and contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']] > contig_df.iloc[check_up][dir_dict[wildcards.direct]['r_end']]):
							df.at[index, f'TRANS_{wildcards.direct}'] = 'SYN-DUP'
							prev_inter = False
							last_trans = False
							last_syn_id = i
						else:
							gap_bed = BedTool.from_dataframe(pd.DataFrame.from_dict({'#chr' : [contig_chr], 'start' : [min(contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']], start_pred)], 'end' : [max(contig_df.iloc[i][dir_dict[wildcards.direct]['r_pos']], start_pred)]}))
							cov_df = gap_bed.coverage(sv_bed).to_dataframe(names=['#chr', 'start', 'end', 'events', 'bases', 'len', 'perc'])
							if cov_df.iloc[0]['perc'] > .90:
								df.at[index, f'TRANS_{wildcards.direct}'] = 'SYNTENIC'
								prev_inter = False
								last_trans = False
								last_syn_id = i
							else:
								df.at[index, f'TRANS_{wildcards.direct}'] = 'TRANSPOSE'
								last_trans = True
								prev_inter = False
								trans_len = df.at[index, 'alen']
		df.to_csv(output.trans, sep='\t', index=False)


rule combine_trans_search:
	input:
		forward = 'tmp/{sample}_transp-for.bed',
		rev = 'tmp/{sample}_transp-rev.bed'
	output:
		trans = 'tmp/{sample}_transp.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		df = pd.merge(pd.read_csv(input.forward, sep='\t'), pd.read_csv(input.rev, sep='\t'))

		df['TRANS'] = df.apply(lambda row: 'SYNTENIC' if 'SYNTENIC' in [row['TRANS_for'], row['TRANS_rev']] else row['TRANS_for'], axis=1)

		df['SVLEN'] = df['alen']
		df['SVTYPE'] = df['TRANS']
		df['ID'] = df['r_chr']+'-'+df['r_pos'].astype('str')+'-'+df['TRANS']+'-'+df['alen'].astype(str)
		df['CONTIG_SEQ'] = df['q_chr']+":"+df['q_pos'].astype(str)+"-"+df['q_end'].astype(str)

		df = df.rename(columns={'r_chr' : '#CHROM', 'r_pos' : 'POS', 'r_end' : 'END', 'q_chr' : 'CONTIG'})

		df.loc[(df['TRANS'] != 'SYNTENIC') & (df['TRANS'] != 'EDGE') & (df['TRANS'] != 'SMALL_CTG')][['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID', 'CONTIG_SEQ']].to_csv(output.trans, sep='\t', index=False)

rule call_from_gaps:
	input:
		orient = 'tmp/{sample}-orient.paf',
		gap = rules.call_gaps.output.gap
	output:
		gap_call = 'tmp/{sample}_gapcall.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.gap, sep='\t', header=None, names=['#CHROM', 'POS', 'END'])

		df['browse'] = df['#CHROM']+":"+df['POS'].astype(str)+"-"+df['END'].astype(str)

		df['GAP_LEN'] = df['END'] - df['POS']

		contig_df = pd.read_csv(input.orient, sep='\t', header=None, names=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_len', 'r_pos', 'r_end', 'nmatch', 'alen', 'mapq', 'id', 'cg'])

		contig_df = contig_df.loc[contig_df['q_len'] > 1000000]

		contig_df['browse'] = contig_df['r_chr']+":"+contig_df['r_pos'].astype(str)+"-"+contig_df['r_end'].astype(str)

		strand_dict_end = {'+' : 'q_pos', '-' : 'q_end'}
		strand_dict_pos = {'+' : 'q_end', '-' : 'q_pos'}

		for index in df.index:
			chrom = df.at[index, '#CHROM']
			gap_start = max(0, df.at[index, 'POS']-50)
			gap_end = df.at[index, 'END']+50
			contig_chr = contig_df.loc[(contig_df['r_chr'] == chrom) & (((contig_df['r_end'] >= gap_start) & (contig_df['r_end'] <= gap_end)) | ((contig_df['r_pos'] >= gap_start) & (contig_df['r_pos'] <= gap_end)))].copy().sort_values('r_pos')
			if len(contig_chr) < 2 or len(contig_chr['q_chr'].unique()) > 1:
				df.at[index, 'GAP_STATUS'] = 'GAP'
				df.at[index, 'CONTIG_LEN'] = 0
			elif len(contig_chr) > 2:
				df.at[index, 'GAP_STATUS'] = 'COMPLEX_ALIGNMENT'
				df.at[index, 'CONTIG_LEN'] = 0
			else:
				q_start = contig_chr.iloc[0]['q_end']
				q_end = contig_chr.iloc[1]['q_pos']
				mid_bases = contig_df.loc[(contig_df['q_end'] <= q_end) & (contig_df['q_pos'] >= q_start)]
				df.at[index, 'CONTIG_LEN'] = q_end - q_start - q_end		
				df.at[index, 'q_end'] = contig_chr.iloc[1][strand_dict_end[contig_chr.iloc[1]['strand']]]
				df.at[index, 'q_pos'] = contig_chr.iloc[0][strand_dict_pos[contig_chr.iloc[0]['strand']]]
				df.at[index, 'q_chr'] = contig_chr.iloc[0]['q_chr']
				if df.at[index, 'CONTIG_LEN'] < df.at[index, 'GAP_LEN']:
					df.at[index, 'GAP_STATUS'] = 'DEL_GAP'
				elif df.at[index, 'CONTIG_LEN'] > df.at[index, 'GAP_LEN']:
					df.at[index, 'GAP_STATUS'] = 'DUP_GAP'
				else:
					df.at[index, 'GAP_STATUS'] = 'GAP_NONSYN'

		df = df.dropna(subset=['q_end', 'q_pos', 'q_chr'])
		df = df.loc[(df['CONTIG_LEN'] >= 30000) | (df['GAP_LEN'] >= 30000)]
		df = df.loc[(df['GAP_STATUS'] != "GAP") & (df['GAP_STATUS'] != "COMPLEX_ALIGNMENT")]
		df['SVLEN'] = np.abs(df['GAP_LEN'] - df['CONTIG_LEN']).astype(int).astype(str)
		df['SVTYPE'] = df['GAP_STATUS']
		df['ID'] = df['#CHROM']+'-'+df['POS'].astype(str)+'-'+df['GAP_STATUS']+'-'+np.abs(df['GAP_LEN'] - df['CONTIG_LEN']).astype(int).astype(str)
		df['CONTIG_SEQ'] = df['q_chr']+":"+df['q_pos'].astype(int).astype(str)+"-"+df['q_end'].astype(int).astype(str)
		df[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID', 'CONTIG_SEQ']].to_csv(output.gap_call, sep='\t', index=False)


rule call_dup:
	input:
		broken = 'tmp/{sample}-broken.paf'
	output:
		dup = 'tmp/{sample}_dup.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.broken, sep='\t', header=None, names=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_len', 'r_pos', 'r_end', 'nmatch', 'alen', 'mapq', 'id', 'cg'])

		df = df.loc[df['q_len'] > 1000000].reset_index(drop=False)

		dup_bed = BedTool.from_dataframe(df[['r_chr', 'r_pos', 'r_end', 'index', 'q_chr', 'q_pos', 'q_end']])
		dup_df = dup_bed.intersect(dup_bed, wa=True, wb=True).to_dataframe(names=['q_chr', 'q_pos', 'q_end', 'orig_index', 'orig_contig', 'orig_pos', 'orig_end', 'dup_chr', 'dup_pos', 'dup_end', 'dup_index', 'dup_contig', 'dup_ctg_pos', 'dup_ctg_end'])

		dup_df = dup_df.loc[(dup_df['orig_index'] != dup_df['dup_index']) & (dup_df['orig_contig'] == dup_df['dup_contig'])].reset_index(drop=True)

		dup_bed = pd.DataFrame()

		dup_df['final_dup_pos'] = dup_df.apply(lambda row: max(row['q_pos'], row['dup_pos']), axis=1)
		dup_df['final_dup_end'] = dup_df.apply(lambda row: min(row['q_end'], row['dup_end']), axis=1)

		dup_df['SVLEN'] = dup_df['final_dup_end'] - dup_df['final_dup_pos']

		# dup_df['browse'] = dup_df['q_chr']+":"+dup_df['final_dup_pos'].astype(str)+"-"+dup_df['final_dup_end'].astype(str)

		dup_df = dup_df.loc[dup_df['SVLEN'] > 30000].drop_duplicates(['q_chr','final_dup_pos', 'final_dup_end'])

		dup_df['SVTYPE'] = 'DUP'

		dup_df['ID'] = dup_df['q_chr']+"-"+dup_df['final_dup_pos'].astype(str)+"-"+dup_df['SVTYPE']+"-"+dup_df['SVLEN'].astype(str)
		dup_df['CONTIG_SEQ'] = dup_df['orig_contig']+":"+dup_df['orig_pos'].astype(str)+"-"+dup_df['orig_end'].astype(str)+","+dup_df['dup_contig']+":"+dup_df['dup_ctg_pos'].astype(str)+"-"+dup_df['dup_ctg_end'].astype(str)

		dup_df = dup_df.rename(columns={'q_chr' : '#CHROM', 'final_dup_pos' : 'POS', 'final_dup_end' : 'END', 'orig_contig' : 'CTG'})

		dup_df[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID', 'CONTIG_SEQ']].sort_values(['#CHROM', 'POS']).to_csv(output.dup, sep='\t', index=False)


rule combine_calls:
	input:
		bed = expand('tmp/{{sample}}_{char}.bed', char=['gapcall', 'transp', 'sv_insdel', 'sv_inv', 'dup'])
	output:
		bed = 'results/raw/{sample}/{sample}_var_all.bed'
	threads: 1
	resources: 
		mem = 8,
		hrs = 24
	run:
		df = pd.concat([pd.read_csv(x, sep='\t', usecols=['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID', 'CONTIG_SEQ']) for x in input.bed ] )
		df['SVTYPE'] = df['ID'].str.split('-', expand=True)[2]
		df['SVLEN'] = df['ID'].str.split('-', expand=True)[3].astype(int)
		df = df.drop_duplicates().sort_values(['#CHROM', 'POS'])
		df.to_csv(output.bed, sep='\t', index=False)



rule var_join:
	input:
		bed = rules.combine_calls.output.bed,
		broken = 'tmp/{sample}-broken.paf'
	output:
		bed = 'results/merged/{sample}/{sample}_var_all.bed'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.bed, sep='\t')
		broken_df = pd.read_csv(input.broken, sep='\t', header=None, usecols=[0,1,2,3,4,5,7,8], names=['q_chr', 'q_len', 'q_pos', 'q_end', 'strand', 'r_chr', 'r_pos', 'r_end'])
		# Only call variants in contigs > 1Mbp
		broken_df = broken_df.loc[broken_df['q_len'] > 1000000]

		merge_pairs = {'DEL' : ['DEL', 'DEL_GAP'], 'DUP' : ['DUP_GAP', 'DUP'], 'DUP_GAP' : ['DUP_GAP', 'DUP'], 'DEL_GAP' : ['DEL', 'DEL_GAP'], 'TRANSPOSE' : ['TRANSPOSE'], 'INTER' : ['INTER'], 'INV' : ['INV'], 'DUP' : ['DUP']}

		df = df.sort_values(['#CHROM', 'POS'])

		gap_df = pd.DataFrame()

		chrom = df.iloc[0]['#CHROM']

		broken_bed = BedTool.from_dataframe(broken_df[['r_chr', 'r_pos', 'r_end']])

		merge_dict = {}


		for i in range(len(df)):
			merge_dict[i] = []
			if i == len(df)-1:
				continue
			elif df.iloc[i+1]['#CHROM'] != chrom:
				chrom = df.iloc[i+1]['#CHROM']
				continue
			else:
				if np.abs(df.iloc[i]['END']+1 - df.iloc[i+1]['POS']-1) < 5000 and df.iloc[i+1]['SVTYPE'] in merge_pairs[df.iloc[i]['SVTYPE']]:
					merge_dict = merge_cat(merge_dict, i)
				else: 
					gap_bed = BedTool.from_dataframe(pd.DataFrame.from_dict({'#CHROM' : [chrom], 'POS' : [max(0,min(df.iloc[i]['END']+1, df.iloc[i+1]['POS']-1))], 'END' : [max(df.iloc[i]['END']+1, df.iloc[i+1]['POS']-1)]}))
					if len(gap_bed.subtract(broken_bed, A=True).to_dataframe()) == 0:
						continue
					elif df.iloc[i+1]['SVTYPE'] in merge_pairs[df.iloc[i]['SVTYPE']]:
						merge_dict = merge_cat(merge_dict, i)

		# gap_df['LEN_GAP'] = gap_df['END'] - gap_df['POS']

		for i in merge_dict:
			merge_dict[i].sort()

		g = nx.DiGraph()
		g.add_nodes_from(merge_dict.keys())
		for k, v in merge_dict.items():
			g.add_edges_from(([(k, t) for t in v]))

		UG = g.to_undirected()

		sub_graphs = sorted(nx.connected_components(UG), key = len, reverse=True)

		merge_df = pd.DataFrame()

		for node in sub_graphs:
			if len(node) == 1:
				merge_df = merge_df.append(df.iloc[list(node)[0]])
			else:
				merges = df.iloc[min(node):max(node)].copy()
				sv_start = min(merges['POS'])
				sv_chrom = merges.iloc[0]['#CHROM']
				sv_end = max(merges['END'])
				merges['SVTYPE'] = merges['SVTYPE'].str.replace('_GAP', '')
				if len(merges['SVTYPE'].unique()) == 1:
					sv_type = merges['SVTYPE'].unique()[0]
					sv_merge_type = merges['SVTYPE'].unique()[0]
				else:
					sv_type = 'COMPLEX'
					sv_merge_type = ";".join(merges['SVTYPE'].unique())
				sv_merge_id = ";".join(merges['ID'])
				sv_len = np.sum(merges['SVLEN'])
				sv_merge_contig = ";".join(merges['CONTIG_SEQ'])
				sv_id = f'{sv_chrom}-{sv_start}-{sv_type}-{sv_len}'
				merge_df = merge_df.append(pd.DataFrame.from_dict({'#CHROM' : [sv_chrom], 'POS' : [sv_start], 'END' : [sv_end], 'SVTYPE' : [sv_type], 'SVTYPE_MERGE' : [sv_merge_type], 'ID' : [sv_id], 'SVLEN' : [sv_len], 'ID_MERGE' : [sv_merge_id], 'CONTIG_SEQ_MERGE' : [sv_merge_contig]}))

		merge_df = merge_df.fillna('SINGLE')
		merge_df['SVTYPE'] = merge_df['SVTYPE'].str.replace('_GAP', '')

		merge_df['SVTYPE_MERGE'] = merge_df.apply(lambda row: row['SVTYPE'] if row['SVTYPE_MERGE'] == 'SINGLE' else row['SVTYPE_MERGE'], axis=1)
		merge_df['ID_MERGE'] = merge_df.apply(lambda row: row['ID'] if row['ID_MERGE'] == 'SINGLE' else row['ID_MERGE'], axis=1)
		merge_df['CONTIG_SEQ_MERGE'] = merge_df.apply(lambda row: row['CONTIG_SEQ'] if row['CONTIG_SEQ_MERGE'] == 'SINGLE' else row['CONTIG_SEQ_MERGE'], axis=1)


		merge_df[['#CHROM','POS', 'END', 'SVTYPE_MERGE', 'SVLEN', 'ID_MERGE', 'CONTIG_SEQ_MERGE']].sort_values(['#CHROM', 'POS']).to_csv(output.bed, index=False, sep='\t')


rule bed9:
	input:
		bed = rules.var_join.output.bed
	output:
		bed = 'results/bed9/{sample}/{sample}_saf.bed9'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.bed, sep='\t')
		color_dict = {'DEL' : '220,0,0', 'DUP' : '0,0,220', 'INTER' : '0,220,0', 'INV' : '220,140,0', 'TRANSPOSE' : '128,0,128'}
		df['COLOR'] = df['SVTYPE_MERGE'].apply(lambda x : color_dict[x])
		df['BED_END'] = df.apply(lambda row: row['END']+row['SVLEN'] if row['SVTYPE_MERGE'] == 'DUP' else row['END'], axis=1)
		df['SCORE'] = '0'
		df['STRAND'] = '+'
		df[['#CHROM', 'POS', 'BED_END', 'ID_MERGE', 'SCORE', 'STRAND', 'SCORE', 'SCORE', 'COLOR']].sort_values(['#CHROM', 'POS']).to_csv(output.bed, sep='\t', header=['#ct','st','en','name','score','strand','tst','ten','color'], index=False)


rule saff_out:
	input:
		paf = 'tmp/{sample}-broken.paf'
	output:
		saf = 'results/saffire/{sample}/{sample}.saf'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	shell:
		'''
		rb stats --paf {input.paf} > {output.saf}
		'''
