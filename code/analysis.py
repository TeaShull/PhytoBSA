from analysis_functions import *

print("""
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Producing plots and identifying putative causal mutations
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
""")


# Establish variables

for key in lines_dict:

	# Read in ems table, calculate metrics, format
	vcftable = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".noknownsnps.table"
	)

	df = pd.read_csv(vcftable, sep="\t")

	## Calculate delta SNP ratio
	df['ratio'] = deltaSNParray(
		df['wt_ref'], 
		df['wt_alt'],
		df['mu_ref'],
		df['mu_alt']
	)

	## Calculate G-statistic to glean more info about potential causal SNPs
	df['G_S'] = gStatistic_Array(
		df['wt_ref'], 
		df['wt_alt'], 
		df['mu_ref'], 
		df['mu_alt']
	)

	## Drop indels, leave only SNPs
	df[
	    df["ref"].apply(lambda x: len(x) == 1)
	    & df["alt"].apply(lambda x: len(x) == 1)
	]

	## Drop NAs. Infrequently will arise while div/0 seemingly
	df.dropna(axis=0, how='any', subset="ratio", inplace=True)

	# lowess Smoothing of ratio and G-statistic by chromosome
	chr_facets=df["chr"].unique()
	df_list = []

	## Set lowess parameters
	lowess = sm.nonparametric.lowess
	lowess_span=0.3

	for i in chr_facets:
	    # Establish chr copies and smooth edges by inverting duplicated data
	    df_chr = df[df['chr']==i].copy()

	    ## establish inverted dataframes

	    positions = df_chr['pos'].to_numpy()
	    deltas=[]
	    # Produce deltas
	    for i, pos in enumerate(positions):
	        if i == 0:
	            deltas.append(pos)
	        if i > 0:
	            deltas.append(positions[i] - positions[i-1])
	    
	    ## Create negative and positive delta extensions
	    deltas_pos_inv = deltas[::-1][-15:-1]
	    deltas_neg_inv = deltas[::-1][1:15]

	    ## Extend deltas, making "mirrored" delta positions at ends of chr
	    deltas_mirrored_ends=[]
	    deltas_mirrored_ends.extend(deltas_pos_inv+deltas+deltas_neg_inv)
	    
	    ## Create new pseudo positions for lowess smoothing 
	    psuedo_pos = []
	    for i, pos in enumerate(deltas_mirrored_ends):
	        if i == 0:
	            psuedo_pos.append(0)
	        if i > 0:
	            psuedo_pos.append(pos+psuedo_pos[i-1])

	    ## Create "mirrored" data from both ends of chr
	    df_chr_inv_neg = df_chr[::-1].iloc[-15:-1]
	    df_chr_inv_pos = df_chr[::-1].iloc[1:15]

	    ## Concat dataframe
	    df_chr_smooth_list = [df_chr_inv_neg, df_chr, df_chr_inv_pos]

	    df_chr = pd.concat(df_chr_smooth_list, ignore_index=False)
	            
	    ## Add psuedo positions to dataframe
	    df_chr['psuedo_pos'] = psuedo_pos

	    ## Fit LOWESS ratio and G-statistic ~ position
	    X=df_chr['psuedo_pos'].values

	    ## Fit ratio 
	    ratio_Y=df_chr['ratio'].values
	    df_chr['ratio_yhat'] = lowess(ratio_Y,X, frac=lowess_span)[:,1]
	   
	    ## Fit G-Statistic
	    G_S_Y = df_chr['G_S'].values
	    df_chr['G_S_yhat'] = lowess(G_S_Y,X, frac=lowess_span)[:,1]
	    
	    ## Produce ratio-scaled G-statistic and fit
	    df_chr['RS_G'] = df_chr['G_S']*df_chr['ratio']
	    
	    RS_G_Y = df_chr['RS_G'].values
	    df_chr['RS_G_yhat'] = lowess(RS_G_Y,X, frac=lowess_span)[:,1]

	    ## remove edge smoothing (mirrored) data
	    df_chr = df_chr[14:-14]
	    df_chr.drop(axis=1, inplace=True, columns='psuedo_pos')
	    df_list.append(df_chr)

	df = pd.concat(df_list)

	# Randomize data 1000 times, calculate metrics to derive cutoffs
	## subset position, wt and mu
	dfShPos = df[['pos']].copy()
	dfShwt = df[['wt_ref', 'wt_alt']].copy()
	dfShmu = df[['mu_ref', 'mu_alt']].copy()

	## Establish lists
	smGstatAll = []
	smRatioAll = []
	RS_GAll = []
	smRS_G_yhatAll = []
	smRatio_yhatAll = []

	## Break the phenotype / read count / position association.
	## Iterated and smooth 1000 times 
	for i in range(1000):
		shuffle(dfShPos, dfShwt, dfShmu)

	# Collect emperical cutoffs from 1000x iter, identify candidates
	G_S_95p = np.percentile(smGstatAll, 95)
	RS_G_95p = np.percentile(RS_GAll, 95)
	RS_G_Y_99p = np.percentile(smRS_G_yhatAll, 99.99)

	df['G_S_001p'] = [1 if (np.isclose(x, G_S_95p) or (x > G_S_95p))
	                else 0 for x in df['G_S']]

	df['RS_G_001p'] = [1 if (np.isclose(x, G_S_95p) or (x > G_S_95p))
	                else 0 for x in df['G_S']]

	df['RS_G_yhat_001p'] = [1 if (np.isclose(x, G_S_95p) or (x > G_S_95p))
	                else 0 for x in df['G_S']]

	df_likely_cands = df.loc[df['RS_G_yhat_001p'] == 1]

	## Use 99.99 percentile cutoff to identify likely candidates
	df_likely_cands_sorted = df_likely_cands.sort_values(
	    by = ['G_S','RS_G_yhat'], 
	    ascending = [False, False], 
	    na_position = 'first'
	    )

	# Save complete table
	finalcompletetablename = os.path.join(
		"./output/" 
		+ key + "/" 
		+ key 
		+ ".complete_table" 
		+ ".tsv"
	)
	
	df.to_csv(finalcompletetablename, sep='\t', index=False)

	## save likely candidate list	
	finallikelycandidatestablename = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".likely_candidates" 
		+ ".tsv"
	)

	df_likely_cands.to_csv(
		finallikelycandidatestablename, 
		sep='\t', 
		index=False
	)

	# Generate Plots
	warnings.filterwarnings( "ignore", module = "plotnine\..*" )
	df['pos_mb'] = df['pos']*0.000001
	## Establish plot names
	GS_plotname = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".GS" 
		+ ".png"
	)
	GS_yhat_plotname = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".GS_yhat" 
		+ ".png"
	)

	RS_GS_plotname = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".RS_GS" 
		+ ".png"
	)
	RS_GS_yhat_plotname = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".RS_GS" 
		+ ".png"
	)

	deltaSNP = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".deltaSNP" 
		+ ".png"
	)
	deltaSNP_yhat_plotname = os.path.join(
		"./output/" 
		+ key 
		+ "/" 
		+ key 
		+ ".deltaSNP_yhat" 
		+ ".png"
	)

	## G-stat plot
	chart = ggplot(df, aes('pos_mb', y='G_S'))
	title = ggtitle("G-stastic")
	points = geom_point(color='goldenrod', size=0.8)
	lines = geom_line(color='blue')
	themes = theme_linedraw()
	facets = facet_grid('. ~ chr', space='free_x', scales='free_x')
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("G-statistic")
	spacing = theme(panel_spacing=0.025)
	plot = (chart 
	    + points 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	)

	plot.save(
		filename = GS_plotname, 
		height=6, 
		width=8, 
		units = 'in', 
		dpi=500
	)

	#G-stat yhat plot
	chart = ggplot(df, aes('pos_mb', y='G_S_yhat'))
	title = ggtitle("Lowess smoothed G-stastic")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Fitted G-statistic")
	plot = (chart 
	    + points
	    + lines
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing

	)
	plot.save(
		filename = GS_yhat_plotname, 
		height=6, 
		width=8, 
		units = 'in', 
		dpi=500
	)
	
	## Ratio-scaled G-Stat plot
	chart = ggplot(df, aes('pos_mb', y='RS_G'))
	title = ggtitle("Ratio-scaled G statistic")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Ratio-scaled G-statistic")

	plot = (chart 
	    + points 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	)
	plot.save(
		filename = RS_GS_plotname, 
		height=6, 
		width=8, 
		units = 'in', 
		dpi=500
	)

	## Fitted Ratio-Scaled G-stat plot
	chart = ggplot(df, aes('pos_mb', y='RS_G_yhat'))
	title = ggtitle("Lowess smoothed ratio-scaled G statistic")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Fitted Ratio-scaled G-statistic")
	cutoff = geom_hline(
		yintercept = RS_G_Y_99p, 
		color='red', 
		linetype="dashed", 
		size=0.3
		)
	
	plot = (chart 
	    + points 
	    + lines 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + cutoff
	    + spacing
	)

	plot.save(
		filename = RS_GS_yhat_plotname, 
		height=6, 
		width=8, 
		units = 'in', 
		dpi=500
	)


	## Delta SNP Ratio plot
	chart = ggplot(df, aes('pos_mb', y='ratio'))
	title = ggtitle("Delta SNP ratio")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Ratio")

	plot = (chart 
	    + points 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	)

	plot.save(
		filename = deltaSNP, 
		height=6, 
		width=8, 
		units = 'in', 
		dpi=500
	)

	## Fitted delta SNP ratio plot
	chart = ggplot(df, aes('pos_mb', y='ratio_yhat'))
	title = ggtitle("Delta SNP ratio")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Ratio")

	plot = (chart 
	    + points
	    + lines
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	)

	plot.save(
		filename = deltaSNP_yhat_plotname, 
		height=6, 
		width=8, 
		units = 'in', 
		dpi=500
	)

	print(f'''
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Results for {key} generated. 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
''')