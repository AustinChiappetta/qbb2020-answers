#!/usr/bin/env python 3

import sys
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt

def D_data_reader(D_values_txt):
	#go line by line through D.txt and parse
	D_file = open(D_values_txt)
	D_list = []
	for line in D_file:
		entry = line.strip('[]\n')
		fields = entry.split(',')
		codon = int(fields[0])
		D = fields[1]
		D = D.strip(' ')
		D = int(D)
		D_list += [[codon, D]]

	D_file.close()

	return(D_list)

def log2_data_reader(log2_values_txt):
	#go line by line through log2_dN_dS.txt and parse
	log2_file = open(log2_values_txt)
	log2_list = []
	for line in log2_file:
		entry = line.strip('[]\n')
		fields = entry.split(',')
		codon = int(fields[0])
		log2val = fields[1]
		log2val = log2val.strip(' ')
		log2val = float(log2val)
		log2_list += [[codon, log2val]]

	log2_file.close()

	counter = 1
	for sublist in log2_list:
		for num in range(0,1):
			if sublist[num] == 'NaN':
				print('Nan here:', sublist, 'counter is:', counter)
			counter += 1

	return(log2_list)

	#for sublist in D_list:
	#	print(sublist)

def calculate_z_scores(D_data):
	#Make pandas df: first column is codon, second column is D value, third column will be Z scores
	D_df = pd.DataFrame(D_data)
	D_df.columns = ['Codon #', 'D value']
	#print(D_df)
	D_array = D_df['D value']
	
	#print(D_df)
	#print(D_array)
	#print(type(D_array))

	#use stats.zscore(array) and add to df
	zscores = stats.zscore(D_array)
	D_df['Z scores'] = zscores
	#print(D_df)
	return(D_df)


def make_log2_df(log2_data):
	log2_df = pd.DataFrame(log2_data)
	log2_df.columns = ['Codon #', 'log2(dN/dS)']
	return(log2_df)

def merge_dfs(D_df, log2_df):
	#D_df.merge(log2_df, on='Codon #')
	#print(log2_df)
	#log2_rearranged_df = pd.DataFrame(log2_df['log2(dN/dS)', 'Codon #'])
	log2_col1 = log2_df['log2(dN/dS)']

	log2_col2 = log2_df['Codon #']
	
	log2_rearranged_df = log2_col1.to_frame()
	log2_rearranged_df['Codon #'] = log2_col2
	#print(log2_rearranged_df)


	#print(log2_df)
	#print(D_df)
	#log2_rearranged_df.merge(D_df, on='Codon #')
	#print(log2_rearranged_df)

	merged_df = pd.merge(log2_rearranged_df, D_df, on='Codon #')
	#print(merged_df)

	significant_zscores_df = merged_df[merged_df['Z scores'] < 0.05]
	insignificant_zscores_df = merged_df[merged_df['Z scores'] >= 0.05]

	#print(significant_zscores_df)
	#print(insignificant_zscores_df)

	return(significant_zscores_df, insignificant_zscores_df)
	
def make_plot(signif_data, insignif_data):
	#2. Plot the codon number vs log2(dN/dS): two scatter plot on the same plot
	# first plot will be dN/dS ratios with a Z score above 0.05
	# second plot will be dN/dS ratios with a Z score below 0.05, colored a different color

	sig_data = signif_data
	insig_data = insignif_data

	#print(sig_data)
	#print(insig_data)

	fig, ax = plt.subplots()
	ax.set_title("log2 ratio of dN/dS across the sequence\nRed points: p < 0.05")
	ax.set_xlabel("Codon #")
	ax.set_ylabel("log2(dN/dS)")
	plt.scatter(sig_data['Codon #'], sig_data['log2(dN/dS)'], color='b', s=10)
	plt.scatter(insig_data['Codon #'], insig_data['log2(dN/dS)'], color='r', s=10)

	fig.show()
	fig.savefig("log2_dN_dS_plot.png")
	plt.close(fig)

	#plt.scatter(X,Y1,color='k')
	#plt.scatter(X,Y2,color='g')
	#plt.show()






if __name__ == "__main__":
	D_values = D_data_reader('D.txt')
	log2_values = log2_data_reader('log2_dN_dS.txt')
	#sig_data, insig_data = calculate_z_scores(D_values, log2_values)
	D_df = calculate_z_scores(D_values)
	log2_df = make_log2_df(log2_values)
	sig_data, insig_data = merge_dfs(D_df, log2_df)
	make_plot(sig_data, insig_data)