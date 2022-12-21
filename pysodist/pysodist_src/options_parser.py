class OptionsParser(): #parses the options hyperparameter file, provides options via the dictionary self.options
	def __init__(self,options_file='options.txt'):
		self.options=dict()
		f=open(options_file,'r')
		file=f.read()
		f.close()
		lines=file.split('\n')
		for line in lines:
			split_line=line.split(' ')
			if len(split_line)>0 and split_line[0]=='->':
				self.options[split_line[1]]=split_line[3]
