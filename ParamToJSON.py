"""
ParamToJSON - crappy script to translate old param-files to JSON-format
no release, no warranty, no license
Autor: Thorolf Schulte (LNM)

Easy method for writing params.txt on linux (ubuntu):
find . -name *.param > params.txt
"""

from string import replace

#settings
FROM_FILE = True
DEBUG = False
SUFFIX_IN = ".param"
SUFFIX_OUT = "_dev.json"

FILES = []
if FROM_FILE:
	paramfile = open("params.txt")
	for line in paramfile:
		arr = line.split(".param")
		FILES.append(arr[0])
	paramfile.close()
else:
	FILES = ["./stokes/stokes"]

#/settings


#some static variables
VAR_START1 = ""
VAR_START2 = ""
LEVEL_START = ""
level = 0
fullComment = ""
inlineComment = ""
inlineCommentOld = ""

for filename in FILES:
	try:
		infile = open(filename + SUFFIX_IN)
	except:
		print("Cannot open " + filename + SUFFIX_IN + ". Skipping.")
		continue

	try:
		outfile = open(filename + SUFFIX_OUT, 'w')
	except:
		print("Cannot open " + filename + SUFFIX_OUT + ". Skipping.")
		continue

	#start file with first level of hierarchy
	outfile.write("{\n")
	level = 1

	for line in infile:

		#comments
		if ("#" in line):
			#inline or full?
			lineArray = line.partition("#")
			if (len(lineArray[0].strip()) == 0):
				fullComment = fullComment + level*"\t" + "//" + lineArray[2]
			elif (len(lineArray[2].strip()) != 0):
				inlineComment = " \t //" + lineArray[2].strip()
			line = lineArray[0]

		#new hierarchy level?
		if ("{" in line):
			if DEBUG: print("NEW LEVEL: \t" + line)
			name = (line.partition("{")[0]).strip()
			outfile.write(LEVEL_START + fullComment + level*"\t" + '"' + name + '":\n' + level*"\t" + "{ " + inlineComment + "\n")
			level = level + 1
			VAR_START1 = ""
			VAR_START2 = ""
			LEVEL_START = ",\n\n"
			fullComment = ""
			inlineComment = ""

		#new variable?
		elif("=" in line):
			if DEBUG: print("VARIABLE: \t" + line)
			var_line = line.partition("=")
			var_str = ""
			#what type the variable is? string needs "...", vector [...]
			try:
				#try number
				var_str = float(var_line[2])
				if ((var_str%1) == 0.0):
					var_str = int(var_str)
				var_str = str(var_str)
				#var_str = (var_line[2]).strip()
			except ValueError:
				#try vector
				var_vec = (var_line[2]).split()
				if len(var_vec) == 3:
					vec_try = True
					for val in var_vec:
						try:
							float(val.strip())
						except ValueError:
							vec_try = False
							break
					if vec_try:
						var_str = var_vec[0] + ", " + var_vec[1] + ", " + var_vec[2]
						var_str = '[' + var_str + ' ]'
			#still empty? String!
			if var_str == "":
				var_str = '"' + (var_line[2]).strip() + '"'

			var_name = (var_line[0]).strip()
			#CHANGED VARIABLES
			if (var_name == "InitialCond" and "poisson" in filename):
				var_name = "RefineSteps"

			#last variable cannot end with "," so it's added for line before
			outfile.write(VAR_START1 + inlineCommentOld + VAR_START2 + fullComment + level*"\t" + '"' + var_name + '":\t\t' + var_str)
			VAR_START1 = ","
			VAR_START2 = "\n"
			fullComment = ""
			inlineCommentOld = ""

		#hierarchy level closed?
		elif ("}" in line):
			if DEBUG: print("CLOSED LEVEL: \t" + line)
			level = level-1
			outfile.write(fullComment + inlineCommentOld + "\n" + level*"\t" + "}")
			fullComment = ""
			inlineCommentOld = ""

		inlineCommentOld = inlineCommentOld + inlineComment
		inlineComment = ""

	outfile.write("\n\n}")
	outfile.close()
	infile.close()
	VAR_START1 = ""
	VAR_START2 = ""
	LEVEL_START = ""
	level = 0
	inlinecomment = ""
	print(filename + " done.")
		
