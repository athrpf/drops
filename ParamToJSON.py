"""
ParamToJSON - crappy script to translate old param-files to JSON-format
no release, no warranty, no license
Autor: Thorolf Schulte (LNM)
"""

#settings
FILES = ["./stokes/stokes"]
SUFFIX_IN = ".param"
SUFFIX_OUT = ".json"
DEBUG = False
#/settings


#some static variables
STARTING_COMMENTS = True
VAR_START = ""
LEVEL_START = ""
level = 0
comment = ""

for filename in FILES:
	try:
		infile = open(filename + SUFFIX_IN)
		outfile = open(filename + SUFFIX_OUT, 'w')
	except:
		if not infile:
			print("Cannot open " + filename + SUFFIX_IN + ". Skipping.")
		if not outfile:
			print("Cannot open " + filename + SUFFIX_OUT + ". Skipping.")
		continue

	#start file with first level of hierarchy
	outfile.write("{\n")
	level = 1

	for line in infile:
		#starting comments in first level
		if STARTING_COMMENTS:
			if (line[0] == "#"):
				comment = comment + line
				continue
			#empty line?
			elif (line.strip() == ""): 
				continue
		
			STARTING_COMMENTS = False
			outfile.write(level*"\t" + '"_comment":\n"' + comment.rstrip() + '",\n\n')
			comment = ""

		#comment?
		if (line[0] == "#"):
			comment = comment + line
			if DEBUG: print("COMMENT: \t" + line)
			continue

		#remove inline comments without saving them
		line = (line.lstrip(" \n")).partition("#")[0]
		#new hierarchy level?
		if ("{" in line):
			if DEBUG: print("NEW LEVEL: \t" + line)
			name = (line.partition("{")[0]).strip()
			outfile.write(LEVEL_START + level*"\t" + '"' + name + '":\n' + level*"\t" + "{\n")
			level = level + 1
			#add comments fetched before
			if (comment != ""):
				comment.rstrip("\n")
				outfile.write(level*"\t" + '"_comment":\n"' + comment.rstrip() + '",\n\n')
				comment = ""
			VAR_START = ""
			LEVEL_START = ",\n\n"
		#new variable?
		elif("=" in line):
			if DEBUG: print("VARIABLE: \t" + line)
			var_line = line.partition("=")
			var_str = ""
			#what type the variable is? string needs "...", vector [...]
			try:
				#try number
				float(var_line[2])
				var_str = (var_line[2]).strip()
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
						for val in var_vec:
							var_str = var_str + " " + val
						var_str = '[' + var_str + ' ]'
			#still empty? String!
			if var_str == "":
				var_str = '"' + (var_line[2]).strip() + '"'

			var_name = (var_line[0]).strip()
			#CHANGED VARIABLES
			if (var_name == "InitialCond" and "poisson" in filename):
				var_name = "RefineSteps"

			#last variable cannot end with "," so it's added for line before
			outfile.write(VAR_START + level*"\t" + '"' + var_name + '":\t\t' + var_str)
			VAR_START = ",\n"
		#hierarchy level closed?
		elif ("}" in line):
			if DEBUG: print("CLOSED LEVEL: \t" + line)
			level = level-1
			outfile.write("\n" + level*"\t" + "}")
	
	outfile.write("\n\n}")
	outfile.close()
	infile.close()
	STARTING_COMMENTS = True
		
