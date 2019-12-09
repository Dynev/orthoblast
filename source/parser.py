import pandas

class Mode:
    def __init__(self, proteins, species, evalue, coverage, draw, outname, pdf):
        self.proteins = proteins
        self.species = species
        self.evalue = evalue
        self.coverage = coverage
        self.draw = draw
        self.outname = outname
        self.pdf = pdf

def getMode(args):
    if args[0] == '':
        return 'Please specify parameters. To run help, type "help"'
    if args[0] == 'help':
        return """Orthoblast Usage Guide:

  Orthoblast is a tool that blasts a set of human proteins vs a set of
  genomes, which are translated, then it counts regions of high homology
  and compares their numbers across genomes.
        
  -p  Proteins: In form of a list, e.g. "amy1 agtr1 bag1" or a path to .xlsx 
      file that has a column labeled "Proteins". This tool will download 
      fastas of specified proteins from uniprot and put them into "proteins"
      folder
        
  -s  Species: For species' names please use a special format:
      <3_first_letters_of_family_name>_<3_first_letters_of_species_name>
      For a bulk run, use a simple list, e.g. "hom_sap het_gla cav_por"
      To do comparisons, use parentheses to specify pairs of results to 
      compare, e.g. "(het_gla cav_por) (hom_sap pan_tro)". Blast databases
      are stored in "genomes" folder. You can make them manually, or the 
      tool will try to download and install the ones which aren't 
      already in the folder

  -e  Evalue for TBLAST (float)
          
  -c  Target percent coverage of query by hit to consider homologous (int)
      If not specified, program will use internal custom criteria
        
  -d  Create a graph for every result to review manually (boolean)
      Default = False (do not create graphs)
        
  -n  Name for a folder inside "results" in program directory
      where the output will be stored
  
  After every run, you can find your results in the folder "results" """
    
    args.append('-')
    proteins = []
    species = []
    pdf = None
    if '-p' in args:
        p = args[args.index('-p') + 1]
        if '.xlsx' in p:
            try:
                pdf = pandas.read_excel(p)
            except IOError:
                return f"Cannot access file {p}. Exiting"
            try:
                proteins = list(pdf['Proteins'])
            except KeyError:
                return f"No column 'Proteins' in {p}. Exiting"
        else:
            for i in range(args.index('-p')+1, len(args)-1):
                if '-' in args[i]:
                    break
                else:
                    proteins.append(args[i])
    else:
        return "Proteins (-p) are required to run. Exiting"
    if not proteins:
        return "Please specify proteins to run. Exiting"
    if '-s' in args:
        s = args[args.index('-s')+1]
        if '(' in s:
            for i in range(args.index('-s')+1, len(args)-1, 2):
                if '-' in args[i]:
                    break
                elif '(' in args[i] and ')' in args[i+1]:
                    species.append([args[i], args[i+1]])
                else:
                    return "Species not properly specified. Exiting"
        else:
            for i in range(args.index('-s')+1, len(args)-1):
                if '-' in args[i]:
                    break
                else:
                    species.append(args[i])
    else:
        return "Species (-s) are required to run. Exiting"
    if not species:
        return "Please specify species to run. Exiting"
    
    evalue = 10
    coverage = None
    draw = False
    outname = None
    if '-e' in args:
        e = args[args.index('-e')+1]
        if '-' not in e:
            try:
                evalue = float(e)
                if evalue <= 0:
                    return "Invalid evalue parameter. Exiting"
            except ValueError:
                return "Cannot convert evalue parameter to float. Exiting"
        else:
            return "Evalue parameter mentioned, but not specified. Exiting"
    if '-c' in args:
        c = args[args.index('-c')+1]
        if '-' not in c:
            try:
                coverage = float(c)
                if coverage < 0 or coverage > 1:
                    return "Invalid coverage parameter. Exiting"
            except ValueError:
                return "Cannot convert coverage parameter to float. Exiting"
        else:
            return "Coverage parameter mentioned, but not specified. Exiting"
    if '-d' in args:
        d = args[args.index('-d')+1]
        if '-' not in d:
            if d == 'True':    
                draw = True
            elif d == 'False':
                draw = False
            else:
                return "Cannot convert draw parameter to boolean. Exiting"
        else:
            return "Draw parameter mentioned, but not specified. Exiting"
    if '-n' in args:
        n = args[args.index('-n')+1]
        if '-' not in n:
            outname = n
        else:
            return "Outname parameter mentioned, but not specified. Exiting"
    
    return Mode(proteins, species, evalue, coverage, draw, outname, pdf)