#! /usr/bin/python -u
__author__  = "Marco Mariotti"
__email__   = "marco.mariotti@crg.eu"
__licence__ = "GPLv3"
__version__ = "0.4"
import sys
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
from ete2_MM import *
import traceback
help_msg="""This program works as a codeml wrapper to run analyses on a set of aligned coding sequences and display the results in the form of a tree.

Many options accept nodes as arguments. These can be specified by their name if they are leaves, or using the first common ancestor principle if they are ancestral ( node1name+node2name ). Multiple nodes can be specified joining them with a comma (  node1name,node2name+node3name,node5name ). 
Some nodes can be specified with keywords:  root, all, all_leaves

## input
-a   CDS alignment file in aligned fasta
-t   phylogenetic tree, with names matching the first words of the fasta headers
-o   output prefix for this project. The file prefix.cdn_ali and prefix.codeml will be created, or loaded if found already present. In this last case, options -a and -t are not compulsory

## analysis with codeml models
-w_change  [nodes]        perform a w_change analysis on a number of nodes defined by arg (no arg->all). This analysis consists in comparing a null, fixed omega model with a w_change model with a LRT test. In the w_change model, all nodes below the tested node (including it) are assumed to have the same omega, which is different from the rest of the tree.
-w_change2 [nodes]        same as w_change, but this tests the change in omega affecting a single node (and not its descendants, as w_change)
-free_w                   compute omegas across all tree using a free ratio model, rich in parameters and long to execute

## sequence analysis 
-KaKs   | -Ks |  -Ka [nodes]  count sites and syn/nonsyn changes to compute the proportion of those fixed among the possible ones, considering all leaves below the considered nodes (no arg->all)
-dKaKs  |-dKs | -dKa [nodes]  compute these indexes comparing each node with its closest ancestor
-sum_dKs | -sum_dKa  [nodes]  compute the sum of the dKs or dKa of all children of this node, parsing until leaves nodes
-dCountS | -dCountA  [nodes]  to show simply the count of changes in respect to the direct ancestor node
-KaKs_with  | ..  [ref:nodes] compute KaKs for target nodes in respect to reference (ancestor). -Ks_with and -Ka_with are also available
-CountS_with | .. [ref:nodes] count changes for target nodes in respect to reference. -CountA_with is also available 
-wsize       affects all sequence analyses; compute values over a sliding window of this size (in codons)
-wstep       define step of sliding windows, if wsize is active

## node graphics       
the option -G is determining what graphical elements are shown in the nodes. If you don't specify any -G options, default elements are added for each analysis run.
Alternatively, you manually add elements, and possibly specify options for each one, in the form -G:analysis:options. Multiple elements are to be added using multiple such elements on the command line.
The options depend on the type of analysis. You can inspect the available options running --help analysis_name. Options are defined in the form   name=value,name2=value2 
Example: 
-G:w_change                       add a w_change graphical element with default options
-G:w_change:circle=2,text_r=1     same, specifying options to the graphical element
If you want to use default elements for most, but not all analysis, use option -G alone to populate with default values the graphical options, and then -G:* options to customize. 

## other graphical options
-scale   arg       source of distances among nodes in the displayed tree. Default is dKs
-title  "arg"      use this as the title displayed (instead of a generic title with the name of this codon alignment)
-ref_style [file]  define style for ref nodes (highlighting of syn/non-syn if -aa or -seq, size, color etc.).  Example of such file (tab separated, with header):
#index size  shape   opacity fgcolor bgcolor ns_pen  s_pen   stp_pen ns_bkg  s_bkg   stp_bkg unk_pen unk_bkg
1      10    square  1.0     #222222 #BABABA #000000 #000000 #000000 #CC2222 #2222CC #CCCC44 #000000 #FFFFFF
-style  [file]     define style for specific leaf nodes. Any number of ete2.nodestyle attributes can be specified. To define attributes of faces of the nodes, use prefix "faces.". Example of such file (tab separated):
node_name  size=8  bgcolor=#DD0000   faces:all.margin_bottom=6
node_name  fgcolor=#DD0000           faces:branch-right.margin_bottom=6

## sequence alignment
-seq         display aligned codon sequences next to leaf nodes, coloring syn and non-syn changes when comparing to -ref nodes
-aa          show aminoacid sequences instead of codons
-ref [nodes] define one or more reference node used for comparison for coloring changes. Only nodes below a ref node are compared and colored. More than a node can be defined: in this case, the colors displayed refers to the closest reference node, allowing a visual time-mapping of the changes. Can use -ref_style (see below) for customization
-small       display small alignment, letters are not readable but you can still have an overview with the colors
-zoom [arg]  the alignment is displayed small, and some region(s) are zoomed. The argument must be in the form    start-end   or   start1-end1,start2-end2 etc.  Start and end are relative to the alignment, both 1-based and included
-fsize [arg] font size for the alignment displayed
-no_ruler    do not display the usual position ruler as an aligned header 

## additional layers of information
-ann    arg  annotation file for displaying domains. Format (tab separated): "ALI start end text row color_line color_bkg color_text" ; all missing values are taken as defaults. You can replace ALI as the first word with any sequence name, if the coordinates are relative to that, and not to the alignment
-plot   arg  plot file, for displaying numerical values associated to each column in form of a plot below the alignment. The format is the following for each line: title TAB python_function [TAB row TAB color].  The python function can use "A" (codon_alignment) and must return a list of scores, one per position of the alignment. Row and plot color can be defined with the optional arguments.

## ancestral sequences   (normally predicted by codeml using null model)
-sankoff           for large alignments, if you just want the visualization: do not compute the null model, derive the ancestral states using sankoff algorithm
-ali_out  [file]   write an fasta alignment file including all nodes and also the ancestral sequences

## misc
-D                 do not open ete enviroment, just compute and print stuff in log
-out               do not open interactive ete enviroment, but save image in the argument file (.pdf or .png extension accepted). -width and -height can be used to control the size in pixels
-P                 open interactive prompt after showing tree
-j    [arg]        do not execute codeml. Instead, write all commands that would be executed to this file, or to prefix.job if no argument is provided. Useful for parallelizing pycodeml
-no_save           do not save object in cdn_ali file at the end of computation. Saving the object will save derived data that will not be recomputed
-no_log            do not save automatically the output of this program as prefix.last_log, where prefix is option -o
-log   arg         write log to this file instead into prefix.last_log
-print_opt         print currently active options
-h OR --help [+]   print this help and exit. With an argument, will display a guide to graphical options for a particular analysis (e.g. w_change) 
"""

command_line_synonyms={'free':'free_w'}
def_opt= {'a':None, 't':None, 'o':None,
'all':0,
'wsize':None, 'wstep':1,
'D':0, 'p':0, 'out':0,  'width':None, 'height':None, 'G':0, 
'ref':'root', 'seq':0, 'aa':0,
'no_save':0, 'sec':0, 'small':0, 'zoom':None, 'clean':0, 
'no_log':0, 'log':0, 
'alpha':0.05,
'sankoff':0, 'scale':'dKs', 'title':None,
'j':0, 'ref_style':0, 'ali_out':0,
'style':0, 
'fsize':10, 'ann':0, 'no_ruler':0, 'plot':0,
'prova':0}

#structure of following hash:   option name (analysis) -> [  function to decorate each node, {options to pass to the function} ]   #a third possible element is a function that given a node, returns True if the decorating function cannot be applied to that node (it is skipped)
###
analysis2decorating_functions={\
'w_change': [EvolPhyloTree.decorate_with_w_change, {'descendants':True}, EvolPhyloTree.is_root] , 
'w_change2':[EvolPhyloTree.decorate_with_w_change, {'descendants':False}, EvolPhyloTree.is_root], 
'free_w':[EvolPhyloTree.decorate_with_free_w, {}, EvolPhyloTree.is_root], 
'KaKs':[ EvolPhyloTree.decorate_with_KaKs_lineage, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['KaKs'], no_leaves=True)   ], 
'Ks':[ EvolPhyloTree.decorate_with_Ks_lineage, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['Ks'], no_leaves=True)   ],  #only for targets of this analysis
'Ka':[ EvolPhyloTree.decorate_with_Ka_lineage, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['Ka'], no_leaves=True)   ], 
'dKaKs':[ EvolPhyloTree.decorate_with_dKaKs, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['dKaKs'], no_root=True)   ], 
'dKs':[ EvolPhyloTree.decorate_with_dKs, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['dKs'], no_root=True)   ],  #only for targets of this analysis
'dKa':[ EvolPhyloTree.decorate_with_dKa, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['dKa'], no_root=True)   ], 
'sum_dKs':[ EvolPhyloTree.decorate_with_sum_dKs, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['sum_dKs'], no_leaves=True)   ],  #only for targets of this analysis
'sum_dKa':[ EvolPhyloTree.decorate_with_sum_dKa, {},  lambda x:not x in A.get_nodes( get_MMlib_var('opt')['sum_dKa'], no_leaves=True)   ], 
'dCountS':[  EvolPhyloTree.decorate_with_dCountS, {}, lambda x:not x in A.get_nodes( get_MMlib_var('opt')['dCountS'], no_root=True)   ],
'dCountA':[  EvolPhyloTree.decorate_with_dCountA, {}, lambda x:not x in A.get_nodes( get_MMlib_var('opt')['dCountA'], no_root=True)   ],
'count_sites':[  lambda x:x , {}, lambda x:False                 ], # [   EvolPhyloTree.decorate_with_count_sites, {}, lambda x:False  ],
'CountA_with':[  lambda x:x , {}, lambda x:False                 ],
'CountS_with':[  lambda x:x , {}, lambda x:False                 ],
'KaKs_with': [ lambda x:x, {},    lambda x:False                 ],
'Ks_with': [ lambda x:x, {},      lambda x:False                 ],
'Ka_with': [ lambda x:x, {},      lambda x:False                 ],
#'positive1': [ lambda x:x, {},      lambda x:False                 ],
#'positive2': [ lambda x:x, {},      lambda x:False                 ],
### TODO: add function to decorate with sites_face
}
for analysis in analysis2decorating_functions: def_opt[analysis]=False

advanced_help={}  
for analysis in analysis2decorating_functions: 
  f=analysis2decorating_functions[analysis][0]
  o=analysis2decorating_functions[analysis][1]
  h= ' ##### Advanced help for: '+analysis+'\nRuns function: '+str(f.__name__) +int(bool(o))*(' with these options (don\'t touch these): '+str(o))+'\n'
  h+='Arguments and defaults: \n\n'
  h+=function_help(f)
  advanced_help[analysis]=h
     
##### styles for sequence coloring when comparing to refnodes, in the order they are defined on command line. overriden by -ref_style
seq_ref_node_styles=[]
seq_ref_node_evolCodonColorSchemes=[]
## 1
seq_ref_node_style1=NodeStyle()
seq_ref_node_style1["size"] = 10;         seq_ref_node_style1["shape"] = "square"
seq_ref_node_style1["fgcolor"]="#222222"; seq_ref_node_style1["bgcolor"]="#BABABA"
seq_ref_node_evolCodonColorSchemes.append(evolCodonColorScheme("#000000","#000000","#000000", '#CC2222', '#2222CC', "#CCCC44", "#000000", "#858585", opacity=1.0))
seq_ref_node_styles.append(seq_ref_node_style1)
## 2
seq_ref_node_style2=NodeStyle()
seq_ref_node_style2["size"] = 10;         seq_ref_node_style2["shape"] = "square"
seq_ref_node_style2["fgcolor"]="#333333"; seq_ref_node_style2["bgcolor"]="#CFCFCF"
seq_ref_node_evolCodonColorSchemes.append(evolCodonColorScheme("#000000","#000000","#000000", '#CC2222', '#2222CC', "#CCCC44", "#000000", "#858585",opacity=0.5))
seq_ref_node_styles.append(seq_ref_node_style2)
##3
seq_ref_node_style3=NodeStyle()
seq_ref_node_style3["size"] = 10;         seq_ref_node_style3["shape"] = "square"
seq_ref_node_style3["fgcolor"]="#444444"; seq_ref_node_style3["bgcolor"]="#E5E5E5"
seq_ref_node_evolCodonColorSchemes.append(evolCodonColorScheme("#000000","#000000","#000000", '#CC2222', '#2222CC', "#CCCC44", "#000000", "#858585", opacity=0.2))
seq_ref_node_styles.append(seq_ref_node_style3)
custom_node_styles={'faces':{}}

p_values_colors={True:'magenta', False:''} # for printing to screen colored values for significant p_values
class notracebackException(Exception):"""when this exception is raised, programs dies with a single message -- less noisy than usual """
class input_error(notracebackException):  """ """

#########################################################
###### start main program function
def main(args={}):
#########################################################
############ loading options

  colored_keywords={'WARNING':'red', 'ERROR': 'red,underscore'};   set_MMlib_var('colored_keywords', colored_keywords)
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms, tolerated_regexp=['G:.+'], strict=True, advanced=advanced_help )
  else:  opt=args
  set_MMlib_var('opt', opt)

  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 

  global A     #### main object
  
  ######################
  #### reading codon_alignment object and options, or die complaining 
  if opt['o']:
    if opt['o'].endswith('.cdn_ali'): opt['o']=opt['o'][:-8]
  if not opt['a'] and not opt['t'] and (not opt['o'] or not is_file(opt['o']+'.cdn_ali')): raise input_error, "ERROR not input files! you must provide either an alignment (-a) and a tree (-t), and/or a previously computed cdn_ali file (-o). Run with -h for help"

  if not opt['a'] and not opt['t']:    #only option -o provided
    service('loading codon_alignment ...  ')    
    A= codon_alignment(path=opt['o'], clean=opt['clean'])
#  elif not opt['t'] or not opt['a']:   A= pass #raise IOError, "ERROR you must provide an an alignment (-a) and a tree (-t)! unless this is not your first run, in which case you can also provide also a  cdn_ali (-o) "
  else: 
     service('loading/making codon_alignment ...  ')
     A=codon_alignment(ali_file=opt['a'], tree_file=opt['t'], path=opt['o'], clean=opt['clean'])
  # options dependencies and default values of activated options
  if opt['aa']:   opt['seq']=1
  if opt['small'] and not opt['zoom']: opt['zoom']=False
  if opt['all']: 
    for analysis in analysis2decorating_functions.keys():      opt[analysis]=1
  if type(opt['w_change'])== int: opt['w_change']='all'
  if type(opt['w_change2'])== int: opt['w_change2']='all'
  if type(opt['KaKs'])== int: opt['KaKs']='all'
  if type(opt['Ka'])== int: opt['Ka']='all'
  if type(opt['Ks'])== int: opt['Ks']='all'
  if type(opt['dKaKs'])== int: opt['dKaKs']='all'
  if type(opt['dKa'])== int: opt['dKa']='all'
  if type(opt['dKs'])== int: opt['dKs']='all'
  if type(opt['sum_dKs'])== int: opt['sum_dKs']='all'
  if type(opt['sum_dKa'])== int: opt['sum_dKa']='all'
  if type(opt['dCountS'])== int: opt['dCountS']='all'
  if type(opt['dCountA'])== int: opt['dCountA']='all'
  #
  fsize_alignment= opt['fsize']
  
  if not opt['no_log']:
    log_filename=A.path+'.last_log'
    if opt['log']: log_filename= opt['log']
    log_file=open(log_filename, 'w')
    set_MMlib_var('log_file', log_file)

  # some noise in log
  write('\n'+ center_str(' pycodeml version '+str(__version__)+' | '+'Date:'+bbash('date')+' ', 80, '-_') +'\n', 1  )
  write('Command line: ')
  if opt['a']: write('-a '+opt['a']+' ')
  if opt['t']: write('-t '+opt['t']+' ')
  if opt['o']: write('-o '+opt['o']+' ')
  for k in def_opt: 
    if not  k in 'ato' and  opt[k]!=def_opt[k]: write('-'+k+' '+str(opt[k])+ ' ')
  write('', 1)  
  if not opt['no_log']:    write('Log file: '+A.path+'.last_log', 1)
  write("\n"+str(A), 1)

  #computing tree if missing
  if not A.tree: 
    write('Tree missing! computing it ...')
    A.build_tree()
    write('done', 1)

  write('\nTree topology: \n'+str(A.tree)+'\n', 1)
  write('-'*80+'\n', 1)

  # leaves style
  if opt['style']:
    global custom_node_styles;         custom_node_styles['faces']={}   ## keeping indexed those nodes for which we want somethign modified in all their faces
    write('Option -style:  reading file '+opt['style'], 1)    
    for line in open(opt['style']):
      if not line.strip() or line.startswith('#'): continue
      splt=line.strip().split('\t');       node_name=splt[0]
      found= A.get_nodes(node_name)
      if len(found)!=1: raise Exception, "-style ERROR node name:  "+node_name+" : "+str(len(found))+' nodes match!'
      st= NodeStyle(); st['size']=0
      for text in splt[1:]:
        attr_name, value= text.split('=')
        if attr_name.startswith('faces:'):
          if not custom_node_styles['faces'].has_key( found[0].get_name() ): custom_node_styles['faces'][found[0].get_name()]={}
          position, attr_name  =attr_name.split(':')[1].split('.')   # e.g. branch-right, or "all"
          if not custom_node_styles['faces'][ found[0].get_name() ].has_key(position): 
            custom_node_styles['faces'][ found[0].get_name() ][ position ] ={}           
          custom_node_styles['faces'][ found[0].get_name() ][ position ][attr_name]=  option_value(value)  
        else:        st[attr_name]=option_value(value) ## converting to appropriate data type

      custom_node_styles[found[0].get_name()]=st   #using found to avoid problems with node name synonyms
      
  #ref style
  if opt['ref_style']:
    write('Option -ref_style:  reading file '+opt['ref_style'], 1)
    global seq_ref_node_evolCodonColorSchemes; seq_ref_node_evolCodonColorSchemes=[]; 
    global seq_ref_node_styles; seq_ref_node_styles=[]
    for line in open( opt['ref_style'] ):
      if not line or line.startswith('#'): continue
      t=line.strip().split('\t')
      seq_ref_node_style1=NodeStyle();      seq_ref_node_style1['size']=int(t[1]);        seq_ref_node_style1['shape']=t[2];  
      seq_ref_node_style1['fgcolor']=t[4];     seq_ref_node_style1['bgcolor']=t[5]
      codon_color_scheme=   evolCodonColorScheme(   t[6], t[7], t[8], t[9], t[10], t[11], t[12], t[13],  opacity= float(t[3]) )
      seq_ref_node_styles.append(seq_ref_node_style1)
      seq_ref_node_evolCodonColorSchemes.append(codon_color_scheme)

  ## list of ref nodes 
  seq_ref_nodes= A.get_nodes(opt['ref'])
  if opt['seq'] and len(seq_ref_nodes)> len(seq_ref_node_styles): raise input_error, "ERROR too many -ref nodes specified! the program implements colors only for "+str(len(seq_ref_node_styles))
  
  ### reading customized graphical options
  graphical_options={}
  any_graphical_options_specified= bool([1 for k in opt if k.startswith('G:')])
  if opt['G'] or not any_graphical_options_specified:
    for analysis in analysis2decorating_functions: 
      if opt[analysis]:
        graphical_options[analysis]=   analysis2decorating_functions[analysis][1] #setting default grpahical options
  for k in opt:
    if k.startswith('G:'):
      try:
        analysis=k.split(':')[1]
        graphical_options[analysis]=   analysis2decorating_functions[analysis][1] #setting default graphical options
        if k.count(':')>1: 
          for expr in k.split(':')[2].split(','):
            option_name=   expr.split('=')[0]
            value=         option_value( expr.split('=')[1] )
            graphical_options[analysis][option_name]=value  #setting customized options
      except: raise input_error, "ERROR syntax error in the graphical option: -"+str(k)+'    See -help'
      if not opt[analysis]: raise input_error, "ERROR you specified a graphical option for an analysis which is not active! Option: "+k+" Analysis: "+analysis
      

  if opt['j']:   ### print job mode, not running
    opt['D']=1; opt['no_save']=1
    filename=    A.path+'.job'
    if opt['j']!=1:   filename= opt['j']
    write('#####   No run mode ( -j )  Instead, all commands will be printed to: '+filename+' ####', 1)
    fileh=open(filename, 'w')
    set_ete2_MM_var('no_run_mode', fileh)
    null_model=A.add_codeml_model(null_codeml_model)

  ############################
  #### ancestral sequences, null model computation if necessary
  else: 
    if opt['sankoff']:
        write('Setting ancestral sequences with the sankoff algorithm', 1)
        A.set_ancestral_sequences('sankoff')
    else:
      if not A.has_model('null'):
        write('Running/loading  codeml model: null', 1)
        null_model=A.add_codeml_model(null_codeml_model)
      else: 
        write('Found codeml model: null', 1)      
        null_model= A.get_model('null')

      write('Setting ancestral sequences from codeml null model', 1)
      A.set_ancestral_sequences('null')
      write('\nW in null model: '+str(round(A.null_model_omega(), 4)), 1)

  if opt['ali_out']:
    write('Option -ali_out: writing alignment of sequences including for ancestral nodes --> '+str(opt['ali_out']), 1)
#    for n in A.ancestral_sequences:
#      print ">"+n
#      print A.ancestral_sequences[n]
    ali=alignment()
    for node in A.tree.traverse(strategy="preorder"):
      if node.is_leaf():        
        ali.add( A.ali.fill_title( node.get_name() ),  A.seq_of(node.get_name()) )
      else:
        if node.is_root(): continue
        ali.add( node.get_name()+" ancestral "+A.ancestral_sequences_source ,  A.ancestral_sequences[ node ] )        
    ali.display(  opt['ali_out']  )

  write('\n'+center_str(' Starting analysis ', 80, '-_'), 1)      
  
#########################
  ### branch analysis
  if opt['w_change']:             ############ 
    write('\n'+center_str(' w_change test (including descendants) ', 80, '_'), 1   )
    targets= A.get_nodes( opt['w_change'], no_root=True)
    for target_node in targets:
      name='w_change_D.'+ target_node.get_name()
      write('Test: '+ljust(name, 80)+' ')
      if not A.has_test(name):        test= A.add_codeml_test(w_change_test, target=target_node, descendants=True) ## running or loading here!        
      else:                           test= A.get_test(name)
      model= A.get_model(name)
      results=test.results()
      if not opt['j']:
        write('W: '+str(round(results[1][1], 4))+' (rest: '+str(round( results[1][0], 4))+') ')
        write(' p-value: '+str(round(test.p_value, 4)), how=p_values_colors[is_significant(test.p_value)])
        if len(targets)>1: 
          test.bonferroni= test.p_value*len(targets)
          write('   Bonferroni: '+str(round(test.bonferroni, 4)), how=p_values_colors[is_significant(test.bonferroni)])
      write('', 1)

  if opt['w_change2']:             ############ 
    write('\n'+center_str(' w_change test (single node) ', 80, '_'), 1   )
    targets= A.get_nodes( opt['w_change2'], no_root=True)
    for target_node in targets:
      name='w_change_S.'+ target_node.get_name()
      write('Test: '+ljust(name, 80)+' ')
      if not A.has_test(name):         test= A.add_codeml_test(w_change_test, target=target_node, descendants=False) ## running or loading here!        
      else:                            test= A.get_test(name)
      model= A.get_model(name)
      results=test.results()
      if not opt['j']:
        write('W: '+str(round(results[1][1], 4))+' (rest: '+str(round( results[1][0], 4))+') ')
        write(' p-value: '+str(round(test.p_value, 4)), how=p_values_colors[is_significant(test.p_value)])
        if len(targets)>1: 
          test.bonferroni= test.p_value*len(targets)
          write('   Bonferroni: '+str(round(test.bonferroni, 4)), how=p_values_colors[is_significant(test.bonferroni)])
      write('', 1)

  if opt['free_w']:             ############ 
    write('\n'+center_str(' free_w ', 80, '_'), 1   )
    name='free_w'
    write('Test: '+ljust(name, 80)+' ')
    if not A.has_test(name):
      test= A.add_codeml_test(free_w_test) ## running or loading here!        
    else: 
      test= A.get_test(name)
    #model= A.get_model(name)
    test.run()
    if not opt['j']:
      write('p-value: '+str(round(test.p_value, 4)), how=p_values_colors[is_significant(test.p_value)])
    write('', 1)
    

#########################    ### NOT  FINISHED!!! ###
  ### site analysis
  if opt['positive1']:             ############ 
    write('\n'+center_str(' positive1 : test#1 for positive selection (M2a vs M1a) ', 80, '_'), 1   )
    name='positive_selection1_sites'
    write('Test: '+ljust(name, 80)+' ')
    if not A.has_test(name):
      test= A.add_codeml_test(positive_selection1_sites_test) ## running or loading here!        
    else: 
      test= A.get_test(name)
    #model= A.get_model(name)
    test.run()
    write('p-value: '+str(round(test.p_value, 4)), 1, how=p_values_colors[is_significant(test.p_value)])
    if is_significant(test.p_value):
      m2a=A.get_model('positive_selection')      
      n_sites_selected_fifty=len(  m2a.sites  )
      write( '\n'+ str(n_sites_selected_fifty) +" sites have Pr(w>1)  > 0.5   Find here below those for which Pr(w>1) is considered significant. The mean w and its deviation are shown ", 1)
      printed_any=False
      for position in sorted( m2a.sites.keys() ):
        site= m2a.sites[position]
        if  1: #is_significant(1.0-site.probability):   DEBUG!!
          write('Codon pos: '+str(position).ljust(4)+'  Pr(w>1): '+str(  round(site.probability, 4)  ).ljust(6)+'  mean w: '+str(  round(site.post_mean, 4)  ).ljust(6)+' +- '+str(  round(site.deviation, 4)  ).ljust(6), 1   )
          printed_any=True
      if not printed_any:  write('(None)', 1)    

  if opt['positive2']:             ############ 
    write('\n'+center_str(' positive2 : test#2 for positive selection (M8 vs M7) ', 80, '_'), 1   )
    name='positive_selection2_sites'
    write('Test: '+ljust(name, 80)+' ')
    if not A.has_test(name):
      test= A.add_codeml_test(positive_selection2_sites_test) ## running or loading here!        
    else: 
      test= A.get_test(name)
    #model= A.get_model(name)
    test.run()
    write('p-value: '+str(round(test.p_value, 4)), 1, how=p_values_colors[is_significant(test.p_value)])
    if is_significant(test.p_value):
      m2a=A.get_model('beta_positive')      
      n_sites_selected_fifty=len(  m2a.sites  )
      write( '\n'+ str(n_sites_selected_fifty) +" sites have Pr(w>1)  > 0.5   Find here below those for which Pr(w>1) is considered significant. The mean w and its deviation are shown ", 1)
      printed_any=False
      for position in sorted( m2a.sites.keys() ):
        site= m2a.sites[position]
        if  is_significant(1.0-site.probability):  
          write('Codon pos: '+str(position).ljust(4)+'  Pr(w>1): '+str(  round(site.probability, 4)  ).ljust(6)+'  mean w: '+str(  round(site.post_mean, 4)  ).ljust(6)+' +- '+str(  round(site.deviation, 4)  ).ljust(6), 1   )
          printed_any=True
      if not printed_any:  write('(None)', 1) 

  ############################################
  ######### analysis without codeml models
  if opt['wsize']: 
    windows=[]
    codon_start=1; window_size=opt['wsize'];  window_step=opt['wstep']
    codon_end=codon_start+window_size-1
    while codon_end<A.length()/3:
      windows.append( [codon_start, codon_end]  )
      codon_start+=window_step
      codon_end=codon_start+window_size-1
  else: windows=[None]

  for window in windows:
    if window is None:    positions=None;     window_name=''
    else:                 positions=[window]; window_name='@{st}-{end}:'.format(st=window[0], end=window[1])

    ## analysis that can't be applied to root
    for analysis_name in [ 'dCountA', 'dCountS', 'dKa', 'dKs', 'dKaKs' ]:
      if opt[analysis_name]:
        write('\n'+(' '+analysis_name+' (pycodeml) ').center(80, '_'), 1   )
        targets= A.get_nodes( opt[analysis_name], no_root=True)            # no root
        for target_node in targets:      
          write( ljust(window_name+target_node.get_name(), 80)+' '+analysis_name+': ')
          value= eval('target_node.'+analysis_name+'(positions=positions)')   
          if analysis_name.startswith('dCount'): write(str(value), 1)
          else:                                  write(str(round(value, 4) ), 1)

    ## analysis that can't be applied to leaves
    for analysis_name in [ 'sum_dKs', 'sum_dKa', 'Ka', 'Ks', 'KaKs' ]:
      if opt[analysis_name]:
        write('\n'+(' '+analysis_name+' (pycodeml) ').center(80, '_'), 1   )
        targets= A.get_nodes( opt[analysis_name], no_leaves=True)           #no leaves
        analysis_name_call=analysis_name
        if not analysis_name.startswith('sum_'): analysis_name_call+='_lineage' ## note: adding "_lineage" to analysis name; KaKs is called KaKs_lineage in ete2_MM
        for target_node in targets:      
          write( ljust(window_name+target_node.get_name(), 80)+' '+analysis_name+' ')
          value= eval('target_node.'+analysis_name_call+'(positions=positions)')   
          write(str(round(value, 4) ), 1)

    ## analysis "with" one or more specific nodes
    for analysis_name in [ 'CountS_with', 'CountA_with', 'Ks_with', 'Ka_with', 'KaKs_with' ]:
      if opt[analysis_name]:
        write('\n'+(' '+analysis_name+' (pycodeml) ').center(80, '_'), 1   )
        list_of_refs_and_targets=     parse_args_analysis_with(  opt[analysis_name], A, analysis_name )  
        for ref_for_this, target_nodes in list_of_refs_and_targets:
          for target_node in target_nodes:
            write( ljust(window_name+target_node.get_name(), 80)+' '+analysis_name+' '+  ref_for_this.get_name()+' : ')
            value= eval('target_node.'+analysis_name+'(ref_for_this, positions=positions)')
            if analysis_name.startswith('Count'): write(str(value), 1)
            else:                                 write(str(round(value, 4) ), 1)        


  ####### scaling tree
  if opt['scale']:
      write('\nSetting branch lengths from: '+str(opt['scale']), 1)
      A.set_branch_length(opt['scale'], silent=True)

  ############################
  #### graphics
  if not opt['D']:
    write('\n'+center_str(' Starting ete2-PyQt4 graphics ', 80, "-_")+'\n', 1 )

    A.tree_style.legend_position=3  # bottom-left
    # adding title
    if not opt['title'] is None:  displayed_title= opt['title']
    else:                         displayed_title= 'Codon alignment: '+base_filename(A.path)
    title_face= faces.TextFace( displayed_title , fsize=16 )
    A.tree_style.title.add_face(title_face, column=1)

    # adding legend with null model kaks
    if A.has_model('null'):
      average_kaks_face=faces.TextFace( 'Null model KaKs: '+str(round( A.null_model_omega()  , 4)), fgcolor=get_color_free_kaks(A.null_model_omega()), fsize=12 )
      A.tree_style.legend.add_face( average_kaks_face,  column=4)    

    nostyle = NodeStyle()  ## default style for nodes
    nostyle["size"] = 0
    for node in A.tree.traverse():
      for analysis in graphical_options:
        decorating_function=        analysis2decorating_functions[analysis][0]
        options_for_this_function=  graphical_options[analysis]
        if len(analysis2decorating_functions[analysis])> 2:  ### checking if this analysis can't be applied to this node  (3 arg in analysis2decorating_functions[analysis] )
          if analysis2decorating_functions[analysis][2](node)          : continue
        try:         decorating_function(node, **options_for_this_function)     ### adding faces for this analysis!
        except:      
          printerr('ERROR with graphical function: '+decorating_function.__name__+' ; arguments passed were node '+node.get_name()+' and options: '+str(options_for_this_function), 1)
          raise

      if node.is_leaf():          
        name_face=faces.TextFace( node.name, fgcolor="#000000", fsize=12)
        name_face.margin_right=10
        node.add_face(name_face, column = node.column_index, position="branch-right")
        node.column_index+=1

      ############ now sequence face!
    if opt['seq']:   #part1
      zoomed_regions_list=None
      if not opt['zoom'] is None and not opt['zoom']: #like with option -small
        zoomed_regions_list=False
      elif opt['zoom']: 
        try:
          zoomed_regions_list=[]
          for region in opt['zoom'].split(','):
            st, end = region.split('-')
            zoomed_regions_list.append( [int(st), int(end)] )
        except: raise notracebackException, "ERROR format for option -zoom was not valid! see -help"

      max_column_index=4
      #counting column index for visualizing aligned faces of sequences in case
      for node in A.tree.traverse():
        if node.is_leaf():          
          if node.column_index>max_column_index: max_column_index=node.column_index        
      ## now attaching seq face
      for node in A.tree.traverse():
        if node.is_leaf():
          node.column_index=max_column_index        
          codon_face=evolCodonsFace(seq=node.sequence(), translate= bool(opt['aa']), zoom= zoomed_regions_list, fsize=fsize_alignment  )
          node.codon_face=codon_face # linking to get it later and color 
          node.add_face(codon_face, column = node.column_index, position="aligned") 
          node.column_index+=1

      ## position ruler
      if not opt['no_ruler']:
        pos_ruler_face=     alignmentRulerFace(  guide=codon_face )     #using last codon face as guide
        A.tree_style.aligned_header.add_face(pos_ruler_face, column=max_column_index )
          
    for node in A.tree.traverse():
      if not custom_node_styles.has_key( node.get_name() ):      node.set_style(nostyle)      
      else:                                      node.set_style(  custom_node_styles[node.get_name()] )
      if custom_node_styles['faces'].has_key( node.get_name() ):
        for position in custom_node_styles['faces'][node.get_name()]:
          if position == 'all': iterate_pos=['branch-right', 'branch-top', 'branch-bottom', 'float', 'float-behind', 'aligned']
          else:                 iterate_pos=[position]
          for pos in iterate_pos:
            fdict = getattr(node.faces, pos)
            for col_index in fdict.keys():
              for attr_name in custom_node_styles['faces'][node.get_name()][position].keys():
                a_faces=fdict[col_index]
                for a_face in a_faces:                 setattr(a_face,  attr_name,  custom_node_styles['faces'][node.get_name()][position][attr_name] )
            #print pos, fdict            
            #raw_input('...')

    if opt['seq']:  #part2
      #now coloring codons according to the changes with the seq_ref_node, but only if the node is under it
      for seq_ref_node_index, seq_ref_node in enumerate(seq_ref_nodes):
        seq_ref_node.set_style(seq_ref_node_styles[seq_ref_node_index])
        for node in seq_ref_node.traverse():
          if node.is_leaf():
            node.codon_face.set_color_codon_changes(seq_ref_node.sequence(), color_scheme= seq_ref_node_evolCodonColorSchemes[seq_ref_node_index] )

      example_codon_face= node.codon_face 


      if opt['ann']:
        write('# Option -ann: reading domain annotation from file: '+str(opt['ann']), 1)
        ann_per_row={}
        for line in open(opt['ann']):
          if not line.strip() or line.startswith('#') : continue
          splt=line.strip().split('\t')
          target= splt[0];     start=int(splt[1]);     end=int(splt[2])
          if target!='ALI':
            start= A.ali.position_in_ali(target, (start-1)*3+1  ) ##nucleotide based, 1 based
            if start is None: start=1
            start= (start+2)/3
            end=   A.ali.position_in_ali(target, (end-1)*3+1)            
            end=   (end+2)/3
          text= splt[3]      if len(splt)>3 else ''
          row = int(splt[4]) if len(splt)>4 else 1
          col_line= splt[5]      if len(splt)>5 else "#000000"
          col_bkg=  splt[6]      if len(splt)>6 else "#FFFFFF"
          col_text= splt[7]      if len(splt)>7 else "#000000"
          annotation = {}
          annotation['start']=start;       annotation['end']=end
          annotation['text']=text;         annotation['col_line']=col_line;
          annotation['col_bkg']=col_bkg;   annotation['col_text']=col_text
          if not ann_per_row.has_key(row): ann_per_row[row]=[]
          ann_per_row[row].append( annotation )
        write('  a total of '+str(  sum(  [ len(ann_per_row[row]) for row in ann_per_row  ] ) ) +' annotations were loaded.', 1  )
        for row in ann_per_row:
          annotation_face = alignmentAnnotationFace(example_codon_face, ann_per_row[row], fsize=fsize_alignment )
          A.tree_style.aligned_foot.add_face( annotation_face,  column=max_column_index)

      if opt['plot']:        
        plot_per_row={}
        write('# Option -plot: reading plot specifics from file: '+str(opt['plot']), 1)
        for line in open(opt['plot']):
          if not line.strip() or line.startswith('#') : continue          
          try: 
            splt=line.strip().split('\t')
            plot={}
            plot['title']=splt[0];  plot['y_values']= eval( splt[1] )
            if len(plot['y_values']) != A.ali.length()/3: 
              raise Exception, "option -plot ERROR with line: "+line.strip()+' : the function must return exactly '+str(A.ali.length()/3) +' values, one per codon alignment position! Instead there are '+str(len(plot['y_values']))+' values'
            row= 1 if len(splt)<=2 else int(splt[2])
            plot['color']= "#000000" if len(splt)<=3 else splt[3]            
            if not plot_per_row.has_key(row): plot_per_row[row]=[]
            plot_per_row[row].append(    plot      )
          except: 
            printerr("option -plot syntax ERROR reading file "+str(opt['plot'])+ " in line: "+line.strip(), 1)
            raise
        write('  a total of '+str(  sum(  [ len(plot_per_row[row]) for row in plot_per_row  ] ) ) +' plots were loaded.', 1  )
        empty_face=faces.Face()
        for row in plot_per_row:
          plot_face =    coloredPlot(  example_codon_face,   plot_per_row[row],  fsize=10, font='Courier' )
          A.tree_style.aligned_foot.add_face( plot_face,  column=max_column_index)
#          title_faces=plot_face.title_faces()
#          for title_face in title_faces:
#            A.tree_style.aligned_foot.add_face( title_face, column=max_column_index)
            

      if opt['prova']:
        p_model=A.get_model('positive_selection')
        sites_list= []
        for p0 in range(A.length()/3):
          if p_model.sites.has_key( p0+1 ):  sites_list.append(   p_model.sites[p0+1]  )
          else:                              sites_list.append(None)
        
        site_face=  selectedSitesFace(sites_list, codon_face= example_codon_face)
        A.tree_style.aligned_foot.add_face( site_face,  column=max_column_index)

    if opt['out']:
      #### rendering graphical output file, png or pdf format are accepeted. using ete2 method
      write('Writing graphical output to file: '+opt['out'], 1)
      A.tree.render(opt['out'], h=opt['height'], w=opt['width'])
    else:
      ####### opening interactive ete
      A.tree.show()

  write('', 1)

  if opt['p']:
    ## opening interactive ipython
    x=A
    interactive_mode(message='The codon_alignment object is loaded into variable: x')()

  if not opt['no_save']:
    A.save()   #it prints stuff within the function
  
  write('\n'+center_str(' pycodeml completed | Date: '+bbash('date')+ ' ', 80, '-_'), 1)

###############
def basic_layout(node):  pass
def parse_args_analysis_with(argument, A, analysis_name):
  """ used to get the arguments for functions such as KaKs_with, countA_with etc"""
  out= []      #   [   [ref1, targets1],  [ref2, targets2] , ... ]
  ref_for_this=None
  if type(argument)==int: raise input_error, "ERROR analysis "+analysis_name+" with wrong syntax! Its argument must start with a reference node followed by a single ':'"
  for comma_sep_piece in argument.split(','):
    if comma_sep_piece.count(':')>1: raise input_error, "ERROR analysis "+analysis_name+" with wrong syntax! Its argument must start with a reference node followed by a single ':'"
    if ':' in comma_sep_piece:          
      ref_for_this, comma_sep_piece = comma_sep_piece.split(':')
      ref_for_this = A.get_nodes(ref_for_this)[0]
    target_nodes = A.get_nodes(comma_sep_piece)
    if ref_for_this is None: raise input_error, "ERROR analysis "+analysis_name+" activated with wrong syntax! Its argument must start with a reference node followed by ':'" 
    out.append(   [ref_for_this, target_nodes]   )
  return out

#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:    pass

  if sys.exc_info()[0]: #an exception was raised. let's print the information using printerr, which puts it in the logfile as well, if any.   
    if issubclass(sys.exc_info()[0], notracebackException):      printerr( sys.exc_info()[1], 1)
    elif issubclass(sys.exc_info()[0], SystemExit):      pass
    else:                                                                  printerr('ERROR '+ traceback.format_exc( sys.exc_info()[2]) , 1)


  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()

