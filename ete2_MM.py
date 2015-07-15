#! /usr/bin/python2.6 -u
from string import *
import sys
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
#sys.path.insert(0, "/software/rg/el6.3/ete2-2.2rev1056/build/lib/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
from ete2 import Tree, PhyloTree, TreeStyle, faces, NodeStyle, treeview
from PyQt4 import QtCore, QtGui
from scipy.stats import chi2
no_run_mode=False

def set_ete2_MM_var(varname, value):  globals()[varname]=value
def get_ete2_MM_var(varname):         return globals()[varname]

def is_significant(pvalue): return pvalue<0.05

n_digits=3

class colorScheme(object):
  """a class to handle color schemes. it is basically just an object with k -> values like an hash """ 
  def __repr__(self):
    o=''
    for k in dir(self):
      if  not k.startswith('__'): o+= str(k)+': '+str(getattr(self, k))+'\n'
    return o

class evolCodonColorScheme(colorScheme):
  """ color scheme for codon mutations: syn or non-syn"""
  def __init__(self, non_syn_pen="#000000", syn_pen="#000000", nonsense_pen="#000000", non_syn_bkg="#FFFFFF", syn_bkg="#FFFFFF", nonsense_bkg="#000000",  unknown_pen="#000000", unknown_bkg="#00FF2F", opacity=1.0):
    self.non_syn_pen=non_syn_pen
    self.syn_pen=syn_pen
    self.non_syn_bkg=non_syn_bkg
    self.syn_bkg=syn_bkg
    self.unknown_pen=unknown_pen
    self.unknown_bkg=unknown_bkg
    self.nonsense_pen=nonsense_pen
    self.nonsense_bkg=nonsense_bkg
    self.opacity=opacity

basicSynNonsynColorScheme=evolCodonColorScheme( "#000000","#000000","#000000", '#DD5555', '#5555DD', "#AAAA44"  )

class kaksColorScheme(colorScheme):
  """ color scheme for  kaks values (purifying, neutral, positive selection)"""
  def __init__(self, purifying, neutral, positive):
    self.purifying=purifying
    self.neutral=neutral
    self.positive=positive

basicKaKsColorScheme=kaksColorScheme( "#0000AA", "#AA0000", "#00AA00"   )
kaks_neutral_boundaries=[0.8, 1.2] ## boundaries INCLUDED!


def kaks_to_color(kaks):
  """ returns a color for any kaks. puryfing is displayed in blue, neutral in red, positive in green. 
  intermediate colors are chosen to have nice gradients"""
  #saturation points
  
  try:
    
    if kaks<0.05:  kaks=0.05
    elif kaks>5.0: kaks=5.0

    if kaks   <= 0.4:              #min  #max   #min    #col min    #col max                    #min is also saturation point, for extreme values 
      color=  color_scale(   (kaks-0.05)/(0.4 - 0.05), '#00006E',  '#0000FF'        )  #dark blue,   full blue
    elif kaks <= 0.75:
      color=  color_scale(   (kaks-0.4)/(0.75 - 0.4),  '#0000FF' ,  '#AAAAFF'       )  #full blue, light blue,   
    elif kaks <= 1.0:
      color=  color_scale(   (kaks-0.75)/(1.0 - 0.75), '#AAAAFF' ,  '#FF0000'       )  #light blue,   full red
    elif kaks <= 2.0: 
      color=  color_scale(   (kaks-1.0)/(2.0 - 1.0),   '#FF0000' ,  '#00BF00'      )   #full red,     green
    else:  # <=5.0
      color=  color_scale(   (kaks-3.0)/(5.0 - 2.0),   '#00BF00' ,  '#00FF00'     )   #green,        full green

    return color
  except: 
    print "\n\n*****", kaks 
    raise


def add_kaks_color_legend(tree_style, position='branch-bottom'):
  pass 
#  for c, value in enumerate(  
#  tree_style.legend.add_face( column=c, position=position )
  

  ### UNFINISHED
class alignmentRulerFace(faces.Face):
  """ """
  def __init__(self, guide): # length, fsize=10 ):
    faces.Face.__init__(self)
    self.guide= guide #fsize=fsize
  def update_pixmap(self):
    # code pasted from SequenceFace  
    font = QtGui.QFont("Courier", self.guide.fsize)
    width = self.guide.pixmap.width() ;     height= self.guide.pixmap.height()
    self.pixmap = QtGui.QPixmap(width,height)
    self.pixmap.fill()
    p = QtGui.QPainter(self.pixmap)
    p.setFont(font)
    letter_brush = QtGui.QBrush(QtGui.QColor(  "#000000"    ))   #bkground
    letter_pen = QtGui.QPen(QtGui.QColor( "#000000"  ))
    p.setPen(letter_pen)
    if not self.guide.zoomed_regions:      mark_list=range( 0,  ( len( self.guide.seq) -2 ), 20  )
    else: 
      mark_list=[]
      for i in range( 0,  ( len( self.guide.seq) -2 ), 10  ):
        if        i%20==0 or  any(   [  i>=zoom_st and i<=zoom_end    for zoom_st, zoom_end in self.guide.zoomed_regions] ):
          mark_list.append( i )    ### keeping every -20 tick; plus, also the -10 but only if it falls in a zoomed region
    mark_list=mark_list[1:]
    for mark in mark_list:
      x_pos1, x_pos2=self.guide.x_of_position(mark)
      y_pos= height/2
      p.drawText(x_pos1, y_pos, str(mark) )
      p.drawLine(x_pos1,  0, x_pos1, height)
    p.end()

class alignmentAnnotationFace(faces.Face):
  """ """
  def __init__(self,  guide, annotations, fsize=10, font='Courier'):
    """ annotations is a list of objects with these attributes:  start, end, text, col_line, col_bkg, col_text """ 
    faces.Face.__init__(self)
    self.guide= guide
    self.annotations=annotations
    self.fsize=fsize; self.font=font
    return 
  def update_pixmap(self):
    width = self.guide.pixmap.width() ;     height= self.guide.pixmap.height()
    self.pixmap = QtGui.QPixmap(width,height)
    self.pixmap.fill()
    p = QtGui.QPainter(self.pixmap)
    for annotation in self.annotations:
      brush = QtGui.QBrush( QtGui.QColor(  annotation['col_bkg'] ))   #bkground
      #drawing rectangle
      x1, z= self.guide.x_of_position(  annotation['start']  ) 
      z, x2= self.guide.x_of_position(  annotation['end']    )   #throwing away z
      p.fillRect(x1,  0,   x2-x1, height, brush)
      pen=    QtGui.QPen( QtGui.QColor(  annotation['col_line'] ) , 1, QtCore.Qt.SolidLine )
      p.setPen(pen)
      p.drawRect(x1, 0, x2-x1, height)
      if annotation['text']:
        font = QtGui.QFont(self.font, self.fsize)
        letter_pen = QtGui.QPen(QtGui.QColor( annotation['col_text']  ))
        p.setPen(letter_pen);         p.setFont(font)
        p.drawText( x1, height*0.6666, annotation['text'] )
    p.end()


class sequenceAnnotation(faces.Face):
  def __init__(self, geneObj, seqlength, color='#000000', fsize=10):
    faces.Face.__init__(self)
    self.fsize=fsize
    self.seqlength=seqlength
    self.color=color

  def update_items(self):
    colored_box=     QtGui.QGraphicsRectItem(0, 0, 0, numbered_box_height)
    colored_box.setBrush(QtGui.QBrush(QtGui.QColor( hexcolor )))
  ### UNFINISHED


class coloredRectangle(faces.Face):
  def __init__(self, color="#FFFFFF", w=20, h=20, line_size=0, line_color="#000000", text='', fsize=10, text_color="#000000", font='Courier'):
    faces.Face.__init__(self)
    self.color=color ;     self.w=w;     self.h=h;     
    self.line_size=line_size;     self.line_color=line_color;     
    self.text=text ;     self.fsize=fsize; self.text_color=text_color
    self.font=font
  def update_pixmap(self):
    # code pasted from SequenceFace  
    self.pixmap = QtGui.QPixmap(self.w, self.h)
    self.pixmap.fill()
    p = QtGui.QPainter(self.pixmap)  
    brush = QtGui.QBrush( QtGui.QColor(  self.color ))   #bkground
    #drawing rectangle
    p.fillRect(0,0, self.w, self.h, brush)      
    if self.line_size:
      pen=    QtGui.QPen( QtGui.QColor(  self.line_color ) , self.line_size, QtCore.Qt.SolidLine )
      p.setPen(pen)
      p.drawRect(0,0, self.w, self.h)        
    if self.text:
      font = QtGui.QFont(self.font, self.fsize)
      p.setFont(font)
      letter_pen = QtGui.QPen(QtGui.QColor( self.text_color  ))
      p.setPen(letter_pen)
      fm = QtGui.QFontMetrics(font)
      #interactive_mode()()
      #print   fm.height()  ,fm.leading(), fm.overlinePos(), fm.underlinePos(), fm.boundingRectChar()
      letter_height = fm.height()
      text_width    = self.fsize*len(self.text)
      x=   (self.w-text_width)    / 2
      y=   (self.h -letter_height)/2 + letter_height
      p.drawText(x, y, self.text)
    p.end()
    
class coloredPlot(faces.Face):
  """ Face to draw a colored plot (line-like) aligned to the sequence face. Guide argument must point to a evol_codon_face """
  def __init__(self, guide,   plots,  fsize=10, font='Courier', height=150, header=20, foot=5, color_bkg="#FFFDC6" ):
    faces.Face.__init__(self)
    self.guide=guide;    self.plots=plots
    self.fsize=fsize;    self.font=font
    max_v=-sys.maxint;   min_v=sys.maxint
    for plot in self.plots:
      max_v= max (   plot['y_values'] + [max_v]    )
      min_v= min (   plot['y_values'] + [min_v]    )
    self.max= max_v;  self.min=min_v
    self.height=height; self.header=header; self.foot=foot
    self.color_bkg=color_bkg

  def update_pixmap(self):
    width = self.guide.pixmap.width() ;     height= self.height ##guide.pixmap.height() *2   ### 
      #vertical separator, above and below the plot
    y_min=height - self.foot;       y_max=self.header
    ## header (with max, then titles)
    ## plot
    ## foot (with min above)
    self.pixmap = QtGui.QPixmap(width,height)
    self.pixmap.fill()
    p = QtGui.QPainter(self.pixmap)    
    font = QtGui.QFont(self.font, self.fsize) #set font
    fm = QtGui.QFontMetrics(font)
    p.setFont(font); 
    brush = QtGui.QBrush( QtGui.QColor(  self.color_bkg ) )   
    p.fillRect(0,0, width, height, brush)  # draw bkground
    p.setPen(    QtGui.QPen( QtGui.QColor(  '#FFFFFF' ) , 1, QtCore.Qt.DashLine )     )
    p.drawLine(0, y_min, width, y_min)   # draw min line
    p.drawLine(0, y_max, width, y_max)   # draw max line

    titles_x=self.fsize*len(str(round(self.max, 3))) + 3 
    for plot in self.plots:
      p.setPen(     QtGui.QPen( QtGui.QColor(  plot['color'] ) , 1, QtCore.Qt.SolidLine )     )
      for pos in range(len(plot['y_values'])):
        #pos 0 based on cds alignment, codon based
        x1, x2 = self.guide.x_of_position(  pos+1  ) 
        new_x  = average( [x1, x2] ) 
        new_y  =  height -  self.foot  - (height-self.header-self.foot) * (plot['y_values'][pos] - self.min) / (self.max - self.min)      #nb 0 means at top
        if not pos==0:           
          p.drawLine(old_x,  old_y, new_x, new_y)   ### DRAWING PLOT
        old_x, old_y = new_x, new_y
      p.setPen(    QtGui.QPen( QtGui.QColor(  '#333333' ) , 1, QtCore.Qt.SolidLine )  )  ## min-max color
      p.drawText(0, y_min, str(round(self.min, 3)))  # draw min value
      p.drawText(0, y_max, str(round(self.max, 3)))  # draw max value
      p.setPen(     QtGui.QPen( QtGui.QColor(  plot['color'] ) , 1, QtCore.Qt.SolidLine )     )      
      p.drawText(titles_x,  self.header-5 , plot['title'])   # draw title   ## 5 = title foot within header
      titles_x+=self.fsize*len(plot['title']) + 3 #spacer titles:3
    p.end()

  def title_faces(self):
    """ Returns faces that serve as titles. Output is a list like:  [   title_face1, ... title_faceN ] with one element for each plot object."""
    out=[]
    for plot in self.plots:      out.append(   faces.TextFace( str(plot['title']) , fgcolor=plot['color'], fsize=self.fsize )   )
    return out

class selectedSitesFace(faces.Face):
  """ to display a bar like graph aligned to the positions positively selected. Should be coupled with an evolCodonsFace, and added as header or footer.
  Must be initialized with a sites_list=   a list with one element per codon, each element is either None or a positively_selected_site instance
  and also: a evolCodonsFace instance, which will be aligned to this face""" 
  def __init__(self,  sites_list=None, codon_face=None, height=200, min_kaks=1.0, max_kaks=None):
    faces.Face.__init__(self)
    self.sites_list=sites_list
    self.codon_face=codon_face
    self.height=height   ### NB the real hight is a little bit more to be sure no circle is partially left out
    self.set_kaks_scale( min_kaks, max_kaks)


  def set_kaks_scale(self, min_kaks, max_kaks=None): 
    """ If the max is set to None, the max value (+0.5) among the ones present will be used"""
    self.min_kaks=min_kaks
    if max_kaks is None:      max_kaks=   0.5 +  max([site.post_mean  for site in self.sites_list if not site is None  ])
    self.max_kaks=max_kaks  

  def update_pixmap(self):

    min_opacity=0.4
    max_opacity= 0.7   #theorical max opacity for computation. when significant, opacity becomes 1 anyway
    vertical_spacer=25
    size_text=20

    total_width= self.codon_face.x_of_position(  len(self.sites_list)  )[1]  #the right boundary of the last position
    total_height=  vertical_spacer+ self.height +vertical_spacer
    self.pixmap = QtGui.QPixmap(total_width,  total_height   )
    self.pixmap.fill()
    p = QtGui.QPainter(self.pixmap)
    font = QtGui.QFont("Courier", size_text)
    fm = QtGui.QFontMetrics(font)
    p.setFont(font); 
    p.setPen(  QtGui.QPen(QtGui.QColor( "#000000" )) )
    bg_brush=QtGui.QBrush(QtGui.QColor(  "#FFFFCC" ));       p.setBrush(bg_brush)
    p.drawRect(  QtCore.QRect(1,vertical_spacer, total_width-2, total_height-2*vertical_spacer )  )

    #bg_internal_brush=QtGui.QBrush(QtGui.QColor(  "#FFFFFF" ))
    #p.fillRect(0, vertical_spacer, total_width, self.height, bg_internal_brush)
    
          
    for p0, site in enumerate(  self.sites_list  ):
      position=p0+1
      if site is None:  continue
      
      x1, x2  = self.codon_face.x_of_position(   position   )

      bg_color= '#CC0000'; pen_color='#CC0000';   pen_width=1

      if is_significant( 1-site.probability ):    
        opacity=    1.0
        pen_color="#000000"
        pen_width=3
      else: 
        opacity =  min_opacity+  (max_opacity-min_opacity )*  ( (site.probability  -  0.5  ) / ( 1.0   -  0.5  ) )  ## getting value normalizing probability between 0.5 and 1.0 p

      letter_brush = QtGui.QBrush(QtGui.QColor(  bg_color ))   #bkground
      letter_pen = QtGui.QPen(QtGui.QColor( pen_color  ))
      letter_pen.setWidth(  pen_width  )
      p.setPen(letter_pen);       p.setBrush(letter_brush);        
      p.setOpacity( opacity )

      y= vertical_spacer+   self.height - self.height * ( site.post_mean - self.min_kaks) / (self.max_kaks - self.min_kaks)

      circle_width=x2-x1+1
      y_top = y-circle_width/2.0
      p.drawEllipse(   QtCore.QRect(  x1, y_top  , circle_width, circle_width  )  )
      p.setOpacity(1)  # Restore opacity 
  

      print 'drawing pos ', p0, '  --> ', y_top, opacity


    #self.probability=probability
    #self.post_mean=post_mean
    #self.deviation= deviation

    p.end()        

class evolCodonsFace(faces.Face):
  """ basically a seq face, thought to display syn and nonsyn changes. initialise with seq=sequence 
  use fsize to modify the size of text 
  use translate if you want to see aminoacids
  use zoom if you want to have an overview alignment, with letters no readable, but colored --> zoom=False ; optionally you can have one or more regions readable: use zoom=[[start, end], [start2, end2] , .. ]; by default (zoom=None) all letters are readable
  
  """
  def __init__(self, seq, fsize=10, translate=False, zoom=None):
    faces.Face.__init__(self)
    self.seq=seq         
    self.fsize= fsize
    self.translate=translate
    if len(self.seq)%3 !=0: raise ValueError, "evolCodonsFace ERROR the sequence must be composed of codons (length multiple of 3)"
    self.codon_colors={}    ## annotates the color for a certain codon position !!! NOTE codon positions are 0 based
    self.type2info = {"s": [2, fsize*2],                        "m": [fsize, fsize*2],                      "*": [2, fsize*2]                  }
    if zoom is None:      self.zoomed_regions= [[1, len(seq)/3]]
    elif not zoom:      self.zoomed_regions= []
    else:               self.zoomed_regions= zoom
    
  def update_pixmap(self):
    # code pasted from SequenceFace  
    font = QtGui.QFont("Courier", self.fsize)
    fm = QtGui.QFontMetrics(font)
    height = fm.leading() + fm.overlinePos() + fm.underlinePos()
    #width  = fm.size(QtCore.Qt.AlignTop, self.seq).width()

    codon_spacer=4
    #computing width and height of whole drawing seq area
    n_letters_zoomed=0
    for start, end  in self.zoomed_regions:      n_letters_zoomed+= end-start+1 
    if self.translate:            w=  n_letters_zoomed * self.type2info['m'][0]
    else:                         w=  n_letters_zoomed * (self.type2info['m'][0]*3     +  codon_spacer )
    w+=    (   len(self.seq)/3 - n_letters_zoomed   ) * self.type2info['s'][0] #### adding unzoomed positions width
    h=0
    if self.zoomed_regions:                h=    self.type2info['m'][1]
    if len(self.seq) - n_letters_zoomed:   h= max([h,  self.type2info['s'][1] ])

    self.pixmap = QtGui.QPixmap(w,h)
    self.pixmap.fill()
    p = QtGui.QPainter(self.pixmap)
    x = 0
    y = height - fm.underlinePos()*2
    p.setFont(font)

    for codon_index in range(len(self.seq)/3):    #0-based!
      codon=self.seq[codon_index*3:codon_index*3+3]
      #gettign color from self.codon_colors
      codon_pen_color="#000000"
      codon_bkg_color="#FFFFFF"
      opacity=1.0
      if self.codon_colors.has_key(codon_index):
        codon_pen_color=self.codon_colors[codon_index][0]
        codon_bkg_color=self.codon_colors[codon_index][1]
        opacity=self.codon_colors[codon_index][2]

      letter_brush = QtGui.QBrush(QtGui.QColor(  codon_bkg_color ))   #bkground
      letter_pen = QtGui.QPen(QtGui.QColor( codon_pen_color  ))
      p.setPen(letter_pen)

      is_zoomed=self.is_pos_zoomed( codon_index+1 )
      #print codon_index+1, is_zoomed

      if is_zoomed:
        letter_width= self.type2info['m'][0]
        letter_height=self.type2info['m'][1]
      else:
        letter_width= self.type2info['s'][0]
        letter_height=self.type2info['s'][1]

      if not self.translate:
        if not is_zoomed:
          p.setOpacity( opacity )
          p.fillRect(x,0, letter_width, letter_height, letter_brush)
          p.setOpacity(1)  # Restore opacity 
          x += letter_width
        else:
          for letter in codon:
            letter = letter.upper()
            #p.fillRect(x,0, width_for_sequence, height,letter_brush)
            p.setOpacity( opacity )
            p.fillRect(x,0, letter_width, letter_height, letter_brush)
            p.setOpacity(1)  # Restore opacity 
            p.drawText(x, y, letter)
            x += letter_width
          if is_zoomed:        x+=codon_spacer
      else:
        letter=transl(codon) ### drawing aminoacids!
        p.setOpacity( opacity )
        p.fillRect(x, 0, letter_width, letter_height, letter_brush)
        p.setOpacity(1)  # Restore opacity 
        if is_zoomed:
          p.drawText(x, y, letter)

        x += letter_width #x += float(width_for_sequence)/ ( len(self.seq)/3)

    p.end()

  def is_pos_zoomed(self, position):
    """returns True is the position ( 1-based, codon-based) is in a zoomed region (or everything is zoomed) """
    if not self.zoomed_regions: return False
    for st, end in self.zoomed_regions:
      if position >= st and position <= end: return True
    return False    

  def set_color_codon_changes(self, ref_node_seq, color_scheme=basicSynNonsynColorScheme, mode=0):
    """This method compares the .sequence in ref_node and the sequence in self.seq   (the sequence in each node.sequence and their face.seq should be the same!) and sets the colors in self.codon_pen_colors accordingly. A different color is chosen for nonsyn and syn changes, and the colors are picked from colorScheme.
    Depending on the argument mode:
    mode -> 0:  skip the positions in which a color is already defined  (don't overwrite)
    mode -> 1:  parses all the positions, replace with a new color in case a change is called (overwrite)
    mode -> 2:  reset all colors in sequence, set colors from scratch   (reset)    
      """
    if len(self.seq)!=len(ref_node_seq): raise Exception, "evolCodonsFace->set_color_codon_changes ERROR self.seq and ref_node.sequence must have the same length! "
    if mode==2: self.codon_colors={}  #resetting codon_pen_colors (deleting all existing annotated pos)
    
    for codon_index in range(len(self.seq)/3):    
      if mode==0 and self.codon_colors.has_key(codon_index): continue  ### skipping colored position if don't overwrite

      codon_self= self.seq[codon_index*3:codon_index*3+3]
      codon_ref = ref_node_seq[codon_index*3:codon_index*3+3]
      if '-' in codon_self and codon_self!='---': raise Exception, "evolCodonsFace->set_color_codon_changes ERROR there are gaps which are not multiple of 3! the alignment must be codon based! "
      if '-' in codon_ref and codon_ref!='---': raise Exception, "evolCodonsFace->set_color_codon_changes ERROR there are gaps which are not multiple of 3! the alignment must be codon based! "
      
      if codon_self == codon_ref and codon_self != '---' :
        continue # don't care, let's go to next codon 
      else:
        aa_self= transl( codon_self )
        aa_ref=  transl( codon_ref)
        if aa_self=='X' or aa_ref=='X' or aa_self=='-' or aa_ref=='-':
          #unknown
          #print self, codon_index, color_scheme.unknown_bkg   ##debug
          self.codon_colors[codon_index]=   [color_scheme.unknown_pen, color_scheme.unknown_bkg, 1.0] ####### setting fixed opacity for unknown color pos
        elif aa_self == aa_ref:
          #syn!
          self.codon_colors[codon_index]=   [color_scheme.syn_pen, color_scheme.syn_bkg, color_scheme.opacity]
        elif aa_self=='*':
          self.codon_colors[codon_index]=   [color_scheme.nonsense_pen, color_scheme.nonsense_bkg, color_scheme.opacity]
        else:
          #non syn
          self.codon_colors[codon_index]=   [color_scheme.non_syn_pen, color_scheme.non_syn_bkg, color_scheme.opacity]


  def x_of_position(self, pos):
    """Returns the (X1, X2) coordinate of  the horizontal boundaries of the box where position pos is drawn. Pos is 1-based, codon-based. This is useful if you want to draw things aligned to it"""
    x=0;       codon_spacer=4
    for index in range(pos-1):   #index is 0 based
      is_zoomed=self.is_pos_zoomed( index+1 )
      #print codon_index+1, is_zoomed
      if is_zoomed:        letter_width= self.type2info['m'][0]
      else:                letter_width= self.type2info['s'][0]
        
      if self.translate:        x+=   letter_width
      else:                     x+=   letter_width*3 + codon_spacer

    last_is_zoomed= self.is_pos_zoomed( pos )
    if last_is_zoomed:   letter_width= self.type2info['m'][0]
    else:                letter_width= self.type2info['s'][0]
    if self.translate:      last_width=letter_width
    else:                   last_width=letter_width*3
  
    return  [x, x+last_width]
  



def get_color_kaks_change( kaks, average_kaks, kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme ):
  """Utility to get the color for a kaks change analysis  """
  if kaks > average_kaks: 
    if kaks> kaks_neutral_boundaries[1]:          text_color=kaks_color_scheme.positive
    else:                                                          text_color=kaks_color_scheme.neutral
  else:
    if kaks> kaks_neutral_boundaries[0]:          text_color=kaks_color_scheme.neutral
    else:                                                          text_color=kaks_color_scheme.purifying
  return text_color    

def get_color_free_kaks( kaks, kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme ):
  """ Utility to get the color for a single kaks"""
  if     kaks < kaks_neutral_boundaries[0]:     text_color=kaks_color_scheme.purifying
  elif   kaks <= kaks_neutral_boundaries[1]:    text_color=kaks_color_scheme.neutral
  else:                                                              text_color=kaks_color_scheme.positive
  return text_color  


############################################################################################################################

class codon_alignment(object):
  """Master class to manage all codeml analsys and coding sequence based evolutionary tests. Instanciate with a codon fasta alignment a tree with the same titles (at least the first word in the alignment, which can have other crap in the headers)
  When called with a ali_file and a tree_file, the codon_alignment is saved into a pickled cdn_ali file, whose path is specified by the path argument when calling codon_alignment(), or by the path of the ali_file without extension.
  If the path.cdn_ali is already existing, python checks if the object you would create is the same of the one is stored there (checksum check). If it is not, it crashes unless you called overwrite=True when calling codon_alignment(). in this case, it deleted previously stored information.
  
  Attributes: ...
  
  Methods:    ...  
  
   """
  #verbose=False
  def __init__(self, ali_file=None, tree_file=None, path=None, overwrite=False, clean=False):#, sp_naming_function=lambda x:x.split()['.'][-1]): 
    if not path is None and path.endswith('cdn_ali'): path= path[:-8]
    if ali_file and not path: path=abspath(join(ali_file.split('.')[:-1], '.'))
    self.path=path
    self.models={} # will store all models run for this codon_alignment
    self.tests={} # will store all tests run for this codon_alignment
    self.full_titles={}  # will store all full titles int he fasta alignment... only short titles are kept in the self.ali
    self.checksum=None
    self.ancestral_sequences=None
    self.ancestral_sequences_source='none'
    self.stop_codon_positions={}   # 0-based, codon-based --> a position means there's a stop codon in some sequence, doens't tell which
    self.derived_data={}
    self.ali=None; self.tree=None
    self.pep_ali=None
    self.tree_style=TreeStyle();     self.tree_style.show_leaf_name = False
    self.layout= None
    self.clean=clean
        
#    self.sp_naming_function=sp_naming_function
    if ali_file: 
      self.load_alignment(ali_file, clean=clean)
      if tree_file:         
        self.load_tree(tree_file)
        if not is_file(self.path+'.cdn_ali'):        self.save()
        else:
          c= checksum( self.tree.write(format=9)+'\n'+str(self.ali),  is_text=True)
          a=codon_alignment(path=path) #loading object from path.cdn_ali
          if a.checksum!=c:  # checksum from this ali and tree are not the same as the ones loaded from cdn_ali file!
            if overwrite:
              printerr('WARNING file checksum inconsistency. Deleting files in '+self.path+'.* , recreating them with the newly provided alignment and tree', 1)
              bash('rm -r  '+self.path+'.cdn_ali '+self.path+'.codeml')
              self.save()             
            else:
             raise Exception, "ERROR loading codon_alignment: the alignment and tree provided are not the same used to build the codon_alignment found in "+self.path+'.cdn_ali. If you wanted to overwrite it, you need to specify overwrite=True when calling codon_alignment()   --- old checksum: "'+a.checksum +'"   -- new checksum: "'+c+'"'
          else:          self.__dict__=a.__dict__  # taking loaded object
    elif path:      self.load_cdn_ali()

  def load_alignment(self, ali_file, silent=False, clean=False):
    self.ali=alignment(ali_file, short_titles=True, check_uniq_titles=True) 
    codons_to_clean={} #codons position , 0 for first , 1 for second and so on
    for title in self.titles(): 
      if len(title)>50: raise Exception, "ERROR the titles in the alignment and tree must be less than 50 chars or codeml will crash!"
      seq = replace( upper(  self.ali.seq_of(title)), 'U', 'T')  
      for codon_index  in range( len(seq)/3):
        codon=seq[codon_index*3:codon_index*3+3]      
        for char in codon:      
          if not char in '-ACGTN': raise Exception, 'load_alignment ERROR unknown character '+char+' in sequence: '+title
        if clean and 'N' in codon:           codons_to_clean[codon_index]=1
      self.ali.set_sequence(title, seq)

    if clean and codons_to_clean:
      if not silent: printerr('WARNING clean option is active -- removing codons in alignment (1-based): '+join(   [ str(i+1) for i in codons_to_clean.keys() ], ', '     ), 1)
      length_before_trim=self.ali.length()
      for title in self.titles(): 
        seq=''
        for pos in range(length_before_trim/3):
          if not pos in codons_to_clean:            seq+=   self.ali.seq_of(title)  [pos*3:pos*3+3]
        self.ali.set_sequence(title,      seq)
    
    #storing full titles
    for line in open(ali_file):
      if line.startswith('>'): self.full_titles[ line.strip().split()[0][1:] ] = line.strip()[1:]
    #checking integrity of alignment
    self.test_codons()
    #annotating stop codons    
    for p in range( self.ali.length()/3 ):
      for title in self.ali.titles():
        codon =   self.ali.seq_of(title)[p*3:p*3+3]
        if codon in ['TAG', 'TGA', 'TAA']: 
          if not silent: printerr('WARNING column in position: '+str(p+1)+' will not be considered in codeml because it contains a stop codon in one or more sequences (one is '+title+')' , 1)
          self.stop_codon_positions[p]=codon    
          break

  def protein_alignment(self):
    """ Returns a lazy-computed protein alignment, reflecting the codon alignemnt"""
    if self.pep_ali is None: 
      self.pep_ali=alignment()
      for t in self.titles():
        self.pep_ali.add(t, transl(self.seq_of(t), get_MMlib_var('opt')['sec'])   )
    return self.pep_ali

  def build_tree(self):
    """ build a protein-based ML tree with salva pipeline. careful: works only you have the files ready """
    tree_f= self.protein_alignment().build_tree(folder=self.path+'.phylogeny/')
    self.load_tree(tree_f)
    try:     bbash('ln -s '+abspath(tree_f)+' '+self.path+'.tree')
    except:  printerr('WARNING I could not link the computed tree:  '+abspath(tree_f)+' to '+self.path+'.tree', 1)

  def load_tree(self, tree_file):
    """ Load a tree from tree_file, checks it, and set it as self.tree"""
    self.tree=EvolPhyloTree(tree_file, parent=self)
    names_in_tree= sorted(  [i.name for i in self.tree.get_leaves() ])
    names_in_ali = sorted( self.titles() )
    if names_in_ali != names_in_tree: raise Exception, "load_tree ERROR the node names do not have a perfect correspondance with the first word of the titles in the alignment! \nTREE:\n"+join(names_in_tree, '\n')+'\n\nALI:\n'+join(names_in_ali, '\n')
    if not self.tree.is_binary():      raise Exception, "load_tree ERROR tree is not binary!"
    if not self.tree.is_rooted():       
      self.tree.set_outgroup( self.tree.get_midpoint_outgroup() )
      for n in self.tree.traverse():        n.parent=self

  def save(self, path=None):
    """ Stores the information of the codon_alignment in a cdn_ali file, as a pickle dumped-file"""
    if path is None: path=self.path 
    if not self.ali or not self.tree: raise Exception, "codon_alignment-> save ERROR can't save! the alignment and the tree must both be loaded!"
    #getting checksums for alignment and tree
    self.checksum= checksum(self.tree.write(format=9)+'\n'+str(self.ali), is_text=True)
    #deleting graphical information for saving -- also cause otherwise it does not work
    ts=self.tree_style
    self.tree_style =None
    for i in self.tree.traverse(): 
      i._speciesFunction=None
      try: del i.codon_face
      except: pass
      try:
        i.faces = treeview.main._FaceAreas()
        i._faces = treeview.main._FaceAreas()      
        i._tmp_faces = None
      except: raise
    #saving pickled obj to file
    fileh=open(self.path+'.cdn_ali', 'w')
    write('Saving pickled content to file: '+self.path+'.cdn_ali', 1)
    pickle.dump(self, fileh)
    fileh.close()
    self.tree_style=ts

  def load_cdn_ali(self):
    """ Load a previously pickled object in the cds_ali file"""
    ffile=self.path+'.cdn_ali'
    if not is_file(ffile):  raise IOError, "ERROR can't find cdn_ali file: "+ffile
    a=pickle.load( open(ffile, 'r') )
    self.__dict__= a.__dict__
    for node in self.tree.traverse():       
      node.parent=self
      node.column_index=0
    self.tree_style=TreeStyle();     self.tree_style.show_leaf_name = False

  def seq_of(self, title):     return self.ali.seq_of(title)
  def titles(self):            return  self.ali.titles()
  def full_titles(self):       return [ self.full_titles[t] for t in self.titles() ]
  def length(self): return self.ali.length()
  __len__ = length

  def summary(self, title=None):
    """ String summary of the object"""
    if title is None: title=' codon_alignment '
    o='---------------'+(' '+title+' ').center(50, '=')+'---------------\n'
    len_tab=len(o)-2 #actually it's length -1
    tmp_s='-path: '+self.path
    o+=(tmp_s).ljust(len_tab)+ '-'*int(len(tmp_s)<len_tab)+  '\n'
    if not self.tree:   o+=('-tree: not loaded').ljust(len_tab)+'\n'
    else:               
      n_leaves=0; n_ancestral=0
      for i in self.tree.traverse():
        if i.is_leaf(): n_leaves+=1
        else:           n_ancestral+=1
      o+=('-tree: '+str(n_leaves)+' leaves').ljust(len_tab)+'-\n'
    if not self.ali:    o+=('-ali: not loaded').ljust(len_tab)+'-\n'
    else:               o+=('-ali: '+str(self.ali.length()/3)+' codons').ljust(len_tab)+'-\n'
    if not self.models: o+=('-models: None').ljust(len_tab)+'-\n'
    else:               
      tmp_s='-models: '+join([self.models[m].name      for m in sorted(self.models.keys())  ] , ', ')
      o+=(tmp_s).ljust(len_tab)+ '-'*int(len(tmp_s)<len_tab)+  '\n'
    try:      o+=('-clean: '+str(self.clean)).ljust(len_tab)+'-\n'
    except:   pass
    o+='-'*(len_tab+1)
    return o         
  __str__ = summary

  def get_nodes(self, expr, no_root=False, no_leaves=False):
    """ returns the target nodes instances, given an expression as the one explained here: 
Nodes can be specified by their name if they are leaves, or using the first common ancestor principle if they are ancestral ( node1name+node2name ). Multiple nodes can be specified joining them with a comma (  node1name,node2name+node3name,node5name ). 
Some nodes can be specified with keywords:  root, all, all_leaves"""
    targets=[]
    if expr=='all':             targets=[i for i in self.tree.traverse() if (not no_root or not i.is_root()) and (not no_leaves or not i.is_leaf() ) ]
    elif expr in ['all_leaves', 'leaves']:    targets=[i for i in self.tree.traverse() if i.is_leaf() ]
    else:
      for it in expr.split(','):
        if it=='root':          targets.append(self.tree)
        elif it.count('+')==0: 
          try:        targets.append(  self.tree&it )   #searching node name provided in command line
          except:     raise ValueError, "ERROR can't find node with name: "+it
        else:
          try: 
            list_of_nodes_names=     it.split('+')
            targets.append( self.tree.get_common_ancestor(*list_of_nodes_names)         )
          except:     raise ValueError, "ERROR can't find common ancestors of nodes with name: "+str(list_of_nodes_names)
    return targets    

  def test_codons(self):
    """ Tests if the sequences are aligned codons. any gap should occur in multiple of 3"""
    l=len( self.ali.seq_of( self.titles()[0]  ) )
    for t in self.titles(): 
      if len(self.seq_of(t))!=l: raise Exception, "ERROR the sequences are not of the same length! "
    for codon_i in range(self.length()/3):
      for t in self.titles():
        codon= self.seq_of(t)[codon_i*3:codon_i*3+3]
        if '-' in codon and codon!='---': raise Exception, "ERROR this is not a codon alignment: all gaps must occur in multiple of 3!"

  def alignment_without_stops(self):
    """ Returns a copy of self.ali with stop codons removed (if any: otherwise, just returns self.ali)"""
    if not self.stop_codon_positions: return self.ali
    a=self.ali.copy()
    for position in sorted( self.stop_codon_positions.keys(), reverse=True):
      for title in self.ali.titles():
        a.set_sequence(title,    a.seq_of(title)[:position*3]+a.seq_of(title)[position*3+3:]   )
    return a 

  def set_ancestral_sequences(self, model_name='null'):
    """ Accepts as good the ancestral_sequences predicted by the specified model. 'sankoff' can be specified to derive ancestral states with this algorithm instead than parsing it from a computed model"""
    if model_name=='sankoff': 
      try:    self.load_ancestral_states_with_sankoff()
      except: self.derive_ancestral_sequences_with_sankoff()
    else:
      if isinstance( model_name, codeml_model ):  model_name=model_name.name
      if not model_name in self.models or not self.models[model_name].ancestral_sequences: raise Exception, "set_ancestral_sequence ERROR the model "+str(model_name)+' is not computed, or the ancestral sequences were not computed!'
      self.ancestral_sequences= self.models[model_name].ancestral_sequences
    #resetting all derived data as it may have used ancestral sequences, in case we changed it
    self.ancestral_sequences_source=model_name

  def derive_ancestral_sequences_with_sankoff(self, silent=False):
    """ Add ancestral states by predicting them with the sankoff algorithm, and save them to .sankoff file.  Do not use directly: use set_ancestral_sequences('sankoff') instead. """
    write('Running sankoff algorithm... ')
    self.ancestral_sequences = sankoff(self.tree, node2seqFn=  EvolPhyloTree.sequence)
    write('ok', 1)
    #resetting all derived data as it may have used ancestral sequences, in case we changed it
    sankoff_file=   self.path+'.sankoff'
    #saving ancestral states to a file 
    sankoff_file_h= open(sankoff_file, 'w')
    for node in self.ancestral_sequences:      print >> sankoff_file_h, node.get_name()+'\t'+self.ancestral_sequences[node]
    sankoff_file_h.close()

  def load_ancestral_states_with_sankoff(self):
    """ Load and set as ancestral states of sequences those previously computed with sankoff """
    sankoff_file=   self.path+'.sankoff'
    sankoff_file_h=open(sankoff_file)
    self.ancestral_sequences={}
    for line in sankoff_file_h:
      splt=line.strip().split('\t')
      self.ancestral_sequences[  self.get_nodes( splt[0]  )[0]  ] = splt[1]
    sankoff_file_h.close()

  def set_branch_length(self, source=None, silent=False):
    """ Accepts as good the branch lengths predicted by the specified model or source.
    If no arguments are provided, the model 'null' is used.
    source valid arguments are loaded models instances or names, or any float-returning method of the nodes: 
     dKs, dKa, sum_dKs, sum_dKa, Ks, Ka, KaKs, ...  """
    if source is None: source='null'
    #if not self.tre: Exception
    is_model= self.models.has_key(source)
    try:     is_model= is_model or self.models.has_key(source.name)
    except:  pass    
    if is_model:
      ## model provided
      model_name=source
      if isinstance( model_name, codeml_model ):  model_name=model_name.name
      if not self.models[model_name].branch_lengths: raise Exception, "set_branch_length ERROR branch lengths in model "+str(model_name)+' were not found!'
      for node in self.tree.traverse():
        node.dist = self.models[model_name].branch_lengths[node]
    else:
      ## method source was provided
      if not source in dir(self.tree):       raise Exception, "set_branch_length ERROR the method "+str(source)+' was not found!'
      for node in self.tree.traverse():
        try:    
          distance= eval('node.'+source+'()')
          node.dist= distance
        except:           ## calling source method for node failed. maybe it's just exceptions, like, dKs for root, but we'll print a warning
          if not silent:  printerr('set_branch_length WARNING failed call function: '+source+' in node: '+node.get_name())
          node.dist=0.0

  def has_derived_data(self, function, title=''):     return self.derived_data.has_key(function.__name__) and self.derived_data[function.__name__].has_key(title)
  def get_derived_data(self, function, title=''):     return self.derived_data[function.__name__][title]
  def save_derived_data(self, function, title='', data_value=None):
    """ Generic entry for lazy computation. Saves a data value generated by a function in self.derived_data[function.__name__][title], where title is generally the name of the node for which it was computed, or an empty string if this doesn't apply"""
    if not     self.derived_data.has_key(     function.__name__): 
      self.derived_data[function.__name__]={}
    self.derived_data[function.__name__][title]=data_value

  def reset_derived_data(self, titles=None, function_names=[]):
    """ reset data containers used for lazy computing, for certain function_names and titles. If any of the two or both are not specified, all data for what you specified is deleted ### NOTE: THIS FUNCTION IS NEVER USED BY INTERNAL METHODS, AS EVERYTHING THAT IS COMPUTED IS USEFUL """
    raise Exception # DEBUG 
    write('Codon alignment: '+self.path+ ': resetting derived data', 1)
    if not function_names:      function_names=self.derived_data.keys()  #not defined -> deleting all derived data
    for function_name in function_names:
      if not self.derived_data.has_key(function_name): continue
      if titles is None: titles=self.derived_data[function_name].keys()
      for title in titles: 
        if not self.derived_data[function_name].has_key(title): continue
        del self.derived_data[function_name][title]
      
  def add_codeml_model(self, model, **args):
    """ add a codeml_model. Use a codeml_model class as argument. arguments that you will call for the model during instanciation can be provided as **args here. Note: if the model is already present, it is not added and the function complaints a bit"""
    m=model(parent=self, **args)
    self.models[m.name]=m
    return m

  def get_model(self, name):
    if not self.has_model(name): raise Exception, "ERROR no model with name: "+name+ " were found for this object! do you remember adding it? "
    return      self.models[name] 

  def has_model(self, name):    return self.models.has_key(name)

  def add_codeml_test(self, test, **args):
    """ add a test. first argument is a codeml_test CLASS! (not an instance). arguments that you would use to instanciate are provided as **args """
    t=test(parent=self, **args)
    self.tests[t.name]=t
    return t

  def get_test(self, name):
    if not self.has_test(name): raise Exception, "ERROR no test with name: "+name+ " were found for this object! do you remember adding it? "
    return      self.tests[name] 
  
  def has_test(self, name):    return self.tests.has_key(name)

  def null_model_omega(self):  return self.get_model('null').w['all']

############################################################################################################################
class positively_selected_site(object):
  """ Simple object with a few attributes"""
  def __init__(self, probability=None, post_mean=None, deviation=None):
    self.probability=probability
    self.post_mean=post_mean
    self.deviation= deviation


class loading_model_error(Exception): 
  """ """

class codeml_model(object):
  """ class for codeml models. They must be instanciated with a codon_alignment as parent, and have methods for running the model  -> run(),  or loading it -> load() 
  initialization: 
  __init__(self, name=None, parent=None, build=True, run=True, lazy=True,  **parameters):
  parent is the codon_alignment instance to which it refers. If build / run are True, the model is built / run . If lazy is True, these operatinos if active are done anyway only if the files are not found already, in which case, it tries to load them
  
  After running or loading, attributes such as these are available:
  -self.lnL =  log likelihood
  -self.np  =   number of parameters in the model
  -self.k   =  kappa (ts/tv)
  -self.w   = omegas  (dN/dS)  ; this is an hash of node -> values to allow complex models. In fixed omega models, a single key is present:   'all':value
  also, the parent codon_alignment is populated with node specific results
  """

  default_parameters={'verbose': ['0', '  1: detailed output, 0: concise output\n'], 'CodonFreq': ['2', ' 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n        ndata  = 10\n'], 'ncatG': ['3', '  # of categories in dG of NSsites models\n'], 'cleandata': ['0', ' remove sites with ambiguity data (1:yes, 0:no)?\n 0: ignore, -1: random, 1: initial, 2: fixed\n'], 'NSsites': ['0', '  0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n  5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\n  10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\n  13:3normal>0\n'], 'fix_omega': ['0', ' 1: omega or omega_1 fixed, 0: estimate\n'], 'clock': ['0', ' 0:no clock, 1:clock; 2:local clock; 3:TipDate\n'], 'seqfile': ['SEQFILE', ' sequence data file name\n'], 'runmode': ['0', '  0: user tree; 1: semi-automatic; 2: automatic\n  3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\n'], 'fix_kappa': ['0', ' 1: kappa fixed, 0: kappa to be estimated\n'], 'fix_alpha': ['1', '  0: estimate gamma shape parameter; 1: fix it at alpha\n'], 'Small_Diff': ['.5e-6', ''], 'method': ['0', ' 0: simultaneous; 1: one branch at a time\n'], 'fix_rho': ['1', ' 0: estimate rho; 1: fix it at rho\n'], 'Malpha': ['0', '  different alphas for genes\n'], 'aaDist': ['0', ' 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n 7:AAClasses\n'], 'RateAncestor': ['0', ' (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n'], 'aaRatefile': ['wag.dat', ' only used for aa seqs with model=empirical(_F)\n dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own\n'], 'outfile': ['codeml.out', ' main result file name\n'], 'icode': ['0', ' 0:universal code; 1:mammalian mt; 2-11:see below\n'], 'rho': ['0.', ' initial or fixed rho,     0:no correlation\n'], 'alpha': ['0.', '  initial or fixed alpha, 0:infinity (constant rate)\n'], 'seqtype': ['1', ' 1:codons; 2:AAs; 3:codons-->AAs\n'], 'omega': ['.4', ' initial or fixed omega, for codons or codon-based AAs\n'], 'getSE': ['0', " 0: don't want them, 1: want S.E.s of estimates\n"], 'noisy': ['9', '  0,1,2,3,9: how much rubbish on the screen\n'], 'Mgene': ['0', ' 0:rates, 1:separate;\n'], 'kappa': ['2', ' initial or fixed kappa\n'], 'treefile': ['TREEFILE', ' tree structure file name\n'], 'model': ['0', ' models for codons:\n 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n models for AAs or codon-translated AAs:\n 0:poisson, 1:proportional,2:Empirical,3:Empirical+F\n 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)\n']}
  parameter_names =['seqfile', 'outfile', 'treefile', 'noisy', 'verbose', 'runmode', 'seqtype', 'CodonFreq', 'clock', 'aaDist', 'aaRatefile', 'model', 'NSsites', 'icode', 'Mgene', 'fix_kappa', 'kappa', 'fix_omega', 'omega', 'fix_alpha', 'alpha', 'Malpha', 'ncatG', 'fix_rho', 'rho', 'getSE', 'RateAncestor', 'Small_Diff', 'cleandata', 'method']  #ordered list to get the ctl file in the same order as the standard ones
  # format: 
  def __init__(self, name=None, parent=None, build=True, run=True, lazy=True,  **parameters):
    self.lnL=None;     self.np=None; self.k=None   # log likelihod, number of parameters,  kappa (trasitions/transversions)
    self.w = None  # hash node2value of omegas (dN/dS) ; in case it is fixed, the key 'all' contains the value.
    self.sites={}  #hash k: position 1-based, codon-based  -> value:   site class, with probability and stuff
    self.ancestral_sequences=None # hash node2sequence #  may be not filled, in many cases we don't compute ancestral states
    self.branch_lengths=None      # hash node2branch_length 
    self.name=name; 
    self.parent=parent  # codon_alignment instance
    #self.target=target  # used in some models, can be a node for example
    self.parameters={}  
    for p in parameters:
      if not p in codeml_model.default_parameters: raise Exception, "ERROR "+str(p)+ ' is not a valid codeml model parameter '
      self.parameters[p]=parameters[p]
    if build: 
      if not lazy or not is_file(self.model_folder()+'codeml.ctl') or  self.ctl_file()+'\n' !=  join( [l for l in open( self.model_folder()+'codeml.ctl')], '' ):
        self.build()
    if run:   
      loaded=False
      if lazy:
        try: 
          self.load() 
          loaded=True
        except loading_model_error: pass
      if not loaded:        
        self.run()

  def model_folder(self, absolute=False):
    """ get path to model folder"""
    codeml_folder=Folder(self.parent.path+'.codeml')
    if absolute: return abspath( Folder(codeml_folder+self.name) )
    else:        return Folder(codeml_folder+self.name)   

  def tree_file(self):
    """ returns the text in newick used when running the model. must override this method for models in which you need flags along the tree"""
    return self.parent.tree.write(format=9)

  def ctl_file(self):
    """ Returns the string for the ctl file associated to this model. It is thought to be put in the model folder (self.parent.path+'.codeml/'+self.name) so some filepaths are relative """
    if not self.parent: raise Exception, "ERROR no parent codon_alignment instance is defined for this model, can't call this method!"
    o=''
    for p in codeml_model.parameter_names:
      if self.parameters.has_key(p):  value=self.parameters[p]
      else:                           value=codeml_model.default_parameters[p][0]
      if   p=='treefile':             value='tree.newick'
      elif p=='seqfile':              value='cds_alignment.faa'
      comment= join( codeml_model.default_parameters[p][1].strip().split('\n'), '\n*')
      o+= p.rjust(16)+' = '+value +' * '+comment+'\n'
    return o

  def build(self):  
    """ Builds the model, meaning that it prepares the files needed for later execution: alignemnt in paml format, tree, codeml.ctl"""  
    if not self.parent: raise Exception, "ERROR no parent codon_alignment instance is defined for this model, can't call this method!"    
    model_folder= self.model_folder()
    write_to_file(  self.parent.alignment_without_stops().codeml_format(),   model_folder+'cds_alignment.faa'  )
    write_to_file(  self.tree_file(),     model_folder+'tree.newick'  )
    write_to_file(  self.ctl_file(),      model_folder+'codeml.ctl'  )

  def run(self):
    """ Runs codeml for this model. Exception is raised is error status ends different than 0. other errors may be raise when loading it"""
    model_folder= self.model_folder()
    cmnd_run= 'cd '+model_folder+'; codeml &> codeml.log'
    if not no_run_mode:     
      b=bash(cmnd_run)  ## running
      if b[0]: raise Exception, "ERROR running model: "+self.name+" on codon_alignment "+self.parent.path       
      self.load()
    else:                   
      write('printing job for: '+self.name+' ')
      print >> no_run_mode, cmnd_run


  def load(self):
    """ must be defined in each subclass of model"""
    raise Exception, "ERROR the load function must be defined for each subclass of codeml_model. It is not in this class: "+str(self.__class__)

  def load_basic(self):
    """ Reads the following things: lnL, np, ancestral_states if they were run"""
    if not is_file(self.model_folder()+'codeml.out' ): raise loading_model_error, "ERROR file not found:  "+self.model_folder()+'codeml.out'
    for line in open(self.model_folder()+'codeml.out'):
      if line.strip().startswith('lnL'):
        self.lnL=   float(line.split()[-2])
        self.np =   int(line.split('np:')[1].split(')')[0].strip())  
      if line.strip().startswith('kappa'):
        self.k  = float(line.strip().split()[-1])
    try: self.load_branch_lenghts()
    except: raise loading_model_error, "ERROR loading branch lengths from file: "+self.model_folder()+'codeml.out'
    if not self.lnL or not self.np or not self.k or not self.branch_lengths: raise loading_model_error, "ERROR lnL or np or kappa or branch_lengths were not found in file:  "+self.model_folder()+'codeml.out'
    try:    self.add_ancestral_sequences() 
    except IOError: pass 

  def add_ancestral_sequences(self):
    """ add ancestral states from the rst file. Add possible stop codon positions using sankoff"""
    self.ancestral_sequences =  add_ancestral_states( self.parent.tree,  self.model_folder()+'rst', add=False )
    # self.ancestral_sequences lacks stop codons positions, if any  --> sankoff to reconstruct!
    if self.parent.stop_codon_positions:
      add_to_ancestral={} #node2seqadd 
      stop_codons_ali = alignment()
      for title in self.parent.titles():
        s=''
        for stop_p in sorted( self.parent.stop_codon_positions.keys() ) :
          s+= self.parent.seq_of(title)[stop_p*3:stop_p*3+3] #codon
        stop_codons_ali.add(title, s)
      stop_codon_anc_states= sankoff(self.parent.tree, node2seqFn=  lambda x: stop_codons_ali.seq_of( x.name )  )
      for node in self.ancestral_sequences:
        seq=self.ancestral_sequences[node]
        for stop_index, stop_p in enumerate(  sorted( self.parent.stop_codon_positions.keys() )  ):
          codon_to_add= stop_codon_anc_states [node] [stop_index*3:stop_index*3+3]
          seq= seq[:stop_p*3] + codon_to_add + seq[stop_p*3:]
        self.ancestral_sequences[node]=seq    

  def load_variable_w(self):
    """ This function is used whenever you already run a viarble omega model and you want those W values to be available.  Values are stored in self.w  as hash node2w (in model class)""" 
    b=bash('grep -A 1 "w ratios as labels for TreeView:" '+self.model_folder()+'codeml.out | tail -1 ')
    if b[0]: raise loading_model_error, "load_variable_w ERROR line: 'w ratios as labels for TreeView' not found in: "+self.model_folder()+'codeml.out'
    t2_text=b[1]
    t2_text=replace( replace(   t2_text,   "'", ''), "#", ':-')  ## making ete2 read the free kaks as if they were distances
    t2=PhyloTree(t2_text)
    self.w={}
    for node1, node2 in mapping_trees(self.parent.tree, t2):
      if node2.dist<0:  #every node apart from root should satisfy this
        self.w[ node1 ] =   -node2.dist

  def load_selected_sites(self):
    """ Reads the BEB output indicating the positively selected sites. NOTE: this is available only for a limited number of models (positive selection and similar, only site models). If no BEB output is found, an exception is raised"""
    if not is_file(self.model_folder()+'codeml.out' ): raise loading_model_error, "ERROR file not found:  "+self.model_folder()+'codeml.out'
    fileh= open(self.model_folder()+'codeml.out')
    line= fileh.readline()
    while not line.startswith(    'Bayes Empirical Bayes (BEB) analysis'   ):      line= fileh.readline()
    if not line: raise loading_model_error, "ERROR Bayes Empirical Bayes (BEB) output not found in : "+self.model_folder()+'codeml.out'
    while line and not line.startswith('The grid'):  #just in case no site is significant
      while line and not line.strip().startswith('Pr(w>1)'):  line= fileh.readline()
      line= fileh.readline(); line= fileh.readline()  
      #now in first line of BEB
      while line and line.strip(): #reading until the first blank line
        splt=line.strip().split()        
        position_1_based=int(splt[0])
        probability=    float(splt[2].strip('*'))       #that w>1
        post_mean  =     float(splt[3])
        deviation  =     float(splt[5])      
        self.sites[position_1_based] = positively_selected_site( probability=probability, post_mean=post_mean, deviation=deviation  )
        line= fileh.readline() 
    fileh.close()

  def load_branch_lenghts(self):
    """ Loads the branch length (number of nucleotide substitutions per codon) from the codeml.out file and stores them in self.branch_lengths. Then they can be accepted in the codon_alignment with set_branch_length functions """
    b=bash('grep -A 4  "tree length = " '+self.model_folder()+'codeml.out  | tail -1')
    if b[0]: raise loading_model_error, "load_branch_lenghts ERROR line: 'tree length = ' not found in: "+self.model_folder()+'codeml.out'
    t2_text=replace(b[1], ' ', '')  ## otherwise lengths are not read!
    t2=PhyloTree(t2_text)
    self.branch_lengths={}   
    for node1, node2 in mapping_trees(self.parent.tree, t2):
      self.branch_lengths[node1] = node2.dist
    
  def omega_for(self, node):
    """ Facilitating retrieval of specific omega values in this model"""
    if not isinstance( node, EvolPhyloTree):
      node=(self.parent.tree)&node
    if 'all' in self.w: return self.w['all']
    else: 
      if not self.w.has_key(node): 
        descr=str(node)
        try: descr+="\n name: "+node.name
        except: pass  
        raise Exception, "codeml_model-> omega_for ERROR can't find node: "+descr
      else: return self.w[node]

  def LRT(self, null_model):
    """ performs a Likelihood Ratio Test of this model against another (null) model. It uses an approximated Chi square distribution to compare the likelihood values and the number of parameters (degrees of freedom) and evaluate if this model is statistically significantly  better than the null model.
    NOTE: one of the LRT assumptions is that the null model is a specific case of the alternative (the self) model, with less parameters. This is not checked by this function!    """
    D= 2*( self.lnL    - null_model.lnL )
    ddf=  self.np - null_model.np 
    if ddf<0: raise Exception, "ERROR the null model cannot have more parameters than this model!"
    p_value = 1 - chi2.cdf(D, ddf)        
    return p_value 


class null_codeml_model(codeml_model):
  """ name='null'
  Standard null model M0. W is fixed along the tree and along the sites. The ancestral states are computed by default with this model, it can be turned off with RateAncestor='0' (See its __init__ function )"""
  def __init__(self,  parent=None, RateAncestor='1'):
    codeml_model.__init__(self, parent=parent, name='null', RateAncestor=RateAncestor) #setting self.name='null'

  def load(self):
    """ Reads the following things: lnL, np, ancestral_states if they were run, and omega (w)"""
    self.load_basic() 
    self.w  = {};
    for line in open(self.model_folder()+'codeml.out'):
      if line.strip().startswith('omega'):
        self.w['all']=   float(line.strip().split()[-1])
    if not self.w: raise loading_model_error, "ERROR line with 'omega' not found in file: "+self.model_folder()+'codeml.out'

#### BRANCH MODELS

class free_w_codeml_model(codeml_model):
  """ name='free_w' ; also called M2.free , it is specified by model = 1 and NSsites = 0 in codeml.ctl. 
    It allows a different w for every branch in the tree. very parameter rich, and time consuming
    """
  def __init__(self,  parent=None):
    name= 'free_w'  
    codeml_model.__init__(self, parent=parent, name=name, model='1', )        

  def load(self):
    self.load_basic() 
    self.load_variable_w()


class many_w_codeml_model(codeml_model):
  """ Generic model for model=2 , Nsite=0 of codeml.  You can specify any combination of different w in different part of the tree. 
  The most important argument when initiating is     targets: this is a hash like: 
    key: node -> value  ; where the value is a numeric flag, starting from 1.  To say that you want the same w in different branches of the tree, you repeat the same numeric flag in multiple nodes.
    
  """
  def __init__(self,  parent=None, targets={}, descendants=False, name=None):
    if name is None: ## you specify a name just for subclasses of this model
      name= 'many_w.'
      for t in sorted( self.targets, key=lambda x:x.get_name() ):
        name+= t.get_name()+'@'+str(self.targets[t])+','
      name=self.name[:-1] 
    if not targets: raise Exception, "ERROR w_change_codeml_model no target nodes defined!"
    self.targets=targets
    codeml_model.__init__(self, parent=parent, name=name,  model='2')
    
  def tree_file(self):
    t=deepcopy(self.parent.tree)
    for n1, n2 in mapping_trees( self.parent.tree, t ):
      if n1 in self.targets:
        n2.add_feature('tag',  str(self.targets[n1] ) )
    return tagged_tree2codeml(t)

  def load(self):
    self.load_basic() 
    self.load_variable_w()

class w_change_codeml_model(many_w_codeml_model):
  """ Specific many_w_codeml_model where there's a single change, either only in a node, or in all branches below a node (descendants=True). """
  def __init__(self,  parent=None, target=None, descendants=False):
    flag_descendants={False:'S', True:'D'}[descendants]
    name='w_change_'+flag_descendants+'.'+target.get_name()
    if not descendants:
      many_w_codeml_model.__init__(self, parent=parent, name=name, targets={target:1})
    else:
      targets={}
      for n in target.traverse():        targets[n]=1
      many_w_codeml_model.__init__(self, parent=parent, name=name,  targets=targets)      
    self.target=target

  def w_target(self):  
    return self.omega_for(self.target)
    
  def w_rest_of_tree(self):
    """ Returns the "external" w in this w_change model, that is to say, the one for the rest of the tree apart from the target (or the target and its descendants) """
    try:   o= self.omega_for(self.target.up)
    except:  #in case the .up is the root, it will not have a omega, we have to look at the other child  
      other_child= [ n for n in self.target.up.get_children() if n!=self.target ] [0]
      o= self.omega_for(other_child)
    return o

#### SITE MODELS
class nearly_neutral_codeml_model(codeml_model):
  """Nearly neutral site model (M1a) described in the codeml manual."""
  def __init__(self,  parent=None):
    name= 'nearly_neutral'  
    codeml_model.__init__(self, parent=parent, name=name, model='0', NSsites='1' )        

  def load(self):
    self.load_basic() 
    #  missing information?    

class positive_selection_codeml_model(codeml_model):
  """Positive selection site model (M2a) described in the codeml manual """
  def __init__(self,  parent=None):
    name= 'positive_selection'  
    codeml_model.__init__(self, parent=parent, name=name, model='0', NSsites='2' )

  def load(self):
    self.load_basic() 
    self.load_selected_sites()

class beta_codeml_model(codeml_model):
  """Beta site model (M7) described in the codeml manual."""
  def __init__(self,  parent=None):
    name= 'beta'  
    codeml_model.__init__(self, parent=parent, name=name, model='0', NSsites='7' )        

  def load(self):
    self.load_basic() 
    #  missing information?    

class beta_positive_codeml_model(codeml_model):
  """Beta + positive selection site model (M8) described in the codeml manual."""
  def __init__(self,  parent=None):
    name= 'beta_positive'  
    codeml_model.__init__(self, parent=parent, name=name, model='0', NSsites='8' )

  def load(self):
    self.load_basic() 
    self.load_selected_sites()

class discrete_codeml_model(codeml_model):
  """discrete site model (M3) described in the codeml manual."""
  def __init__(self,  parent=None):
    name= 'discrete'  
    codeml_model.__init__(self, parent=parent, name=name, model='0', NSsites='3' )

  def load(self):
    self.load_basic() 
    #  missing information?    





############################################################################################################################

class codeml_test(object):
  """ routine to add right models and test them. typically you want to instanciate from codon_alignment  -> .add_codeml_test()  which gives you a test object from which you can call the results() function of each codeml_test . This varies, but it generally is a list of things, with the first object being a p_value of the test
  """
  def __init__(self, parent=None, run=True, **args):
    self.parent=parent
    self.name= None          # --> depends on the subclass of test and arguments (e.g. target)
    for a in args:
      exec('self.'+a+'=args[a]')
      self.p_value=None
    if run: self.run()

  def run(self): pass ## -> define this in subclasses
  def results(self): 
    if self.p_value is None: self.run()
    return [self.p_value]  #could returns something more  in subclass specific results function, so we say it always retursn a list with p_value as first value
    
class w_change_test(codeml_test):
  """ Initialise with a target=node, and descendants=(T/F)
  self.result() ->    [  self.p_value,  self.w_change  ]   where self.w_change is    [w_rest_of_tree, w_target_node ]
  """

  def __init__(self, parent=None, target=None, descendants=False, run=True):
    codeml_test.__init__(self, parent=parent, target=target, descendants=descendants, run=False)
    self.name=self.get_name()
    self.w_change=[]
    if run: self.run()

  def get_name(self):
    flag_descendants={False:'S', True:'D'}[self.descendants]
    return 'w_change_'+flag_descendants+'.'+self.target.get_name()

  def run(self):
    if not self.parent.has_model('null'): self.parent.add_codeml_model(null_codeml_model)
    w_change_model_name=self.get_name()
    if not self.parent.has_model(w_change_model_name): self.parent.add_codeml_model( w_change_codeml_model, target=self.target, descendants=self.descendants  )
    if no_run_mode: return
    m0=self.parent.get_model('null')
    m1=self.parent.get_model(w_change_model_name)
    self.p_value= m1.LRT(m0)
    w_rest_of_tree=m1.w_rest_of_tree()
    self.w_change= [ w_rest_of_tree     ,         m1.omega_for(self.target) ]
    
  def results(self):
    codeml_test.results(self)
    return [  self.p_value,  self.w_change  ]


class free_w_test(codeml_test):
  """ LRT test between null and free_w model"""
  def __init__(self, parent=None, run=True):
    codeml_test.__init__(self, parent=parent, run=False)
    self.name='free_w'
    if run: self.run()

  def run(self):
    if not self.parent.has_model('null'): self.parent.add_codeml_model(null_codeml_model)    
    if not self.parent.has_model('free_w'): self.parent.add_codeml_model( free_w_codeml_model)
    if no_run_mode: return
    m0=self.parent.get_model('null')
    m1=self.parent.get_model('free_w')
    self.p_value= m1.LRT(m0)


class positive_selection1_sites_test(codeml_test):
  """ LRT test between M1a (nearly neutral) and M2a (positive selection) model"""
  def __init__(self, parent=None, run=True):
    codeml_test.__init__(self, parent=parent, run=False)
    self.name='positive_selection1_sites'
    if run: self.run()

  def run(self):
    if not self.parent.has_model('nearly_neutral'): self.parent.add_codeml_model(nearly_neutral_codeml_model)    
    if not self.parent.has_model('positive_selection'): self.parent.add_codeml_model( positive_selection_codeml_model)
    if no_run_mode: return
    m1=self.parent.get_model('nearly_neutral')
    m2=self.parent.get_model('positive_selection')
    self.p_value= m2.LRT(m1)

class positive_selection2_sites_test(codeml_test):
  """ LRT test between M7 (beta) and M8 (beta+positive selection) model"""
  def __init__(self, parent=None, run=True):
    codeml_test.__init__(self, parent=parent, run=False)
    self.name='positive_selection2_sites'
    if run: self.run()    
    
  def run(self):
    if not self.parent.has_model('beta'): self.parent.add_codeml_model(beta_codeml_model)    
    if not self.parent.has_model('beta_positive'): self.parent.add_codeml_model( beta_positive_codeml_model)
    if no_run_mode: return
    m1=self.parent.get_model('beta')
    m2=self.parent.get_model('beta_positive')
    self.p_value= m2.LRT(m1)
    
class variable_w_sites_test(codeml_test):
  """ LRT test between M3 (discrete) and M0 (null -- one kaks) model"""
  def __init__(self, parent=None, run=True):
    codeml_test.__init__(self, parent=parent, run=False)
    self.name='variable_w_sites'
    if run: self.run()    
    
  def run(self):
    if not self.parent.has_model('null'): self.parent.add_codeml_model(null_codeml_model)    
    if not self.parent.has_model('discrete'): self.parent.add_codeml_model(discrete_codeml_model)    
    if no_run_mode: return
    m0=self.parent.get_model('null')
    m1=self.parent.get_model('discrete')
    self.p_value= m1.LRT(m0)
        
############################################################################################################################

class EvolPhyloTree(PhyloTree):
  """ Ete2 PhyloTree subclass to manage codon_alignment.tree 
  .parent is a fundamental attribute and points to the relevant codon_alignment object. nothing works if this is not defined.
  most useful methods:
  .sequence() -> get sequence (for non-leaf nodes, must be accepted from some model or computed)
  .get_name() -> returns the name of the node. can be .name for leaf, or a uniq name computed from leaves names for non-leaf nodes
  .Ka_with(n) .Ks_with(n) .KaKs_with(n)  -> compute dN dS dN/dS of this node compared to a refnode which is conceptually an ancestor (sites are counted on the refnode)
  .dKa() .dKs() .dKaKs()      -> compute dN dS dN/dS of this node compared to its parent (first direct ancestor)
  .Ka_lineage() .Ks_lineage() .KaKs_lineage()  -> compute dN dS dN/dS of this node considering with all unique changes observed in leaves under this node
  """

  def __init__(self, newick=None, parent=None):     ##########3 note has to add .parent=codon_alignment isntance in every node!!!
#    if parent: 
      #sp_naming_function=self.parent.sp_naming_function
#      PhyloTree.__init__(self, newick=newick, sp_naming_function=sp_naming_function)
#    else:       
    PhyloTree.__init__(self, newick=newick) 
    self.column_index=0
    self.parent=parent
    for n in self.traverse():      n.parent=parent

  def get_nodes(self, expr, no_root=False, no_leaves=False):
    """ returns the target nodes instances, given an expression as the one explained here: 
Nodes can be specified by their name if they are leaves, or using the first common ancestor principle if they are ancestral ( node1name+node2name ). Multiple nodes can be specified joining them with a comma (  node1name,node2name+node3name,node5name ). 
Some nodes can be specified with keywords:  root, all, all_leaves"""
    targets=[]
    if expr=='all':             targets=[i for i in self.traverse() if (not no_root or not i.is_root()) and (not no_leaves or not i.is_leaf() ) ]
    elif expr in ['all_leaves', 'leaves']:    targets=[i for i in self.traverse() if i.is_leaf() ]
    else:
      for it in expr.split(','):
        if it=='root':          targets.append(self)
        elif it.count('+')==0: 
          try:        targets.append(  self&it )   #searching node name provided in command line
          except:     raise ValueError, "ERROR can't find node with name: "+it
        else:
          try: 
            list_of_nodes_names=     it.split('+')
            targets.append( self.get_common_ancestor(*list_of_nodes_names)         )
          except:     raise ValueError, "ERROR can't find common ancestors of nodes with name: "+str(list_of_nodes_names)
    return targets    

  def show(self, **args):
    if not 'layout' in args: 
      if self.parent.layout is None: args['layout']=lambda x:x # null function
      else:                          args['layout']=self.parent.layout
    if not 'tree_style' in args: args['tree_style']=self.parent.tree_style
    PhyloTree.show(self, **args)
      
  def render(self, fileout, **args):
    if not 'layout' in args: 
      if self.parent.layout is None: args['layout']=lambda x:x # null function
      else:                          args['layout']=self.parent.layout
    if not 'tree_style' in args: args['tree_style']=self.parent.tree_style
    PhyloTree.render(self, fileout,  **args)

  def get_name(self):
    if self.is_leaf(): return self.name
    elif self.is_root(): return 'root'
    else:   
      first_n = sorted (  self.get_children()[0].get_leaves()  , key= lambda x:x.name )[0]
      second_n= sorted (  self.get_children()[1].get_leaves()  , key= lambda x:x.name )[0]
      l= sorted( [first_n.name, second_n.name] )
      return join(l, '+')

  def is_rooted(self):    return len(self.get_children())==2
  def is_binary(self):    
    """ checks all nodes apart from root"""
    for i in self.traverse():
      if not i.is_root() and not i.is_leaf(): 
        if len(i.get_children())!=2: 
          return False
    return True

  def sequence(self): 
    """ Returns the sequence for this node. If it is ancestral, raise an Exception unless the ancestral_sequences were computed and accepted"""
    if self.is_leaf(): return self.parent.seq_of( self.name  )
    else:             
      if not self.parent.ancestral_sequences: raise Exception, "node: "+self.get_name()+' -> sequence() ERROR ancestral_sequences are not defined! call method codon_alignment.set_ancestral_sequences on your favourite model'
      return self.parent.ancestral_sequences[self]
  
  def codons(self):
    """ Iterate through the codons of this node. Returns list of [  codon_index, codon   ] with codon_index 0-based """
    out=[]
    for i in range(  len(self.sequence())/3  ) :       out.append( i, self.sequence()[i*3:i*3+3] )
    return out

  def count_sites(self, silent=True, positions=None):
    """  See MMlib function count_sites for details
   NOTE:  nonsense (stop) mutations are ALWAYS differentiated from nonsyn mutations HERE. 
   the function returns [ nonSyn, Syn, NonSense,  CpG_nonSyn, CpG_syn, CpG_nonsense ]
   argument "positions" can be used to partition certain columns of the alignment. Argument should be list of tuples, each one with start-end (1-based, codon-based)
   """
    if not positions is None:
      subseq=self.subseq(positions)
      return count_sites(nogap(subseq), silent=silent, split_nonsense=True )
    else:
      if not self.parent.has_derived_data(self.count_sites, self.get_name()):
        data= count_sites(nogap(self.sequence()), silent=silent, split_nonsense=True )     
        self.parent.save_derived_data(  self.count_sites , self.get_name(),   data   )
      return self.parent.get_derived_data(self.count_sites, self.get_name()    )

  def subseq(self, positions):
    """ Argument should be list of tuples, each one with start-end (1-based, codon-based). Gap-containing sequence is returned    """
    seq=self.sequence()
    subseq=''
    for st, end in positions:  subseq+=  seq [ (st-1)*3: (end)*3 ]
    #print self.get_name(), str(positions).ljust(20), seq, subseq
    return subseq        
   
  def count_changes_with(self, node, silent=True, split_nonsense=True, positions=None):
    """   See MMlib function count_changes for details
 nonsense (stop) mutations are differentiated from nonsyn mutations. the function  returns [ nonSyn, Syn, NonSense,  CpG_nonSyn, CpG_syn, CpG_nonsense ] """
    if positions is None:      
      title= self.get_name()+':'+node.get_name()
      if not self.is_leaf() or not node.is_leaf():  title+=':AS='+self.parent.ancestral_sequences_source
      title+=':NS='+str(int(split_nonsense))
      if not self.parent.has_derived_data(self.count_changes_with,  title):
        seq2=node.sequence()      
        data= count_changes(  self.sequence(), seq2,  silent=silent, split_nonsense=split_nonsense )# [ nonSyn, Syn, NonSense,  CpG_nonSyn, CpG_syn, CpG_nonsense ]
        self.parent.save_derived_data( self.count_changes_with , title,   data )
      return self.parent.get_derived_data(self.count_changes_with, title  )
    else: 
      seq1 = self.subseq(positions=positions)
      seq2 = node.subseq(positions=positions)
      return count_changes(  seq1, seq2,  silent=silent, split_nonsense=split_nonsense )        ## some position limits specified; not saving in this case

 
  def CountS_with(self, node):
    """ returns the number of synonymous changes between the sequence at the self node and the sequence at the "node" argument """
    counted_changes= self.count_changes_with( node ) 
    return counted_changes[1]
  def CountA_with(self, node):
    """ returns the number of synonymous changes between the sequence at the self node and the sequence at the "node" argument """
    counted_changes= self.count_changes_with( node ) 
    return counted_changes[0]

   #### Methods to compare a node with another
  def Ks_with(self, refnode, positions=None):
    """ conceptually, refnode is a node above self -- but it works with any node. we compute the proportion of fixed syn mutation over the possible ones (computed on refnode)    """
    if not positions is None:
      counted_sites_ref = refnode.count_sites(positions=positions)
      counted_changes   = self.count_changes_with(refnode, positions=positions)
      return float(counted_changes [1])/counted_sites_ref[1]
    else:
      title=self.get_name()+':'+refnode.get_name()
      if not self.is_leaf() or not refnode.is_leaf():  title+=':AS='+self.parent.ancestral_sequences_source
      if not self.parent.has_derived_data(self.Ks_with,  title ):
        counted_sites_ref = refnode.count_sites()
        counted_changes   = self.count_changes_with(refnode)
        data= float(counted_changes [1])/counted_sites_ref[1]
        self.parent.save_derived_data( self.Ks_with , title,   data   )
      return self.parent.get_derived_data(self.Ks_with , title  )    

  def Ka_with(self, refnode, positions=None):
    """proportion of fixed nonsyn mutation over the possible ones . NOTE: nonsense mutations are not counted! """
    if not positions is None:
      counted_sites_ref = refnode.count_sites(positions=positions)
      counted_changes   = self.count_changes_with(refnode, positions=positions)
      return float(counted_changes [0]+counted_changes [2])/(counted_sites_ref[0]+counted_sites_ref[2])
    else:
      title=self.get_name()+':'+refnode.get_name()
      if not self.is_leaf() or not refnode.is_leaf():  title+=':AS='+self.parent.ancestral_sequences_source
      if not self.parent.has_derived_data(self.Ka_with,  title ):
        counted_sites_ref = refnode.count_sites()
        counted_changes   = self.count_changes_with(refnode)
        data= float(counted_changes [0]+counted_changes [2])/(counted_sites_ref[0]+counted_sites_ref[2])
        self.parent.save_derived_data( self.Ka_with , title,   data   )
      return self.parent.get_derived_data(self.Ka_with , title  )
    
  def KaKs_with(self, refnode, positions=None):
    """ compute w - omega - dN/dS - KaPS - KaKs, however you want to call it -- against a reference node. Returns 999.0 if infite (Ks=0)"""
    if not positions is None:
      counted_sites_ref = refnode.count_sites(positions=positions)
      counted_changes   = self.count_changes_with(refnode, positions=positions)
      Ka= float(counted_changes [0]+counted_changes [2])/(counted_sites_ref[0]+counted_sites_ref[2])
      Ks= float(counted_changes [1])/counted_sites_ref[1]
      if Ka==0.0:    data= 0.0
      elif Ks==0.0:  data= 999.0
      else:          data= Ka/Ks
      #print data,
      return data
    else:
      title=self.get_name()+':'+refnode.get_name()
      if not self.is_leaf() or not refnode.is_leaf():  title+=':AS='+self.parent.ancestral_sequences_source
      if not self.parent.has_derived_data(self.KaKs_with,  title ):
        counted_sites_ref = refnode.count_sites()
        counted_changes   = self.count_changes_with(refnode)
        Ka= float(counted_changes [0]+counted_changes [2])/(counted_sites_ref[0]+counted_sites_ref[2])
        Ks= float(counted_changes [1])/counted_sites_ref[1]
        if Ka==0.0:    data= 0.0
        elif Ks==0.0:  data= 999.0
        else:          data= Ka/Ks
        self.parent.save_derived_data( self.KaKs_with , title,   data   )
      return self.parent.get_derived_data(self.KaKs_with , title  )

  #### Methods to compare a node with its parent, that is to say, compute a sortof first derivative indexes
  def dKs(self):
    """ compute Ks with the first ancestor of this node. equivalent to     self.Ks_with(self.up)    """
    if not self.up: raise Exception, "ERROR can't compute dKs for the root node!"
    return self.Ks_with(self.up)

  def dKa(self):
    """ compute Ka with the first ancestor of this node. equivalent to     self.Ka_with(self.up)    """
    if not self.up: raise Exception, "ERROR can't compute dKa for the root node!"
    return self.Ka_with(self.up)

  def dKaKs(self):
    """ compute KaPS with the first ancestor of this node. equivalent to     self.Ks_with(self.up)    """
    if not self.up: raise Exception, "ERROR can't compute dKaKs for the root node!"
    return self.KaKs_with(self.up)
    
  def dCountS(self):
    """ returns the number of synonymous changes with the first ancestor of this node. """
    if not self.up: raise Exception, "ERROR can't compute dCountS for the root node!"    
    counted= self.count_changes_with( self.up ) #     nonSyn, Syn, 
    return counted[1] 

  def dCountA(self):
    """ returns the number of synonymous changes with the first ancestor of this node. """
    if not self.up: raise Exception, "ERROR can't compute dCountA for the root node!"    
    counted= self.count_changes_with( self.up ) #     nonSyn, Syn, 
    return counted[0] 
  dCountN=dCountA

  ### methods to compare this node with all its leafs
  def Ks_lineage(self, non_leaves=False):
    """ compute Ks comparing the sequence in self with all its leaves, counting each change only once.  non_leaves are also considered if non_leaves==True"""
    if self.is_leaf(): raise Exception, "ERROR can't compute Ks_lineage for leaf node!"
    title=self.get_name()+':NL='+str(int(non_leaves))
    if non_leaves: title+=':AS='+self.parent.ancestral_sequences_source
    if not self.parent.has_derived_data(self.Ks_lineage , title  ):
      counted_sites = self.count_sites()    
      other_cds=[ n.sequence() for n in self.traverse() if n.is_leaf() or non_leaves ]    
      counted_uniq_changes=count_unique_changes( self.sequence(), other_cds, split_nonsense=True )
      data= float( counted_uniq_changes[1] )/counted_sites[1]
      self.parent.save_derived_data( self.Ks_lineage , title, data)
    return self.parent.get_derived_data(self.Ks_lineage , title)

  def Ka_lineage(self, non_leaves=False):
    """ compute Ka comparing the sequence in self with all its leaves, counting each change only once.  non_leaves are also considered if non_leaves==True"""
    if self.is_leaf(): raise Exception, "ERROR can't compute Ka_lineage for leaf node!"
    title=self.get_name()+':NL='+str(int(non_leaves))
    if non_leaves: title+=':AS='+self.parent.ancestral_sequences_source
    if not self.parent.has_derived_data(self.Ka_lineage , title  ):
      counted_sites = self.count_sites()    
      other_cds=[ n.sequence() for n in self.traverse() if n.is_leaf() or non_leaves ]    
      counted_uniq_changes=count_unique_changes( self.sequence(), other_cds, split_nonsense=True )
      data= float( counted_uniq_changes[0]+counted_uniq_changes[2] )/(counted_sites[0]+counted_sites[2])
      self.parent.save_derived_data( self.Ka_lineage , title,   data   )
    return self.parent.get_derived_data(self.Ka_lineage , title)

  def KaKs_lineage(self, non_leaves=False):
    """ compute KaKs comparing the sequence in self with all its leaves, counting each change only once.  non_leaves are also considered if non_leaves==True
        Returns 999.0 if infinite (Ks=0)    """
    if self.is_leaf(): raise Exception, "ERROR can't compute KaKs_lineage for leaf node!"
    title=self.get_name()+':NL='+str(int(non_leaves))
    if non_leaves: title+=':AS='+self.parent.ancestral_sequences_source
    if not self.parent.has_derived_data(self.KaKs_lineage , title):
      counted_sites = self.count_sites()    
      other_cds=[ n.sequence() for n in self.traverse() if n.is_leaf() or non_leaves ]    
      counted_uniq_changes=count_unique_changes( self.sequence(), other_cds, split_nonsense=True )
      Ka= float( counted_uniq_changes[0]+counted_uniq_changes[2] )/(counted_sites[0]+counted_sites[2])
      Ks= float( counted_uniq_changes[1] )/counted_sites[1]
      if Ka==0.0: data=0.0  ## --> useful only if Ks is 0.0 as well
      elif Ks==0.0: data=999.0
      else: data= Ka/Ks    
      self.parent.save_derived_data( self.KaKs_lineage , title, data)
    return self.parent.get_derived_data(self.KaKs_lineage , title)

  def sum_dKs(self):
    """ compute the sum of all dKs of children until leaves of this non-leaf node"""
    if self.is_leaf(): raise Exception, "ERROR can't compute sum_dKs for leaf node!"
    title=self.get_name()+':AS='+self.parent.ancestral_sequences_source
    if not self.parent.has_derived_data(self.sum_dKs , title):    
      data=sum( [node.dKs()    for node in self.traverse() if node != self ] )
      self.parent.save_derived_data( self.sum_dKs , title, data)      
    return self.parent.get_derived_data(self.sum_dKs , title)      

  def sum_dKa(self):
    """ compute the sum of all dKa of children until leaves of this non-leaf node"""
    if self.is_leaf(): raise Exception, "ERROR can't compute sum_dKa for leaf node!"
    title=self.get_name()+':AS='+self.parent.ancestral_sequences_source
    if not self.parent.has_derived_data(self.sum_dKa , title):    
      data=sum( [node.dKa()    for node in self.traverse() if node != self ] )
      self.parent.save_derived_data( self.sum_dKa , title, data)      
    return self.parent.get_derived_data(self.sum_dKa, title)      

  def count_unique_changes(self, non_leaves=False):
    """ compute the uniq changes of this sequence comparing with the leaves below this node
    returns: [ nonSyn, Syn, NonSense,  CpG_nonSyn, CpG_syn, CpG_nonsense ]
    see MMlib.count_unique_changes
    """
    if self.is_leaf(): raise Exception, "ERROR can't compute uniq_changes for leaf node... nothing to compare with!"
    title=self.get_name()+':NL='+str(int(non_leaves))
    if non_leaves: title+=':AS='+self.parent.ancestral_sequences_source    
    if not self.parent.has_derived_data(self.count_unique_changes , title):
      other_cds=[ n.sequence() for n in self.traverse() if n.is_leaf() or non_leaves ]    
      data= count_unique_changes( self.sequence(), other_cds, split_nonsense=True  )
      self.parent.save_derived_data( self.count_unique_changes , title, data)
    return self.parent.get_derived_data(self.count_unique_changes , title)

  #################### decorating function: to add faces to nodes

  def decorate_with_text(self, text, fsize=8, color='#000000', pos='branch-top', add_column=True, margin_left=2, margin_right=2):
    """ Add a face on the node with specified attributes"""
    t_face=faces.TextFace( str(text) , fgcolor=color, fsize=fsize)         
    t_face.margin_left=margin_left;   t_face.margin_right=margin_right
    self.add_face(t_face, column = self.column_index, position=pos)
    if add_column: self.column_index+=1

  """def decorate_with_count_sites(self, text, fsize=8, color='#000000', pos='branch-top', add_column=True, margin_left=2, margin_right=2):
    "" Add a face on the node with specified attributes""

    t_face=faces.TextFace( str(text) , fgcolor=color, fsize=fsize)         
    t_face.margin_left=margin_left;   t_face.margin_right=margin_right
    self.add_face(t_face, column = self.column_index, position=pos)
    if add_column: self.column_index+=1"""

  def decorate_with_w_change(self, descendants=False, circle=1, text_w=1, text_p=0, text_r=0,  pos_w='branch-top', pos_p='branch-bottom', pos_r='branch-bottom',  kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme, fsize=8):
    """Add a face representing the kaks change analysis. A model named 'w_change_?.'+self.get_name() must be defined in self.models (? is D if descendants==True, or S if not). If p_value is checked, a self.tests[same_name] must also be present. A circle and some text may be drawn, depending on the arguments:
  circle: 0 -> don't draw     
          1 -> draw if significant
          2->  draw as sphere if not significant, filled if significant      
          3 -> draw as sphere always
text_w, text_r, text_p:            #text_w is the numeric value of omega in the target node, text_r is the numeric value of omega in the rest of the tree,   text_p is about p_value
          0 -> don't put it             
          1 -> put it if significant                                          
          2 -> put it always
Position of the texts depends on variables pos_w, pos_r, pos_p
The color depends on whether the kaks for the branch is greater or minor than the kaks for the rest of tree, following the kaks_color_scheme (  if kaks branch is > than the average, and < than maximum value called neutral, but even < than the minimum value called neutral, it is shown as neutral since that is the direction of the change)
  kaks_neutral_boundaries defines what is considered neutral (just for color choosing purposes)"""  
  
    flag_descendants={False:'S', True:'D'}[descendants]
    model_name='w_change_'+flag_descendants+'.'+self.get_name()
    if not self.parent.has_model(model_name):     raise Exception, "decorate_with_w_change ERROR model not found: "+model_name
    model=self.parent.get_model(model_name)
    if (circle in [1, 2] or text_w==1 or text_p !=0):
      if not self.parent.has_test(model_name):  raise Exception, "decorate_with_w_change ERROR test not found: "+model_name
      test=self.parent.get_test(model_name)
    omega=     model.omega_for(self)
    omega_rest=model.w_rest_of_tree()
    #null_omega=self.parent.null_model_omega()
    color=get_color_kaks_change(omega, omega_rest, kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme)    
    ### text_w for omega in the rest of tree
    if text_r==2 or (text_r==1 and is_significant( test.p_value )):
      color_rest=get_color_free_kaks(omega_rest, kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme)    
      self.decorate_with_text(   str(round(omega_rest, n_digits)), color=color_rest, fsize=fsize, pos=pos_r)
    ### text_w for omega in this node
    if text_w==2 or (text_w==1 and is_significant( test.p_value )):
      self.decorate_with_text(   str(round(omega, n_digits)), color=color, fsize=fsize, pos=pos_w)
    ### text_p for omega in this node
    if text_p==2 or (text_p==1 and is_significant( test.p_value )):
      self.decorate_with_text(   str(round(test.p_value, n_digits)), color='#000000', fsize=fsize, pos=pos_p)
    ### circle representation of w_change
    if circle>2 or (circle==1 and is_significant( test.p_value )):
      non_significant_opacity=0.2;     min_significant_opacity=0.5;     max_opacity=1.0
      max_circle_size=20; factor_size=4
      if omega > omega_rest:     size= factor_size *    float(omega)/omega_rest ## >=1
      else:                      size= factor_size *  omega_rest/float(omega) ## >=1
      if size >= max_circle_size: size=max_circle_size
      if circle==3: 
        opacity=max_opacity
        style="sphere"
      elif not is_significant(test.p_value):   
        opacity=non_significant_opacity
        style="sphere"
      else:                     
        opacity=       (test.p_value) / (  get_MMlib_var('opt')['alpha'] )  * ( max_opacity - min_significant_opacity ) + min_significant_opacity
        style="circle"
      diff_kaks_face=faces.CircleFace(size, color, style=style)
      diff_kaks_face.opacity=opacity
      self.add_face(diff_kaks_face, column = self.column_index, position="branch-right")    
      self.column_index+=1

  def decorate_with_w(self, w, circle=1, text_w=0, pos_w="branch-top", kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme, fsize=8,  style='sphere', color=0):
    """ This function is adds face(s) representing a W -- omega -- KaKs -- KaKs. A colored circle and text can be added depending on options:
    circle = 0 --> dont put it          circle = 1   --> put circle   ## its style (sphere|circle) is defined by style argument of this function 
    text_w = 0 --> dont put it          text_w = 1   --> put text 
    position of text depends on pos_w argument
    The color of the text_w and circle is determined by its class, defined by the kaks_neutral_boundaries ( if <, puryfing, < <, neutral,  > positive) and by the kaks_color_scheme """

    min_circle_size=4;  max_circle_size=20;  opacity=1.0;  kaks_positive_saturation=5
    if not color:
      color=kaks_to_color(w) 
    #color=get_color_free_kaks(w,   kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme)
    if text_w: ### adding text
      self.decorate_with_text(   str(round(w, n_digits)), color=color, fsize=fsize, pos=pos_w)      
    if circle: ### adding circle
      if w< kaks_neutral_boundaries[0]:          
        #puryfing
        size= ((kaks_neutral_boundaries[0]-w)*(max_circle_size - min_circle_size)/kaks_neutral_boundaries[0] )+min_circle_size
      elif w<= kaks_neutral_boundaries[1]:          
        #netrual
        size= ( (w - kaks_neutral_boundaries[0] )    *(max_circle_size - min_circle_size)/( kaks_neutral_boundaries[1]- kaks_neutral_boundaries[0]) )+min_circle_size 
      else:    
        #positive
        size= ( (w - kaks_neutral_boundaries[1] )    *(max_circle_size - min_circle_size)/( kaks_positive_saturation- kaks_neutral_boundaries[1]) )+min_circle_size 
        if size > max_circle_size: size=max_circle_size

      size=10 ############# debug
      free_kaks_face=faces.CircleFace(size, color, style=style)
      free_kaks_face.opacity=opacity
      self.add_face(free_kaks_face, column = self.column_index, position="branch-right")    
      self.column_index+=1

  def decorate_with_free_w(self, circle=1, text_w=0, pos_w="branch-top",  kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme, fsize=8,  style='sphere'):
    """ Adds a face with the W value for this node computed with a free_w codeml model. It uses the function decorate_with_w, here its help:\n""" ##see below 
    if not self.parent.has_model('free_w'):     raise Exception, "decorate_with_free_w ERROR model not found: free_w"
    model=self.parent.get_model('free_w')
    omega=     model.omega_for(self)
    self.decorate_with_w(omega, circle=circle, text_w=text_w,  pos_w=pos_w,  kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=kaks_color_scheme, fsize=fsize,  style=style)
  decorate_with_free_w.__doc__+=function_help(decorate_with_w)

  #### Ks_lineage, Ka_lineage, KaKs_lineage
  def decorate_with_KaKs_lineage(self, circle=0, text_w=1,  pos_w="branch-bottom", kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme, fsize=8,  style='sphere'):
    """ Adds a face with W value computed comparing the sequence in this node with all the leaves below this node.  It cannot be applied to leaves. It uses the function decorate_with_w, here its help;\n"""  
    self.decorate_with_w( self.KaKs_lineage() , circle=circle, text_w=text_w, pos_w=pos_w, kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=kaks_color_scheme, fsize=fsize,  style=style)
  decorate_with_KaKs_lineage.__doc__+=function_help(decorate_with_w)    

  def decorate_with_Ks_lineage(self, fsize=8, color='#007070', pos='branch-top'):
    """ Adds a face with Ks_lineage computed comparing the sequence in this node with all the leaves below this node.  It cannot be applied to leaves. 
Its text is added to the node. It uses the following function: \n"""  
    self.decorate_with_text(   str(round(self.Ks_lineage(), n_digits)), color=color, fsize=fsize, pos=pos)
  decorate_with_Ks_lineage.__doc__+=function_help(decorate_with_text)      
      
  def decorate_with_Ka_lineage(self, fsize=8, color='#CC33FF', pos='branch-top'):
    """ Adds a face with Ka_lineage computed comparing the sequence in this node with all the leaves below this node.  It cannot be applied to leaves. 
Its text is added to the node. It uses the following function:\n"""  
    self.decorate_with_text(   str(round(self.Ka_lineage(), n_digits)), color=color, fsize=fsize, pos=pos)
  decorate_with_Ka_lineage.__doc__+=function_help(decorate_with_text)      
      
  #### dKs, dKa, dKaKs
  def decorate_with_dKaKs(self, circle=0, text_w=1,  pos="branch-top", kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=basicKaKsColorScheme, fsize=8,  style='sphere', color=0):
    """ Adds a face with W value computed comparing the sequence in this node with the one just above it.  It cannot be applied to root. It uses the function decorate_with_w, here its help;\n"""  
    self.decorate_with_w( self.dKaKs() , circle=circle, text_w=text_w,  pos_w=pos, kaks_neutral_boundaries=kaks_neutral_boundaries, kaks_color_scheme=kaks_color_scheme, fsize=fsize,  style=style, color=color)
  decorate_with_dKaKs.__doc__+=function_help(decorate_with_w)    

  def decorate_with_dKs(self, fsize=8, color='#007070', pos='branch-bottom'):
    """ Adds a face with dKs computed comparing the sequence in this node with the one just above it. It cannot be applied to root. 
Its text is added to the node. It uses the following function: \n"""  
    self.decorate_with_text(   str(round(self.dKs(), n_digits)), color=color, fsize=fsize, pos=pos)
  decorate_with_dKs.__doc__+=function_help(decorate_with_text)      

  def decorate_with_dKa(self, fsize=8, color='#CC33FF', pos='branch-bottom'):
    """ Adds a face with Ka_lineage computed comparing the sequence in this node with the one just above it. It cannot be applied to root. 
Its text is added to the node. It uses the following function:\n"""  
    self.decorate_with_text(   str(round(self.dKa(), n_digits)), color=color, fsize=fsize, pos=pos)
  decorate_with_dKa.__doc__+=function_help(decorate_with_text)    

  def decorate_with_sum_dKs(self, fsize=8, color='#009999', pos='branch-bottom'):
    """ Adds a face with the sum of dKs of all descendants of this node. It cannot be applied to leaves. 
Its text is added to the node. It uses the following function: \n"""  
    self.decorate_with_text(   str(round(self.sum_dKs(), n_digits)), color=color, fsize=fsize, pos=pos)
  decorate_with_sum_dKs.__doc__+=function_help(decorate_with_text)      

  def decorate_with_sum_dKa(self, fsize=8, color='#DD44FF', pos='branch-bottom'):
    """ Adds a face with the sum of dKa of all descendants of this node. It cannot be applied to leaves. 
Its text is added to the node. It uses the following function: \n"""  
    self.decorate_with_text(   str(round(self.sum_dKa(), n_digits)), color=color, fsize=fsize, pos=pos)
  decorate_with_sum_dKa.__doc__+=function_help(decorate_with_text)      
  
  def decorate_with_dCountS(self, fsize=8, color='#2222DD', pos='branch-bottom'):
    """ Add a face with the count of syn changes of this node compared to its direct ancestor. It uses the following function: """
    self.decorate_with_text(   str(self.dCountS()), color=color, fsize=fsize, pos=pos)    
  decorate_with_dCountS.__doc__+=function_help(decorate_with_text)          
        
  def decorate_with_dCountA(self, fsize=8, color='#DD2222', pos='branch-bottom'):
    """ Add a face with the count of nonsyn changes of this node compared to its direct ancestor. It uses the following function: """
    self.decorate_with_text(   str(self.dCountA()), color=color, fsize=fsize, pos=pos)    
  decorate_with_dCountA.__doc__+=function_help(decorate_with_text)          
        

def add_colored_table(tree, table_hash, column, value_limits=None, colors=["#BB0000", "#00BB00"], size=[50, 50], ordered_fields=None, special_values={}):
  """Given a ete2 tree, it uses the function tree.add_face to add a table made of colored rectangles, with colors representing numeric values (or even label values -- see special_values below). 
  The faces are added with column indexes starting from the argument column. For each value, the corresponding color is computed, by linear interpolation of the minimum - maximum value, given the min and max colors. If argument value_limits is not specified, the min and max values are taken from the data. If colors is not specified, they are red (min) and green (max).
  You can also use three colors for coloring, with the third item in the provided list which is used as a midpoint of the first two
  Size determines the size of the rectangles, in [width, height] -- default: [50, 50].
  The variable ordered_fields, if provided, is a list of the column titles in the order they will appear left ot right (otherwise, alphabetical order is used).
  the hash special_values can be provided to draw a specific color for a certain data value, without any interpolation or anything:   {value -> color }
  The input data is table_hash, which is a hash of hashes:  {  column_name  ->  { node or node_name ->  value } }

  probably, after calling this you want to add the column titles to your tree_style object, with something like this:
  for title in ordered_fields:
      title_face=TextFace(title, fsize=12)
      title_face.hz_align = 1
      tree_style.aligned_header.add_face(title_face, column=len(families_order))  
  """
  if ordered_fields is None: 
    ordered_fields= sorted(table_hash.keys())
  if value_limits is None: 
    min_v=sys.maxint
    max_v=-sys.maxint
    for k in table_hash.keys():
      for n in table_hash[k]:
        value=table_hash[k][n]
        if type(value) in [int, long, float]:
          if value < min_v: min_v= value
          if value > max_v: max_v= value        

  else:     min_v, max_v =     value_limits
  for node in tree.traverse():
    if node.is_leaf():
      for k_index, k in enumerate(ordered_fields):
        #parsing table columns
        if table_hash[k].has_key(node):          value=table_hash[k][node]
        elif table_hash[k].has_key(node.name):   value=table_hash[k][node.name]
        else:                                    value=None

        if not value is None:
          if not special_values.has_key(value):
            if not type(value) in [int, long, float]: raise Exception, "ERROR can't interpolate value of type "+str(type(value)) + ' '+str(value)
            normalized= (value-min_v)/float(max_v - min_v)
            if len(colors)==2:
              color=color_scale(normalized,  colors[0], colors[1])
            elif len(colors)==3:
              color=color_scale_midpoint(normalized,  colors[0], colors[2], colors[1])
            else: raise Exception, "wrong length for colors: it can be just 2 or 3. it is: "+str(colors)
          else:
            color=special_values[value]
          #print node, normalized, colors[0], colors[1], color
          f= coloredRectangle( color, w=size[0], h=size[1], line_size=0, line_color="#000000")
          node.add_face(f, column=column+ k_index,  position='aligned' )
