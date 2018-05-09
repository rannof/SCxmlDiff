#!/usr/bin/env python
# by Ran Novitsky Nof (ran.nof@gmail.com), 2015 @ BSL

# Compare two Seiscomp3 xml files
#
#
#
# ***********************************************************************************
# *    Copyright (C) by Ran Novitsky Nof                                            *
# *                                                                                 *
# *    SCxmlDiff.py is free software: you can redistribute it and/or modify         *
# *    it under the terms of the GNU Lesser General Public License as published by  *
# *    the Free Software Foundation, either version 3 of the License, or            *
# *    (at your option) any later version.                                          *
# *                                                                                 *
# *    This program is distributed in the hope that it will be useful,              *
# *    but WITHOUT ANY WARRANTY; without even the implied warranty of               *
# *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
# *    GNU Lesser General Public License for more details.                          *
# *                                                                                 *
# *    You should have received a copy of the GNU Lesser General Public License     *
# *    along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
# ***********************************************************************************


import argparse,sys,os,re,datetime,tempfile
import matplotlib as mpl
import datetime
from PyQt4.Qt import QLabel
mpl.use('QT4Agg')
from matplotlib.backend_bases import NavigationToolbar2, Event
from matplotlib.backends.backend_qt4agg import(
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from matplotlib import cm
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import uic
from pylab import *
from xml.dom import minidom
import pyproj

parser = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description='''SCxmlDiff - Compare two Seiscomp xml files''',
         epilog='''Created by Ran Novitsky Nof (ran.nof@gmail.com), 2015 @ BSL''')
parser.add_argument('-i',metavar='file',nargs='+',help='input E2 result files (xml or log)',type=str)
parser.add_argument('-r',metavar='file',nargs='+',help='input reference files (xml or csv)\nxml - Seiscomp3 xml file\ncsv - comma separated value (GII DB forat [epiid,ml,typ,lat,lon,abs_ot_d,abs_ot_t]. ',type=str)
parser.add_argument('-b',metavar='bounds',nargs=4,help='Region bounding box (west east south north)',type=float,default=(-180.0,180.0,-90,90))
#variables
AGENCYID='GII'
AUTHOR='E2'
FONTSIZE=8
GRIDON=False # grid on or off [True | False]
VERBOSE=False # printout message
geo = pyproj.Geod(ellps='WGS84')


EvaluationStatus = {
   0: "preliminary",
   1: "confirmed",
   2: "reviewed",
   3: "final",
   4: "rejected",
   5: "reported"
  }

HEADERLABELS = ['EvtID','OT','M','Lat','Lon','Status','refID','Evaluation']

def concatenateSC3XML(xmls):
  import seiscomp3.IO as scio
  import seiscomp3.DataModel as scdatamodel
  import seiscomp3.Core as sccore
  'concatenates Seiscomp3 xml files. Minimal error checks.'
  ar = scio.XMLArchive()
  xml = None
  for x in xmls:
    if not ar.open(x): continue# validate a valid SC3 xml file
    ar.close()
    if not xml: # get the first file and append all data into it
      xml = minidom.parse(x) # parse xml file
      eparams = xml.getElementsByTagName('EventParameters')[0] # get the eventparameter element
    else: # since we have a main xml, we add to it
      [eparams.appendChild(c) for c in minidom.parse(x).getElementsByTagName('EventParameters')[0].childNodes if not c.nodeName=='#text'] # append to main eventparameters
  tmp = tempfile.NamedTemporaryFile('w',dir='.') # create a temp file
  f = open(tmp.name,'w') # open a new handler so we can close it
  f.write(xml.toxml()) # write to temp file
  f.close() # close the file so buffer is flushed
  if not ar.open(tmp.name): return # try to open the temp file as scdatamodel xml
  eparams = scdatamodel.EventParameters.Cast(ar.readObject()) # get event parameters
  tmp.close() # remove temp file
  return eparams # scdatamodel eventparameters

def concatenateCSV(CSVs):
  giidbtype = dtype([('epiid','|S50'),('ml',float),('typ',int),('lat',float),('lon',float),('abs_ot_d','|S50'),('abs_ot_t','|S50')])
  rettype = dtype([('EID','|S50'),('ot','|S50'),('lat',float),('lon',float),('mag',float)])
  ret = array([],dtype=rettype)
  for f in CSVs:
    try:
      fdata = loadtxt(f,delimiter=',',dtype=giidbtype,skiprows=1)
    except ValueError:
      raise ValueError('%s is not a valid scv file'%(f))
    ret = append(ret,array([(v[0],'T'.join([v[-2],v[-1]])+'Z',v[3],v[4],v[1]) for v in fdata],dtype=ret.dtype))
  return ret

def eparams2seisEvents(eparams,starttime=None,endtime=None,minmag=-10.0,maxmag=10.0,bbox=None,preferred='p',includeunassoc=False):
  '''
  Converts a seiscomp3 eventParameters object to seisEvent object.
  starttime and endtime should be a datatime value for begin/end of return values
  minmag,maxmag - magnitude range
  bbox should be a lllon lllat urlon urlat values of locations bounding box
  preferred should be 'p': preferred origin ; 'f': first origin ; 'l': last origin
  includeunassoc - should we output also unassociated origins (default - False)
  will return a dictionary of SeisEvent objects, where id is event id or a yymmddhhmmss.s value
  '''
  if not preferred in ['p','f','l']:
    if VERBOSE: print >> sys.stderr,"preferred origin should be 'p': preferred origin ; 'f': first origin ; 'l': last origin"
    return {}
  # get events origins
  events = [eparams.event(i) for i in xrange(eparams.eventCount())]
  if preferred=='p': # preferred origin
    origins = [eparams.findOrigin(event.preferredOriginID()) for event in events]
  elif preferred=='f': # first origin
    origins = [eparams.findOrigin(event.originReference(0).originID()) for event in events]
  elif preferred=='l': # last origin
    origins = [eparams.findOrigin(event.originReference(event.originReferenceCount()-1).originID()) for event in events]
  D = {}
  for origin,event in zip(origins,events):
    try:
      stat = EvaluationStatus[origin.evaluationStatus()]
    except:
      stat = EvaluationStatus[0]
    D[event.publicID()]=SeisEvent(event.publicID(),origin.time().value().iso(),
                                  origin.latitude().value(),origin.longitude().value(),
                                  origin.magnitude(0).magnitude().value(),evalstat=stat)
  if includeunassoc: # unassociated origins
    referrencedOriginIds = sum([[event.originReference(i).originID() for i in xrange(event.originReferenceCount())] for event in events])
    origins = [eparams.origin(i) for i in xrange(eparams.originCount()) if not eparams.origin(i).publicID() in referrencedOriginIds]
    [D.update({origin.time().value().iso().iso().replace('-','').replace('T','').replace(':','').split('.')[0]:
               SeisEvent(event.publicID(),origin.time().value().iso(),origin.latitude().value(),origin.longitude().value(),origin.magnitude(0).magnitude().value())}) for origin in origins]
  return D

def db2seisEvents(DBDATA,starttime=None,endtime=None,minmag=-10.0,maxmag=10.0,bbox=None,preferred='p',includeunassoc=False):
  '''
  Converts a db array to seisEvent object.
  db dtype should be: dtype([('EID','|S50'),('ot','|S50'),('lat',float),('lon',float),('mag',float)]). See concatenateCSV for example
  starttime and endtime should be a datatime value for begin/end of return values
  minmag,maxmag - magnitude range
  bbox should be a lllon lllat urlon urlat values of locations bounding box
  preferred is ignored
  includeunassoc - should we output also unassociated origins (default - False)
  will return a dictionary of SeisEvent objects
  '''
  #convert DBDATA to SeisEvent dictionary
  D = dict([(r[0],SeisEvent(*r,evalstat='reported')) for r in DBDATA])
  return D

class SeisEvent(object):
  'event object'
  def __init__(self,EID,ot,lat,lon,mag,magType='M',evalstat='preliminary'):
    self.id = EID
    self.lat = lat
    self.lon = lon
    self.mag = mag
    self.magType = magType
    self.SCevaluationStatus = evalstat # evaluation status.
    self.status = '' # color code. False event: 'r' ; Missed event: 'FFA500'; True event: 'g'
    self.ot = datetime.datetime.strptime(ot,'%Y-%m-%dT%H:%M:%S.%fZ')
    self.line = None
    self.refid = ''
  def __sub__(self,other):
    return geo.inv(self.lon,self.lat,other.lon,other.lat)
  def __repr__(self):
    return "%s %.6f %.6f %.2f %s"%(self.id,self.lat,self.lon,self.mag,self.ot)

class ConnectedLineEdit(QLineEdit):
  'QLineEdit class for connected limits (like max/min magnitude)'
  def __init__(self,val,minval,maxval,deci=2):
    QLineEdit.__init__(self)
    self.val = float(val)
    self.setText(val)
    self.minval=minval
    self.maxval=maxval
    self.validator1 = QDoubleValidator()
    self.validator1.setRange(minval,maxval,deci)
    #self.setValidator(validator)
    self.editingFinished.connect(self.validate)
    self.setToolTip('%0.2lf <= Value >= %0.2lf'%(minval,maxval))
    self.connectedLineEdit = None
    self.connectedLineEditVal = None
    self.connectedVal = None
    self.limDirection = None
    self.connected = False
  def setConnected(self,widget,param):
    assert param in ['top','bottom']
    self.connectedLineEdit = widget
    self.connectedLineEditParam = param # should be top or bottom
    self.connected=True
  def validate(self):
    if not self.validator1.validate(self.text(),2)==(2,2): # reset value if not valid
      self.setText(str(self.val))
      return 0
    if self.connected: # update connected widget limit or reset
      param = self.connectedLineEditParam
      widget = self.connectedLineEdit
      if param=='top':
        o = '<='
      else:
        o = '>='
      if eval(str(widget.text())+o+str(self.text())):
        self.val = float(self.text())
        eval('widget.validator1.set'+param.capitalize()+'('+str(self.val)+')')
        minval,maxval = (widget.validator1.bottom(),widget.validator1.top())
        widget.setToolTip("%0.2lf <= Value >= %0.2lf"%(minval,maxval))
        self.emit(SIGNAL('ok'))
        return 1
      else:
        self.setText(str(self.val))
        return 0
    self.val = float(self.text())
    self.emit(SIGNAL('ok'))
    return 1

# map area dialog
class bboxForm(QDialog):
  def __init__(self,parent=None):
    QDialog.__init__(self,parent=parent)
    self.setWindowTitle('SCxmlDiff - Set area rectangle')
    vbox = QVBoxLayout(self)
    w = QWidget()
    grid = QGridLayout(w)
    vbox.addWidget(w)
    westlabel = QLabel('West')
    W = ConnectedLineEdit('-180',-180,180)
    grid.addWidget(westlabel,2,1)
    grid.addWidget(W,2,2)
    eastlabel = QLabel('East')
    E = ConnectedLineEdit('180',-180,180)
    grid.addWidget(eastlabel,2,5)
    grid.addWidget(E,2,6)
    northlabel = QLabel('North')
    N = ConnectedLineEdit('90',-90,90)
    grid.addWidget(northlabel,1,3)
    grid.addWidget(N,1,4)
    southlabel = QLabel('South')
    S = ConnectedLineEdit('-90',-90,90)
    grid.addWidget(southlabel,3,3)
    grid.addWidget(S,3,4)
    self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel,Qt.Horizontal, self)
    vbox.addWidget(self.buttons)
    self.buttons.accepted.connect(self.accept)
    self.buttons.rejected.connect(self.reject)
    self.W = W
    self.E = E
    self.S = S
    self.N = N
  def setLims(self,w,e,s,n):
    self.W.setText(str(w))
    self.E.setText(str(e))
    self.S.setText(str(s))
    self.N.setText(str(n))
  def getLims(self):
    w = self.W.text().toDouble()[0]
    e = self.E.text().toDouble()[0]
    s = self.S.text().toDouble()[0]
    n = self.N.text().toDouble()[0]
    return w,e,s,n
  def validate(self):
    w,e,s,n = self.getLims()
    if w>=e and s>=n:
      QMessageBox.warning(self,'SCxmlDiff - Error','West value should be lower than East value.\nSouth value should be lower than North value.')
      return 0
    elif w>=e:
      QMessageBox.warning(self,'SCxmlDiff - Error','West value should be lower than East value.')
      return 0
    elif s>=n:
      QMessageBox.warning(self,'SCxmlDiff - Error','South value should be lower than North value.')
      return 0
    else:
      return 1

class AppForm(QMainWindow):
  'QT form for application'
  def __init__(self, splash,args,parent=None):
    self.starttime=None # start time of data
    self.endtime=None # end time of data
    self.bbox=args.b # geographical bounding box [w,e,s,n]
    self.preferredOrigin='p' # prefferred should be 'p': preferred origin ; 'f': first origin ; 'l': last origin
    self.includeunassoc=False # should we look also at unassociated origins (default - False)
    self.minmag=0 # minimum magnitude
    self.maxmag=10.0 # maximum magnitude
    self.deltaT=60.0 # absolute time difference in seconds to associate events
    self.deltaR=100.0 # absolute distance difference in km to associate events
    self._datadict = {} # dictionary of input events
    self._refdict = {} # dictionary of reference events
    self.data2reflut = {} # input to ref look up table
    self.ref2datalut = {} # ref to input look up table
    self.Bbox = bboxForm()
    self.Bbox.setLims(*self.bbox)
    splash.showMessage('Initializing...',Qt.AlignCenter)
    QApplication.processEvents()
    QMainWindow.__init__(self, parent)
    self.setWindowTitle('SC3 XML Diff')
    self.args = args # save command line arguments
    self.meter = mpl.lines.Line2D([],[],color='r') # line for measuring distances along canvas
    splash.showMessage('Creating application...',Qt.AlignCenter)
    QApplication.processEvents()
    self.create_menu() # create the app top menu
    self.create_status_bar() # add a status bar
    self.create_main_frame() # create the main frame.
    self.create_toolbox() # create a toolbox
    if args.i:
      splash.showMessage('Reading input data...',Qt.AlignCenter)
      QApplication.processEvents()
      self.get_input_file(args.i)
    if args.r:
      self.get_ref_file(args.r)
      splash.showMessage('Reading reference data...',Qt.AlignCenter)
      QApplication.processEvents()
    self.init_connections() # initialize signal connections
    splash.showMessage('Populating tables...',Qt.AlignCenter)
    QApplication.processEvents()
    self.updatetables()

  def init_connections(self):
    '''Connect signals to functions.
       Using signals to run functions from subprocesses.'''
    self.Bbox.accepted.connect(self.onBboxAccepted) # connect bbox to dialog ok button
    self.connect(self.minmagLine, SIGNAL('ok'),self.updateMagLims)
    self.connect(self.maxmagLine, SIGNAL('ok'),self.updateMagLims)
    self.connect(self.deltaTLine, SIGNAL('ok'),self.updateDeltaT)
    self.connect(self.deltaRLine, SIGNAL('ok'),self.updateDeltaR)
    self.starttimeLine.dateTimeChanged.connect(self.updatestarttime)
    self.endtimeLine.dateTimeChanged.connect(self.updateendtime)
    self.includeunassocLine.stateChanged.connect(self.updateUnassoc)
    self.preferredOriginLine.currentIndexChanged[str].connect(self.updatePreferredOrigin)
#    self.canvas.mpl_connect('button_press_event',self.on_click) # connect click on canvas
#    self.canvas.mpl_connect('button_release_event',self.on_unclick) # connect mouse button release on canvas
#    self.canvas.mpl_connect('motion_notify_event',self.on_move) # connect mouse motion on canvas

  def create_main_frame(self):
    'Create the main frame of application'
    # create widget
    self.main_frame = QWidget()
    # create a general layout
    vbox = QVBoxLayout()
    # create a side by side layout for tables
    self.tb = self.addToolBar('Tools')
    # tool bar can be moved around
    self.tb.setMovable(True)
    self.tb.setFloatable(True)
    hsplit = QSplitter(Qt.Horizontal)
    self.headerlabels = HEADERLABELS
    ######################################
    # create a layout for reference data #
    ######################################
    vboxref = QVBoxLayout()
    # add a title
    label = QLabel('Reference Data')
    # create the reference table
    self.reftable = QTableWidget(0,len(self.headerlabels))
    # Add header to table
    self.reftable.setHorizontalHeaderLabels(self.headerlabels)
    #self.reftable.setSortingEnabled(True)
    # settings for drag&drop
    self.reftable.setDragEnabled(True)
    self.reftable.setDragDropOverwriteMode(False)
    self.reftable.setDragDropMode(3) # can drag and drop
    # connect click events to row selection
    self.reftable.cellClicked.connect(self.reftable_cellClicked)
    # can select only rows
    self.reftable.setSelectionBehavior(QAbstractItemView.SelectRows)
    # no edits please
    #self.reftable.setEditTriggers(QAbstractItemView.NoEditTriggers)
    # No auto scrolling
    self.reftable.setAutoScroll(False)
    # connect drop events to adding events
    self.reftable.dropEvent = self.reftable_dropEvent
    # add label to ref layout
    vboxref.addWidget(label)
    # add table to ref layout
    vboxref.addWidget(self.reftable)
    ##################################
    # create a layout for input data #
    ##################################
    vboxdata = QVBoxLayout()
    # add a title
    label = QLabel('Input data')
    # create the input table
    self.datatable = QTableWidget(0,len(self.headerlabels))
    # Add header to table
    self.datatable.setHorizontalHeaderLabels(self.headerlabels)
    #self.datatable.setSortingEnabled(True)
    # settings for drag&drop
    self.datatable.setDragEnabled(True)
    self.datatable.setDragDropOverwriteMode(False)
    self.datatable.setDragDropMode(3) # can drag and drop
    # connect click events to row selection
    self.datatable.cellClicked.connect(self.datatable_cellClicked)
    # can select only rows
    self.datatable.setSelectionBehavior(QAbstractItemView.SelectRows)
    # no edits please
    #self.datatable.setEditTriggers(QAbstractItemView.NoEditTriggers)
    # No auto scrolling
    self.datatable.setAutoScroll(False)
    # connect drop events to ID association
    self.datatable.dropEvent = self.datatable_dropEvent
    # add label to SC3 layout
    vboxdata.addWidget(label)
    # add table to SC3 layout
    vboxdata.addWidget(self.datatable)
    # Add layouts to side by side layout
    w = QWidget()
    w.setLayout(vboxref)
    hsplit.addWidget(w)
    w = QWidget()
    w.setLayout(vboxdata)
    hsplit.addWidget(w)
    # add side by side to general layout
    vbox.insertWidget(1,hsplit)
    # add layout to widget
    self.main_frame.setLayout(vbox)
    # set widget to be central widget
    self.setCentralWidget(self.main_frame)

  def reftable_cellClicked(self,row,col):
    'any click on table will select the whole row'
    # select all row
    self.reftable.selectRow(row)
    # clear any selection in data table
    self.datatable.clearSelection()

  def datatable_cellClicked(self,row,col):
    'any click on table will select the whole row'
    # select all row
    self.datatable.selectRow(row)
    # clear any selection in ref table
    self.reftable.clearSelection()

  def datatable_dropEvent(self,event):
    'Associate ref event with input event if ref row is dragged and dropped on input row'
    if not event.source()==self.reftable: return # make sure row is from reference table
    if not self.datatable.itemAt(event.pos()): return # make sure we dropped reference row on an input row
    # get the referrence ID
    refID = str(event.source().selectedItems()[0].text())
    # get the ID item on input table
    item = self.datatable.item(self.datatable.itemAt(event.pos()).row(),0)
    inputID = str(item.text())
    # confirm with user to associate
    if refID in self.data2reflut or inputID in self.data2reflut:
      self.statusBar().showMessage('%s and/or %s are already associated. Aborting.'%(refID,inputID) , 2000)
      return
    ok = QMessageBox.question(self,'Associate Event',
                           'Are you sure you want to associate\n%s to %s?'%(refID,inputID),
                           'No','Yes')
    if ok:
      # update the lut
      self.ref2datalut[refID]=inputID
      self.data2reflut[inputID]=refID
      self.updateStatus()
      # report to status bar
      self.statusBar().showMessage('Associated %s with %s'%(refID,inputID) , 2000)

  def reftable_dropEvent(self,event):
    'Add new events to reference table if input row is dragged and dropped on ref table'
    if not event.source()==self.datatable: return # make sure its an input row
    # get row items
    items = event.source().selectedItems()
    # get input ID
    ID = str(items[0].text())
    # make sure ID is associated already
    if ID in self.data2reflut:
      self.statusBar().showMessage('Event %s already associated with %s'%(ID,self.data2reflut[ID]) , 2000)
      return
    # disable sorting
    self.datatable.setSortingEnabled(False)
    self.reftable.setSortingEnabled(False)
    #update lut
    self.data2reflut[ID]=ID
    self.ref2datalut[ID]=ID
    # add a new row to reference table
    i = self.reftable.rowCount()
    self.reftable.insertRow(i)
    # copy line from input table
    for j,item in enumerate(items):
      newitem = QTableWidgetItem(item.text())
      self.reftable.setItem(i,j,newitem)
    # update status
    self.updateStatus()
    # sort ref table
    self.datatable.setSortingEnabled(True)
    self.reftable.setSortingEnabled(True)
    # report to statusbar
    self.statusBar().showMessage('Added Event %s'%(ID) , 2000)

  def delete_row(self):
    'Delete a selected row from table'
    # get selected row(s) on SC3 table
    items = [i for i in self.datatable.selectedItems()+self.reftable.selectedItems() if i.column()==0]
    for item in items:
      # get selected row number
      row = item.row()
      # get event ID
      ID = str(item.text())
      # remove the row from table
      item.tableWidget().removeRow(row)
      # remove lut records
      [self.data2reflut.pop(k) for k,v in self.data2reflut.items() if k==ID or v==ID]
      [self.ref2datalut.pop(k) for k,v in self.ref2datalut.items() if k==ID or v==ID]
      self.updateStatus()
      # report to statusbar
      self.statusBar().showMessage('Removed Event %s from table'%ID , 2000)
    if len(items)>1: self.statusBar().showMessage('%d Events removed from table'%len(items) , 2000)

  def clear_table(self,table):
    'Clear a table'
    for i in range(table.rowCount())[::-1]:
      table.removeRow(i)

  def clear_tables(self):
    [self.clear_table(t) for t in [self.datatable,self.reftable]]

  def create_tooltip_widget(self):
    'creates tooltip of figure elements'
    self.TTWidget = QLabel() # take a QLabel
    self.TTWidget.setFrameShape(QFrame.StyledPanel) # add a frame
    self.TTWidget.setWindowFlags(Qt.ToolTip) # make window look like a tooltip
    self.TTWidget.setAttribute(Qt.WA_TransparentForMouseEvents) # mouse events can't affet it.
    self.TTWidget.hide() # hide for now. see self.on_move function on how to use.

  def get_input_file(self,filesurl=None):
    'open a dialog to get a file name'
    if not filesurl: filesurl = QFileDialog.getOpenFileNames(self, 'Open Input data file(s)',filter='*.xml') # get the file name
    self.ref2datalut={} # init look up table
    self.data2reflut={} # init look up table
    if len(filesurl):
      xmls = [str(f) for f in filesurl]
      Dataep = concatenateSC3XML(xmls)
      if not Dataep:
        self.statusBar().showMessage('No Events in file.', 2000)
        return
      self._datadict = eparams2seisEvents(Dataep,starttime=self.starttime,endtime=self.endtime,minmag=self.minmag,maxmag=self.maxmag,bbox=self.bbox,preferred=self.preferredOrigin,includeunassoc=self.includeunassoc)
      self.statusBar().showMessage('Loaded %d events'%len(self._datadict), 2000)
    self.updatetables()

  def get_ref_file(self,filesurl=None,refFileType=None):
    'open a dialog to get a file name'
    if not filesurl:
      filesurl,refFileType = QFileDialog.getOpenFileNamesAndFilter(self, 'Open Reference data file(s)',filter='XML Files (*.xml);;CSV Files (*.csv)') # get the file name
    else:
      refFileType =  os.path.splitext(filesurl[0])[1].upper()
    self.ref2datalut={}# init look up table
    self.data2reflut={}# init look up table
    if 'XML' in refFileType:
      xmls = [str(f) for f in filesurl]
      Refep = concatenateSC3XML(xmls)
      if not Refep:
        self.statusBar().showMessage('No Events in file.', 2000)
        return
      self._refdict = eparams2seisEvents(Refep,starttime=self.starttime,endtime=self.endtime,minmag=self.minmag,maxmag=self.maxmag,bbox=self.bbox,preferred=self.preferredOrigin,includeunassoc=self.includeunassoc)
      self.statusBar().showMessage('Loaded %d events'%len(self._refdict), 2000)
      return self.updatetables()
    if 'CSV' in refFileType:
      CSVs = [str(f) for f in filesurl]
      Refdb = concatenateCSV(CSVs)
      if not len(Refdb):
        self.statusBar().showMessage('No Events in file.', 2000)
        return self.updatetables()
      self._refdict = db2seisEvents(Refdb,starttime=self.starttime,endtime=self.endtime,minmag=self.minmag,maxmag=self.maxmag,bbox=self.bbox,preferred=self.preferredOrigin,includeunassoc=self.includeunassoc)
      self.statusBar().showMessage('Loaded %d events'%len(self._refdict), 2000)
      return self.updatetables()
    self.statusBar().showMessage("Can't load reference file(s)", 3000)

  def filterEvents(self,D):
    'filter events by time,location and magnitude limits'
    [D.pop(i) for i in [k for k in D if D[k].mag>self.maxmag or D[k].mag<self.minmag]]
    if self.starttime: [D.pop(i) for i in [k for k in D if D[k].ot<self.starttime]]
    if self.endtime: [D.pop(i) for i in [k for k in D if D[k].ot>self.endtime]]
    [D.pop(i) for i in [k for k in D if D[k].lon>self.bbox[1] or D[k].lon<self.bbox[0]]]
    [D.pop(i) for i in [k for k in D if D[k].lat>self.bbox[3] or D[k].lat<self.bbox[2]]]
    return D

  def associate(self):
    'associate events'
    self.ref2datalut = {}
    self.data2reflut = {}
    # filter unwanted events
    self.refdict = {}
    self.datadict = {}
    self.refdict.update(self._refdict)
    self.datadict.update(self._datadict)
    self.refdict = self.filterEvents(self.refdict)
    self.datadict= self.filterEvents(self.datadict)
    for ref in self.refdict.keys(): # get potentials associated events for reference. (time diff,seisevent obj) list
      potentials = [(abs(self.datadict[k].ot-self.refdict[ref].ot).total_seconds(),k) for k in self.datadict if
                     abs(self.datadict[k].ot-self.refdict[ref].ot).total_seconds()<=self.deltaT
                     and
                     (self.datadict[k]-self.refdict[ref])[-1]/1000.0<=self.deltaR]

      potentials.sort() # sort according to time difference
      if potentials:
        self.ref2datalut[ref] = potentials[0][1] # get first one.
        self.data2reflut[potentials[0][1]] = ref

  def recolor_table(self,table):
    ' add colors to lines. Green - True event, Red - false event, orange - missed event'
    # mark missed events on table
    stat = HEADERLABELS.index('Status')
    for i in xrange(table.rowCount()):
      table.selectRow(i)
      for item in table.selectedItems():
        f = item.foreground()
        if str(table.item(i,stat).text())=='M':
          f.setColor(QColor('#FFA500'))
        elif str(table.item(i,stat).text())=='F':
          f.setColor(QColor('red'))
        elif str(table.item(i,stat).text())=='T':
          f.setColor(QColor('k'))
        item.setForeground(f)
    table.clearSelection()


  def recolor_tables(self):
    [self.recolor_table(table) for table in [self.datatable,self.reftable]]

  def build_tables(self):
    'populate the gui tables from dictionaries'
    # populate ref table
    self.reftable.setSortingEnabled(False)
    for i,e in enumerate(self.refdict):
      # add a row
      self.reftable.insertRow(i)
      # create an ID item
      item = QTableWidgetItem(e)
      # add ID to row
      self.reftable.setItem(i,0,item)
      # insert event date to row
      for j,v in enumerate([self.refdict[e].__dict__[k] for k in ['ot','mag','lat','lon','status','refid','SCevaluationStatus']]):
        item = QTableWidgetItem(str(v))
        self.reftable.setItem(i,j+1,item)
    # sort the table
    self.reftable.setSortingEnabled(True)
    # adjust columns size
    self.reftable.resizeColumnsToContents()
    # populate data table
    self.datatable.setSortingEnabled(False)
    for i,e in enumerate(self.datadict):
      # add a row
      self.datatable.insertRow(i)
      # create an ID cell
      item = QTableWidgetItem(e)
      # add ID to row
      self.datatable.setItem(i,0,item)
      # insert event date to row
      for j,v in enumerate([self.datadict[e].__dict__[k] for k in ['ot','mag','lat','lon','status','refid','SCevaluationStatus']]):
        item = QTableWidgetItem(str(v))
        self.datatable.setItem(i,j+1,item)
    # sort the table
    self.datatable.setSortingEnabled(True)
    # adjust columns size
    self.datatable.resizeColumnsToContents()
    # update status
    self.updateStatus()

  def updateStatus(self):
    'update status of event (T,M,F)'
    # referrece table
    stat = HEADERLABELS.index('Status')
    refID = HEADERLABELS.index('refID')
    self.reftable.setSortingEnabled(False)
    for i in xrange(self.reftable.rowCount()):
      ID = str(self.reftable.item(i,0).text())
      if ID in self.ref2datalut:
        self.reftable.item(i,stat).setText('T')
        self.reftable.item(i,refID).setText(self.ref2datalut[ID])
      else:
        self.reftable.item(i,stat).setText('M')
        self.reftable.item(i,refID).setText('')
    self.reftable.setSortingEnabled(True)
    self.datatable.setSortingEnabled(False)
    for i in xrange(self.datatable.rowCount()):
      ID = str(self.datatable.item(i,0).text())
      if ID in self.data2reflut:
        self.datatable.item(i,stat).setText('T')
        self.datatable.item(i,refID).setText(self.data2reflut[ID])
      else:
        self.datatable.item(i,stat).setText('F')
        self.datatable.item(i,refID).setText('')
    self.datatable.setSortingEnabled(True)
    # repaint missed False and true events
    self.recolor_tables()

  def updatetables(self):
    self.associate()
    self.clear_tables()
    self.build_tables()
    self.statusBar().showMessage('Tables updated',2000)

  def updateMagLims(self):
    self.minmag = float(self.minmagLine.text())
    self.maxmag = float(self.maxmagLine.text())
    self.statusBar().showMessage('Magnitude limits updated',2000)

  def updateDeltaT(self):
    self.deltaT = float(self.deltaTLine.text())
    self.statusBar().showMessage('Association time limit updated',2000)

  def updateDeltaR(self):
    self.deltaR = float(self.deltaRLine.text())
    self.statusBar().showMessage('Association distance limit updated',2000)

  def updatestarttime(self,starttime):
    self.starttime = starttime.toPyDateTime()
    self.statusBar().showMessage('Start time limit updated',2000)

  def updateendtime(self,endtime):
    self.endtime = endtime.toPyDateTime()
    self.statusBar().showMessage('End time limit updated',2000)

  def updateUnassoc(self,stat):
    self.includeunassoc = bool(stat)
    self.statusBar().showMessage('Including unassociated origins set to %s'%self.includeunassoc,2000)

  def updatePreferredOrigin(self,pref):
    self.preferredOrigin = pref[0]
    self.statusBar().showMessage('Preferred Origins set to %s'%pref,2000)

  def evaluate(self):
    i = self.headerlabels.index('Status')
    table = self.datatable
    T = len([table.item(j,i) for j in xrange(table.rowCount()) if str(table.item(j,i).text())=='T'])
    F = len([table.item(j,i) for j in xrange(table.rowCount()) if str(table.item(j,i).text())=='F'])
    table = self.reftable
    M = len([table.item(j,i) for j in xrange(table.rowCount()) if str(table.item(j,i).text())=='M'])
    self.statusBar().showMessage('T: %d ; M: %d ; F: %d ; Score: %.1f%% for %d events'%(T,M,F,100.0*T/sum([T,F,M]),sum([T,F,M])), 5000)
    return T,F,M,100.0*T/sum([T,F,M])

  def saveAs_figure(self):
    pass

  def create_menu(self):
    'Creates main menu'
    # Populate the menubar:
    # Add File submenu
    self.file_menu = self.menuBar().addMenu("&File")
    # load input data
    load_input_action = self.create_action("Load &Input",
            shortcut="Ctrl+I", slot=self.get_input_file,
            icon='document-open',tip="Load input data from a file(s)")
    # load reference data
    load_ref_action = self.create_action("Load &Reference",
            shortcut="Ctrl+R", slot=self.get_ref_file,
            icon='document-open',tip="Load reference data from a file(s)")
    # Save As...
    saveAs_action = self.create_action("S&ave As...",
            shortcut="Shift+S", slot=self.saveAs_figure,
            icon='filesaveas',tip="Save the figure")
    # Quit
    quit_action = self.create_action("&Quit", slot=self.close,
            icon='system-shutdown',shortcut="Ctrl+Q", tip="Close the application")
    # populate the file submenu
    self.add_actions(self.file_menu,
            (load_input_action,load_ref_action, saveAs_action, None, quit_action))
    # Add Edit submenu
    self.Edit_menu = self.menuBar().addMenu("&Edit")
        # delete an event
    delete_action   = self.create_action("&Delete",
            slot=self.delete_row,shortcut="Del",
            tip="Delete selected row from table")
    # populate tools submenu
    self.add_actions(self.Edit_menu,[delete_action])
    # Add help submenu
    self.help_menu = self.menuBar().addMenu("&Help")
    # Help
    help_action = self.create_action("&Help",
            shortcut='F1', slot=self.on_help,
            icon='Help',tip='help')
    # About
    about_action = self.create_action("&About",
            shortcut='F2', slot=self.on_about,
            tip='About This Application')
    # About QT
    aboutQt_action = self.create_action("&About QT",
            shortcut='F3', slot=self.on_aboutQt,
            tip='About QT')
    # License
    license_action = self.create_action("&License",
            shortcut='F4', slot=self.on_license,
            tip='Application License')
    # Populate help submenu
    self.add_actions(self.help_menu, (help_action,None,about_action,aboutQt_action,license_action))

  def create_toolbox(self):
    w = QWidget()
    l = QGridLayout(w)
    self.refreshbutton = self.create_pushButton('Refresh', slot=self.updatetables, shortcut=None, icon=None, tip='Refresh Tables')
    self.scorebutton = self.create_pushButton('Score', slot=self.evaluate, shortcut=None, icon=None, tip='Evaluate Score')
    l.addWidget(self.refreshbutton,0,0)
    l.addWidget(self.scorebutton,1,0)
    self.tb.addWidget(w)
    w = QWidget()
    l = QGridLayout(w)
    self.minmagLine = ConnectedLineEdit(str(self.minmag),0,self.maxmag)
    self.maxmagLine = ConnectedLineEdit(str(self.maxmag),self.minmag,10)
    self.minmagLine.setConnected(self.maxmagLine, 'bottom')
    self.maxmagLine.setConnected(self.minmagLine, 'top')
    l.addWidget(QLabel('Magnitude Limits'),0,1)
    l.addWidget(QLabel('Min'),1,0)
    l.addWidget(QLabel('Max'),2,0)
    l.addWidget(self.minmagLine,1,1)
    l.addWidget(self.maxmagLine,2,1)
    self.tb.addWidget(w)
    w = QWidget()
    l = QGridLayout(w)
    self.deltaTLine = ConnectedLineEdit(str(self.deltaT),0,18000,3)
    self.deltaRLine = ConnectedLineEdit(str(self.deltaR),0,1000,3)
    l.addWidget(QLabel('Association Limits'),0,1)
    l.addWidget(QLabel('Time (sec)'),1,0)
    l.addWidget(QLabel('Range (km)'),2,0)
    l.addWidget(self.deltaTLine,1,1)
    l.addWidget(self.deltaRLine,2,1)
    self.tb.addWidget(w)
    w = QWidget()
    v = QVBoxLayout(w)
    ww = QWidget()
    v.addWidget(ww)
    l = QGridLayout(ww)
    W,E,S,N = self.bbox
    self.EAST = QLabel(str(E))
    self.WEST = QLabel(str(W))
    self.NORTH = QLabel(str(N))
    self.SOUTH = QLabel(str(S))
    self.bboxbutton = self.create_pushButton('Region', slot=self.Bbox.show, shortcut=None, icon=None, tip='Define a geographic region')
    l.addWidget(self.EAST,1,2,1,1,Qt.Alignment(Qt.AlignLeft))
    l.addWidget(self.WEST,1,0,1,1,Qt.Alignment(Qt.AlignRight))
    l.addWidget(self.NORTH,0,1,1,1,Qt.Alignment(Qt.AlignHCenter))
    l.addWidget(self.SOUTH,2,1,1,1,Qt.Alignment(Qt.AlignHCenter))
    v.addWidget(self.bboxbutton)
    self.tb.addWidget(w)
    w = QWidget()
    l = QGridLayout(w)
    self.starttimeLine = QDateTimeEdit()
    self.starttimeLine.setDisplayFormat('yyyy-MM-dd HH:mm:ss')
    self.starttimeLine.setCalendarPopup(True)
    self.endtimeLine = QDateTimeEdit()
    self.endtimeLine.setDateTime(QDateTime.currentDateTimeUtc())
    self.endtimeLine.setDisplayFormat('yyyy-MM-dd HH:mm:ss')
    self.endtimeLine.setCalendarPopup(True)
    l.addWidget(QLabel('Start:'),0,0)
    l.addWidget(self.starttimeLine,0,1)
    l.addWidget(QLabel('End:'),1,0)
    l.addWidget(self.endtimeLine,1,1)
    self.tb.addWidget(w)
    w = QWidget()
    l = QGridLayout(w)
    self.preferredOriginLine = QComboBox()
    self.preferredOriginLine.addItems(['preferred','first','last'])
    self.includeunassocLine = QCheckBox('Unassociated')
    self.includeunassocLine.setCheckState(self.includeunassoc)
    l.addWidget(QLabel('Origins:'),0,0)
    l.addWidget(self.preferredOriginLine,1,0)
    l.addWidget(self.includeunassocLine,2,0)
    self.tb.addWidget(w)

  def onBboxAccepted(self):
    '''get region limits from BboxForm.
       fires when Bbox is accepted'''
    if self.Bbox.validate(): # make sure limits are acceptable
      self.bbox = self.Bbox.getLims() # get limits
      w,e,s,n = self.bbox
      self.WEST.setText(str(w))
      self.EAST.setText(str(e))
      self.NORTH.setText(str(n))
      self.SOUTH.setText(str(s))

  def add_actions(self, target, actions):
    'Utility function for menu creation'
    for action in actions:
      if action is None:
        target.addSeparator()
      else:
        target.addAction(action)

  def create_action(  self, text, slot=None, shortcut=None,
                      icon=None, tip=None, checkable=False,
                      signal="triggered()"):
    'Utility function for menu actions creation'
    action = QAction(text, self)
    action.setIconVisibleInMenu(True)
    if icon is not None:
      i = QIcon.fromTheme(icon,QIcon(":/%s.png" % icon))
      action.setIcon(i)
    if shortcut is not None:
      action.setShortcut(shortcut)
    if tip is not None:
      action.setToolTip(tip)
      action.setStatusTip(tip)
    if slot is not None:
      self.connect(action, SIGNAL(signal), slot)
    if checkable:
      action.setCheckable(True)
    return action

  def create_pushButton(self,text,toolbar=None, slot=None, shortcut=None, icon=None, tip=None):
    'Utility function for button creation'
    # create the button
    button = QPushButton(text,self)
    # populate properties
    if slot:
      # connect a function
      button.clicked.connect(slot)
    if icon:
      # add icon
      i = QIcon.fromTheme(icon,QIcon(":/%s.png" % icon))
      button.setIcon(i)
      button.setIconSize(QSize(24,24))
    if shortcut:
      # set the shortcut
      button.setShortcut(shortcut)
    if tip:
      # add tooltip and status tip
      button.setToolTip(tip)
      button.setStatusTip(tip)
    if toolbar:
      # add the button to a toolbar (or any widget)
      toolbar.addWidget(button)
    return button

  def create_status_bar(self):
    'Add a status bar'
    # set default message
    self.status_text = QLabel("Ready")
    self.connstatLabel = QLabel()
    self.statusBar().addWidget(self.status_text, 1)
    self.statusBar().addPermanentWidget(self.connstatLabel)


  def on_about(self):
    'show a messagebox about the application'
    msg = "<p align='center'><big>SCxmlDiff</big><br><br> \
    Compare two Seiscomp3 XML files<br><br> \
    <small>Created<br> \
    by<br> \
    Ran Novitsky Nof @ BSL, 2015</small><br><br>\
    <a href='http://ran.rnof.info/'>http://ran.rnof.info</a><p>"
    QMessageBox.about(self,"About", msg.strip())

  def on_aboutQt(self):
    'show a messagebox about QT'
    QMessageBox.aboutQt(self,'')

  def on_license(self):
    'GPL licanse message'
    msg = "<p><b>This</b> is a free software; you can redistribute it and/or modify it under the \
terms of the GNU General Public License as published by the Free Software \
Foundation; either version 3 of the License, or (at your option) any later \
version.</p>\
<p><b>This application</b> is distributed in the hope that it will be useful, but WITHOUT ANY \
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR \
A PARTICULAR PURPOSE.  See the GNU General Public License for more details.</p> \
<p>You should have received a copy of the GNU General Public License along with \
this application; if not, see <a href='http://www.gnu.org/licenses/'>http://www.gnu.org/licenses/</a>.</p>"
    QMessageBox.about(self,"Application Licanse", msg.strip())

  def on_help(self):
    'Show help on a message window. Uses the argparse help'
    msg = '<pre>'+parser.format_help()#.replace('\n','</p><p>')
    msg = msg.replace('ran.nof@gmail.com',"<a href='mailto:ran.nof@gmail.com'>ran.nof@gmail.com</a>").strip()+'<\pre>'
    QMessageBox.about(self,"Help", msg)

  def message(self,msg,title='Error'):
    'a simple message window'
    QMessageBox.about(self,title,msg)

def main(args):
  # create the application
  app = QApplication(sys.argv)
  #splash
  splash_pix = QPixmap('splash.png')
  splash_pix.scaledToHeight(50)
  splash = QSplashScreen(splash_pix)
  splash.setMask(splash_pix.mask())
  splash.show()
  app.processEvents()
  # populate the QT4 form
  appwin = AppForm(splash,args)
  appwin.show()
  splash.finish(appwin)
  # run the application
  sys.exit(app.exec_())

if __name__=="__main__":
  # parse the arguments
  args = parser.parse_args(sys.argv[1:])
  mpl.rcParams['font.size']=float(FONTSIZE)
  main(args)



w ='''
def getParams(xmlFile,Eid=None):
  ar = scio.XMLArchive() # seiscomp xml creator
  ar.setFormattedOutput(True) # output formatted xml file
  if not ar.open(xmlFile): sys.exit("Can't parse %s"%xmlFile)
  eparams = scdatamodel.EventParameters.Cast(ar.readObject())
  if Eid:
    events = [eparams.findEvent(Eid)]
  else:
    events = [eparams.event(i) for i in xrange(eparams.eventCount())]
  origins = []
  [origins.extend([eparams.findOrigin(e.originReference(i).originID()) for i in xrange(e.originReferenceCount())]) for e in events]
  mags = [origin.magnitude(origin.magnitudeCount()-1).magnitude().value() for origin in origins if origin]
  lats = [origin.latitude().value() for origin in origins if origin]
  lons = [origin.longitude().value() for origin in origins if origin]
  nSs  = [origin.quality().usedStationCount() for origin in origins if origin]
  nTs  = [origin.quality().usedPhaseCount() for origin in origins if origin]
  ots  = [origin.time().value().toDouble() for origin in origins if origin]
  return mags,lats,lons,nSs,nTs,ots

def plotParams(mags1,lats1,lons1,nSs1,nTs1,ots1,t1,mags2,lats2,lons2,nSs2,nTs2,ots2,t2,baseName):
  fig = figure(figsize=(8,6))
  ax = gca()
  title('Event M=%0.1lf ot=%s lat=%0.3lf lon=%0.3lf'%(mags1[0],datetime.datetime.fromtimestamp(ots1[0]).isoformat(),lats1[0],lons1[0]),size=14,va='top',position=[0.5,1.07])
  xlabel('Solution #',size=14)
  l1 = plot(mags1,'kx-',lw=3,label=t1)[0]
  l2 = plot(mags2,'x--',c='gray',lw=3,label=t2)[0]
  legend(frameon=False)
  ylabel('Magnitude',size=14)
  savefig(baseName+'_mags.jpg')
  l1.set_ydata(lats1)
  l2.set_ydata(lats2)
  mn,mx = min(lats1+lats2),max(lats1+lats2)
  d = (mx-mn)/5.0
  ylim(mn-d,mx+d)
  legend(frameon=False)
  ylabel('Latitude',size=14)
  savefig(baseName+'_lats.jpg')
  l1.set_ydata(lons1)
  l2.set_ydata(lons2)
  mn,mx = min(lons1+lons2),max(lons1+lons2)
  d = (mx-mn)/5.0
  ylim(mn-d,mx+d)
  legend(frameon=False)
  ylabel('Longitude',size=14)
  savefig(baseName+'_lons.jpg')
  l1.set_ydata(nTs1)
  l2.set_ydata(nTs2)
  mn,mx = min(nTs1+nTs2),max(nTs1+nTs2)
  d = (mx-mn)/5.0
  ylim(mn-d,mx+d)
  legend(frameon=False)
  ylabel('Number of Triggers',size=14)
  savefig(baseName+'_nTs.jpg')
  l1.set_ydata(nSs1)
  l2.set_ydata(nSs2)
  mn,mx = min(nSs1+nSs2),max(nSs1+nSs2)
  d = (mx-mn)/5.0
  ylim(mn-d,mx+d)
  legend(frameon=False)
  ylabel('Number of Stations',size=14)
  savefig(baseName+'_nSs.jpg')
  dot1 = array(ots1)-ots1[0]
  dot2 = array(ots2)-ots2[0]
  l1.set_ydata(dot1)
  l2.set_ydata(dot2)
  mn,mx = min(append(dot1,dot2)),max(append(dot1,dot2))
  d = (mx-mn)/5.0
  ylim(mn-d,mx+d)
  legend(frameon=False)
  ylabel('Origin Time (Relative to %s'%datetime.datetime.fromtimestamp(ots1[0]).isoformat().split('T')[1][:-3],size=14)
  savefig(baseName+'_ots.jpg')

def main(args):
  mags1,lats1,lons1,nSs1,nTs1,ots1 = getParams(args.xml1,args.E1)
  mags2,lats2,lons2,nSs2,nTs2,ots2 = getParams(args.xml2,args.E2)
  plotParams(mags1,lats1,lons1,nSs1,nTs1,ots1,args.t1,mags2,lats2,lons2,nSs2,nTs2,ots2,args.t2,args.o)


if __name__=='__main__':
  # parse the arguments
  args = parser.parse_args(sys.argv[1:])
  main(args)
'''