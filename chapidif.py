
import os

import sys
import webbrowser as browser
#from functools import partial
from pathlib import Path
from PyQt6.QtCore import QRect
from PyQt6.QtGui import  QDoubleValidator,QIntValidator, QAction, QGuiApplication, QBrush,  QCloseEvent,QColor
from PyQt6.QtWidgets import (QApplication,QLabel, QMainWindow,
    QMenu, QPushButton, QTabWidget,
    QWidget, QGroupBox, QLineEdit,
    QTableWidget, QRadioButton, QTableWidgetItem,
    QCheckBox,QButtonGroup ,QFrame  
)

### MODIF: IMPORTANT MODIFICATION POUR ECRIRE DES FLOAT ET ACTIVER LES TEXTBOX QUI NE FONCTIONNAIENT PAS ###
from PyQt6.QtCore import QLocale
QLocale.setDefault(QLocale(QLocale.Language.English, QLocale.Country.UnitedStates))
#####################################################################################################

from calculate import calculate
from run_and_plot import run_and_plot
from read_write_files import read_write_files
class myWindow(QMainWindow):
    """Main Window."""
    
    def __init__(self, parent=None):
        
        """Initializer."""
        super().__init__(parent)
        self.setWindowTitle("Chapidif, PyQt6, nanobind version (Jan. 2026)")
        self.MainWindowWidth=1130
        self.MainWindowHeight=800
        app.setStyleSheet("""
            QGroupBox {
                border: 2px solid blue;
                border-radius: 5px;
                margin-top: 1ex; /* leave space for the title */
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top center; /* centers all group box titles */
                padding: 0 3px;
                color: red;
            }
        """)
        screensize = QGuiApplication.primaryScreen().geometry()
        self.setGeometry(QRect(screensize.width()-self.MainWindowWidth,20,self.MainWindowWidth,self.MainWindowHeight))
        self.centralWidget = QWidget()      

        self.setCentralWidget(self.centralWidget)
      #  self.centralWidget.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
        self.centralWidget.setStyleSheet("font-family: Arial; font-size: 13px;")
        self.calc = calculate(self)
        self.runplot = run_and_plot(self)
        self.readwrite = read_write_files(self)
        self.runplot.screen_width=screensize.width()
        self.runplot.screen_height=screensize.height()
        self.current_dir = str(Path.home())
        menuBar =self.menuBar( )
     
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        self.ReadDFAction = QAction( "&Read DF from file", self)
        fileMenu.addAction(self.ReadDFAction)
        self.ReadDFAction.triggered.connect (self.readwrite.read_df)
        
        helpMenu = QMenu("&Help", self)
        menuBar.addMenu(helpMenu)
        self.HelpAction = QAction( "&Open documentation", self)
        helpMenu.addAction(self.HelpAction)
        self.HelpAction.triggered.connect (self.Help)
        self.CheckForUpdates_Action = QAction( "&Check for updates", self)
        helpMenu.addAction(self.CheckForUpdates_Action)
        self.CheckForUpdates_Action.triggered.connect (self.check_for_update)
        
        
        self.WriteDFAction = QAction("&Write DF to file",self)
        fileMenu.addAction(self.WriteDFAction)
        self.WriteDFAction.triggered.connect(self.readwrite.write_df)
        fileMenu.addSeparator()
        self.SaveCalcAction = QAction("&Save calculation to file",self)
        fileMenu.addAction(self.SaveCalcAction)
        self.SaveCalcAction.triggered.connect(self.readwrite.save_calculation)
        

        self.tabWidget = QTabWidget(parent=self.centralWidget)
        self.tabWidget.setEnabled(True)
        self.tabWidget.setGeometry(QRect(0, 0, self.MainWindowWidth, self.MainWindowHeight-70))
        self.tabWidget.setObjectName("tabWidget")
        self.tab_1 = QWidget()
        self.tab_1.setObjectName("tab1")
        self.tabWidget.addTab(self.tab_1, "define dielectric function")
        self.tab_2 = QWidget()
        self.tab_2.setObjectName("tab2")
        self.tabWidget.addTab(self.tab_2, "perform calculation")
        self.tab_3 = QWidget()
        self.tab_3.setObjectName("tab3")
        self.tabWidget.addTab(self.tab_3, "edit defaults")
  #====================================================================    
        self.MainWindowIO()
  #=======tab1 contents================================         
        
     
        self.calc.InitDF()   
        self.OscillatorTable_init() 
        self.GOSTable_init() 
        self.BelkacemTable_init()       
        self.KanekoTable_init()
        self.ModelChoice_init()
        self.addchi_elf_init()
        self.RPA_choice_init()
        self.Dispersion_init()
        self.Relativity_init()
        self.init_DF_comment()
        self.init_target_properties_box()
      
        self.calc.get_df_properties() 
        self.input_target_properties()
        self.zero_DF_Btn()
     #   self.update_all_tables()
#tab 2 content
        self.calculate_which_property()
        self.PlotOptions()
        self.REELSOptions()
        self.comparison_data_box()
        self.Projectile_target_Interaction()
        self.transform_df()
#tab3 content
        self.calculation_options()
        self.plot_options()
        
        self.update_all_tables()
        
    def OscillatorTable_init(self):
        OscBox = QGroupBox(parent=self.tab_1,title="Oscillator Parameters:")
        OscBox.setGeometry(QRect(4,4, 315, 530))
        self.OscTable = QTableWidget(self.calc.maxOscillators,5,parent = OscBox)
        self.OscTable.setGeometry(QRect(7, 30, 300, 485))
        self.OscTable.setColumnWidth(0,65)
        self.OscTable.setColumnWidth(1,50)
        self.OscTable.setColumnWidth(2,50)
        self.OscTable.setColumnWidth(3,35)
        self.OscTable.setColumnWidth(4,40)
     
        for i in range(self.calc.maxOscillators):
            self.OscTable.setRowHeight(i,15)
        self.OscTable.cellChanged[int,int].connect(self.OscTableChangedAction)
       

    def OscTableChangedAction(self, x, y):  # changes only one element
        if self.Updating: return
        a= self.OscTable.item(x,y).text()
    
        try:
            b=float(a)
        except Exception:
            self.OscTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))
            return 
        if y == 0: self.calc.Amps[x] = b
        elif y == 1:self.calc.Omegas[x] = b
        elif y == 2: self.calc.Gammas[x] = b
        elif y == 3:
            if not( self.calc.DFmodel== "Mermin" or self.calc.DFmodel== "Vlasov"):
                self.calc.Alphas[x] = b
        elif y == 4: self.calc.Us[x] = b
        
        if (self.calc.DFmodel== "Tauc" or self.calc.DFmodel== "TL_an"):
            
            if self.calc.Omegas[x] < self.calc.Us[x]:
                self.UpdateStatus(f"ERROR, osc. {x+1}: gap U larger than Eₓ")
                return
            if self.calc.Gammas[x] > 2*self.calc.Omegas[x]:
                self.UpdateStatus(f"ERROR, osc. {x+1}: width Cx larger than 2*E_x")
                return
        self.OscTable.item(x,y).setForeground(QBrush(QColor(0, 0, 0))) 
        self.calc.get_df_properties()
        self.Update_target_properties()   
        self.updating=False
        
            
    def Update_OscTable(self):  #rewrites all elements
        self.Updating=True
        
        for x in range (self.calc.maxOscillators):
            self.OscTable.setItem(x,0,QTableWidgetItem( self.MyFormat1(self.calc.Amps[x])))
            self.OscTable.setItem(x,1,QTableWidgetItem( self.MyFormat2(self.calc.Omegas[x])))
            self.OscTable.setItem(x,2,QTableWidgetItem( self.MyFormat2(self.calc.Gammas[x])))
            if  not( self.calc.DFmodel== "Mermin" or self.calc.DFmodel== "Vlasov"):
                self.OscTable.setItem(x,3,QTableWidgetItem(str(self.calc.Alphas[x])))
            else:
                self.OscTable.setItem(x,3,QTableWidgetItem("-"))
            self.OscTable.setItem(x,4,QTableWidgetItem(self.MyFormat2(self.calc.Us[x])))
        self.Updating=False 
       
        self.additional_input.setGeometry(QRect(10,140,50,18))  
        self.additional_input_value.setGeometry(QRect(60,140,70,18)) 
   
#========================oscillator model related=================================            
    def ModelChoice_init(self):
        self.ModelChoiceBox = QGroupBox(parent=self.tab_1, title="Oscillator Model:")
        self.ModelChoiceBox.setGeometry(QRect(4,537, 315,162))
        DFchoice = ("Extended Drude (+ U)","Drude-Lindhard (+ U)","Mermin (+ U)", "Vlasov","Tauc-Lorentz", "TL-analytic ","Tauc-Mermin", "Forouhi-Bloomer", "Brendel-Bormann", "Orosco-Coimbra") #tuples here (round brackets, as there is no need for them to change
        self.DFchoice_short = ("Drude","DL","Mermin",  "Vlasov","Tauc","TL_an","TL_Mermin","FB","BB","OC")
        self.rb_modelchoice = []
        self.additional_input= QLabel("", parent=self.ModelChoiceBox) 
        self.additional_param=1.0
        self.additional_input_value = QLineEdit(str(self.additional_param),parent= self.ModelChoiceBox ) 
      #  self.additional_input_value.setVisible(False)
        self.additional_input_value.editingFinished.connect(self.Update_additional_param)
             
   
        for i in range(len(self.DFchoice_short)):
            rb=QRadioButton(DFchoice[i], parent= self.ModelChoiceBox )
            if i < 4:
                rb.setGeometry(QRect(10,21*i+25, 165,25))
            else:
                rb.setGeometry(QRect(185,21*(i-4)+25, 165,25))
            rb.clicked.connect(self.rb_modelchoice_clicked)
            self.rb_modelchoice.append(rb)
        self.Update_OscTable()  

    def Update_additional_param(self):
        if self.calc.DFmodel == "Drude" or self.calc.DFmodel == "Tauc":
            self.calc.eps_bkg = float( self.additional_input_value.text())
        elif self.calc.DFmodel == "TL_an":    
            self.calc.a_TL_an =  float(self.additional_input_value.text())
        elif self.calc.DFmodel == "FB":    
            self.calc.n_infty = float( self.additional_input_value.text())
            

    def rb_modelchoice_clicked(self):
        for i in range (len(self.DFchoice_short)):
            if self.rb_modelchoice[i].isChecked():
                self.calc.DFmodel=self.DFchoice_short[i]
                self.update_header_osc_table()
        self.Update_OscTable()  
        self.calc.get_df_properties()  
        self.Update_target_properties()       
                                
    def set_modelchoice(self):
        try:
            i= self.DFchoice_short.index(self.calc.DFmodel)
        except ValueError:
            self.UpdateStatus("model not implemented:"+ self.calc.DFmodel)
        self.rb_modelchoice[i].setChecked(True)  
        self.update_header_osc_table() 
        
    
    def update_header_osc_table(self):
        if self.calc.DFmodel == "Drude":
            self.myHeaderLabels=["Aₓ\n(eV²)","ωₓ\n(eV)","Γₓ\n(eV)","αₓ","U\n(eV)"]
            self.additional_input.setText("ϵ(∞)=")
            self.additional_input_value.setText(str(self.calc.eps_bkg))
            self.additional_input_value.setValidator(QDoubleValidator(1, 100, 6))
            self.additional_input_value.setVisible(True)
           # self.Change_U.setText("U Factor")
            self.changeDFButtonText = "Drude to DL"
        elif self.calc.DFmodel == "DL":
            self.myHeaderLabels=["Aₓ","ωₓ\n(eV)","Γₓ\n(eV)","αₓ","U\n(eV)"]
            self.additional_input.setText("")
            self.additional_input_value.setVisible(False)
            self.changeDFButtonText = "DL to Drude"
           
        elif self.calc.DFmodel == "Mermin":
            self.myHeaderLabels=["Aₓ","ωₓ\n(eV)","Γₓ\n(eV)","-","U\n(eV)"] 
            self.additional_input.setText("")
            self.additional_input_value.setVisible(False)
            self.changeDFButtonText = "Mermin to Drude"
          
        elif self.calc.DFmodel == "Tauc" :
            self.myHeaderLabels=["Aₓ\n(eV)","Eₓ\n(eV)","Cₓ","αₓ","Egapₓ\n(eV)"] 
            self.additional_input.setText("ϵ(∞)=")
            self.additional_input_value.setText(str(self.calc.eps_bkg))
            self.additional_input_value.setValidator(QDoubleValidator(1, 100, 6))
            self.additional_input_value.setVisible(True)
           
        elif self.calc.DFmodel == "TL_an":
            self.myHeaderLabels=["Aₓ\n(eV)","Eₓ\n(eV)","Cₓ","αₓ","Egapₓ\n(eV)"]  
            self.additional_input.setText("A (eV)=")   
            self.additional_input_value.setText(str(self.calc.a_TL_an))
            self.additional_input_value.setValidator(QDoubleValidator(1e-6, 10, 6))
            self.additional_input_value.setVisible(True) 
            
        elif self.calc.DFmodel == "TL_Mermin":
            self.myHeaderLabels=["Aₓ\n(eV)","Eₓ\n(eV)","Cₓ","αₓ","Egapₓ\n(eV)"]  
            self.additional_input.setText("ϵ(∞)=")
            self.additional_input_value.setText(str(self.calc.eps_bkg))
            self.additional_input_value.setValidator(QDoubleValidator(1, 100, 6))
            self.additional_input_value.setVisible(True)
           
        elif self.calc.DFmodel == "Vlasov":
            self.myHeaderLabels=["Aₓ\n","Qₓ\n(a.u.)","Γₓ\n(eV)","-","U \n(eV)"]  
            self.additional_input.setText("")   
            self.additional_input_value.setVisible(False)
          
        elif self.calc.DFmodel == "FB": 
            self.myHeaderLabels=["Aₓ\n","Bₓ\n(eV)","Cₓ\n(eV²)","αₓ","Egapₓ\n(eV)"] 
            self.additional_input.setText("n(∞)=") 
            self.additional_input_value.setText(str(self.calc.n_infty))
            self.additional_input_value.setValidator(QDoubleValidator(1, 100, 6))
            self.additional_input_value.setVisible(True)
           
        elif self.calc.DFmodel == "BB" or self.calc.DFmodel == "OC": 
            self.myHeaderLabels=["Aₓ\n (eV²)","ωₓ\n(eV)","Γₓ\n(eV)","σₓ","U\n(eV)"]
            self.additional_input.setText("")
            self.additional_input_value.setVisible(False) 
           
        self.OscTable.setHorizontalHeaderLabels(self.myHeaderLabels)  
        if self.calc.DFmodel=="DL" or self.calc.DFmodel== "Mermin":
            # self.Change_U.setText("change U")
            # self.Change_U.setVisible(True) 
            self.ulabel.setVisible(True)
            self.U_inp.setVisible(True)
            self.DF_to_Kaneko.setVisible(True)
            self.l_inp.setVisible(True)
            self.Q_inp.setVisible(True)
            self.au_label.setVisible(True)
            self.Q_label.setVisible(True)
            self.l_label.setVisible(True)
            self.Change_DF_model.setVisible(True)
            if self.calc.DFmodel=="DL" :
                self.Change_DF_model.setText("DL to Drude")
            else:
                self.Change_DF_model.setText("Memin to Drude")
                
            
        else:
          #  self.Change_U.setVisible(False) 
            self.ulabel.setVisible(False)
            self.U_inp.setVisible(False)
            self.DF_to_Kaneko.setVisible(False)
            self.l_inp.setVisible(False)
            self.Q_inp.setVisible(False)
            self.au_label.setVisible(False)
            self.Q_label.setVisible(False)
            self.l_label.setVisible(False)
            self.Change_DF_model.setVisible(False)
        if self.calc.DFmodel=="Drude":
            self.ulabel.setVisible(True)
            self.U_inp.setVisible(True)
            self.Change_DF_model.setText("Drude to DL")
            self.Change_DF_model.setVisible(True)
          
        if self.calc.DFmodel == "Tauc"  or self.calc.DFmodel == "TL_an":
            self.dispersion_based_on.setVisible(True)
            self.rb_DispersionChoice_based_on_button[0].setVisible(True)       
            self.rb_DispersionChoice_based_on_button[1].setVisible(True)
            self.rb_DispersionChoice[0].setVisible(False)   
            self.rb_DispersionChoice[1].setVisible(False) 
            self.extendeddispersionlabel.setVisible(False)      
        else:
            self.dispersion_based_on.setVisible(False) 
            self.rb_DispersionChoice_based_on_button[0].setVisible(False)       
            self.rb_DispersionChoice_based_on_button[1].setVisible(False)  
            self.rb_DispersionChoice[0].setVisible(True)   
            self.rb_DispersionChoice[1].setVisible(True)   
            self.extendeddispersionlabel.setVisible(True)   
        if self.calc.DFmodel == "Tauc"  or self.calc.DFmodel == "TL_an" or self.calc.DFmodel == "Drude" or self.calc.DFmodel == "DL" :   
            self.rb_Gamma_Dispersion_button[0].setVisible(True) 
            self.rb_Gamma_Dispersion_button[1].setVisible(True) 
            self.disperionlabel.setVisible(True)
        else:
            self.rb_Gamma_Dispersion_button[0].setVisible(False) 
            self.rb_Gamma_Dispersion_button[1].setVisible(False)  
            self.disperionlabel.setVisible(False)   
               
            
#==============================end Oscillator I/O==begin GOS I/0==========================================================================            
            
    def GOSTable_init(self):
        self.GOSBox = QGroupBox(parent=self.tab_1,title="Hydrogenic GOS Parameters:")#\n (per unit cell):")
        self.GOSBox.setGeometry(QRect(325,4, 235, 363))
        
        self.GOSTable = QTableWidget(self.calc.maxGOS,4,parent = self.GOSBox)
        self.GOSTable.setGeometry(QRect(10, 22, 210, 220))
        self.GOSTable.setHorizontalHeaderLabels(("Nₓ","Uₓ(eV)","10*n+l","Z"))
        
        self.GOSTable.setColumnWidth(0,32)
        self.GOSTable.setColumnWidth(1,56)
        self.GOSTable.setColumnWidth(2,61)
        self.GOSTable.setColumnWidth(3,33)
        for i in range(self.calc.maxGOS):
            self.GOSTable.setRowHeight(i,15)
        self.GOSTable.cellChanged[int,int].connect(self.GOSTableChangedAction) 
        self.GOS_options_init()                   
            
    def Update_GOSTable(self):  
        self.Updating=True
        for x in range (self.calc.maxGOS):
            self.GOSTable.setItem(x,0,QTableWidgetItem(str(self.calc.ConcGOS[x])))
            self.GOSTable.setItem(x,1,QTableWidgetItem(str(self.calc.EdgeGOS[x])))
            self.GOSTable.setItem(x,2,QTableWidgetItem(str(self.calc.nlGOS[x])))
            self.GOSTable.setItem(x,3,QTableWidgetItem(str(self.calc.ZGOS[x])))
        self.Updating= False    

 
    def GOSTableChangedAction(self, x, y):  
        if self.Updating: return
        myinput= self.GOSTable.item(x,y).text()
        try:
            a=float(myinput)
        except Exception:
            # nb, you can only have one setForeground color per funtion call, else race condition
            self.GOSTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))
            return
        if y == 0: self.calc.ConcGOS[x] = a
        elif y == 1: self.calc.EdgeGOS[x] = a
        elif y == 2:
            try:
                a=int(myinput)
                n = int(a/10)
                L = a-10*n 
                self.calc.nlGOS[x]= a
                if (L >= n) or (a <= 0) or (n > 3):
                    self.GOSTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))
                    return
                self.calc.nlGOS[x]= a    
            except Exception:
                self.GOSTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))
                return   
        elif y == 3:
            try:
                a=int(myinput)
                if( a <= 0) or (a>100):
                    self.GOSTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))  
                    return          
                self.calc.ZGOS[x]= a
            except Exception:
                self.GOSTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))
                return    
        self.GOSTable.item(x,y).setForeground(QBrush(QColor(0, 0, 0)))            
        self.calc.get_df_properties()
        self.Update_target_properties()   
                       
    def GOS_options_init(self):
        ypos=250
        QLabel("<font color='red'>GOS options:",parent=self.GOSBox).setGeometry(QRect(5,ypos, 300,18))
        
        self.PrecisionCheckbox = QCheckBox('precise (slow)', self.GOSBox)
        self.PrecisionCheckbox.setChecked(False)
        self.PrecisionCheckbox.setGeometry(QRect(117,ypos, 110,18))
        self.PrecisionCheckbox.stateChanged.connect(self.PrecisionCheckbox_state_changed)
        
        self.rb_GOSrescaling = []
        rb =QRadioButton("rescaling off", parent= self.GOSBox )
        rb.setGeometry(QRect(5,ypos+12,100,40))
        rb.clicked.connect(self.rb_GOSrescaling_clicked)
        self.rb_GOSrescaling.append(rb)
        rb =QRadioButton("rescaling on", parent= self.GOSBox  )
        rb.setGeometry(QRect(120,ypos+12,100,40))
        rb.clicked.connect(self.rb_GOSrescaling_clicked)
        self.rb_GOSrescaling.append(rb) 
      
        QLabel("Max ω density correction (eV):",parent=self.GOSBox ).setGeometry(QRect(5,ypos+45, 210,18))
        self.GOSDesityCor = QLineEdit(str(self.calc.maxEnergyDensityEffect),parent=self.GOSBox)
        self.GOSDesityCor.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.GOSDesityCor.setGeometry(QRect(185,ypos+45, 40,18))
        self.GOSDesityCor.editingFinished.connect(self.get_maxEnergyDensityEffect_from_lineEdit)  
         
    def  rb_GOSrescaling_clicked(self):
        for i in range (2):
            if self.rb_GOSrescaling[i].isChecked():
                 self.calc.ApplySumRuleToGOS=i  
                 
    def set_rb_GOSrescaling(self,i): self.rb_GOSrescaling[i].setChecked(True) 

    def get_maxEnergyDensityEffect_from_lineEdit(self):self.calc.maxEnergyDensityEffect= float(self.GOSDesityCor.text()) 
    def set_maxEnergyDensityEffect(self,a):self.GOSDesityCor.setText(str(a)) 
    def PrecisionCheckbox_state_changed(self):self.calc.PreciseDenityCor = self.PrecisionCheckbox.isChecked()
            
    def set_massunitcell(self,a):
        self.molweight.setText(str(a))     
    def set_TargetDensity(self,a):
        self.densinp.setText(str(a)) 

                    
#==============================end GOS I/O==begin Belkacem I/0==========================================================================            
            
    def BelkacemTable_init(self):
        BelkacemBox = QGroupBox(parent=self.tab_1, title="Belkacem Parameters:")
        BelkacemBox.setGeometry(QRect(325,365, 235, 170))
        self.BelkacemTable = QTableWidget(self.calc.maxBelkacem,3,parent = BelkacemBox)
        self.BelkacemTable.setGeometry(QRect(10, 20, 210,145))
        self.BelkacemTable.setHorizontalHeaderLabels(("Nₓ","ωₓ (eV)","Γₓ (eV)"))
        self.BelkacemTable.setColumnWidth(0,60)
        self.BelkacemTable.setColumnWidth(1,60)
        self.BelkacemTable.setColumnWidth(2,60)
        for i in range(self.calc.maxBelkacem):
            self.BelkacemTable.setRowHeight(i,15)
        self.BelkacemTable.cellChanged[int,int].connect(self.BelkacemTableChangedAction)
            
   
    def Update_BelkacemTable(self): 
        self.Updating = True 
        for x in range (self.calc.maxBelkacem):
            self.BelkacemTable.setItem(x,0,QTableWidgetItem(str(self.calc.Conc_Belkacem [x])))
            self.BelkacemTable.setItem(x,1,QTableWidgetItem(str(self.calc.w_Belkacem[x])))
            self.BelkacemTable.setItem(x,2,QTableWidgetItem(str(self.calc.gamma_Belkacem[x])))
        self.Updating= False    
        self.calc.get_df_properties()
        self.Update_target_properties()
                 
    def BelkacemTableChangedAction(self, x, y): 
        if self.Updating: return 
        myinput = self.BelkacemTable.item(x,y).text()
        try:
            a=float(myinput)
        except Exception:
            self.BelkacemTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))
            return
        if y == 0: self.calc.Conc_Belkacem[x] = a
        elif y == 1: self.calc.w_Belkacem[x] = a
        elif y == 2: self.calc.gamma_Belkacem[x] = a 
        self.BelkacemTable.item(x,y).setForeground(QBrush(QColor(0, 0, 0)))  
        self.calc.get_df_properties()
        self.Update_target_properties()           
#==============================end Belkacem I/O==begin Kaneko I/0==========================================================================
    def KanekoTable_init(self):
        KanekoBox = QGroupBox(parent=self.tab_1, title="Kaneko Parameters:" )
        KanekoBox.setGeometry(QRect(565,4, 318, 530))
        self.KanekoTable = QTableWidget(self.calc.maxKaneko,6,parent = KanekoBox)
        self.KanekoTable.setGeometry(QRect(7, 30, 300, 455))
        self.KanekoTable.setHorizontalHeaderLabels(("Nₓ","Qₓ\n(au)","Γₓ\n(eV)","lₓ","Uₓ\n(eV)","γₓ"))
        for i in  range (6):
            self.KanekoTable.setColumnWidth(i,40)
        for i in range(self.calc.maxKaneko):
            self.KanekoTable.setRowHeight(i,15)
        self.KanekoTable.cellChanged[int,int].connect(self.KanekoChangedAction)
        
        
        QLabel("Kaneko option:",parent=KanekoBox ).setGeometry(QRect(5,490, 125,20))
        KanekoChoice = { 0: "Modified Kaneko", 1: "Original kaneko"}
        self.rb_Kaneko_Choice = []
        for i in sorted(KanekoChoice.keys()):
            label = "%s" % (KanekoChoice[i])
            rb=QRadioButton(label, parent=  KanekoBox )
            rb.setGeometry(QRect(5+145*i,505, 140,25))
            rb.clicked.connect(self.rb_Kaneko_Choice_clicked)
            self.rb_Kaneko_Choice.append(rb)
        self.set_Kaneko_Choice (self.calc.Kaneko_choice)
       
    def Update_KanekoTable(self):  
        self.Updating = True
        for x in range (self.calc.maxKaneko):
            self.KanekoTable.setItem(x,0,QTableWidgetItem(str(self.calc.N_Kaneko[x])))
            self.KanekoTable.setItem(x,1,QTableWidgetItem(str(self.calc.Q_Kaneko[x])))
            self.KanekoTable.setItem(x,2,QTableWidgetItem(str(self.calc.width_Kaneko[x])))
            self.KanekoTable.setItem(x,3,QTableWidgetItem(str(self.calc.l_Kaneko[x])))
            self.KanekoTable.setItem(x,4,QTableWidgetItem(str(self.calc.Edge_Kaneko[x])))      
            self.KanekoTable.setItem(x,5,QTableWidgetItem(str(self.calc.gamma_Kaneko[x]))) 
        self.Updating = False             
    
    def KanekoChangedAction(self, x, y): 
        if self.Updating: return 
  
        myinput= self.KanekoTable.item(x,y).text()
        try:
            a=float(myinput)
        except Exception:
            self.KanekoTable.item(x,y).setForeground(QBrush(QColor(255, 0, 0)))
            return
        if y == 0: self.calc.N_Kaneko[x] = a
        elif y == 1: self.calc.Q_Kaneko[x] = a
        elif y == 2: self.calc.width_Kaneko[x] = a 
        elif y == 3: self.calc.l_Kaneko[x] = round(a)  
        elif y == 4: self.calc.Edge_Kaneko[x]= a   
        elif y == 5: self.calc.gamma_Kaneko[x] = a  
        self.KanekoTable.item(x,y).setForeground(QBrush(QColor(0, 0, 0)))
        self.calc.get_df_properties()
        self.Update_target_properties()  
        
        
    def rb_Kaneko_Choice_clicked(self):
        for i in range (2):
            if  self.rb_Kaneko_Choice[i].isChecked(): self.calc.Kaneko_choice= i    #    calc.Kaneko_choice=0 means modified
        
    def set_Kaneko_Choice (self,i): self.rb_Kaneko_Choice[i].setChecked(True)             
              
 #=====================================End Kaneko I/O====begin eps model choice=====================
              
    def update_all_tables(self):
     
        self.calc.get_df_properties()
        self.Update_target_properties()
        self.Update_OscTable() 
        self.Update_KanekoTable() 
        self.Update_BelkacemTable()
        self.Update_GOSTable()
        #now update also the other boxes
        self.remainder_df()
   
     
    def remainder_df(self):
        self.set_modelchoice()  
        self.set_RelativityChoice(self.calc. Dispersion_relativistic)
        self.set_rb_GOSrescaling(self.calc.ApplySumRuleToGOS)  
        self.set_DispersionChoice(self.calc.Dispersion_choice) 
        self.set_DispersionChoice_based_on_button(self.calc.delayed_dispersion)
        self.set_Gamma_dispersion(self.calc.Add_Doppler_Width)
        self.set_addChi_ELF(self.calc.AddELF) 
        self.set_RPAchoice(self.calc.Merminize)   
        self.set_maxEnergyDensityEffect(self.calc.maxEnergyDensityEffect)
        self.set_massunitcell(self.calc.massunitcell)
        self.set_TargetDensity(self.calc.specificweight) 
        self.show_DFcomment()
 
                          
#=============================== add chi or elf======================
    def addchi_elf_init(self):
        addChiElfBox = QGroupBox(parent=self.tab_1, title="all contributions:")
        addChiElfBox.setGeometry(QRect(895,653, 155,45))
        AddChiAddElf = {0: "Add chi", 1: "add ELF"}

        self.rb_addChi_ELF = []
        for i in sorted(AddChiAddElf.keys()):
            label = "%s" % (AddChiAddElf[i])
            rb=QRadioButton(label, parent= addChiElfBox )
            rb.setGeometry(QRect(5+75*i,20, 65,25))
            rb.clicked.connect(self.rb_AddChiAddElf_clicked)
            self.rb_addChi_ELF.append(rb)
        
    def rb_AddChiAddElf_clicked(self):
        for i in range (2):
            if  self.rb_addChi_ELF[i].isChecked(): self.calc.AddELF=i        
        
    def set_addChi_ELF (self,i): self.rb_addChi_ELF[i].setChecked(True)       
#==============RPA options=============================
    def RPA_choice_init(self):
        RPAChoiceBox = QGroupBox(parent=self.tab_1, title="RPA option:")
        RPAChoiceBox.setGeometry(QRect(325,540, 140,105))
        QLabel("(Lindhard, Kaneko)",parent=RPAChoiceBox ).setGeometry(QRect(5,15, 130,20))
        RPAchoice=["Plain","Mermin-corrected","Direct"]
        self.rb_RPAchoice = []
        for i in range(3):
            rb=QRadioButton(RPAchoice[i], parent= RPAChoiceBox )
            if i == 0:  rb.setGeometry(QRect(10,35, 100,20))
            elif i == 1:  rb.setGeometry(QRect(10,55, 150,20))
            elif i == 2:  rb.setGeometry(QRect(10,75, 70,20))
            rb.clicked.connect(self.rb_RPAchoice_clicked)
            self.rb_RPAchoice.append(rb)
       
  
    def rb_RPAchoice_clicked(self):
        for i in range (3):
            if self.rb_RPAchoice[i].isChecked():self.calc.Merminize = i #0 means Lindhard, 1 means Merminize, 2 Direct approach=i
        if   self.calc.Merminize == 0: self.rb_modelchoice[2].setText("Lindhard or LL")
        elif self.calc.Merminize == 1: self.rb_modelchoice[2].setText("Mermin or MLL")
        elif self.calc.Merminize == 2: self.rb_modelchoice[2].setText("RPA direct")
     
                                
    def set_RPAchoice(self,i):self.rb_RPAchoice[i].setChecked(True)
        
#===================================== Dispersion frame   
    def Dispersion_init(self):
        # use QButtonGroup
        DispersionBox = QGroupBox(parent=self.tab_1, title="Dispersion options:")
        DispersionBox.setGeometry(QRect(630,540, 250,105))
        quadratic_full_group=QButtonGroup(DispersionBox)
        DispersionChoice = {0: "ω(q) = ω(0) + αQ", 1:""}
        
        self.rb_DispersionChoice = []
        for i in sorted(DispersionChoice.keys()):
            label = "%s" % (DispersionChoice[i])
            rb=QRadioButton(label, parent= DispersionBox ) 
            if i == 0:  rb.setGeometry(QRect(5,25, 150,20))  
            if i == 1:  rb.setGeometry(QRect(5,53, 24,20))
            rb.clicked.connect(self.rb_DispersionChoice_clicked)
            quadratic_full_group.addButton(rb)
            self.rb_DispersionChoice.append(rb) 
        self.extendeddispersionlabel=QLabel("ω(q)² = ω(0)² + 2/3 α(v<sub>f</sub>)² Q + Q²",parent=DispersionBox)
        self.extendeddispersionlabel.setGeometry(20,53,200,20)     
        
        Dispersion_basis_group=QButtonGroup(DispersionBox)    
        self.dispersion_based_on=QLabel("dispersion based on:",parent=DispersionBox)   
        self.dispersion_based_on.setGeometry(QRect(5,20, 125,20)) 
        self.dispersion_based_on.setVisible(False)
        self.rb_DispersionChoice_based_on_button =[]
        Dispersion_based_on_list = {0: "Eₓ", 1: "√(Eₓ² − Egapₓ²)"}
        for i in sorted(Dispersion_based_on_list.keys()):
            label = "%s" % (Dispersion_based_on_list[i])
            rb=QRadioButton(label, parent= DispersionBox ) 
            if i == 0:  rb.setGeometry(QRect(5,40, 60,20))  
            if i == 1:  rb.setGeometry(QRect(55,40, 120,20))
             
            rb.clicked.connect(self.rb_Dispersion_based_on_clicked)
            Dispersion_basis_group.addButton(rb)
            self.rb_DispersionChoice_based_on_button.append(rb)
            self.rb_DispersionChoice_based_on_button[i].setVisible(False) 
        self.drawHLine(80,DispersionBox)
    
        Gamma_width_group=QButtonGroup(DispersionBox)  
        self.rb_Gamma_Dispersion_button =[]
        Gamma_Disersion_list = {0: "Γ(q)=Γ(0)", 1: ""}
        self.disperionlabel = QLabel("Γ(q)=√(Γ(0)² + (q*k<sub>f</sub>)²)",parent=DispersionBox)
        self.disperionlabel.setGeometry(QRect(108,79,160,25)) 
        self.disperionlabel.setVisible(False)
        for i in sorted(Gamma_Disersion_list.keys()):
            label = "%s" % (Gamma_Disersion_list[i])
            rb=QRadioButton(label, parent= DispersionBox ) 
            if i == 0:  rb.setGeometry(QRect(5,80,80,20))  
            if i == 1:  rb.setGeometry(QRect(90,80, 150,20))
             
            rb.clicked.connect(self.rb_Gamma_Dispersion_button_clicked)
            Gamma_width_group.addButton(rb)
            self.rb_Gamma_Dispersion_button.append(rb)
            self.rb_Gamma_Dispersion_button[i].setVisible(False) 
            
    def rb_DispersionChoice_clicked(self):
        for i in range (2):
            if self.rb_DispersionChoice[i].isChecked(): self.calc.Dispersion_choice = i #0 simple quadratic, 1 full dispersion
            
    def set_DispersionChoice(self,i): self.rb_DispersionChoice[i].setChecked(True) 
            
    def rb_Gamma_Dispersion_button_clicked(self):
        for i in range (2):
            if self.rb_Gamma_Dispersion_button[i].isChecked(): self.calc.Add_Doppler_Width = i #0  constant width, 1 add Doppler width  prop. to q.K_f
       
    def set_Gamma_dispersion(self,i):  self.rb_Gamma_Dispersion_button[i].setChecked(True)
            
    def rb_Dispersion_based_on_clicked(self):
        for i in range (2):
            if self.rb_DispersionChoice_based_on_button[i].isChecked(): 
                self.calc.delayed_dispersion = i  # 0 means dispersion relative to E_1, 1 relatve to sqrt(E_r^2-gap^2). used in Tuac and Tauc_an only 
                
    def set_DispersionChoice_based_on_button(self,i): self.rb_DispersionChoice_based_on_button[i].setChecked(True) 
        
    def Relativity_init(self):
        RelativityBox = QGroupBox(parent=self.tab_1, title="Recoil Energy Q:")
        RelativityBox.setGeometry(QRect(475,540,140,105))
        RelativityChoice = {0: "Q = q²/2", 1: "Q = √(c²q²+m²c⁴)\n       - mc²"}
        self.rb_RelativityChoice = []
        for i in sorted(RelativityChoice.keys()):
            label = "%s" % (RelativityChoice[i])
            rb=QRadioButton(label, parent= RelativityBox ) 
            if i == 0:  rb.setGeometry(QRect( 3,35, 80,20))  
            if i == 1:  rb.setGeometry(QRect(3,60, 165,45))
            rb.clicked.connect(self.rb_RelativityChoice_clicked)
            self.rb_RelativityChoice.append(rb) 
       
        
    def rb_RelativityChoice_clicked(self):
        for i in range (2):
            if self.rb_RelativityChoice [i].isChecked(): self.calc.Dispersion_choice = i #0 non-relativistic, 1 relativistic
                
    def set_RelativityChoice(self,i): self.rb_RelativityChoice[i].setChecked(True)  
        
    def zero_DF_Btn(self):
        self.zeroBtn = QPushButton("zero df", self.tab_1)  
        self.zeroBtn.setGeometry(QRect(1060,665, 50,25))
        self.zeroBtn.clicked.connect(self.calc.ZeroDF)  
       
        
                                 
#===================================== Target property box 
    def init_target_properties_box(self):    
        ypos=26
        dy=26
        yheight=20
        self.calc.get_df_properties()
        self.TargetPropertiesBox = QGroupBox(parent=self.tab_1, title="Target Properties:")
        self.TargetPropertiesBox.setGeometry(QRect(895,77, 220,350)) 
        self.DensityLabel=QLabel("",parent=self.TargetPropertiesBox )
        self.DensityLabel.setGeometry(QRect(5,ypos, 190,2*yheight))
        ypos+=2*dy
        self.wp_1e_Label=QLabel("",parent=self.TargetPropertiesBox )
        self.wp_1e_Label.setGeometry(QRect(5,ypos, 190,yheight))
        ypos+=dy
        self.Osc_Label=QLabel("",parent=self.TargetPropertiesBox )
        self.Osc_Label.setGeometry(QRect(5,ypos, 190,2*yheight))
        ypos+=2*dy
        self.SumAi_Label=QLabel("",parent=self.TargetPropertiesBox )
        self.SumAi_Label.setGeometry(QRect(5,ypos-10, 190,yheight))
        ypos+=dy
        self.GOS_Label=QLabel("",parent=self.TargetPropertiesBox )
        self.GOS_Label.setGeometry(QRect(5,ypos, 190,2*yheight))
        ypos+=2*dy
        self.Belkacem_Label=QLabel("",parent=self.TargetPropertiesBox )
        self.Belkacem_Label.setGeometry(QRect(5,ypos, 190,2*yheight))
        ypos+=2*dy
        self.Kaneko_Label=QLabel("",parent=self.TargetPropertiesBox )
        self.Kaneko_Label.setGeometry(QRect(5,ypos, 190,3*yheight))
    
    def input_target_properties(self):
        TargetBox = QGroupBox(parent=self.tab_1, title="Target:")
        TargetBox.setGeometry(QRect(895,4, 220,71)) 
        QLabel("density (gr/cm³):",parent=TargetBox ).setGeometry(QRect(5,20, 125,20))
        self.densinp= QLineEdit(str(self.calc.specificweight),parent=TargetBox)
        self.densinp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.densinp.setGeometry(QRect(130,20, 65,18))
        self.densinp.editingFinished.connect(self.UpdateDensity)
        
        QLabel("molar weight (gr):",parent=TargetBox ).setGeometry(QRect(5,45, 125,20))
        self.molweight= QLineEdit(str(self.calc.massunitcell),parent=TargetBox)
        self.molweight.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.molweight.setGeometry(QRect(130,45, 65,18))
        self.molweight.editingFinished.connect(self.Updatemolweight)
        
    def UpdateDensity(self):
        self.calc.specificweight=float(self.densinp.text())
        self.calc.get_df_properties()
        self.Update_target_properties()  
        
    def Updatemolweight(self):
        self.calc.massunitcell=float(self.molweight.text())
        self.calc.get_df_properties()
        self.Update_target_properties()    
        
    def Update_target_properties(self): 
        self.DensityLabel.setText(self.calc.UnitCellDensityText) 
        self.wp_1e_Label.setText(self.calc.w_p_1e_per_uc)
        self.Osc_Label.setText(self.calc.DF_prop_text)
        if self.calc.SumAi <= 1.0:
            self.SumAi_Label.setStyleSheet("QLabel { color : black; }")
        else:
            self.SumAi_Label.setStyleSheet("QLabel {  color : red; }")
        self.SumAi_Label.setText(self.calc.SumAiText)
        self.GOS_Label.setText(self.calc.GOStext)
        self.Belkacem_Label.setText(self.calc.Belkacemtext)
        self.Kaneko_Label.setText(self.calc.Kanekotext)
        
    def init_DF_comment(self):
        self.DF_commentBox = QGroupBox(parent=self.tab_1, title="DF comment ")
        self.DF_commentBox.setGeometry(QRect(325,650, 562,47)) 
        self.DF_filename = QLabel(self.readwrite.fname,parent=self.DF_commentBox )
        self.DF_filename.setGeometry(QRect(90,2, 400,20))
        self.DF_comment_input= QLineEdit(str(self.readwrite.DF_comment),parent=self.DF_commentBox)
        self.DF_comment_input.setGeometry(QRect(5,20, 540,18))
        self.DF_comment_input.editingFinished.connect(self.update_DFcomment)
        
    def show_DFcomment(self):
        self.DF_filename.setText(self.readwrite.fname)
        self.DF_comment_input.setText(self.readwrite.DF_comment)
        
    def update_DFcomment(self):self.readwrite.DF_comment=self.DF_comment_input.text()
    
    def transform_df(self):
        xpos=5
        ypos=0
        deltay=25
        boxwidth=220
        self.Transform_DF_Box = QGroupBox(parent=self.tab_1,title="DF transformations")
        self.Transform_DF_Box.setGeometry(QRect(895,430,boxwidth,220))  
        ypos+= deltay-5
        self.Change_DF_model = QPushButton( self.Transform_DF_Box) 
        
        self.PENN_Transfrom = QPushButton("'Penn'", self.Transform_DF_Box)  
        self.PENN_Transfrom.setGeometry(QRect(xpos,ypos, 50,23))
        self.PENN_Transfrom.clicked.connect(self.calc.Penn_from_ELF)   

        self.Change_DF_model.setText("Memin to Drude")
        self.Change_DF_model.setGeometry(QRect(xpos+100,ypos, 110,22))
        self.Change_DF_model.clicked.connect(self.calc.convert_DF)  
        self.Change_DF_model.setVisible(True)
        
        
        
        ypos+=  deltay
        self.drawHLine(ypos,self.Transform_DF_Box)
        ypos+= +5
        QLabel("using Oscillator" ,parent=self.Transform_DF_Box).setGeometry(QRect(xpos,ypos,100,18))  
        self.first_osc_used = QLineEdit(str(self.calc.First_oscillator_transform+1),parent= self.Transform_DF_Box ) 
        self.first_osc_used.setValidator(QIntValidator(1,self.calc.maxOscillators))
        self.first_osc_used.setGeometry(QRect(100,ypos, 25,18))  
        self.first_osc_used.editingFinished.connect(self.Update_first_osc_used)
        QLabel("to" ,parent=self.Transform_DF_Box).setGeometry(QRect(xpos+130,ypos,120,18))  
        self.last_osc_used = QLineEdit(str(self.calc.Last_oscillator_transform+1),parent= self.Transform_DF_Box ) 
        self.last_osc_used.setValidator(QIntValidator(1,self.calc.maxOscillators))
        self.last_osc_used.setGeometry(QRect(160,ypos, 25,18))  
        self.last_osc_used.editingFinished.connect(self.Update_last_osc_used)
        
        ypos+= deltay
        self.DF_from_OOS = QPushButton("DF from OOS", self.Transform_DF_Box)  
        self.DF_from_OOS.setGeometry(QRect(xpos,ypos, 100,22))
        self.DF_from_OOS.clicked.connect(self.readwrite.DF_from_OOS)    
        
   
      
        
        ypos +=deltay
        self.drawHLine(ypos,self.Transform_DF_Box)
        ypos+= +5
       
        self.Change_U = QPushButton( self.Transform_DF_Box)  
        self.Change_U.setText("change U")
        self.Change_U.setGeometry(QRect(xpos,ypos, 90,22))
        self.Change_U.clicked.connect(self.calc.Change_U)  
        ypos+=  deltay
        self.ulabel=QLabel("U factor (0< x <0.99)",parent=self.Transform_DF_Box)
        self.ulabel.setGeometry(QRect(xpos,ypos,140,18)) 
        self.U_inp = QLineEdit(str(self.calc.U_factor),parent= self.Transform_DF_Box ) 
        self.U_inp.setValidator(QDoubleValidator(0.0,1.0e8,6))
        self.U_inp.setGeometry(QRect(130,ypos, 30,18))  
        self.U_inp.editingFinished.connect(self.Update_U_factor)  

        ypos += deltay
        self.drawHLine(ypos,self.Transform_DF_Box)
        ypos+= +5
        self.DF_to_Kaneko = QPushButton("DL/Mermin to Kaneko", self.Transform_DF_Box)  
        self.DF_to_Kaneko.setGeometry(QRect(xpos,ypos, 150,25))
        self.DF_to_Kaneko.clicked.connect(self.calc.convert_DL_to_Kaneko)  
        ypos+=deltay 
        self.l_label=QLabel("l = ",parent=self.Transform_DF_Box)
        self.l_label.setGeometry(QRect(xpos,ypos,50,18)) 
        self.l_inp = QLineEdit(str(self.calc.l_Kaneko_transform),parent= self.Transform_DF_Box ) 
        self.l_inp.setValidator(QIntValidator(0,3))
        self.l_inp.setGeometry(QRect(35,ypos, 20,18))  
        self.l_inp.editingFinished.connect(self.Update_l_Kaneko_transform)  
        self.Q_label=QLabel(", Q = ",parent=self.Transform_DF_Box)
        self.Q_label.setGeometry(QRect(xpos+55,ypos,50,18)) 
        self.Q_inp = QLineEdit(str(self.calc.Q_Kaneko_transform),parent= self.Transform_DF_Box ) 
        self.Q_inp.setValidator(QDoubleValidator(0.0,1000.0,6))
        self.Q_inp.setGeometry(QRect(95,ypos, 40,18))  
        self.Q_inp.editingFinished.connect(self.Update_Q_Kaneko_transform)  
        self.au_label=QLabel("a.u.",parent=self.Transform_DF_Box)
        self.au_label.setGeometry(QRect(xpos+135,ypos,50,18)) 
        
        
  
    def Update_first_osc_used(self):self.calc.First_oscillator_transform = int(self.first_osc_used.text())-1  # first oscillator is at index 0
    def Update_last_osc_used(self):self.calc.Last_oscillator_transform = int(self.last_osc_used.text())-1  # first oscillator is at index 0
   
    def Update_l_Kaneko_transform(self):self.calc.l_Kaneko_transform = int(self.l_inp.text())
    def Update_Q_Kaneko_transform(self):self.calc.Q_Kaneko_transform = float(self.Q_inp.text())     
    
    def Update_U_factor(self):
        self.calc.U_factor = float(self.U_inp.text())  
        
                 
#start tab2========================================== calculations tab======================
#-main start calculationbox
    def calculate_which_property(self):
        self.CalculationChoiseBox = QGroupBox(parent=self.tab_2,title="Calculatiobn Choice:")
        self.CalculationChoiseBox.setGeometry(QRect(4,4,430,695))
        self.CalcChoice=[]
        self.yCoord=[]
        self.CalcChoice.append({"desc": "<font color='red'>Dielectric Function:", "subroutine_called":"label", "pos": 7})
        self.CalcChoice.append({"desc": "Re[ϵ(ω)],Im[ϵ(ω)]", "subroutine_called":"eps_w", "pos": 0})
        self.CalcChoice.append({"desc": "Re[1/ϵ(ω)], Im[-1/ϵ(ω)]", "subroutine_called":"one_over_eps_w", "pos": 0})
        self.CalcChoice.append({"desc": "Re[ϵ(q)],Im[ϵ(q)]", "subroutine_called":"eps_q", "pos": 0})
        self.CalcChoice.append({"desc": "Re[1/ϵ(q)], Im[-1/ϵ(q)]", "subroutine_called":"one_over_eps_q", "pos": 0})
        self.CalcChoice.append({"desc": "n, k", "subroutine_called":"n_k", "pos": 0})
        self.CalcChoice.append({"desc": "<font color='red'>Kramers-Kronig tests:", "subroutine_called":"label", "pos": 7})

        self.CalcChoice.append({"desc": "Re[ϵ(ω)],Im[ϵ(ω)]", "subroutine_called":"eps_w_kk", "pos": 0})
        self.CalcChoice.append({"desc": "Re[1/ϵ(ω)], Im[-1/ϵ(ω)]", "subroutine_called":"one_over_eps_w_kk", "pos": 0})
        self.CalcChoice.append({"desc": "n, k", "subroutine_called":"n_k_kk", "pos": 0})
        
        self.CalcChoice.append({"desc": "<font color='red'>Others:", "subroutine_called":"label", "pos": 7})
        self.CalcChoice.append({"desc": "color plot  Im[-1/ϵ(q,ω)]", "subroutine_called":"eq_plot", "pos": 0})
        self.CalcChoice.append({"desc": "Oscillator strength", "subroutine_called":"oscillator_strength", "pos": 0})
        self.CalcChoice.append({"desc": "Sum rules", "subroutine_called":"sum_rules", "pos": 0})
        self.CalcChoice.append({"desc": "Inertial sum rules", "subroutine_called":"inertial_rules", "pos": 0})
        self.CalcChoice.append({"desc": "Mean excitation energy", "subroutine_called":"Mean_Excitation_Energy", "pos": 0})
        self.CalcChoice.append({"desc": "S(q,ω) and its sum rule", "subroutine_called":"S_k_omega_rule", "pos": 0})
        self.CalcChoice.append({"desc": "Compton profile at q=            a.u.", "subroutine_called":"Compton", "pos": 0})
        self.CalcChoice.append({"desc": "Fresnel refl. coef., ℏω=          eV", "subroutine_called":"Fresnel_at_E", "pos": 0})
        self.CalcChoice.append({"desc": "Fresnel refl. coef., ϕ=             °", "subroutine_called":"Fresnel_at_angle", "pos": 0})
        self.CalcChoice.append({"desc": "X-ray mass absorption coef.", "subroutine_called":"xray_absorption", "pos": 0})
        self.CalcChoice.append({"desc": "<font color='red'>Local density approx.:", "subroutine_called":"label", "pos": 5})
        self.CalcChoice.append({"desc": "pseudo charge-density", "subroutine_called":"pseudo_charge_density", "pos": 0})
        self.CalcChoice.append({"desc": "radial pseudo charge dens.", "subroutine_called":"radial_charge_density", "pos": 0})
        self.CalcChoice.append({"desc": "ωₚ(r), λ(r)", "subroutine_called":"radial_w_p_lambda", "pos": 0})
        
        self.CalcChoice.append({"desc": "<font color='red'>Quantities involving", "subroutine_called":"label", "pos": 12})
        self.CalcChoice.append({"desc": "<font color='red'>q integration", "subroutine_called":"label", "pos": 7})
        self.CalcChoice.append({"desc": "diimfp", "subroutine_called":"diimfp", "pos": 0})
        self.CalcChoice.append({"desc": "Partial diimfp", "subroutine_called":"partial_DIIMFP", "pos": 0})
        self.CalcChoice.append({"desc": "Partial stopping", "subroutine_called":"partial_stopping", "pos": 0})
        self.CalcChoice.append({"desc": "λ(E₀), dE/dx(E₀),  dE²/dx(E₀)", "subroutine_called":"IMFP_stop_strag", "pos": 0})
        self.CalcChoice.append({"desc": "L (stopping number)", "subroutine_called":"shell_effect_all", "pos": 0})
        self.CalcChoice.append({"desc": "CSDA range", "subroutine_called":"calc_projectile_range", "pos": 0})
        self.CalcChoice.append({"desc": "energy deposition vs depth", "subroutine_called":"Energy_Deposition_Depth", "pos": 0})
        
        self.CalcChoice.append({"desc": "<font color='red'>Theta plots", "subroutine_called":"label", "pos": 7})
        self.CalcChoice.append({"desc": "dσ/(dωdΩ) at ω=                 eV)", "subroutine_called":"DDCS_at_omega", "pos": 0})
        self.CalcChoice.append({"desc": "dσ/(dωdΩ) at θ=                 mrad)", "subroutine_called":"DDCS_at_theta", "pos": 0})
        self.CalcChoice.append({"desc": "dcs  dσ/(dΩ)", "subroutine_called":"dcs_plot", "pos": 0})
        self.CalcChoice.append({"desc": "dσ/(dωdΩ), no retardation", "subroutine_called":"dcs_omega_eq_plot", "pos": 0})
        self.CalcChoice.append({"desc": "dσ/(dωdΩ), incl. retardation", "subroutine_called":"Cerenkov", "pos": 0})
        self.CalcChoice.append({"desc": "retarded - non-retarded", "subroutine_called":"difference_due_to_Cerenkov", "pos": 0})
        
        self.CalcChoice.append({"desc": "<font color='red'>Electron specific:", "subroutine_called":"label", "pos": 5})
        self.CalcChoice.append({"desc": "SELF, DSEP", "subroutine_called":"SELF_DSEP", "pos": 0})
        self.CalcChoice.append({"desc": "color plot SELF", "subroutine_called":"self_eq_plot", "pos": 0})
        self.CalcChoice.append({"desc": "(R)EELS spectrum", "subroutine_called":"REELS", "pos": 0})        
        xpos=5
        self.x_2nd_column=217
        ypos=25
        deltay=25
        self.rb_calculation_Choice = []
        for i in range(len(self.CalcChoice)):
            if(self.CalcChoice[i]["subroutine_called"] == "label"):
                QLabel(self.CalcChoice[i]["desc"],parent=self.CalculationChoiseBox).setGeometry(QRect(xpos+20,ypos+self.CalcChoice[i]["pos"], 200,18))
            else:    
                rb=QRadioButton(self.CalcChoice[i]["desc"], parent= self.CalculationChoiseBox)
                rb.subroutine_called=self.CalcChoice[i]["subroutine_called"]
                rb.toggled.connect(self.onClicked)
                if i== 2:rb.setChecked(True)  
                rb.setGeometry(QRect(xpos,ypos, 208,20))
                self.rb_calculation_Choice.append(rb)
            self.yCoord.append(ypos)
            if self.CalcChoice[i]["subroutine_called"] == "IMFP_stop_strag": 
                self.index_imfp_calc=len(self.rb_calculation_Choice)-1
            if self.CalcChoice[i]["subroutine_called"] == "radial_w_p_lambda": # start second column
                xpos=self.x_2nd_column
                ypos=25
            else:    
                ypos+=deltay
        self.LineEdit_Q_Compton() 
        self.LineEdit_E_Fresnel() 
        self.LineEdit_phi_ellipsometry()     
        self.LineEdit_ddcs_at_w() 
        self.LineEdit_ddcs_at_theta() 
             
    def onClicked(self,checked):
        if checked:
            self.plotchoice= self.sender().subroutine_called
       
    def LineEdit_Q_Compton(self):
        #get index
        index= self.CalcChoice.index(next(filter(lambda n: n.get('subroutine_called') == "Compton", self.CalcChoice)))
        self.q_compton_inp = QLineEdit(str(self.calc.q_Compton),parent= self.CalculationChoiseBox) 
        self.q_compton_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.q_compton_inp.setGeometry(QRect(150,self.yCoord[index], 35,18))  
        self.q_compton_inp.editingFinished.connect(self.UpdateQ_Compton) 
          
    def UpdateQ_Compton(self): self.calc.q_Compton=float(self.q_compton_inp.text())      

    def LineEdit_E_Fresnel(self):
        index= self.CalcChoice.index(next(filter(lambda n: n.get('subroutine_called') == "Fresnel_at_E", self.CalcChoice)))
        self.E_Fresnel_inp = QLineEdit(str(self.calc.E_Fresnel),parent= self.CalculationChoiseBox) 
        self.E_Fresnel_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.E_Fresnel_inp.setGeometry(QRect(158,self.yCoord[index], 30,18))  
        self.E_Fresnel_inp.editingFinished.connect(self.UpdateE_Fresnel)        
   
    def UpdateE_Fresnel(self): self.calc.E_Fresnel=float(self.E_Fresnel_inp.text()) 
     
    def LineEdit_phi_ellipsometry(self):
        index= self.CalcChoice.index(next(filter(lambda n: n.get('subroutine_called') == "Fresnel_at_angle", self.CalcChoice)))
        self.phi_ellipsometry_inp = QLineEdit(str(self.calc.phi_ellipsometry),parent= self.CalculationChoiseBox) 
        self.phi_ellipsometry_inp.setValidator(QDoubleValidator(0, 90.0, 6))
        self.phi_ellipsometry_inp.setGeometry(QRect(158,self.yCoord[index], 30,18))  
        self.phi_ellipsometry_inp.editingFinished.connect(self.Update_phi_ellipsometry)        
   
    def Update_phi_ellipsometry(self): self.calc.phi_ellipsometry=float(self.phi_ellipsometry_inp.text())  
           
    def LineEdit_ddcs_at_w(self):
        index= self.CalcChoice.index(next(filter(lambda n: n.get('subroutine_called') == "DDCS_at_omega", self.CalcChoice)))
        self.ddcs_at_w_inp = QLineEdit(str(self.calc.omega_ddcs),parent= self.CalculationChoiseBox) 
        self.ddcs_at_w_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.ddcs_at_w_inp.setGeometry(QRect(self.x_2nd_column+120,self.yCoord[index]+4, 40,18))  
        self.ddcs_at_w_inp.editingFinished.connect(self.Update_ddcs_at_w)   
        
    def Update_ddcs_at_w(self):       
        self.calc.omega_ddcs=float(self.ddcs_at_w_inp.text())
        
    def LineEdit_ddcs_at_theta(self):
        index= self.CalcChoice.index(next(filter(lambda n: n.get('subroutine_called') == "DDCS_at_theta", self.CalcChoice)))
        self.ddcs_at_theta_inp = QLineEdit(str(self.calc.theta_ddcs),parent= self.CalculationChoiseBox) 
        self.ddcs_at_theta_inp.setValidator(QDoubleValidator(0.0, 3141, 6))
        self.ddcs_at_theta_inp.setGeometry(QRect(self.x_2nd_column+120,self.yCoord[index]+4, 40,18))  
        self.ddcs_at_theta_inp.editingFinished.connect(self.Update_ddcs_at_theta)   
        
    def Update_ddcs_at_theta(self): self.calc.theta_ddcs=float(self.ddcs_at_theta_inp.text())

#-----------------------------plot option boxes-------------------------------------
    def  PlotOptions(self):
        xpos=5
        ypos=17
        deltay=25
        self.PlotOptionBox = QGroupBox(parent=self.tab_2,title="Plot Options:")
        self.PlotOptionBox.setGeometry(QRect(440,4,320,305))
        QLabel("<font color='red'>Energy plots:",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 300,18))
        QLabel("at q=               a.u.",parent=self.PlotOptionBox).setGeometry(QRect(xpos+130,ypos, 300,18))        
        self.q_Eplot_inp = QLineEdit(str(self.calc.q),parent= self.PlotOptionBox ) 
        self.q_Eplot_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.q_Eplot_inp.setGeometry(QRect(xpos+130+32,ypos, 45,17))  
        self.q_Eplot_inp.editingFinished.connect(self.Update_q_Eplot)
           
        ypos+=deltay
        QLabel("from ",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 50,18))
        QLabel("eV to ",parent=self.PlotOptionBox).setGeometry(QRect(xpos+85,ypos, 50,18))
        QLabel("eV,  step",parent=self.PlotOptionBox).setGeometry(QRect(xpos+170,ypos, 50,18))
        QLabel("eV",parent=self.PlotOptionBox).setGeometry(QRect(xpos+270,ypos, 50,18))
        
        self.LowerElim_inp = QLineEdit(str(self.calc.LowerELimit),parent= self.PlotOptionBox ) 
        self.LowerElim_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.LowerElim_inp.setGeometry(QRect(xpos+32,ypos, 45,18))  
        self.LowerElim_inp.editingFinished.connect(self.UpdateLowerElim)   
        
        self.UpperElim_inp = QLineEdit(str(self.calc.UpperELimit),parent= self.PlotOptionBox ) 
        self.UpperElim_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.UpperElim_inp.setGeometry(QRect(xpos+125,ypos, 45,18))  
        self.UpperElim_inp.editingFinished.connect(self.UpdateUpperElim)  
        
        self.StepSize_inp = QLineEdit(str(self.calc.Stepsize),parent= self.PlotOptionBox ) 
        self.StepSize_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.StepSize_inp.setGeometry(QRect(xpos+225,ypos, 40,18))  
        self.StepSize_inp.editingFinished.connect(self.UpdateStepSize)
        
        ypos+=deltay   
        self.drawHLine(ypos, self.PlotOptionBox) 
        ypos+=5
        QLabel("<font color='red'>Momentum plots:",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 300,18))
        QLabel("at ω=               eV",parent=self.PlotOptionBox).setGeometry(QRect(xpos+130,ypos, 300,18))
        self.omega_qplot_inp = QLineEdit(str(self.calc.Energy_qplot ),parent= self.PlotOptionBox ) 
        self.omega_qplot_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.omega_qplot_inp.setGeometry(QRect(xpos+162,ypos, 45,18))  
        self.omega_qplot_inp.editingFinished.connect(self.Update_omega_qplot)  
        
        ypos+=deltay
        QLabel("from 0 a.u. to q= ",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 100,18))
        QLabel("a.u., step             a.u.",parent=self.PlotOptionBox).setGeometry(QRect(xpos+155,ypos,100,18))
        QLabel("a.u.",parent=self.PlotOptionBox).setGeometry(QRect(xpos+265,ypos,100,18))
        self.Upperqlim_inp = QLineEdit(str(self.calc.UpperqLimit),parent= self.PlotOptionBox ) 
        self.Upperqlim_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.Upperqlim_inp.setGeometry(QRect(xpos+105,ypos, 40,18))  
        self.Upperqlim_inp.editingFinished.connect(self.UpdateUpperqlim) 
          
        self.qStepSize_inp = QLineEdit(str(self.calc.Stepsize_qplot ),parent= self.PlotOptionBox ) 
        self.qStepSize_inp.setValidator(QDoubleValidator(0.0, 1e9, 6))
        self.qStepSize_inp.setGeometry(QRect(xpos+220,ypos, 40,18))  
        self.qStepSize_inp.editingFinished.connect(self.UpdateqStepSize)  
        
        ypos+=deltay-2
        self.drawHLine(ypos, self.PlotOptionBox) 
        ypos+=3
        QLabel("<font color='red'>Theta plots:",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 300,18))
        ypos+=deltay
        QLabel("from 0 to θ= ",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 100,18))
        self.theta_max_inp = QLineEdit(str(self.calc.theta_max),parent= self.PlotOptionBox) 
        self.theta_max_inp.setValidator(QDoubleValidator(0.0, 3.1415*1000, 6))
        self.theta_max_inp.setGeometry(QRect(xpos+80,ypos, 40,18))  
        self.theta_max_inp.editingFinished.connect(self.Update_theta_max)   
        QLabel("mrad in               steps",parent=self.PlotOptionBox).setGeometry(QRect(xpos+125,ypos, 200,18))
        self.NThetaStep_inp = QLineEdit(str(self.calc.NThetaStep),parent= self.PlotOptionBox) 
        self.NThetaStep_inp.setValidator(QIntValidator(0, 1000000))
        self.NThetaStep_inp.setGeometry(QRect(xpos+180,ypos, 40,18))  
        self.NThetaStep_inp.editingFinished.connect(self.Update_NThetaStep)   
            
       
        ypos+=deltay-2
        self.drawHLine(ypos,self.PlotOptionBox)
        ypos+=+3
        
        QLabel("Energy scale:",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 310,18))
        self.rb_energyScale = []
        rb =QRadioButton("eV", parent= self.PlotOptionBox )
        rb.setGeometry(QRect(xpos+100, ypos,100,20))
        rb.clicked.connect(self.rb_energyScale_clicked)
        self.rb_energyScale.append(rb)
        
        rb =QRadioButton("nm", parent= self.PlotOptionBox )
        rb.setGeometry(QRect(xpos+170, ypos,50,20))
        rb.clicked.connect(self.rb_energyScale_clicked)
        self.rb_energyScale.append(rb)
        
        rb =QRadioButton("cm⁻¹", parent= self.PlotOptionBox )
        rb.setGeometry(QRect(xpos+240, ypos,50,20))
        rb.clicked.connect(self.rb_energyScale_clicked)
        self.rb_energyScale.append(rb)
        
        
        self.rb_energyScale[self.calc.Energy_Scale_choice].setChecked(True) 
        
        ypos+=deltay
        self.LogX_Checkbox = QCheckBox('Log scale x', self.PlotOptionBox )
        self.LogX_Checkbox.setChecked(self.calc.LogX)
        self.LogX_Checkbox.setGeometry(QRect(xpos,ypos, 110,18))
        self.LogX_Checkbox.stateChanged.connect(self.LogX_Checkbox_state_changed)
        
        self.LogY_Checkbox = QCheckBox('Log scale y', self.PlotOptionBox )
        self.LogY_Checkbox.setChecked(self.calc.LogY)
        self.LogY_Checkbox.setGeometry(QRect(xpos+150,ypos, 110,18))
        self.LogY_Checkbox.stateChanged.connect(self.LogY_Checkbox_state_changed)
        
        ypos+=deltay
        self.LogXY_Checkbox = QCheckBox('Log scale color plots, range 10ˣ,      x=', self.PlotOptionBox )
        self.LogXY_Checkbox.setChecked(self.calc.LogXY)
        self.LogXY_Checkbox.setGeometry(QRect(xpos,ypos, 270,18))
        self.LogXY_Checkbox.stateChanged.connect(self.LogXY_Checkbox_state_changed)
        
        self.LogRange_inp = QLineEdit(str(self.calc.log_range),parent= self.PlotOptionBox ) 
        self.LogRange_inp.setValidator(QDoubleValidator(0.0, 100, 6))
        self.LogRange_inp.setGeometry(QRect(xpos+260,ypos, 20,18))  
        self.LogRange_inp.editingFinished.connect(self.UpdateLogRange)
         
        ypos+=deltay  
        QLabel("Max. color plot (autoscale=0)",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 310,18))
        self.Max_eq_inp = QLineEdit(str(self.calc.max_eq),parent= self.PlotOptionBox ) 
        self.Max_eq_inp.setValidator(QDoubleValidator(0, 1e9, 6))
        self.Max_eq_inp.setGeometry(QRect(xpos+220,ypos, 60,18))  
        self.Max_eq_inp.editingFinished.connect(self.UpdateMax_eq)
        
        ypos+=deltay 
        QLabel("ε or χ:",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos, 40,18))
        self.rb_epsilon_chi = []
        rb =QRadioButton("ε", parent= self.PlotOptionBox )
        rb.setGeometry(QRect(xpos+45, ypos,100,20))
        rb.clicked.connect(self.rb_epsilon_chi_clicked)
        self.rb_epsilon_chi.append(rb)
        rb =QRadioButton("χ", parent= self.PlotOptionBox )
        rb.setGeometry(QRect(xpos+80, ypos,50,20))
        rb.clicked.connect(self.rb_epsilon_chi_clicked)
        self.rb_epsilon_chi.append(rb)
        self.rb_epsilon_chi[self.calc.epsilon_chi_choice].setChecked(True) 

        
    def  UpdateLowerElim(self): self.calc.LowerELimit=float(self.LowerElim_inp.text())      
    def  UpdateUpperElim(self): self.calc.UpperELimit=float(self.UpperElim_inp.text()) 
       
    def  UpdateStepSize(self):  self.calc.Stepsize=float(self.StepSize_inp.text())         
    def  Update_q_Eplot(self):  self.calc.q=float(self.q_Eplot_inp.text())  
    def  UpdateUpperqlim(self): self.calc.UpperqLimit=float(self.Upperqlim_inp.text())    
    def  Update_theta_max(self): self.calc.theta_max=float(self.theta_max_inp.text()) 
    def  Update_NThetaStep(self): self.calc.NThetaStep=int(self.NThetaStep_inp.text())
    def  UpdateqStepSize(self): self.calc.Stepsize_qplot=float(self.qStepSize_inp.text())    
    def  Update_omega_qplot(self): self.calc.Energy_qplot=float(self.omega_qplot_inp.text()) 
    def  LogX_Checkbox_state_changed(self): self.calc.LogX = self.LogX_Checkbox.isChecked()
    def  LogY_Checkbox_state_changed(self): self.calc.LogY = self.LogY_Checkbox.isChecked()
    def  LogXY_Checkbox_state_changed(self): self.calc.LogXY = self.LogXY_Checkbox.isChecked()
    def  UpdateLogRange(self):  self.calc.log_range=round(float(self.LogRange_inp.text())) 
    def  UpdateMax_eq(self):  self.calc.max_eq=float(self.Max_eq_inp.text())          
    
    def  rb_energyScale_clicked(self):
        for i in range (3): 
            if self.rb_energyScale[i].isChecked(): self.calc.Energy_Scale_choice = i 
    def  rb_epsilon_chi_clicked(self):
        for i in range (2): 
            if self.rb_epsilon_chi[i].isChecked(): self.calc.epsilon_chi_choice = i         
                         
    def REELSOptions(self):
        xpos=5
        ypos=17
        deltay=25
        self.REELSOptionBox = QGroupBox(parent=self.tab_2,title="(R)EELS parameters:")
        self.REELSOptionBox.setGeometry(QRect(440,315,320,220))  
        QLabel("Energy res. (eV):",parent=self.REELSOptionBox).setGeometry(QRect(xpos,ypos, 1200,18))      
        self.Eres_inp = QLineEdit(str(self.calc.Eres),parent= self.REELSOptionBox ) 
        self.Eres_inp.setValidator(QDoubleValidator(0, 1e9, 6))
        self.Eres_inp.setGeometry(QRect(xpos+110,ypos, 30,18))  
        self.Eres_inp.editingFinished.connect(self.UpdateEres)
        
        ypos+=deltay
        QLabel("θ₀:",parent=self.REELSOptionBox).setGeometry(QRect(xpos,ypos, 50,18))
        QLabel("θ1:",parent=self.REELSOptionBox).setGeometry(QRect(xpos+100,ypos, 50,18))
        
        self.Theta_in_inp = QLineEdit(str(self.calc.thetaIn),parent= self.REELSOptionBox ) 
        self.Theta_in_inp.setValidator(QDoubleValidator(-90, 90, 6))
        self.Theta_in_inp.setGeometry(QRect(30,ypos, 40,18))  
        self.Theta_in_inp.editingFinished.connect(self.UpdateTheta_in)
        
        self.Theta_out_outp = QLineEdit(str(self.calc.thetaOut),parent= self.REELSOptionBox ) 
        self.Theta_out_outp.setValidator(QDoubleValidator(-90, 90, 6))
        self.Theta_out_outp.setGeometry(QRect(130,ypos, 40,18))  
        self.Theta_out_outp.editingFinished.connect(self.UpdateTheta_out)
        
        ypos+=deltay
        QLabel("Fraction DIIMFP in energy loss range:",parent=self.REELSOptionBox).setGeometry(QRect(xpos,ypos,250,18))
        self.Fraction_DIIMFP_inp = QLineEdit(str(self.calc.fraction_DIIMFP),parent= self.REELSOptionBox ) 
        self.Fraction_DIIMFP_inp.setValidator(QDoubleValidator(0, 1.0, 6))
        self.Fraction_DIIMFP_inp.setGeometry(QRect(250,ypos, 40,18))  
        self.Fraction_DIIMFP_inp.editingFinished.connect(self.UpdateFraction_DIIMFP)
        
        ypos+=deltay
        QLabel("Surface plasmon scaling factor:",parent=self.REELSOptionBox).setGeometry(QRect(xpos,ypos,250,18))
        self.surf_pl_scaling_inp = QLineEdit(str(self.calc.surf_ex_factor),parent= self.REELSOptionBox ) 
        self.surf_pl_scaling_inp.setValidator(QDoubleValidator(0, 1.0, 6))
        self.surf_pl_scaling_inp.setGeometry(QRect(250,ypos, 40,18))  
        self.surf_pl_scaling_inp.editingFinished.connect(self.Update_surf_pl_scaling)
        
        # ypos+=deltay-10
        # self_group = QButtonGroup(self.REELSOptionBox )
        # self.rb_DSEP_Choice = []
        # rb = QRadioButton("DSEP from total eps", parent= self.REELSOptionBox )
        # rb.setGeometry(QRect(xpos, ypos,140,20))
        # rb.clicked.connect(self.rb_DSEP_Choice_clicked)
        # self.rb_DSEP_Choice.append(rb)
        # self_group.addButton(self.rb_DSEP_Choice[0])
        # rb =QRadioButton("DSEP per osc.", parent= self.REELSOptionBox )
        # rb.setGeometry(QRect(xpos+180, ypos,200,20))
        # rb.clicked.connect(self.rb_DSEP_Choice_clicked)
        # self.rb_DSEP_Choice.append(rb)
        # self.rb_DSEP_Choice[self.calc.DSEP_choice].setChecked(True) 
        # self_group.addButton(self.rb_DSEP_Choice[1])
        ypos+=deltay
        self.drawHLine(ypos,self.REELSOptionBox)
        ypos+=5
        REELS_EELS_group = QButtonGroup(self.REELSOptionBox )
        self.rb_Reflection_Choice = []
        rb =QRadioButton("", parent= self.REELSOptionBox )
        rb.setGeometry(QRect(xpos,ypos-3,200,20))
        rb.clicked.connect(self.rb_Reflection_Choice_clicked)
        self.rb_Reflection_Choice.append(rb)
        REELS_EELS_group.addButton(self.rb_Reflection_Choice[0])
        
        QLabel("<font color='red'>Reflection",parent=self.REELSOptionBox).setGeometry(QRect(xpos+20,ypos,250,18))
        QLabel("Part.Int(i)=1+c₁ i + c₂ i² +c₃ i³",parent=self.REELSOptionBox).setGeometry(QRect(120,ypos, 1200,18))      
      
        ypos+=deltay
        QLabel("c₁ =",parent=self.REELSOptionBox).setGeometry(QRect(xpos,ypos,30,18))
        QLabel("c₂ =",parent=self.REELSOptionBox).setGeometry(QRect(xpos+80,ypos,30,18))
        QLabel("c₃ =",parent=self.REELSOptionBox).setGeometry(QRect(xpos+160,ypos,33,18))
        
        self.c1_inp = QLineEdit(str(self.calc.coef1),parent= self.REELSOptionBox ) 
        self.c1_inp.setValidator(QDoubleValidator(-1,1, 6))
        self.c1_inp.setGeometry(QRect(30,ypos, 40,18))  
        self.c1_inp.editingFinished.connect(self.UpdateC1_inp)
        
        self.c2_inp = QLineEdit(str(self.calc.coef2),parent= self.REELSOptionBox ) 
        self.c2_inp.setValidator(QDoubleValidator(-1,1, 6))
        self.c2_inp.setGeometry(QRect(112,ypos, 40,18))  
        self.c2_inp.editingFinished.connect(self.UpdateC2_inp)
        
        self.c3_inp = QLineEdit(str(self.calc.coef3),parent= self.REELSOptionBox ) 
        self.c3_inp.setValidator(QDoubleValidator(-1,1, 6))
        self.c3_inp.setGeometry(QRect(194,ypos, 40,18))  
        self.c3_inp.editingFinished.connect(self.UpdateC3_inp)
        
        ypos+=deltay
        self.drawHLine(ypos,self.REELSOptionBox)
        ypos +=5
        QLabel("<font color='red'>Transmission",parent=self.REELSOptionBox).setGeometry(QRect(xpos+20,ypos,250,18))
        rb2=QRadioButton("", parent= self.REELSOptionBox )
        rb2.setGeometry(QRect(xpos, ypos,200,20))
        rb2.clicked.connect(self.rb_Reflection_Choice_clicked)
        self.rb_Reflection_Choice.append(rb2)
        REELS_EELS_group.addButton(self.rb_Reflection_Choice[1])
        self.rb_Reflection_Choice[self.calc.EELS].setChecked(True) 
        
        
        QLabel("Thickness (Å):",parent=self.REELSOptionBox).setGeometry(QRect(xpos+150,ypos,250,18))
        self.Thickness_inp = QLineEdit(str(self.calc.EELS_thickness),parent= self.REELSOptionBox ) 
        self.Thickness_inp.setValidator(QDoubleValidator(0,1e9, 6))
        self.Thickness_inp.setGeometry(QRect(250,ypos, 40,18))  
        self.Thickness_inp.editingFinished.connect(self.Update_EELS_thickness)
      
    def  UpdateEres(self):  self.calc.Eres=float(self.Eres_inp.text())   
    def  UpdateTheta_in(self):  self.calc.thetaIn=float(self.Theta_in_inp.text())   
    def  UpdateTheta_out(self):  self.calc.thetaOut=float(self.Theta_out_inp.text()) 
    def  UpdateFraction_DIIMFP(self):  self.calc.fraction_DIIMFP=float(self.Fraction_DIIMFP_inp.text())      
    def  Update_surf_pl_scaling(self):   self.calc.surf_ex_factor=float(self.surf_pl_scaling_inp.text()) 
    def  UpdateC1_inp(self):  self.calc.coef1=float(self.c1_inp.text())   
    def  UpdateC2_inp(self):  self.calc.coef2=float(self.c2_inp.text())  
    def  UpdateC3_inp(self):  self.calc.coef3=float(self.c3_inp.text())    
    def  Update_EELS_thickness(self): self.calc.EELS_thickness=float(self.Thickness_inp.text())
        
      
    # def  rb_DSEP_Choice_clicked(self):
        # for i in range (2): 
            # if self.rb_DSEP_Choice[i].isChecked(): self.calc.DSEP_Choice = i 
            
    def  rb_Reflection_Choice_clicked(self):
        for i in range (2): 
            if self.rb_Reflection_Choice[i].isChecked(): self.calc.EELS = i    
            
    def comparison_data_box(self):   
        xpos=5
        ypos=25
        deltay=26
        self.ComparisonBox = QGroupBox(parent=self.tab_2,title="Compare with Literature")   
        self.ComparisonBox.setGeometry(QRect(440,540,320,159)) 
        self.CompPlotBtn = QPushButton("Load comp. data",parent=self.ComparisonBox ) 
        self.CompPlotBtn.setGeometry(xpos,ypos,120,25)
        self.CompPlotBtn.clicked.connect(self.readwrite.load_comp_data) 
        
        self.PlotCompDatabox = QCheckBox('plot comp. data', self.ComparisonBox)
        self.PlotCompDatabox.setChecked(False)
        self.PlotCompDatabox.setGeometry(QRect(xpos+150,ypos+2, 150,18))
        self.PlotCompDatabox.stateChanged.connect(self.PlotCompDatabox_state_changed)
        
        ypos+=deltay+10
        QLabel("legend comp. data:" ,parent=self.ComparisonBox).setGeometry(QRect(xpos,ypos,250,18))
        self.Legend_inp = QLineEdit(self.runplot.LiteratureDescription,parent= self.ComparisonBox ) 
        self.Legend_inp.setGeometry(QRect(122,ypos, 160,18))  
        self.Legend_inp.editingFinished.connect(self.Update_Legend_inp)
        
        ypos+=deltay+10
        QLabel("Plot with:",parent=self.ComparisonBox).setGeometry(QRect(xpos,ypos,125,18))
        LabelList=["IMFP", "Stopping", "Straggling"]
        self.rb_3Plotchoice = []
        for i in range(3):
            rb=QRadioButton(LabelList[i], parent= self.ComparisonBox )
            rb.setGeometry(QRect(60+xpos+80*i,ypos, 100,20))
            rb.clicked.connect(self.rb_3Plotchoice_clicked)
            self.rb_3Plotchoice.append(rb)
        self.rb_3Plotchoice[self.runplot.comp_option_choice].setChecked(True) 
        
    def PlotCompDatabox_state_changed(self):self.runplot.Overplot=  self.PlotCompDatabox.isChecked()                              
    def Update_Legend_inp(self):self.runplot.LiteratureDescription=self.Legend_inp.text()
    def rb_3Plotchoice_clicked(self):
         for i in range (3):
            if self.rb_3Plotchoice[i].isChecked():self.runplot.comp_option_choice = i    
            
            
    def  Projectile_target_Interaction(self):
        xpos=5
        ypos=27
        deltay=25
        self.ProjectileBox = QGroupBox(parent=self.tab_2,title="Projectile-Target interactions")
        self.ProjectileBox.setGeometry(QRect(770,4,350,440))   
        
        QLabel("Projectile:",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,125,18))
        self.rb_ProjectileChoice = []
        LabelList=["e⁻", "H⁺"] 
        self.particlegroup= QButtonGroup(self.ProjectileBox)
        for i in range(2):
            rb=QRadioButton(LabelList[i], parent= self.ProjectileBox )
            rb.setGeometry(QRect(100+xpos+80*i,ypos,40,20))
            rb.clicked.connect(self.rb_ProjectileChoice_clicked)
            self.rb_ProjectileChoice.append(rb)
            self.particlegroup.addButton(rb)
            
        if self.calc.particle=="electron":
            self.rb_ProjectileChoice[0].setChecked(True)
        else:   
            self.rb_ProjectileChoice[1].setChecked(True)
            
        ypos += deltay
        QLabel("E₀ (keV):",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,125,18))
        self.E0_inp = QLineEdit(str(self.calc.E0),parent= self.ProjectileBox ) 
        self.E0_inp.setValidator(QDoubleValidator(0,1e9, 6))
        self.E0_inp.setGeometry(QRect(70,ypos, 50,18))  
        self.E0_inp.editingFinished.connect(self.Update_E0)
        self.E0_velocityLabel=QLabel(self.calc.myvelocitytext,parent=self.ProjectileBox)
        self.E0_velocityLabel.setGeometry(QRect(xpos+130,ypos,150,18))
        
        ypos +=deltay
        self.drawHLine(ypos,self.ProjectileBox)
            
        ypos += 5
        self.MottCheckbox = QCheckBox("H⁺: incl. Mott correction", self.ProjectileBox)
        self.MottCheckbox.setChecked(self.calc.MottCorrection)
        self.MottCheckbox.setGeometry(QRect(xpos,ypos,200,18))
        self.MottCheckbox.stateChanged.connect(self.MottCheckbox_state_changed)
        
        ypos +=deltay
        self.RadLosCheckbox = QCheckBox("e⁻: incl. estimates radiative losses", self.ProjectileBox)
        self.RadLosCheckbox.setChecked(self.calc.RadiativeLosses)
        self.RadLosCheckbox.setGeometry(QRect(xpos,ypos,300,18))
        self.RadLosCheckbox.stateChanged.connect(self.RadLosCheckbox_state_changed)
        
        ypos += deltay
        self.rb_ExchangeChoice = []
        exchangegroup=QButtonGroup(self.ProjectileBox)
        QLabel("e⁻:",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,175,18))
        rb=QRadioButton("no exchange", parent= self.ProjectileBox )
        rb.setGeometry(QRect(xpos+30,ypos,150,20))
        rb.clicked.connect(self.rb_ExchangeChoice_clicked)
        self.rb_ExchangeChoice.append(rb)
        exchangegroup.addButton(rb)
        rb=QRadioButton("SBethe approach", parent= self.ProjectileBox )
        rb.setGeometry(QRect(xpos+180,ypos,150,20))
        rb.clicked.connect(self.rb_ExchangeChoice_clicked)
        self.rb_ExchangeChoice.append(rb)
        exchangegroup.addButton(rb)
        ypos += deltay
       
        rb=QRadioButton("Ashley approach with B.E. (eV):", parent= self.ProjectileBox )
        rb.setGeometry(QRect(xpos+30,ypos,230,20))
        rb.clicked.connect(self.rb_ExchangeChoice_clicked)
        self.rb_ExchangeChoice.append(rb)
        exchangegroup.addButton(rb)
        self.rb_ExchangeChoice[0].setChecked(True)
        
        self.calc.ExchangeCorrection = False
        self.Core_Level_BE_inp = QLineEdit(str(self.calc.BE_for_exchange),parent= self.ProjectileBox ) 
        self.Core_Level_BE_inp.setValidator(QDoubleValidator(0,1e9, 6))
        self.Core_Level_BE_inp.setGeometry(QRect(270,ypos, 50,18))  
        self.Core_Level_BE_inp.editingFinished.connect(self.Update_BE_for_exchange)
        ypos +=deltay
        self.drawHLine(ypos,self.ProjectileBox)        
        ypos += 5            
        QLabel("DIIMFP integration:",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,125,18))
        self.rb_DIIMFP_Integration = []
        self.integrationQualityGroup = QButtonGroup(self.ProjectileBox)
        LabelList=["rough", "medium","fine"] 
        for i in range(3):
            rb=QRadioButton(LabelList[i], parent= self.ProjectileBox )
            rb.setGeometry(QRect(130+i*80,ypos,150,20))
            rb.clicked.connect(self.rb_DIIMFP_Integration_clicked)
            self.rb_DIIMFP_Integration.append(rb)
            self.integrationQualityGroup.addButton(rb)
        self.rb_DIIMFP_Integration[self.calc.Stopping_calc_quality].setChecked(True)
       
        
        ypos +=deltay
        QLabel("Plot:",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,30,18))
        self.rb_IMFPPlotChoice = []
        self.imfpgroup =  QButtonGroup(self.ProjectileBox)
        LabelList=["imfp", "cross section per unit cell"] 
        for i in range(2):
            rb=QRadioButton(LabelList[i], parent= self.ProjectileBox )
            rb.setGeometry(QRect(60+i*80,ypos,200,20))
            rb.clicked.connect(self.getIMFPPlotChoice)
            self.rb_IMFPPlotChoice.append(rb)
            self.imfpgroup.addButton(rb)
        self.rb_IMFPPlotChoice[0].setChecked(self.runplot.plot_imfp)
    #    self.rb_IMFPPlotChoice[1].setChecked( not self.runplot.plot_imfp)
      
        ypos +=deltay
        QLabel("IMFP, Stop., Strag. versus:",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,200,18))
        self.x_axis_unitgroup = QButtonGroup(self.ProjectileBox) 
        self.rb_x_axis = []
        LabelList=["velocity","Energy"]
        for i in range(2):
            rb=QRadioButton(LabelList[i], parent= self.ProjectileBox )
            rb.setGeometry(QRect(180+i*80,ypos,200,20))
            rb.clicked.connect(self.getEnergyVelocityChoice)
            self.rb_x_axis.append(rb)
            self.x_axis_unitgroup.addButton(rb)
        self.rb_x_axis[0].setChecked(not self.runplot.x_axis_keV)
        self.rb_x_axis[1].setChecked(self.runplot.x_axis_keV) 
        
        ypos +=deltay
        QLabel("Stopping Units:",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,200,18))

        self.StoppingUnitsGroup = QButtonGroup(self.ProjectileBox) 
        self.rb_su = []
        LabelList=["eV / Å","eV / (10¹⁵ atoms / cm²)","MeV / (mg / cm²)"]
        for i in range(3):
            rb=QRadioButton(LabelList[i], parent= self.ProjectileBox )
            rb.setGeometry(QRect(i*90,ypos,200,20))
            rb.clicked.connect(self.setStoppingUnitChoice)
            self.rb_su.append(rb)
            self.StoppingUnitsGroup.addButton(rb)
            if i == self.runplot.StoppingUnits: self.rb_su[i].setChecked(True)
        self.rb_su[0].setGeometry(QRect(190,ypos,300,20))  
        ypos +=deltay
        self.rb_su[1].setGeometry(QRect(xpos,ypos,350,20)) 
        self.rb_su[2].setGeometry(QRect(190,ypos,350,20)) 

        ypos += deltay
        self.drawHLine(ypos,self.ProjectileBox)
         
        ypos += 5   
        self.Update_target_properties()   
        QLabel("Stopping Conversions:",parent=self.ProjectileBox).setGeometry(QRect(xpos,ypos,200,18))
        ypos += deltay
        self.conversiontext1= QLabel(self.calc.myconversiontext1,parent=self.ProjectileBox)
        self.conversiontext1.setGeometry(QRect(xpos,ypos,350,18))
        ypos += deltay
        self.conversiontext2= QLabel(self.calc.myconversiontext2,parent=self.ProjectileBox)
        self.conversiontext2.setGeometry(QRect(xpos,ypos,350,18))
        
        ypos +=deltay
        self.drawHLine(ypos,self.ProjectileBox)
            
        ypos += +5 
        self.Approximationsbox = QCheckBox("Plot approximations.        ωₚ  for TPP (eV):", self.ProjectileBox)
        self.Approximationsbox.setGeometry(QRect(xpos,ypos,290,18))
        self.Approximationsbox.stateChanged.connect(self.Approximationsbox_state_changed) 
        
        self.wp_TPPinp = QLineEdit(str(self.calc.w_p_TPP),parent= self.ProjectileBox ) 
        self.wp_TPPinp.setValidator(QDoubleValidator(0,1e9, 6))
        self.wp_TPPinp.setGeometry(QRect(290,ypos, 50,18))  
        self.wp_TPPinp.editingFinished.connect(self.Update_wp_TPP)
        
        
    def rb_ProjectileChoice_clicked(self):
        if self.rb_ProjectileChoice[0].isChecked(): 
            self.calc.particle="electron"
        else:
            self.calc.particle="proton"
        self.calc.my_updateProjectileEnergy()  
        self.E0_velocityLabel.setText(self.calc.myvelocitytext)
              
    def MottCheckbox_state_changed(self): self.calc.MottCorrection  = self.MottCheckbox.isChecked()
    def Approximationsbox_state_changed(self): 
        self.calc.Approximations  = self.Approximationsbox.isChecked()
   
    def RadLosCheckbox_state_changed(self): self.calc.RadiativeLosses  = self.RadLosCheckbox.isChecked()
    def rb_ExchangeChoice_clicked(self):
        self.calc.ExchangeCorrection = not self.rb_ExchangeChoice[0].isChecked()
        self.calc.Exchange_as_in_SBethe = self.rb_ExchangeChoice[1].isChecked()
    def Update_BE_for_exchange(self): self.calc.BE_for_exchange=float(self.Core_Level_BE_inp.text())
  
    def Update_wp_TPP(self): self.calc.w_p_TPP=float(self.wp_TPPinp.text())  
    def rb_DIIMFP_Integration_clicked(self):
        for i in range (3):
            if self.rb_DIIMFP_Integration[i].isChecked():self.calc.Stopping_calc_quality = i
            
    def getIMFPPlotChoice(self): 
        self.runplot.plot_imfp= self.rb_IMFPPlotChoice[0].isChecked()
        self.Update_calcLabel()
        
    def Update_calcLabel(self):    
        if self.runplot.plot_imfp:
            if self.runplot.x_axis_keV:
                self.rb_calculation_Choice[self.index_imfp_calc].setText("λ(E₀), dE/dx(E₀),  dE²/dx(E₀)")
            else: 
                self.rb_calculation_Choice[self.index_imfp_calc].setText("λ(v), dE/dx(v),  dE²/dx(v)")   
        else: 
            if self.runplot.x_axis_keV:   
                self.rb_calculation_Choice[self.index_imfp_calc].setText("σ(E₀), dE/dx(E₀),  dE²/dx(E₀)")
            else:
                self.rb_calculation_Choice[self.index_imfp_calc].setText("σ(v), dE/dx(v),  dE²/dx(v)")
                    
            
    def getEnergyVelocityChoice(self):
        self.runplot.x_axis_keV=self.rb_x_axis[1].isChecked()
        self.Update_calcLabel()
        
    def setStoppingUnitChoice(self):
        for i in range(3):
            if  self.rb_su[i].isChecked():  self.runplot.StoppingUnits = i 
            
    def  Update_E0(self):  
        self.calc.E0=float(self.E0_inp.text())   
        self.calc.my_updateProjectileEnergy()  
        self.E0_velocityLabel.setText(self.calc.myvelocitytext)      

#start tab3========================================== calculations tab======================            
    def calculation_options(self):
        xpos=5
        ypos=30
        deltay=25
        self.CalculationOptionBox = QGroupBox(parent=self.tab_3,title="Calculation options")
        self.CalculationOptionBox.setGeometry(QRect(4,4,430,695)) 
              
        QLabel("maximum number of ω values considered ",parent=self.CalculationOptionBox).setGeometry(QRect(xpos,ypos,300,18))
        self.NMaxEnergyStepsInput = QLineEdit(str(self.calc.MaxNPoints),parent= self.CalculationOptionBox )  
        self.NMaxEnergyStepsInput.setGeometry(QRect(xpos+320,ypos, 45,18))
        self.NMaxEnergyStepsInput.setValidator(QIntValidator(2,999999))
        self.NMaxEnergyStepsInput.editingFinished.connect(self.MaxNPoints_changed) 
        
        ypos+=deltay
        QLabel("Number of projectile energies calulated for stopping curves",parent=self.CalculationOptionBox).setGeometry(QRect(xpos,ypos,300,18))
        self.NStoppingInput = QLineEdit(str(self.calc.NStopping),parent= self.CalculationOptionBox )  
        self.NStoppingInput.setGeometry(QRect(xpos+320,ypos, 45,18))
        self.NStoppingInput.setValidator(QIntValidator(1,1000))
        self.NStoppingInput.editingFinished.connect(self.NStopping_changed) 
        
        ypos+=deltay
        QLabel("Stopping: Energy increment factor (>1)",parent=self.CalculationOptionBox).setGeometry(QRect(xpos,ypos,300,18))
        self.EIncrementInput = QLineEdit(str(self.calc.IncrFactor),parent= self.CalculationOptionBox )  
        self.EIncrementInput.setGeometry(QRect(xpos+320,ypos, 45,18))
        self.EIncrementInput.setValidator(QDoubleValidator(1.0,3.0,6))
        self.EIncrementInput.editingFinished.connect(self.EIncrement_changed) 
        
        ypos+=deltay
        QLabel("First energy stopping curve e⁻ (keV):",parent=self.CalculationOptionBox).setGeometry(QRect(xpos,ypos,300,18))
        self.FirstEnergyElectronInput = QLineEdit(str(self.calc.first_electron_energy),parent= self.CalculationOptionBox )  
        self.FirstEnergyElectronInput.setGeometry(QRect(xpos+320,ypos, 75,18))
        self.FirstEnergyElectronInput.setValidator(QDoubleValidator(0.0,1e9,6))
        self.FirstEnergyElectronInput.editingFinished.connect(self.FirstEnergyElectron_changed)    
        
        ypos+=deltay
        QLabel("First energy stopping curve H⁺ (keV):",parent=self.CalculationOptionBox).setGeometry(QRect(xpos,ypos,300,18))
        self.FirstEnergyProtonInput = QLineEdit(str(self.calc.first_proton_energy),parent= self.CalculationOptionBox )  
        self.FirstEnergyProtonInput.setGeometry(QRect(xpos+320,ypos, 75,18))
        self.FirstEnergyProtonInput.setValidator(QDoubleValidator(0.0,1e9,6))
        self.FirstEnergyProtonInput.editingFinished.connect(self.FirstEnergyProton_changed)   
        
        ypos+=4*deltay
        QLabel("Replace RPA calculation by Drude-Lindhard for q < Cₜ*ω (a.u.)",parent=self.CalculationOptionBox).setGeometry(QRect(xpos,ypos,415,18)) 
        ypos+=deltay
        QLabel("Cₜ=",parent=self.CalculationOptionBox).setGeometry(QRect(xpos+320,ypos,345,18)) 
        self.c_transitionInput = QLineEdit(str(self.calc.c_transition),parent= self.CalculationOptionBox )  
        self.c_transitionInput.setGeometry(QRect(xpos+345,ypos, 75,18))
        self.c_transitionInput.setValidator(QDoubleValidator(0.0,0.1,6))
        self.c_transitionInput.editingFinished.connect(self.c_transition_changed)   
    def MaxNPoints_changed(self):self.calc.MaxNPoints=int(self.NMaxEnergyStepsInput.text())
    def NStopping_changed(self):self.calc.NStopping=int(self.NStoppingInput.text())   
    def EIncrement_changed(self):self.calc.IncrFactor=float(self.EIncrementInput.text())   
    def FirstEnergyElectron_changed(self):self.calc.first_electron_energy=float(self.FirstEnergyElectronInput.text()) 
    def FirstEnergyProton_changed(self):  self.calc.first_proton_energy=float(self.FirstEnergyProtonInput.text())     
    def c_transition_changed(self):   self.calc.c_transition=float(self.c_transitionInput.text())     
       
    def plot_options(self):
        xpos=5
        ypos=30
        deltay=25
        self.PlotOptionBox = QGroupBox(parent=self.tab_3,title="Plot options")
        self.PlotOptionBox.setGeometry(QRect(445,4,430,696)) 
        # self.LatexModeCheckbox = QCheckBox('use Latex mode', self.PlotOptionBox )  # discarded, too many problems for little gain
        # self.LatexModeCheckbox.setChecked(self.runplot.LaTeXlike)
        # self.LatexModeCheckbox.setGeometry(QRect(xpos,ypos, 150,18))
        # self.LatexModeCheckbox.stateChanged.connect(self.LatexModeCheckbox_state_changed) 
        #   def LatexModeCheckbox_state_changed(self): self.runplot.LaTeXlike = self.LatexModeCheckbox.isChecked() 
        
        ypos+=deltay
        QLabel("default file format figures:",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,200,18)) 
        self.rb_figOutputFormat = []
        rb =QRadioButton("pdf", parent= self.PlotOptionBox )
        rb.setGeometry(QRect(xpos+180, ypos,100,20))
        rb.clicked.connect(self.rb_figOutputFormat_clicked)
        self.rb_figOutputFormat.append(rb)
        rb =QRadioButton("png", parent= self.PlotOptionBox )
        rb.setGeometry(QRect(xpos+230, ypos,50,20))
        rb.clicked.connect(self.rb_figOutputFormat_clicked)
        self.rb_figOutputFormat.append(rb)
        if self.runplot.fileformat== "pdf":
            self.rb_figOutputFormat[0].setChecked(True)
        else:    
            self.rb_figOutputFormat[1].setChecked(True)
            
        ypos +=deltay
        QLabel("figure width (cm)",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,100,18)) 
        self.FigureWidth_Inp = QLineEdit(str(self.runplot.figurewidth),parent= self.PlotOptionBox ) 
        self.FigureWidth_Inp.setValidator(QDoubleValidator(0.0,100.0,6))
        self.FigureWidth_Inp.setGeometry(QRect(xpos+120,ypos, 40,18))  
        self.FigureWidth_Inp.editingFinished.connect(self.Update_FigureWidth_Inp)  
        
        ypos +=deltay
        QLabel("figure height (cm)",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,100,18)) 
        self.FigureHeight_Inp = QLineEdit(str(self.runplot.figureheight),parent= self.PlotOptionBox ) 
        self.FigureHeight_Inp.setValidator(QDoubleValidator(0.0,100.0,6))
        self.FigureHeight_Inp.setGeometry(QRect(xpos+120,ypos, 40,18))  
        self.FigureHeight_Inp.editingFinished.connect(self.Update_FigureHeight_Inp)  
       
        ypos +=deltay
        QLabel("Title font size",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,120,18)) 
        self.Title_Fontsize_Inp = QLineEdit(str(self.runplot.title_fontsize),parent= self.PlotOptionBox ) 
        self.Title_Fontsize_Inp.setValidator(QIntValidator(0,50))
        self.Title_Fontsize_Inp.setGeometry(QRect(xpos+120,ypos, 40,18))  
        self.Title_Fontsize_Inp.editingFinished.connect(self.Update_Title_Fontsize_Inp)  
        
        ypos +=deltay
        QLabel("Legend font size",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,120,18)) 
        self.Legend_Fontsize_Inp = QLineEdit(str(self.runplot.legendfontsize),parent= self.PlotOptionBox ) 
        self.Legend_Fontsize_Inp.setValidator(QIntValidator(0,50))
        self.Legend_Fontsize_Inp.setGeometry(QRect(xpos+120,ypos, 40,18))  
        self.Legend_Fontsize_Inp.editingFinished.connect(self.Update_Legend_Fontsize_Inp)  
        
        ypos +=deltay
        QLabel("x tick Label font size",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,120,18)) 
        self.xLabel_Fontsize_Inp = QLineEdit(str(self.runplot.xlabelsize ),parent= self.PlotOptionBox ) 
        self.xLabel_Fontsize_Inp.setValidator(QIntValidator(0,50))
        self.xLabel_Fontsize_Inp.setGeometry(QRect(xpos+120,ypos, 40,18))  
        self.xLabel_Fontsize_Inp.editingFinished.connect(self.Update_xLabel_Fontsize_Inp)  
        
        ypos +=deltay
        QLabel("y tick Label font size",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,120,18)) 
        self.yLabel_Fontsize_Inp = QLineEdit(str(self.runplot.ylabelsize ),parent= self.PlotOptionBox ) 
        self.yLabel_Fontsize_Inp.setValidator(QIntValidator(0,50))
        self.yLabel_Fontsize_Inp.setGeometry(QRect(xpos+120,ypos, 40,18))  
        self.yLabel_Fontsize_Inp.editingFinished.connect(self.Update_yLabel_Fontsize_Inp)  
       
        ypos +=deltay
        QLabel("axis Label font size",parent=self.PlotOptionBox).setGeometry(QRect(xpos,ypos,120,18)) 
        self.axisLabel_Fontsize_Inp = QLineEdit(str(self.runplot.axislabelsize ),parent= self.PlotOptionBox ) 
        self.axisLabel_Fontsize_Inp.setValidator(QIntValidator(0,50))
        self.axisLabel_Fontsize_Inp.setGeometry(QRect(xpos+120,ypos, 40,18))  
        self.axisLabel_Fontsize_Inp.editingFinished.connect(self.Update_axisLabel_Fontsize_Inp)          

        ypos +=deltay
        text="""from matplotlib documentation, this seems buggy but useful
        
**********************************************************************
* INTERACTIVE KEYMAPS                                                    *
**********************************************************************
Event keys to interact with figures/plots via keyboard.
See https://matplotlib.org/stable/users/explain/interactive.html 
fullscreen:\t f, ctrl+f    
home: \t\t h, r, home         
back:\t\t left arrow, c
forward:\t\t right arrow, v
pan:\t\t p                   
zoom:\t\t o                  
save:\t\t s 
closse figure:\t ctrl+w
grid:\t\t g                  
grid_minor:\t G            
log/linear yscale:\t l               
log/linear xscale\t k, L             
ctrl+c, cmd+c:\t copy figure to clipboard     """
        QLabel(text,parent=self.PlotOptionBox).setGeometry(QRect(xpos,390,430,310)) 
     
        
        

    def Update_FigureWidth_Inp(self):self.runplot.figurewidth= float(self.FigureWidth_Inp.text())   
    def Update_FigureHeight_Inp(self):self.runplot.figureheight= float(self.FigureHeight_Inp.text())  
    def Update_Title_Fontsize_Inp(self):self.runplot.title_fontsize= int(self.Title_Fontsize_Inp.text()) 
    def Update_Legend_Fontsize_Inp(self):self.runplot.legendfontsize = int(self.Legend_Fontsize_Inp.text())      
    def Update_xLabel_Fontsize_Inp(self):self.runplot.xlabelsize = int(self.xLabel_Fontsize_Inp.text())  
    def Update_yLabel_Fontsize_Inp(self):self.runplot.ylabelsize = int(self.yLabel_Fontsize_Inp.text())  
    def Update_axisLabel_Fontsize_Inp(self):self.runplot.axislabelsize = int(self.axisLabel_Fontsize_Inp.text()) 
    def rb_figOutputFormat_clicked(self):
        if      self.rb_figOutputFormat[0].isChecked():
            self.runplot.fileformat= "pdf"
        else:
            self.runplot.fileformat= "png"            
            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def MainWindowIO(self):
        self.CalcPlotBtn = QPushButton("calculate and plot",parent=self.centralWidget) 
        self.CalcPlotBtn.setStyleSheet("border: 2px solid red;")
        self.CalcPlotBtn.setObjectName("cpbtn") 
        self.CalcPlotBtn.setGeometry(self.MainWindowWidth-330,self.MainWindowHeight-60,130,30)  
        
        self.CalcPlotBtn.clicked.connect(self.runplot.start_calc_and_plot)  

        self.ReplotBtn = QPushButton("replot",parent=self.centralWidget)     
        self.ReplotBtn.setGeometry(self.MainWindowWidth-180,self.MainWindowHeight-60,50,30)  
        self.ReplotBtn.clicked.connect(self.runplot.replot)  

        self.ClosePlotsBtn = QPushButton("close all plots",parent=self.centralWidget)     
        self.ClosePlotsBtn.setGeometry(self.MainWindowWidth-110,self.MainWindowHeight-60,100,30)  
        self.ClosePlotsBtn.clicked.connect(self.runplot.Close_all_plots)  

        QLabel("Calculation Description:",parent=self.centralWidget).setGeometry(QRect(5,self.MainWindowHeight-68,140,18))
        self.description_inp = QLineEdit(self.runplot.CalcDescription,parent= self.centralWidget) 
        self.description_inp.setGeometry(QRect(150,self.MainWindowHeight-68, 350,18))  
        self.description_inp.editingFinished.connect(self.UpdateDescription_inp)
        
        QLabel("Calculation Status:",parent=self.centralWidget).setGeometry(QRect(5,self.MainWindowHeight-46,120,18))

        self.DisplayStatus= QLabel("",parent=self.centralWidget)
        self.DisplayStatus.setStyleSheet("QLabel { background-color : black; color : white; }")
        self.DisplayStatus.setGeometry(QRect(150,self.MainWindowHeight-46,350,18))
        
    def UpdateCalcPlotBtn(self, newtext):
        self.CalcPlotBtn.setText(newtext) 
        self.CalcPlotBtn.repaint()
        
        
    def UpdateStatus(self,text):
       
        self.DisplayStatus.setText(text)
        text_lower_case = text.lower()
        x=text_lower_case.find("error")
        if x >= 0:
            self.DisplayStatus.setStyleSheet("QLabel { background-color : red; color : blue; }")
        else:
            self.DisplayStatus.setStyleSheet("QLabel { background-color : black; color : white; }")
            
        self.DisplayStatus.repaint()
      
    
    def UpdateDescription_inp(self):self.runplot.CalcDescription=self.description_inp.text()    
  
 
                
    def Help(self):
        a=os.path.join(os.getcwd(), "chapidif_manual_plus_background.pdf")
        browser.open_new("file:" + a)     
         
    def check_for_update(self):
        browser.open_new("https://github.com/MaartenVos/Chapidif")    
      
              
    def MyFormat1(self, myfloat):  # custom formatter to use table space best
        try:
            a=abs(myfloat)
        except Exception:
            return    
        if a == 0: mystring = "0.0"
        elif a < 0.01: mystring =  f"{myfloat:.3e}".replace("e-0", "e-")
        elif a < 0.1:  mystring =  f"{myfloat:.5f}".replace("0.", ".")
        elif a < 1.0:  mystring =  f"{myfloat:.4f}"
        elif a < 10.0: mystring =  f"{myfloat:.3f}"
        elif a < 100.0: mystring =  f"{myfloat:.2f}"
        elif a < 1000.0: mystring =  f"{myfloat:.1f}"
        else: mystring =  f"{myfloat:.0f}"
        return mystring
            
    def MyFormat2(self, myfloat):  # custom formatter to use table space best
        a=abs(myfloat) 
        if a == 0: mystring = "0.0"
        elif a < 10.0: mystring =  f"{myfloat:.3f}"
        elif a < 100.0: mystring =  f"{myfloat:.2f}"
        elif a < 1000.0: mystring =  f"{myfloat:.1f}"
        else: mystring =  f"{myfloat:.0f}"
        return mystring            
        
        
    def closeEvent(self, event: QCloseEvent):
        self.runplot.Close_all_plots()
        event.accept()  # Allow the window to close
    
    def drawHLine(self, start_point_y,Box):
        width=Box.width()
        frame = QFrame(Box)
        frame.setFrameShape(QFrame.Shape.HLine)
        frame.setStyleSheet("QFrame { border: 1px solid blue; }")
        offset=5
        frame.setGeometry(offset,start_point_y, width-(2*offset),1)
    
                 
if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    app.setStyle('Fusion') # was 'Windows'
    win = myWindow()
    win.show()
    # Run the event loop
    sys.exit(app.exec())
  
