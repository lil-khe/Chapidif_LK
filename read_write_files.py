import re  # for regular expression to extract numbers from string

import os

#from PyQt6.QtWidgets import QFileDialog   ### MODIF 11/09 LK: I think it is not needed for our use, and it triggers an error in Windows compilation
   

import numpy as np
#from PyQt6 import QtWidgets,QtCore
from pathlib import Path
try:
    from icecream import ic

    ic.configureOutput(includeContext=True)
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

   # ================================readfile

class read_write_files:
   # def __init__(self, *args, **kwargs):
    def __init__(self,parent):
        self.MyChapApp=parent
        self.calc = self.MyChapApp.calc
        self.runplot  = self.MyChapApp.runplot 
        self.current_dir = str(Path.home())
        self.DF_comment="" # comment for DF file
        self.fname=""      #header in the interface of the comment
        self.runplot.xCompArray=[]
        self.runplot.yCompArray=[]   
        self.Vuesz_mode=True
        if self.Vuesz_mode:  # output file optimesd for importing in Vuesz
            self.Vuesz_descriptor="descriptor `"
        else:
            self.Vuesz_descriptor="`"     
        
    def read_df(self): 
        self.calc.ZeroDF()
        self.MyChapApp.update_all_tables()
        file_dialog = QFileDialog(self.MyChapApp)
        file_dialog.setWindowTitle("Open DF File")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setViewMode(QFileDialog.ViewMode.Detail)
        file_dialog.setNameFilter("Chapidif DF Files (*.cpd *.dat);;Dat Files (*.dat);;All Files (*)")
        if (file_dialog.exec()):
            fileNames = file_dialog.selectedFiles()
        try:
            fileameonly = os.path.basename(fileNames[0])
        except Exception:
            self.MyChapApp.UpdateStatus("no file")
            return 
        self.fname =fileNames[0]
        self.current_dir=os.path.dirname(str(fileNames[0]))
        self.MyChapApp.fname="Comment for/from DF file: " + fileameonly  #header of comment used in interface
        try:
            myfile = open(fileNames[0], "r")
        except Exception:
            text="could not open file:"+fileNames[0]
            self.MyChapApp.UpdateStatus(text)
            return    
        myfile.readline()  # skip line
        self.calc.AddELF=1
        line2 = myfile.readline()
        if "Extended" in line2:
            self.calc.DFmodel='Drude'
            self.calc.AddELF = 0
        elif "Drude Lindhard" in line2:
            self.calc.DFmodel='DL'
        elif " Mermin " in line2 or  "plain Lindhard" in line2:  # keep space around Mermin to avoid being triggered by (Mermin-renormalisation) in the Kaneko case
            self.calc.DFmodel='Mermin'
        elif "analytical" in line2:  # put this before traditional Tauc-lorentz one
            self.calc.DFmodel='TL_an'
        elif "Tauc-Lorentz" in line2:
            self.calc.DFmodel='Tauc'
        elif "Vlasov" in line2:
            self.calc.DFmodel='Vlasov'   
        elif "Forouhi-Bloomer" in line2:
            self.calc.DFmodel='FB' 
        if "(non-relativistic)" in line2: 
            self.calc.Dispersion_relativistic = 0
        if "(relativistic)" in line2: 
            self.calc.Dispersion_relativistic = 1
        if "quadratic dispersion" in line2:    
            self.calc.Dispersion_choice = 0
        if "full dispersion" in line2:    
            self.calc.Dispersion_choice = 1
        if "delayed dispersion" in line2:
            self.calc.delayed_dispersion = 1
        else:
            self.calc.delayed_dispersion = 0        
        if "add Doppler width" in line2:
            self.calc.Add_Doppler_Width = 1
        else:
            self.calc.Add_Doppler_Width = 0     
        if "add Chi" in line2:
            self.calc.AddELF = 0
        elif "add ELF" in line2:
            self.calc.AddELF = 1
        if "(Mermin-renormalisation)"in line2:
            self.calc.Merminize = 1
        elif "Direct" in line2:
            self.calc.Merminize = 2    
        elif self.calc.DFmodel =='Mermin':
            self.calc.Merminize = 0
      
        
        commentline = myfile.readline()  # this is line 3
        self.DF_comment=commentline.strip()
        if "#" in self.DF_comment : self.DF_comment  = self.DF_comment.replace("#", "")
            
        if self.calc.DFmodel == "Drude" or self.calc.DFmodel  == "Tauc"  or self.calc.DFmodel  == "TL_an" or self.calc.DFmodel  == "FB":  
            line4 = myfile.readline()
            result = re.findall(r"\d*\.?\d+", line4)  # from import re
            if  self.calc.DFmodel  == "TL_an":
                self.calc.a_TL_an = float(result[0])
            elif self.calc.DFmodel  == "FB": 
                self.calc.n_infty = float(result[0])
            else:    
                self.calc.eps_bkg = float(result[0])

        line4 = myfile.readline()
        if line4.find("density") >= 0:
            elements = line4.split(",")
            self.calc.specificweight = float(elements[1])
            self.calc.massunitcell = float(elements[3])
        myfile.readline()  # skip line

        i = 0
        igos = 0
        iKaneko = 0
        iBelkacem = 0
        while 1:
            line = myfile.readline()
            if "#" in line: line = line.replace("#", "")
            
            nums = line.split()
            if len(nums) < 3:
                break
            if not line:
                break
            if "Amplitude" in line:  
                continue
            if "N elec" in line:
                continue
                    
            elif "GOS" in line:
                self.calc.ConcGOS[igos]= float(nums[1])
                self.calc.EdgeGOS[igos] = float(nums[2])
                self.calc.nlGOS[igos] = int(nums[3])
                self.calc.ZGOS[igos] = int(nums[4])
                try:
                    self.calc.maxEnergyDensityEffect = float( nums[10])
                except Exception:
                    self.MyChapApp.UpdateStatus("no maximum energy density effect given")    
                igos = igos + 1

                if "rescaling on" in line:
                    self.calc.ApplySumRuleToGOS = 1
                elif  "rescaling off" in line:  
                    self.calc.ApplySumRuleToGOS = 0
                    
            elif "Kaneko" in line:
                self.calc.N_Kaneko[iKaneko] = float(nums[1])
                self.calc.Q_Kaneko[iKaneko] = float(nums[2])
                self.calc.width_Kaneko[iKaneko] = float(nums[3])
                self.calc.Edge_Kaneko[iKaneko] = float(nums[4])
                
                self.calc.l_Kaneko[iKaneko] =  int(nums[5])
                self.calc.gamma_Kaneko[iKaneko] = float(nums[6])
                iKaneko = iKaneko + 1
                if "modified Kaneko" in line:
                    self.calc.Kaneko_choice = 0
                elif "original Kaneko" in line:  
                    self.calc.Kaneko_choice = 1 
                        
            elif "Belkacem" in line:
                self.calc.Conc_Belkacem[iBelkacem] = float(nums[1])
                self.calc.w_Belkacem[iBelkacem] = float(nums[2])
                self.calc.gamma_Belkacem[iBelkacem] = float(nums[3])
                iBelkacem = iBelkacem + 1
                
            else:
                self.calc.Amps[i]= float(nums[0])
                self.calc.Gammas[i] = float(nums[1])
                self.calc.Omegas[i] = float(nums[2])
                self.calc.Us[i] = 0.0
                self.calc.Alphas[i] = 1.0
                if len(nums) == 4:
                    if not (self.calc.DFmodel == "Drude" or self.calc.DFmodel == "DL"):  #old files have no U for these models 1 and 2 but new ones have!
                        self.calc.Us[i] = float(nums[3])
                    else:
                        self.calc.Alphas[i]=(float(nums[3]))
                        self.calc.Us[i] = 0.0
                elif len(nums) == 5:
                    self.calc.Alphas[i] = float(nums[3])  
                    self.calc.Us[i] = float(nums[4])
            i = i + 1
        myfile.close()
        self.MyChapApp.update_all_tables()
       
 # ================================read oos file
        
    def DF_from_OOS(self):
        #Assumes files in the format the .mat files produced by Salvat's "cbethe"
        #self.calc.init()
        self.calc.DF_prop_text=""
        self.MyChapApp.fname="Comment for/from DF file:"
        file_dialog = QFileDialog(self.MyChapApp)
        file_dialog.setWindowTitle("Open ,MAT File")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setViewMode(QFileDialog.ViewMode.Detail)
        file_dialog.setNameFilter("df files (*.mat *.*)")
        if (file_dialog.exec()):
            fileNames = file_dialog.selectedFiles()
        try:
            fileameonly = os.path.basename(fileNames[0])
        except Exception:
            return 2

        self.fname ="filename: " + fileNames[0]

        self.current_dir=os.path.dirname(str(fileNames[0]))
        
        self.MyChapApp.fname="Comment for/from DF file: " + fileameonly
        self.MyChapApp.runplot.initcalc()
        

        with open(fileNames[0], 'r') as fp:
            for FileLength, line in enumerate(fp):
                pass
        try:
            myfile = open(fileNames[0], "r")
        except OSError:
            self.MyChapApp.UpdateStatus("could not open file "+ fileNames[0] )
            return 
        Header = True
        HeaderLength = 0
    
        while  Header:  
            line = myfile.readline()
            HeaderLength += 1
            if line.find("Molecular weight") >= 0:
               elements = line.split(" ")
               self.calc.massunitcell=float(elements[1])
            if line.find("Mass density") >= 0:
               elements = line.split(" ")
               self.calc.specificweight= float(elements[1])
            if line.find("OOS (1/eV") >= 0:
                 Header = False 
        self.calc.N_OOS=FileLength-HeaderLength-1  # there is an empty line at end .mat file. so my file length is off by 1
        for i in range(self.calc.N_OOS):
            line = myfile.readline()
            elements = line.split(" ")
            self.calc.OOSEnergy[i]= float(elements[1])
            self.calc.OOS[i]=float(elements[3])
        myfile.close()
        self.calc.DFChoice=2  
        self.calc.DFmodel="DL"     
        self.calc.DL_from_OOS()
        self.MyChapApp.update_all_tables()
        return
        
 #========create output files=======================================================================================       
    def write_df_only(self):
        self.write_df()
        
        
        
    
        #self.MyChapApp.fname= "Comment for/from DF file:"

    def save_calculation(self):  
        if self.write_df():  #write dielectric function (can be parsed back into chapidif) 
            self.write_calculation_result() #write calculation details and results (cancurrently not be parsed back into chapidif) 
            if self.Vuesz_mode:
                self.write_VueszFile()
          
 
                

    # ==========write current configuration and optional calculation result

    def write_df(self):
     
        file_dialog = QFileDialog(self.MyChapApp)
        file_dialog.setWindowTitle("Save DF parameters")
        file_dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
        file_dialog.setNameFilter("Chapidif DF Files (*.cpd);;Dat Files (*.dat);;All Files (*)")
        file_dialog.setDefaultSuffix("cpd")  # Set the default extension to .cpd
        file_dialog.exec()
        try:
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                filenameList = file_dialog.selectedFiles()
                self.OutputFile=filenameList[0]
           
            self.myVeuszFile= os.path.splitext(self.OutputFile)[0]+'.vsz' 
            self.current_dir=os.path.dirname(str(self.OutputFile))
        except Exception:
            return  False 
        self.now_write()
        return True
        
    def now_write(self):
           # file written so compatible with Sesinipac
        # line one description of file type
        # line 2 model used: Mermin drude etc
        # line 3 comment line
        # line 4 background dielectric constant
        # line 5 amp gamma omega alpha u
        # following lines oscillators
        # lines starting with GOS GOS parameters (not read by sesinipac)
        # file_filter = "Chapidif Files (*.cpd);;Dat Files (*.dat)"
        # self.filename=QFileDialog.getSaveFileName(self.MyChapApp, "Save DF to file", self.current_dir, file_filter)
        with open(self.OutputFile, "w", encoding="utf-8") as file:
            file.write("#input file for Chapidif\n")
            model_used="#"
            if self.calc.DFmodel == "Drude":
                model_used += "Dielectric function model: Extended Drude"
            elif self.calc.DFmodel == "DL":
                model_used += "Dielectric function model: Drude Lindhard"
            elif self.calc.DFmodel == "Mermin" and  self.calc.Merminize == 0:
                model_used += "Dielectric function model: plain Lindhard or LL"
            elif self.calc.DFmodel == "Mermin" and  self.calc.Merminize == 1:
                model_used += "Dielectric function model: Mermin MLL"
            elif self.calc.DFmodel == "Mermin" and  self.calc.Merminize== 2:
                model_used += "Dielectric function model: Lindhard-direct"
            elif self.calc.DFmodel == "Tauc" :
                model_used += "Dielectric function model: Tauc-Lorentz" 
            elif self.calc.DFmodel == "TL_an" :
                model_used += "Dielectric function model: Tauc-Lorentz analytical"   
            elif self.calc.DFmodel == "FB" :
                model_used += "Dielectric function model: Forouhi-Bloomer"   
            if self.calc.AddELF == 0 or self.calc.DFmodel == "Tauc" or self.calc.DFmodel == "TL_an" :
                model_used += ", add Chi"
            else:
                model_used += ", add ELF"
            if self.calc.Merminize == 0 and  (self.calc.DFmodel == "Mermin" or self.calc.DFmodel == "Vlasov") :    
                model_used +=  ", (without Mermin-renormalisation)"
            elif self.calc.Merminize == 1 and (self.calc.DFmodel == "Mermin" or self.calc.DFmodel == "Vlasov") :    
                model_used +=  ", (Mermin-renormalisation)"
            elif  (self.calc.DFmodel == "Mermin" or self.calc.DFmodel == "Vlasov") :     
                model_used +=  ", (Direct)"
            if self.calc.Dispersion_choice== 0:
                model_used +=  ", quadratic dispersion"
            else:
                model_used += ", full dispersion"
            if self.calc.Dispersion_relativistic == 0:
                model_used +=  ", (non-relativistic)"
            else:
                model_used +=  ", (relativistic)"     
            if self.calc.delayed_dispersion == 1 and (self.calc.DFmodel == "Tauc"  or self.calc.DFmodel == "TL_an"): 
                model_used +=  ", (delayed dispersion)"       
            if self.calc.Add_Doppler_Width == 1 and (self.calc.DFmodel == "Tauc"  or self.calc.DFmodel == "TL_an" or self.calc.DFmodel == "DL"):  
                model_used +=  ", (add Doppler width)"    
            file.write(model_used+"\n")    
            file.write("#"+self.DF_comment + "\n")
            if self.calc.DFmodel  == "Drude" or self.calc.DFmodel == "Tauc" :
                file.write("#Background dielectric constant:" + str(self.calc.eps_bkg) + "\n")  #questionable variable should basically be 1
            elif self.calc.DFmodel  == "TL_an": 
                file.write("#a TL- analytic model: " + str(self.calc.a_TL_an) + "\n")  
            elif self.calc.DFmodel  == "FB": 
                file.write("#n at infinity: " + str(self.calc.n_infty) + "\n")  #questionable variable should basically be 1
            file.write("#density," + str(self.calc.specificweight)
                + ", unit cell mass," + str(self.calc.massunitcell) + "\n")
            header=self.MyChapApp.myHeaderLabels
            header[1]  ,header[2]= header[2],header[1]  # interchange two elements, as we write gamma before omega  
            text='\t'.join(header)
            text = text.replace("\n", "")
            text = "#"+text.replace("-\t", "")
            print(text)
            file.write(text+"\n")
        
            for i in range(self.calc.maxOscillators):
                if self.calc.Amps[i] != 0.0:
                    if self.calc.DFmodel == "Mermin" or  self.calc.DFmodel == "Vlasov":
                        file.write("#"+ str(self.calc.Amps[i]) + " \t"
                            + str(self.calc.Gammas[i]) + " \t"
                            + str(self.calc.Omegas[i]) + " \t"
                            + str(self.calc.Us[i]) + "\n" )
                    else:
                        file.write("#" + str(self.calc.Amps[i]) + " \t"
                            + str(self.calc.Gammas[i]) + " \t"
                            + str(self.calc.Omegas[i]) + " \t"
                            + str(self.calc.Alphas[i])  + "\t"
                            + str(self.calc.Us[i]) + "\n" )        
                i = i + 1
           
              
            file.write("#GOS N elec \tedge(eV) \t 10 * n + l\t  atomic number Z \n")                               
            for i in range(self.calc.maxGOS):
                if  self.calc.ConcGOS[i] > 0.0:
                    GOSText = "#GOS \t" + str(self.calc.ConcGOS[i]) + "\t" + \
                        str(self.calc.EdgeGOS[i]) + " \t" + str(round(self.calc.nlGOS[i])) + " \t" + str(round(self.calc.ZGOS[i])) + "\t"
                    GOSText+= " density effect included up  to: "+ str(self.calc.maxEnergyDensityEffect) + " (eV),\t"
                   
                    if self.calc.ApplySumRuleToGOS== 1:
                        GOSText+= "rescaling on \n"
                    else:
                        GOSText+= "rescaling off\n"
                    file.write( GOSText)
                i = i + 1
            file.write("#Belkacem N elec \tomega(eV) \t Gamma \n")    
            for i in range(self.calc.maxBelkacem):
                if self.calc.Conc_Belkacem[i] > 0.0:
                    file.write( "#Belkacem \t"
                        + str(self.calc.Conc_Belkacem[i]) + " \t"
                        + str(self.calc.w_Belkacem[i]) + " \t"
                        + str(self.calc.gamma_Belkacem[i]) + "\n")
                i = i + 1
            file.write("#Kaneko N elec \t Q (a.u.) \t width (eV)\t gap U (eV)\t ang. mom l \t gamma Arista \n")    
            for i in range(self.calc.maxKaneko):
                if self.calc.N_Kaneko[i] > 0.0:
                    KanekoText = "#Kaneko \t" + str(self.calc.N_Kaneko[i]) + " \t" \
                        + str(self.calc.Q_Kaneko[i]) + " \t" +str(self.calc.width_Kaneko[i]) + "\t" \
                        + str(self.calc.Edge_Kaneko[i]) + " \t" \
                        + str(round(self.calc.l_Kaneko[i])) + " \t" + str(self.calc.gamma_Kaneko[i]) + "\t"
                    if self.calc.Kaneko_choice == 0:
                        KanekoText += " modified Kaneko\n"
                    else:
                        KanekoText += " original Kaneko\n"        
                    file.write( KanekoText)
                       
                i = i + 1
            file.write("#\n#==========end description dielectric function==========\n")
            
            
    def write_calculation_result(self):  
        cmd=self.MyChapApp.plotchoice
        if cmd == "eq_plot" or cmd == "self_eq_plot" or cmd == "dcs_omega_eq_plot":
            self.write_colorplot_files() 
            return
        if cmd == "REELS":
            self.add_REELS_details()    
        with open(self.OutputFile, "a", encoding="utf-8") as file:
           
            
            file.write("# calculation executed:" + cmd + "\n")
            projectile =self.runplot.particle_LaTeX
           
            if projectile == "H⁺":
                projectileoptions= "H⁺"
                if self.calc.MottCorrection==1:
                    projectileoptions += ": incl. Mott correction"
                else:
                    projectileoptions += ": no Mott correction"    
            elif projectile == "e⁻": 
                projectileoptions= "e⁻"
                if self.calc.ExchangeCorrection:
                    if self.calc.Exchange_as_in_SBethe:
                        projectileoptions+=": exchange as in SBethe"
                    else:
                        projectileoptions+=": exchange as in Ashley, BE= {} eV".format(self.calc.BE_for_exchange) 
                else:
                    projectileoptions+=": no exchange"            
                         
            file.write("# "+projectileoptions)

            file.write("\n#==========start calculation  output section==========\n")
          
        if cmd == "partial_DIIMFP" or cmd == "partial_stopping": 
            self.write_partial_calc_results()
            return
        elif cmd == "IMFP_stop_strag" and self.calc.Approximations:
            self.write_IMFP_stop_strag_plus_approx()
            return
       
           
        elif self.runplot.Text_label[3] != "--":    
            with open(self.OutputFile, "a", encoding="utf-8") as file: 

                file.write(self.Vuesz_descriptor + self.runplot.xlabel + "`\t`" + self.runplot.Text_label[0] +"`\t`" + self.runplot.Text_label[1] + "`\t`" + self.runplot.Text_label[2] + "`\t`" + self.runplot.Text_label[3] + "`\n") 
              
                for i in range(len(self.calc.Result1)):
                    a = "{:.4e} \t{:.4E} \t{:.4E} \t{:.4E}\t{:.4E}\n".format(
                        self.calc.x_axis[i], self.calc.Result1[i],self.calc.Result2[i],self.calc.Result3[i],self.calc.Result4[i])
                    file.write(a)    
        elif self.runplot.Text_label[2] != "--":  
            with open(self.OutputFile, "a", encoding="utf-8") as file: 

                file.write(self.Vuesz_descriptor + self.runplot.xlabel + "`\t`" + self.runplot.Text_label[0] +"`\t`" + self.runplot.Text_label[1] + "`\t`" + self.runplot.Text_label[2] + "`\n") 
                if cmd=="IMFP_stop_strag":
                    for i in range(len(self.calc.CurvesEnergy)):
                        a = "{:.4e} \t{:.4E} \t{:.4E} \t{:.4E}\n".format( self.calc.CurvesEnergy[i],self.calc.IMFPEnergy[i],
                            self.calc.StoppingEnergy[i]*self.runplot.StoppingFactor, self.calc.StragglingEnergy[i])
                        file.write(a)  
                else:    
                    for i in range(len(self.calc.Result1)):
                        a = "{:.4e} \t{:.4E} \t{:.4E} \t{:.4E}\n".format(
                            self.calc.x_axis[i], self.calc.Result1[i],self.calc.Result2[i],self.calc.Result3[i])
                        file.write(a)      
        elif self.runplot.Text_label[1] != "--": 
            with open(self.OutputFile, "a", encoding="utf-8") as file: 

                file.write(self.Vuesz_descriptor + self.runplot.xlabel + "`\t`" + self.runplot.Text_label[0] + "`\t`" + self.runplot.Text_label[1] + "`\n") 
                
                for i in range(len(self.calc.Result1)):
                    a = "{:.3f} \t{:.4E} \t{:.4E}\n".format(
                        self.calc.x_axis[i], self.calc.Result1[i], self.calc.Result2[i])
                    file.write(a)  
            return
        elif self.runplot.Text_label[0] != "--":     
   
            with open(self.OutputFile, "a", encoding="utf-8") as file: 

                file.write(self.Vuesz_descriptor + self.runplot.xlabel + "`\t`" + self.runplot.Text_label[0] + "`\n") 
                
                for i in range(len(self.calc.Result1)):
                    a = "{:.4e} \t{:.4E}\n".format(
                        self.calc.x_axis[i], self.calc.Result1[i]
                    )
                    file.write(a)
        else:
            print("no results")    
           
    def write_VueszFile(self):
        # first write the DF parameters as comment. These may disappear if you resave the file in veusz itself
        with open(self.myVeuszFile, "w",  encoding="utf-8") as file:
            file.write("SetCompatLevel(1)\n")
            file.write("AddImportPath('{}')\n".format(self.current_dir))
            file.write("ImportFile(\'{}\', '', ignoretext=True, linked=True)\n".format(self.OutputFile))
            file.write("Set('colorTheme', 'default-latest')\n")
            file.write("Set('StyleSheet/axis-function/autoRange', 'next-tick')\n")
            file.write("Add('page', name='page1', autoadd=False)\n")
            file.write("To('page1')\n")
            file.write("Add('graph', name='graph1', autoadd=False)\n")
            file.write("To('graph1')\n")
            file.write("Add('axis', name='x', autoadd=False)\n")
            file.write("To('x')\n")
            file.write("Set('label','"+self.runplot.xlabel+"')\n")
            file.write("To('..')\n")
            file.write("Add('axis', name='y', autoadd=False)\n")
            file.write("To('y')\n")
            file.write("Set('direction', 'vertical')\n")
            file.write("To('..')\n")
            file.write("Add('xy', name='xy1', autoadd=False)\n")
            file.write("To('xy1')\n")
            file.write("Set('xData', '"+self.runplot.xlabel+"')\n")
            file.write("Set('yData', '"+self.runplot.Text_label[0]+"')\n")
            file.write("Set('key','"+ self.runplot.Text_label[0]+"')\n")
            file.write("To('..')\n")
            if self.runplot.Text_label[1] != "--":
                file.write("Add('xy', name='xy2', autoadd=False)\n")
                file.write("To('xy2')\n")
                file.write("Set('xData', '"+self.runplot.xlabel+"')\n")
                file.write("Set('yData', '"+self.runplot.Text_label[1]+"')\n")
                file.write("Set('key','"+ self.runplot.Text_label[1]+"')\n")
                file.write("To('..')\n")
            if self.runplot.Text_label[2] != "--":
                file.write("Add('xy', name='xy3', autoadd=False)\n")
                file.write("To('xy3')\n")
                file.write("Set('xData', '"+self.runplot.xlabel+"')\n")
                file.write("Set('yData', '"+self.runplot.Text_label[2]+"')\n")
                file.write("Set('key','"+ self.runplot.Text_label[2]+"')\n")
                file.write("To('..')\n")    
            if self.runplot.Text_label[3] != "--":
                file.write("Add('xy', name='xy4', autoadd=False)\n")
                file.write("To('xy4')\n")
                file.write("Set('xData', '"+self.runplot.xlabel+"')\n")
                file.write("Set('yData', '"+self.runplot.Text_label[3]+"')\n")
                file.write("Set('key','"+ self.runplot.Text_label[3]+"')\n")
                file.write("To('..')\n")      
            file.write("Add('key', name='key1', autoadd=False)\n")
            file.write("To('..')\n")
            file.write("To('..')\n")
            file.write("\n")
            file.write("\n")
            file.write("\n")       
                    
                    
    def write_colorplot_files(self):
        cmd=self.MyChapApp.plotchoice
        Nqstep=int(self.calc.UpperqLimit / self.calc.Stepsize_qplot)
        with open(self.OutputFile, "a", encoding="utf-8") as file:
            file.write("#"+self.MyChapApp.plotchoice) 
            file.write("\n#==========start calculation  output section==========\n")
            if cmd == "dcs_omega_eq_plot":
               stepsize_b =  self.calc.thetamax.value / Nqstep /1000  # in mrad  
               file.write("theta [mrad] \t omega (eV) \t DDCS\n")                
            elif cmd =="eq_plot" :
               stepsize_b =  self.calc.Stepsize_qplot
               file.write("q [A.U.] \t omega (eV) \t ELF\n")
            else:
               stepsize_b =  self.calc.Stepsize_qplot
               file.write("q [A.U.] \t omega (eV) \t SELF\n")    
            current_b = 0.5 *stepsize_b  # initial value
                
            for i in range( Nqstep):
               for j in range(self.calc.NPoints):
                   current_E = self.calc.CenterFirstBin + j * self.calc.Stepsize
                   a = "{:.3f} \t{:.4E} \t{:.4E}\n".format(
                       current_b, current_E, self.calc.my_image[j, i]
                   )
                   file.write(a)
               current_b += stepsize_b
        # end Simon's patch 06/12/2020
        # also write matrix for imagej
        self.save_in_anotherformat()
                     
    def save_in_anotherformat(self):
        print("in sub")
        file_dialog = QFileDialog(self.MyChapApp)
        file_dialog.setWindowTitle("Save image as matrix")
        file_dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
        file_dialog.setNameFilter("matrix Files (*.eq);;All Files (*.*)")
        file_dialog.setDefaultSuffix("eq")  # Set the default extension to .eq
        file_dialog.exec()
        try:
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                self.filename = file_dialog.selectedFiles()
             
            myfilename=  selected_files[0]
            self.current_dir=os.path.dirname(str(myfilename))
        except Exception:
            return   
        if myfilename != '':  # i.e. dialog was not cancelled
            with open(myfilename, "w", encoding="utf-8") as file:
                file.seek(0)  # go to the beginning in case the file existed already
                np.savetxt(file, self.calc.my_image, delimiter=",", fmt="%.4e")
    
    def add_REELS_details(self): 
        with open(self.OutputFile, "a", encoding="utf-8") as file:    
            file.write("#projectile energy: " + str(self.calc.E0) + "keV\n")
            file.write("#Energy res: " + str(self.calc.Eres) + "eV\n")
            file.write( "#Theta in: "
                + str(self.calc.thetaIn) + ", Theta Out: "
                + str(self.calc.thetaOut) + "\n"
            )
            file.write( "#C1: " + str(self.calc.coef1)
                + ", C2: "     + str(self.calc.coef2)
                + ", C3: "     + str(self.calc.coef3) + "\n"
            )
            file.write( "#surface excitation adjustment factor "
                + str(self.calc.surf_ex_factor) + "\n" )
            file.write("#Fraction DIIMFP in calculated energy range "
                + str(self.calc.fraction_DIIMFP) + "\n" )
            if self.calc.DSEP_choice == 0:
                file.write("#Dsep from global eps\n")
            else:
                file.write("#Dsep per oscillator\n")  
                
                        
    def write_IMFP_stop_strag_plus_approx(self):
         with open(self.OutputFile, "a", encoding="utf-8") as file:    
            file.write(
                "#equivalent parameter single oscillator for IMFP: C = "
                + f"{self.calc.C0 :.5f}"
                + ", W_p = "
                + f"{self.calc.I0 :.1f}"
                + "eV\n"
            )
            file.write(
                "#equivalent parameter single oscillator for stopping: C = "
                + f"{self.calc.C1:.5f}"
                + ", W_p = "
                + f"{self.calc.MIE :.1f}"
                + "eV\n"
            )
            if np.amax(self.calc.DL_StragglingEnergy) > 0.0:
                file.write(
                    "equivalent parameter single oscillator for straggling: C ="
                    + f"{self.calc.C2:.5f}"
                    + ", W_p = "
                    + f"{self.calc.Istraggling:.1f}"
                    + "eV\n"
                )
            else:
                file.write("#No meaningful single oscillator for straggling\n")
            file.write(
                "#cross section per atom (angstrom^2) is (1.0/imfp)/"
                + f"{self.calc.UnitCellDensity:.3f}"
                + "\n"
            )

            file.write(
                self.runplot.xlabel
                + "\t"
                + self.runplot.Text_label[0]
                + "\t"
                + self.runplot.Text_label[1]
                + "\t"
                + self.runplot.Text_label[2])
             
            file.write(    
                "\tIMFP Bethe Disp. \tStop.  Bethe Disp.  \tStrag.  Bethe Disp.  "
                + "\tDLIMFP average_comp \tDL stopping average_comp \tDL straggling average comp "
                + "\tDL IMFP sum osc. \tDL stopping sum osc. \tDL strag. sum osc."
                + "\tDL IMFP from ELF\tDL stopping from ELF \tDL strag. from ELF"
                + "\t IMFP tpp2m\n")
            
            

            for i in range(self.calc.NStopping):
                a = ("{:.4E}\t"*17).format(
                    self.calc.x_axis[i],
                    self.calc.Result1[i],
                    self.calc.Result2[i],
                    self.calc.Result3[i],
                    self.calc.BetheIMFPEnergy[i],
                    self.calc.BetheStoppingEnergy[i]*self.runplot.StoppingFactor,
                    self.calc.BetheStragglingEnergy[i],
                    self.calc.DL_IMFPaverage_Energy[i],
                    self.calc.DL_StoppingEnergy[i]*self.runplot.StoppingFactor,
                    self.calc.DL_StragglingEnergy[i],
                    self.calc.DL_IMFP_sum_Energy[i],
                    self.calc.DL_Stopping_sum_Energy[i]*self.runplot.StoppingFactor,
                    self.calc.DL_Straggling_sum_Energy[i],
                    self.calc.DL_IMFP_from_ELF[i],
                    self.calc.DL_Stopping_from_ELF[i]*self.runplot.StoppingFactor,
                    self.calc.DL_Straggling_from_ELF[i],
                    self.calc.TPP_IMFPEnergy[i]
                )
                file.write(a+"\n")
           
                
          
                
    def write_partial_calc_results(self):
        with open(self.OutputFile, "a", encoding="utf-8") as file:
            a = self.runplot.xlabel
            for j in range(10):
                        a += "\t"+self.runplot.my10labels[j]
            file.write(a+"\n")  
            for i in range(self.calc.NPoints):
                a ="{:.3f}".format(self.calc.x_axis[i])
                for j in range(10):
                    a += "\t{:.4E}".format(self.calc.partialresults[i,j])
                a += "\n"
                file.write(a)            
        
    # def save_in_EQBetheformat(self):  # not sure where this routine is useful for (currently not called)

        # filename = filedialog.asksaveasfilename(
            # initialdir=os.getcwd,
            # title="write file containing matrix for imagej",
            # filetypes=(("dat files", "*.dat"), ("all files", "*.*")),
        # )
        # with open(filename, "w", encoding="utf-8") as file:
            # file.seek(0)  # go to the beginning in case the file existed already
            # np.savetxt(file, self.my_image, delimiter=",", fmt="%.4e")
            
    def load_comp_data(self):
        file_dialog = QFileDialog(self.MyChapApp)
        file_dialog.setWindowTitle("Open File with other theory or experiment")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setViewMode(QFileDialog.ViewMode.Detail)
        file_dialog.setNameFilter("any file (*.*)")
        if (file_dialog.exec()):
            fileNames = file_dialog.selectedFiles()
        try:
            os.path.basename(fileNames[0])
        except Exception:
            return 2
        self.runplot.xCompArray=[]
        self.runplot.yCompArray=[]          
        with open(fileNames[0], "r", encoding="utf-8") as file:
            line = file.readline()  #first line always taken to be comment
            while True:
                line = file.readline()
                if not line:
                    break
                if line[0] == "#":
                    print("comment line:",line)
                else:    
                    line.replace(","," ")  # change a comma delimeted file to a space delimited.
                    line.replace("\t"," ") ## change a tab delimeted file to a space delimited.
                    nums = line.split()  
                    if len(nums) < 2:
                        break
                    self.runplot.xCompArray.append(float(nums[0]))
                    self.runplot.yCompArray.append(float(float(nums[1])))
                   
                    

             
