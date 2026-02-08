

import matplotlib
import numpy as np
import constants as cnst
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
matplotlib.use('QtAgg')
plt.ion()




class run_and_plot:
   # def __init__(self, *args, **kwargs):
    def __init__(self,parent):
        self.MyChapApp=parent
        self.calc = self.MyChapApp.calc
        self.CalcDescription = ""
        self.LiteratureDescription ="lit. data"
        self.DisplayGrid = False
        self.CompDataLine = True
 
        self.legendfontsize = 11
        self.axislabelsize = 12
        self.linewidth =2
        self.markersize = 4
        self.xlabelsize = 12  # was 12
        self.ylabelsize = 12
        self.title_fontsize = 12    
        self.fileformat = "pdf"
               
        self.figurewidth= 15 # in cm
        self.figureheight = 12
        self.legendposition= 0
 
        self.max_eq=0.0
        self.Overplot = False
        self.comp_option_choice = 0 # 0 IMFP 1 stoppping 2 straggling
        self.plot_imfp = True
        self.x_axis_keV = True  # keV else velocity
        self.StoppingUnits = 0 # 0: eV/Å 1:eV / (10¹⁵ atoms/cm²)  2:MeV / (mg/cm²)
        self.update_plot_settings()
        self.px_per_cm=100
        self.start_x, self.start_y, self.dx, self.dy = (50, 0, int(self.figurewidth*self.px_per_cm),  int(self.figureheight*self.px_per_cm)) 
   
        

   
#===================================================================
    def start_calc_and_plot(self):
        self.MyChapApp.UpdateStatus("") # reset message window
        errorcode = self.initcalc() 
        self.MyChapApp.UpdateCalcPlotBtn("calculating.....")
        if errorcode == 0:
            routine = getattr(self, self.MyChapApp.plotchoice)
            routine()
        self.MyChapApp.UpdateCalcPlotBtn("calculate and plot") 
        
    def replot(self):
        try:
            self.PlotDescription
        except Exception:
            print("no calculation done yet, nothing to replot")
            self.calc.ErrorMessage ="no calculation done yet, nothing to replot"
            return
        self.update_plot_settings() 
        plotchoice=self.MyChapApp.plotchoice
        if plotchoice == "IMFP_stop_strag":
            self.plot_result3()
        elif plotchoice == "eq_plot"  or  plotchoice == "self_eq_plot"  or plotchoice == "dcs_omega_eq_plot"\
            or  plotchoice == "dcs_omega_eq_plot" or  plotchoice =="Cerenkov" or  plotchoice == "difference_due_to_Cerenkov":
                self.calc.scale_image()
                self.colorplot_result() 
        elif plotchoice == "partial_DIIMFP" or plotchoice == "partial_stopping":   
            self.plot_result10()
        else:
            self.plot_graph()
                      
#-----------------------------------------------------------------------    
    def initcalc(self):
        self.calc.ErrorMessage = "no errors"
        if self.calc.particle == "proton":
            self.particle_LaTeX ="H⁺"
        else:
            self.particle_LaTeX ="e⁻"  
  
        self.ylabel = ""
        self.ylabel2 =""
        self.PlotDescription = ""
        self.LaTex_label = ["","","",""]
        self.Text_label =  ["","","",""]
        self.second_y=[False,False,False,False]
       
        self.calc.my_updateProjectileEnergy()
        
        self.update_plot_settings()
        error_code=self.calc.initParArray()
        return error_code
        
        
    def update_plot_settings(self): 

        theme = {'axes.grid':  self.DisplayGrid,
             'grid.linestyle': '--',
             'legend.framealpha': 1,
             'legend.facecolor': 'white',
             'legend.shadow': False,
             'legend.fontsize':self.legendfontsize,
             'legend.title_fontsize': 14,
             'xtick.labelsize': self.xlabelsize,
             'ytick.labelsize': self.ylabelsize,
             'axes.labelsize': self.axislabelsize,
             'lines.linewidth': self.linewidth,
             'lines.markersize': self.markersize,
             'figure.dpi': 100,
             'savefig.format':self.fileformat,
             'axes.titlesize':self.title_fontsize,
             'figure.autolayout':True,
             'legend.frameon':False,
             'axes.formatter.use_mathtext' : True,
             'axes.formatter.limits': [-3, 4] 
             }
         
        plt.rcParams.update(theme) 

        self.nl="\n "  #here python makes the newline so no raw    
        plt.rcParams['mathtext.fontset'] = 'stix'
        self.cm=1.0/2.54  
        self.px_per_cm= 100*self.cm
        plt.rcParams["figure.figsize"] = [self.figurewidth,self.figureheight] 
 
        self.dx = int(self.figurewidth*self.px_per_cm)  # in pixels
        self.dy = int(self.figureheight*self.px_per_cm)
   
        if self.CompDataLine:
            self.overplotline='dashed'
        else:
            self.overplotline='none'

    # ================================== start calculations=======================================
    def eps_w(self):
        self.calc.eps1eps2()
        self.xlabel =  "ω (eV)"
        self.omega_xaxis = True
        if self.calc.epsilon_chi_choice == 0:
            self.Text_label[0] = "Re[ε(ω,q=%s)]" %str(self.calc.q)
            self.Text_label[1] = "Im[ε(ω,q=%s)]" %str(self.calc.q)
        else:
            self.calc.Result1 -= 1
            self.Text_label[0] = "Re[χ(ω,q=%s)]" %str(self.calc.q)
            self.Text_label[1] = "Im[χ(ω,q=%s)]" %str(self.calc.q)
        self.PlotDescription =  ""
        self.plot_graph()
        
    def one_over_eps_w(self):  
        self.calc.oneovereps1eps2()
        self.xlabel = "ω (eV)"
        if self.calc.epsilon_chi_choice == 0:
            self.Text_label[0] = "Re[1/ε(ω,q=%s)]" %str(self.calc.q)
            self.Text_label[1] = "Im[-1/ε(ω,q=%s)]" %str(self.calc.q)
            
        else:
            for i in range(self.calc.NPoints): self.calc.Result1[i] -= 1
            self.Text_label[0] = "Re[1/χ(ω,q=%s)]" %str(self.calc.q)
            self.Text_label[1] = "Im[-1/χ(ω,q=%s)]" %str(self.calc.q)
        self.PlotDescription = ""
        self.plot_graph() 
        
    def eps_q(self): 
        self.calc.eps1eps2_q()
        self.xlabel = "q (a.u.)"
        if self.calc.epsilon_chi_choice == 0:
            # self.LaTex_label[0] =  r"$\mathrm{Re} \left[\epsilon (ω=%s, q)\right]$"%str(self.calc.Energy_qplot) 
            # self.LaTex_label[1] =  r"$\mathrm{Im} \left[\epsilon (ω=%s, q)\right]$"%str(self.calc.Energy_qplot) 
            self.Text_label[0] = "Re[ε(ω=%s,q)]" %str(self.calc.Energy_qplot)
            self.Text_label[1] = "Im[ε(ω=%s,q)]" %str(self.calc.Energy_qplot)
        else:
            for i in range(self.calc.NPoints): self.calc.Result1[i] -= 1
            # self.LaTex_label[0] =  r"$\mathrm{Re} \left[\chi (ω=%s, q)\right]$"%str(self.calc.Energy_qplot) 
            # self.LaTex_label[1] =  r"$\mathrm{Im} \left[\chi (ω=%s, q)\right]$"%str(self.calc.Energy_qplot) 
            self.Text_label[0] = "Re[χ(ω=%s,q)]" %str(self.calc.Energy_qplot)
            self.Text_label[1] = "Im[χ(ω=%s,q)]" %str(self.calc.Energy_qplot)
            
        self.ylabel = ""
        self.PlotDescription = ""
        self.plot_graph()      
        
    def one_over_eps_q(self):   
        self.calc.oneovereps1eps2_q()
        self.xlabel = "q (a.u.)"
        if self.calc.epsilon_chi_choice == 0:
            # self.LaTex_label[0] = r"$\mathrm{Re} \left[\frac{1}{\epsilon (ω=%s, q)}\right]$" %str(self.calc.Energy_qplot)
            # self.LaTex_label[1] = r"$\mathrm{Im}\left[\frac{-1}{\epsilon (ω=%s, q)}\right]$" %str(self.calc.Energy_qplot)
            self.Text_label[0] = "Re[1/ε(ω=%s,q)]" %str(self.calc.Energy_qplot)
            self.Text_label[1] = "Im[-1/ε(ω=%s,q)]" %str(self.calc.Energy_qplot)
        else:   
            for i in range(self.calc.NPoints): self.calc.Result1[i] -= 1
            # self.LaTex_label[0] = r"$\mathrm{Re} \left[\frac{1}{\epsilon (ω=%s, q)}\right]-1$" %str(self.calc.Energy_qplot)
            # self.LaTex_label[1] = r"$\mathrm{Im}\left[\frac{-1}{\epsilon (ω=%s, q)}\right]$" %str(self.calc.Energy_qplot) 
            self.Text_label[0] = "Re[χ(ω=%s,q)]" %str(self.calc.Energy_qplot)
            self.Text_label[1] = "Im[χ(ω=%s,q)]" %str(self.calc.Energy_qplot)
            
        self.ylabel = ""
        self.PlotDescription = ""
        self.plot_graph() 
        
    def eps_w_kk(self):   
        error=self.calc.eps_kk_test()
        if error != 0: return
        self.xlabel = "ω (eV)"
        self.LaTex_label[0] = r"$\mathrm{Re} \left[\epsilon (ω, q= %s)\right]$" %str(self.calc.q)
        self.LaTex_label[1] = r"$\mathrm{Im} \left[\epsilon (ω, q= %s)\right]$" %str(self.calc.q)
        self.LaTex_label[2] = r"$1 + \frac{2}{π} {\cal P}\int_0^{%s}\,\, \frac{ω'\, {\rm Im} \left[\epsilon (ω',\, q)\right]}{(ω')^2-ω^2}  dω' $"%str(self.calc.UpperELimit)
        self.LaTex_label[3] = r"$ -\frac{2 ω}{π}  {\cal P}\int_0^{%s}\,\, \frac{ {\rm Re} \left[\epsilon (ω',\, q)\right]-1}{(ω')^2-ω^2}  dω'$"%str(self.calc.UpperELimit)
        self.Text_label[0] = "Re[ε(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[1] = "Im[ε(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[2]= "KK-transform Im[ε(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[3]= "KK-transform Re[ε(ω,q=%s)]" %str(self.calc.q)
        self.PlotDescription = ""
        self.plot_graph()  
        
    def one_over_eps_w_kk(self):   
        error=self.calc.one_over_eps_kk_test()
        if error != 0: return
        self.xlabel = "ω (eV)"
        self.LaTex_label[0] =r"$\mathrm{Re} \left[\frac{1}{\epsilon (ω, q= %s)}\right]$" %str(self.calc.q)
        self.LaTex_label[1] = r"$\mathrm{Im}\left[\frac{-1}{\epsilon (ω, q= %s)}\right]$" %str(self.calc.q)
        self.LaTex_label[2] = r"$1 + \frac{2}{π} {\cal P}\int_0^{%s}\,\, \frac{ω'\, {\rm Im} \left[\frac{1}{\epsilon (ω',\, q)}\right]}{(ω')^2-ω^2}\,  dω' $"%str(self.calc.UpperELimit)
        self.LaTex_label[3] = r"$ -\frac{2 ω}{π}  {\cal P}\int_0^{%s} \,\, \frac{ {\rm Re} \left[\frac{1}{\epsilon (ω',\, q)}\right]-1}{(ω')^2-ω^2}\,  dω'$"%str(self.calc.UpperELimit)
        self.Text_label[0] = "Re[1/ε(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[1] = "Im[-1/ε(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[2]= "KK-transform Im[1/ε(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[3]= "KK-transform Re[-1/ε(ω,q=%s)]" %str(self.calc.q)
        self.PlotDescription = ""
        self.plot_graph()   
        
    def eq_plot(self):  
        self.xlabel = "--"
        self.calc.bulk_eq = True
        if self.calc.LogXY:
            self.PlotDescription = r"$\log\left({\rm{Im}\left[ \frac{-1}{\epsilon(ω,k)}\right]}\right)$"
        else:
            self.PlotDescription = r"$\rm{Im}\left[ \frac{-1}{\epsilon(ω,q)}\right]$"
        self.calc.colorplot_lossfunction()  
        self.calc.scale_image()  
        self.colorplot_result()      
        
    def n_k_kk(self):   
        error=self.calc.n_and_k_kk_test()
        if error != 0: return
        self.xlabel = "ω (eV)"
        self.LaTex_label[0] = r"$n(ω,q=%s)$"%str(self.calc.q)
        self.LaTex_label[1] =  r"$k(ω,q=%s)$"%str(self.calc.q)
        self.LaTex_label[2] = r"$1 + \frac{2}{π} {\cal P}\int_0^{%s}\,\, \frac{ω'\, k (ω',\, q)}{(ω')^2-ω^2}\,  dω' $"%str(self.calc.UpperELimit)
        self.LaTex_label[3] = r"$ -\frac{2 ω}{π}  {\cal P}\int_0^{%s} \,\, \frac{ n(ω',\, q)-1}{(ω')^2-ω^2}\,  dω'$"%str(self.calc.UpperELimit)
        self.Text_label[0] = "n(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[1] = "k(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[2]= "KK-transform k(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[3]= "KK-transform n(ω,q=%s)]" %str(self.calc.q)
        self.PlotDescription = ""
        self.plot_graph()   
 
            
    def n_k(self): 
        self.calc.n_and_k_from_eps1_eps2()
        self.xlabel = "ω (eV)"
        # self.LaTex_label[0] = r"$n(ω,q=%s)$"%str(self.calc.q)
        # self.LaTex_label[1] = r"$k(ω,q=%s)$"%str(self.calc.q)
        self.Text_label[0] = "n(ω,q=%s)]" %str(self.calc.q)
        self.Text_label[1] = "k(ω,q=%s)]" %str(self.calc.q)
        self.PlotDescription = ""
        self.plot_graph()   
        
    def sum_rules(self):  
        self.calc.sum_rules()
        self.xlabel = "ω (eV)"
        self.LaTex_label[0] = r"$\frac{2}{π Ω_p^2}\int_0^ω ω \, {\rm Im} [\epsilon (ω, q=%s)] dω$"%str(self.calc.q)
        self.LaTex_label[1] = r"$\frac{2}{π Ω_p^2}\int_0^ω ω \, {\rm Im} [-1/\epsilon (ω, q=%s)] dω$"%str(self.calc.q)
        self.LaTex_label[2] = r"$\frac{4}{π Ω_p^2} \int_0^ω ω \, k(ω,q=%s) dω$"%str(self.calc.q)
        self.LaTex_label[3] = r"$\frac{2}{π }\int_0^ω \frac{1}{ω}$" + \
            r"$\, {\rm Im} [-1/\epsilon (ω, q=%s)] dω$"%str(self.calc.q)
        self.Text_label[0] ="F sum rule q = %s"  %str(self.calc.q)  
        self.Text_label[1] ="Bethe sum rule q = %s"  %str(self.calc.q)  
        self.Text_label[2] ="sum rule for k, q = %s"  %str(self.calc.q)  
        self.Text_label[3] ="KK (perfect screening) sum rule , q = %s"  %str(self.calc.q)  
        self.ylabel = "electrons per unit cell "
        self.ylabel2 = "KK sum rule"
        self.PlotDescription = ""
        self.second_y[3] = True
        self.plot_graph()  
        
        
    def inertial_rules(self):  
        self.calc.inertial_rules()
        self.xlabel = "ω (eV)"
        self.LaTex_label[0] = r"$\int_0^ω  (Re \left[\frac{1}{\epsilon (ω')}\right]-1) dω' $"
        self.LaTex_label[1] = r"$\int_0^ω  (n(ω')-1) dω' $"
        self.LaTex_label[2] = r"$\lim_{ω \to 0}\,\,  \frac{ω}{4 π}\,  \epsilon_2 (ω,q) $" 
        self.LaTex_label[3] = r"$\frac{-1}{2 π^2}\int_0^ω  (\epsilon_1(ω')-1) dω' $"
        self.Text_label[0] ="Inertia sum rule  Re[1/ε(ω)]-1, q = %s"  %str(self.calc.q)  
        self.Text_label[1] ="Inertia sum rule  n(ω)-1, q = %s"  %str(self.calc.q)  
        self.Text_label[2] ="limit at %s "  %str(self.calc.q)  
        self.Text_label[3] ="Inertia sum rule  Re[ε(ω)]-1, q = %s"  %str(self.calc.q)  
       
 
        self.ylabel = r"inertial sum rule (eV$^{-1})$ "
        self.PlotDescription = "evaluated at q={} a.u.".format(self.calc.q)
        self.plot_graph()      
        
    def diimfp(self):    
        self.calc.calcDIIMFP()
        self.xlabel = self.calc.particle + " energy loss ω (eV)"
        self.ylabel =  r"DIIMFP (eV$^{-1}$ Å$^{-1}$)"
        if self.calc.E0 > 1: self.Text_label[0] = "DIIMFP, E₀={:.1f} keV".format( self.calc.E0)
        else:  self.Text_label[0] = "DIIMFP, E₀={:.1f} eV".format(1000* self.calc.E0)
        self.PlotDescription = r"IMFP: {:.3g} Å".format(float(self.calc.StoppingResultArray[0]))
        #self.PlotDescription = r"IMFP: {:.1g}$ \rm{\AA}$".format(float(self.calc.StoppingResultArray[0]))
        self.PlotDescription +=r", stopping:  {:.3g} eV/Å".format(float(self.calc.StoppingResultArray[1]))
        self.PlotDescription +=  r", straggling: {:.3g} eV$^2$/Å".format(float(self.calc.StoppingResultArray[2]))
        self.plot_graph()
        
    def IMFP_stop_strag(self): 
        
        self.calc.calccurves(False)
 

    #    self.prepare_for_writing_plotting()
        self.plot_result3() 
        
    def prepare_for_writing_plotting(self):
        for i in range(len(self.calc.x_LinearApprox_lowE)-1):
            self.calc.x_LinearApprox_lowE[i+1]=self.calc.CurvesEnergy[i]
        if self.plot_imfp:
            self.calc.Result1=self.calc.IMFPEnergy
        else:
            self.calc.Result1 =self.calc.CrosssectionEnergy
        self.calc.Result2 = self.calc.StoppingEnergy * self.StoppingFactor
        self.calc.Result3 = self.calc.StragglingEnergy       
        
 
        
    def self_eq_plot(self):  #self_ means surface energy loss function here!
        self.xlabel = "--"
        self.calc.bulk_eq = False  
        if self.calc.LogXY:
            self.PlotDescription = r"$\log \rm{Im} \frac{(\epsilon-1)^2}{\epsilon (\epsilon+1)}$"
        else:
            self.PlotDescription = r"$\rm{Im} \frac{(\epsilon-1)^2}{\epsilon (\epsilon+1)}$"
        self.calc.colorplot_lossfunction()
        self.calc.scale_image()    
        self.colorplot_result()  
        
    def Mean_Excitation_Energy(self):   
        self.calc.Mean_Excitation_Energy()
        self.xlabel = "ω (eV)"
        self.LaTex_label[0] =  r"$I_0(ω)= \frac{\int_0^ω \,\log{ω'}\, {\rm Im}" + \
            r"[-1/\epsilon (ω',0)] dω'}{\int_0^ω {\rm Im} [-1/\epsilon (ω',0)] dω'}$" + self.nl
        a= r"($I_0$({:.0f}) =  {:.2f} eV, $C_0=$ {:.2E})".format(self.calc.UpperELimit,  self.calc.I0, self.calc.C0)
        self.LaTex_label[0] += a
        
        self.LaTex_label[1] =  r"$I_1(ω)= \frac{\int_0^ω \, ω'\log{ω'}\, {\rm Im}" + \
            r"[-1/\epsilon (ω',0)] dω'}{\int_0^ω   ω'"+ \
            r"{\rm Im} [-1/\epsilon (ω',0)] dω'}  $" + self.nl
        a= r"($I_1$({:.0f}) =  {:.2f} eV, $C_1=$ {:.2E})".format(self.calc.UpperELimit,  self.calc.MIE, self.calc.C1)
        self.LaTex_label[1] += a
        self.Text_label[0]="mean ionization  energy for IMFP"
        self.Text_label[1]="mean ionization energy for stopping" 
        self.ylabel = "Mean Ionization Energy (eV)"
        self.PlotDescription = ""
        self.plot_graph() 
    
    def SELF_DSEP(self):
        self.calc.surfaceloss()
        self.xlabel =self.calc.particle + " energy loss(eV)"
        self.Text_label[0] = "surf Loss Func"
        self.Text_label[1] =  "DSEP, θ={:.1f}⁰, ".format(self.calc.thetaIn)+"integr.prob.:{:.3g}".format(
            self.calc.SurfExProb) 
        self.PlotDescription = "SELF evaluated at q={} a.u.".format( self.calc.q)
        self.ylabel  = "Surf. loss function "
        self.ylabel2 = "DSEP at E₀={:.1f}".format(self.calc.E0)+"keV"
        self.second_y[1] = True
        self.plot_graph()
        
    def Compton(self):  
        self.calc.CalcCompton()
        self.xlabel = "k (a.u.)"
        self.Text_label[0] =  "Compton profile \nq=%s a.u."%str(self.calc.q_Compton) 
        self.PlotDescription = ""
        self.ylabel = "electrons per a.u. per atom"
        self.plot_graph()  
     
    def Fresnel_at_E(self):
        self.calc.Fresnel_at_E()   
        self.xlabel = "angle (deg.)"
        self.Text_label[0] = "Rₛ"
        self.Text_label[1] = "Rₚ"
        self.second_y[2] = True
        self.Text_label[2] = "Ψ"
        self.Text_label[3] = "Δ"
        self.PlotDescription = "Fresnel coefficients ℏω= {:.2f} eV ({:.3g} μm)".format(self.calc.E_Fresnel,1.2389841/self.calc.E_Fresnel)
        self.ylabel = "Reflectance Rₛ, Rₚ"
        self.ylabel2 = "Ψ, Δ (⁰)"
        self.plot_graph()  
    
    def Fresnel_at_angle(self):
        self.calc.Fresnel_at_angle()   
        self.xlabel = "ω (eV)"
        self.Text_label[0] = "Rₛ"
        self.Text_label[1] = "Rₚ"
        self.second_y[2] = True
        self.Text_label[2] = "Ψ"
        self.Text_label[3] = "Δ"
        self.PlotDescription = "Fresnel coeff. at ϕ={:.2f}°".format(self.calc.phi_ellipsometry)
        self.ylabel = "Reflectance Rₛ, Rₚ"
        self.ylabel2 = "Δ, Ψ (⁰)"  
        self.plot_graph()    
        
        
    def xray_absorption(self):
        print("x-ray abs.")  
        self.calc.xray_absorption() 
        self.xlabel = "photon energy (eV)" 
        self.Text_label[0] = "μ (cm²/g)"
        self.PlotDescription = "Photon absorption coefficient"
        self.plot_graph()
        
    def partial_DIIMFP(self):
        self.calc.weight = 0.0
        self.calc.PartDIIMFP()
        self.ylabel =  "DIIMFP (eV⁻¹Å⁻¹)"
        self.xlabel = "ω (eV)"
        self.PlotDescription = "partial DIIMFP, {}, E₀={:.1f}keV".format(self.particle_LaTeX,self.calc.E0)
        self.plot_result10() 
        
    def partial_stopping(self):
        self.calc.weight = 1.0
        self.calc.PartDIIMFP()
        self.xlabel = "ω (eV)"
        self.PlotDescription = r"partial stopping, {}, E₀={:.1f} keV".format(self.particle_LaTeX,self.calc.E0)
        self.plot_result10()
        
    def oscillator_strength(self):
        self.calc.Calc_Os_Strength()
        self.xlabel = "ω (eV)"
        self.Text_label[0] = "Osc. Strength from Im[(ϵ(q,ω)]"
        
        self.PlotDescription = "evaluated at q={} a.u.".format( self.calc.q)
        self.ylabel = "eV⁻¹"
        self.plot_graph()
        
        
    def S_k_omega_rule(self):
        self.calc.dyn_struct_factor()
        self.xlabel = "ω (eV)"
        self.Text_label[0]= "S(k,ω)"
        self.LaTex_label[1] = r"$\int_0^ω S(k,ω) ω dω$"
        self.Text_label[1]= "sum rule S(k,ω)"
        self.PlotDescription = "evaluated at q={} a.u.".format( self.calc.q)
        self.ylabel = "S(k,ω), eV⁻¹"
        self.ylabel2  = r"sum rule"
        self.second_y[1] = True
        self.plot_graph()        
        
    def pseudo_charge_density(self):
        self.calc.PseudoChargeDensity()
        self.xlabel = "fraction of unit cell"
        self.Text_label[0] = "minimum charge density"+self.nl+r"from ELF up to " + str(self.calc.UpperELimit) + " eV"
        self.PlotDescription = ""
        self.ylabel = "pseudocharge density (e⁻/Å³)"
        self.plot_graph()
        
    def radial_charge_density(self):   
        self.calc.PseudoChargeDensity()
        self.calc.ConvertToRadialPseudoChargeDensity()
        self.xlabel = "r (Å)"
        self.Text_label[0] = r" pseudo charge density (e⁻/Å³)"+self.nl+r"from ELF up to %s"%str(self.calc.UpperELimit)+" eV"\
           +self.nl+"MT Radius %6.2f Å"%self.calc.MT_radius
        self.PlotDescription = ""
        self.ylabel ="pseudocharge density (e⁻/Å³)"
        self.plot_graph()    
         
    def radial_w_p_lambda(self):
        self.calc.PseudoChargeDensity()
        self.calc.ConvertToRadialPseudoChargeDensity()
        self.calc.stopping_IMFP_w_p_versus_r()
        self.PlotDescription = ( self.particle_LaTeX + ", E₀={:.1f} keV, ".format(self.calc.E0)+ "from ELF up to " + str(self.calc.UpperELimit)+" eV")
        self.xlabel = "r (Å)"
        self.Text_label[0] = "ω_p (eV)"
        self.Text_label[1] = "IMFP λ (Å)"
        self.Text_label[2] = "stopping (eV/Å)"
        self.Text_label[3] = "straggling (eV²/Å)"
        self.plot_graph()    
        
    def shell_effect_all(self):
        self.calc.shell_effect(False)
        
        if self.x_axis_keV:
            self.calc.x_axis = self.calc.CurvesEnergy
            self.xlabel = self.calc.particle + " energy (keV)"
        else:
            self.calc.x_axis =  self.calc.CurvesVelocity  
            self.xlabel = self.calc.particle +" velocity (a.u.)"  
        self.PlotDescription = "shell effect (I = {:.3f} eV)".format(self.calc.MIE)
        
        self.Text_label[0] = "A = L Diel. Func."
        self.LaTex_label[1] = r"B =  $\max(\ln{\frac{2v²}{I}}+ \ln{γ²}- β² + 0.5f(γ), 0)$"
        self.Text_label[1] = "B (Bethe limit)"
        self.Text_label[2] = "shell effect: A - B (≈ -C/Z₂)"
        self.Text_label[3] =  "0.5f(γ)" 
        self.plot_graph()
        
    def calc_projectile_range(self):
        self.calc.projectile_range()
        self.PlotDescription = self.particle_LaTeX+ " range (CSDA approx.)"
        self.Text_label[0] = "range (Å)"
        self.xlabel = self.calc.particle + " energy (keV)"
        self.ylabel=""
        self.plot_graph()
        
    def Energy_Deposition_Depth(self):
        self.calc.Energy_Depth_Dist()
        self.PlotDescription = self.particle_LaTeX+", "+ str(self.calc.E0)+" keV" 
        self.Text_label[0] = "energy deposition (eV/Å), no straggling"
        self.Text_label[1] = "energy deposition (eV/Å), incl. straggling"
        self.xlabel = self.calc.particle + " depth (Å)"
        self.ylabel="eV/Å"
        self.plot_graph()
            
            
        
    # def shell_effect_soft(self):
        # self.calc.shell_effect(True)
        # if self.x_axis_keV:
            # self.xlabel = self.calc.particle + " energy (keV)"
        # else:
            # self.xlabel = self.calc.particle + " velocity (a.u.)"
            # for i in range (self.calc.NStopping):
                # self.calc.x_axis[i]= self.calc.CurvesVelocity[i]

        # self.PlotDescription = "shell effect, soft col., ($I =$ {:.3f} eV)".format(self.calc.MIE)
        # self.LaTex_label[0] = "L Bethe - L calc"
        # self.LaTex_label[1] = "L Bethe rel. - L calc"
        # self.plot_graph()
        
    def DDCS_at_omega(self): 
        self.calc.DDCS_at_omega()
        self.PlotDescription = r"$ \frac{dσ}{dω dΩ}$"+" {:.0f} keV,{} ".format(self.calc.E0, self.particle_LaTeX)
        self.ylabel = "DDCS  Å²/sr/eV, per unit cell"
        self.xlabel = "θ (mrad)"
        self.Text_label[0]= "ω =  {:.1f} eV, no retardation".format(self.calc.omega_ddcs)
        self.Text_label[1]= "ω =  {:.1f} eV, with retardation".format(self.calc.omega_ddcs)
        self.plot_graph() 
        
    def DDCS_at_theta(self):  
        self.xlabel = "ω (eV)"
        self.ylabel = "DDCS per unit cell  (Å²/sr/eV)" 
        self.Text_label[0] = "θ =  {:.3g} mrad, no retardation".format(self.calc.theta_ddcs)
        self.Text_label[1] = "θ =  {:.3g} mrad, incl. retardation".format(self.calc.theta_ddcs)
        self.PlotDescription = (self.particle_LaTeX + " {:.0f} keV,".format(self.calc.E0)) 
        self.calc.DDCS_at_theta()
        self.plot_graph()         
        
    def dcs_omega_eq_plot(self):
        self.calc.colorplot_ddcs()
        self.xlabel = "--"

        if not self.calc.LogXY:
            self.PlotDescription = (
                r"$ \frac{dσ}{dω dΩ}$" + "  ( Å²/sr/eV, per unit cell )  {:.1f} keV,".format(self.calc.E0)
                +self.particle_LaTeX
            )
        else:
            self.PlotDescription = (
               r"$ \log (\frac{dσ}{dω dΩ})$" + "  ( Å²/sr/eV, per unit cell ) {:.1f} keV,".format(self.calc.E0)
                + self.particle_LaTeX
            )
        self.calc.scale_image()    
        self.colorplot_result()
        
    def Cerenkov(self):
        self.calc.colorplot_Cerenkov()
        self.xlabel = "--"
 
        if not self.calc.LogXY:
            self.PlotDescription = (
                r"$ \frac{dσ}{dω dΩ}$"+" ( Å²/sr/eV, per unit cell ) {:.1f} keV,".format(self.calc.E0)
                +self.particle_LaTeX
            )
        else:
            self.PlotDescription = (
               r"$ \log (\frac{dσ}{dω dΩ})$"+" ( Å²/sr/eV, per unit cell ) {:.1f} keV,".format(self.calc.E0)
                + self.particle_LaTeX
            )
        self.calc.scale_image()    
        self.colorplot_result()
        
    def  difference_due_to_Cerenkov(self):
        self.calc.difference_due_to_Cerenkov()
        self.xlabel = "--"
        if not self.calc.LogXY:
            self.PlotDescription = ( self.particle_LaTeX +
                r", $ \frac{dσ}{dω dΩ}$ Cerenkov  $\frac{dE}{dx}=$" + "  {:.3g} ".format(self.calc.stopping_due_to_photons ) + r" eV/$\rm{\AA}$")
        else:
            self.PlotDescription = (self.particle_LaTeX +
               r"$, \log (\frac{dσ}{dω dΩ})$ Cerenkov  $\frac{dE}{dx}=$" + "  {:.3g} ".format(self.calc.stopping_due_to_photons ) + r" eV/$\rm{\AA}$")
        self.calc.scale_image()    
        self.colorplot_result()   
            
        

        
    def dcs_plot(self):
        self.calc.DCS()
        self.xlabel = "θ (mrad)"
        self.ylabel =   r"$  (\frac{dσ}{ dΩ})$  ( Å²/sr, per unit cell )" 
        self.Text_label[0] =  "DCS "+self.particle_LaTeX+ ", {:.0f} keV,".format(self.calc.E0)
        self.Text_label[1] =  "DCS incl. retardation"
        self.Text_label[2] = "Rutherford" 
        self.PlotDescription = "σ (this θ range, ω<{:.0f}eV): {:.2e}Å²".format(self.calc.UpperELimit,self.calc.sigma_from_DCS )
        self.plot_graph() 
        
       
        
  
        
    def REELS(self):
        if self.calc.particle == "proton":
            return
        self.calc.REELS_spectrum()
        self.xlabel = "ω (eV)"
        self.Text_label[0] = "Intensity"
        self.PlotDescription = "REELS Spectrum, {:.1f} keV,".format(self.calc.E0)
        self.plot_graph()
            
# ==========================plotting==========================================

    def plot_graph(self):
        x = self.calc.x_axis
        Local_xlabel= self.xlabel
        if "ω (eV)" in Local_xlabel:
            if self.calc.Energy_Scale_choice == 1:
                 Local_xlabel = "nm"
                 x= 1239.84/x
            elif  self.calc.Energy_Scale_choice == 2:
                 Local_xlabel = "cm⁻¹" 
                 x=8065.6*x
        if " energy (keV)" in Local_xlabel:  
            if self.x_axis_keV:
                if self.calc.particle== "electron":
                    Local_xlabel = "electron energy (keV)"
                else:
                    x = x/1000.0
                    Local_xlabel = "proton energy (MeV)"
            else:
                Local_xlabel = self.calc.particle + " velocity (a.u.)"

              
 
        fig, axs = plt.subplots(1, 1)
        if(self.title_fontsize> 0): axs.set_title(self.CalcDescription+" "+self.PlotDescription, pad=10)
        axs.tick_params(direction="in", which="both", right=1, top=1)
        axs.ticklabel_format(axis="both", style="sci", scilimits=(-3,4))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.11)
        plt.xlabel( Local_xlabel)
        plt.ylabel(self.ylabel)
        for Iplot in range(4):
             if self.LaTex_label[Iplot] =="": self.LaTex_label[Iplot] = self.Text_label[Iplot]
             
        if self.Text_label[0] != "": 
            plt.plot(x, self.calc.Result1, color= '#1f77b4',label=self.LaTex_label[0])
        if self.second_y[1]:
             plt.legend(loc = 'center right')
             ax2 = axs.twinx() 
             ax2.set_ylabel(self.ylabel2, rotation=270,labelpad=15)  
        if self.Text_label[1] != "": 
            plt.plot(x, self.calc.Result2, linestyle="dashed", color='#ff7f0e',label=self.LaTex_label[1])
        if self.second_y[2]:
             plt.legend(loc = 'center right')
             ax2 = axs.twinx() 
             ax2.set_ylabel(self.ylabel2, rotation=270,labelpad=15)  
        if self.Text_label[2] != "": 
            plt.plot(x, self.calc.Result3, linestyle="dotted", color='k',label=self.LaTex_label[2])  
        if self.second_y[3]:
             plt.legend(loc = 'center right')
             ax2 = axs.twinx() 
             ax2.set_ylabel(self.ylabel2, rotation=270,labelpad=15)  
           
        if self.Text_label[3] != "": 
            plt.plot(x, self.calc.Result4, linestyle="dashdot",linewidth=1, color= '#9467bd',label=self.LaTex_label[3])        
        if (self.Overplot == 1) and (len(self.xCompArray) > 0):
            if not self.calc.LogX:
                minimum= min(x)
                maximum= max(x)
                extrabit=(maximum-minimum)/50.0
                plt.xlim((minimum-extrabit,maximum+extrabit))
            plt.plot(self.xCompArray, self.yCompArray, label=self.LiteratureDescription, color="red",marker='.', linestyle=self.overplotline)

        if(self.calc.LogX == 1):
            plt.xscale("log")
        if (self.calc.LogY == 1):
            plt.yscale("log")
         

        if self.ylabel2 != "":
            plt.ylabel(self.ylabel2, rotation=270)
        else:
            plt.legend(loc = 'center right')
        if(self.legendfontsize > 0):plt.legend(loc='best')
        self.plotshow(plt)
    

 
    def plot_result3(self):
        if self.x_axis_keV:
            x = self.calc.CurvesEnergy
            self.xlabel = self.calc.particle + " energy (keV)"
        else:
            x= self.calc.CurvesVelocity  
            self.xlabel = self.calc.particle +" velocity (a.u.)"  
        if self.plot_imfp:
            Result1=self.calc.IMFPEnergy
            self.Text_label[0] = "IMFP (Å)"
        else:
            Result1 =self.calc.CrosssectionEnergy   
            self.Text_label[0] = "σ (Å²)"  
        if self.StoppingUnits==0:
            self.Text_label[1] = "stopping (eV/Å)"
            self.StoppingFactor=1.0
        elif self.StoppingUnits==1:
            self.Text_label[1] = "stopping (eV / (1E¹⁵atoms/cm²))"
            self.StoppingFactor=self.calc.massunitcell/(self.calc.specificweight*6.022)
        elif self.StoppingUnits==2:
            self.Text_label[1] = "stopping (MeV/(mg/cm²))"    
            self.StoppingFactor=0.1/self.calc.specificweight
        self.Text_label[2] = "straggling (eV²/Å)"    
        plt.rcParams["figure.figsize"] = [self.figurewidth,2*self.figureheight] 
        plt.rcParams["figure.autolayout"] = False
        self.dy=int(2*self.figureheight*self.px_per_cm)  # this one really controls the size
        fig, axs = plt.subplots(3, 1, sharex=True)

        fig.subplots_adjust(hspace=0)
        plt.figtext(0.02, 0.95, self.CalcDescription, fontsize=self.legendfontsize)

        if(self.title_fontsize> 0): 
            plt.figtext(0.5, 0.90, self.PlotDescription, horizontalalignment="center", fontsize=self.title_fontsize)
        
        for i in range(len(self.calc.x_LinearApprox_lowE)-1):
            self.calc.x_LinearApprox_lowE[i+1]=x[i]
                
        axs[0].plot(x, Result1, label= self.Text_label[0])
        if self.plot_imfp:
            axs[0].set_ylim(top=1.2 * self.calc.IMFPEnergy[self.calc.NStopping - 1])
            if self.calc.particle == "electron" and np.amax(self.calc.TPP_IMFPEnergy) > 0.0:
                axs[0].plot( x, self.calc.TPP_IMFPEnergy, color="purple",
                    label="TPP-2m ω_p="+ f"{self.calc.w_p_TPP:.2f}" + "eV\n\ρ="
                    + f"{self.calc.specificweight:.2f}" + " g/cm³",linestyle="dotted")
   

            if not self.calc.Approximations:
                axs[0].plot(x, self.calc.DL_IMFPaverage_Energy,  color="orchid",
                    label="Average DL IMFP: Cₐᵥ" + f"{self.calc.C0:.3f}" + ", I'ₐᵥ=" + f"{self.calc.I0 :.1f}"
                    + "eV", linestyle="-.", )

                axs[0].plot( x,self.calc.DL_IMFP_sum_Energy,color="black",
                    label="summed DL oscillators",linestyle="dotted")
                axs[0].plot( x, self.calc.DL_IMFP_from_ELF,color="limegreen",
                    label="DL IMFP from ELF", linestyle="dotted" )
                if np.amax(self.calc.TPP_IMFPEnergy) > 0.0:
                    axs[0].plot( x, self.calc.TPP_IMFPEnergy, color="purple",
                        label="TPP-2m ω_p=" + self.calc.Omegas[0]+ "eV\n$ ρ=$" + f"{self.calc.specificweight:.2f}"
                        + " g/cm³",linestyle="dotted")
        
     
        mylabel=r"$\frac{4π}{v²}$ N Z(L₀ + 0.5*F(γ))" + f"\n(I={self.calc.MIE:.1f}eV)"
        axs[1].plot(x, self.calc.StoppingEnergy*self.StoppingFactor, label=self.LaTex_label[1])
        axs[1].plot(x, self.calc.BetheStoppingEnergy_Salvat*self.StoppingFactor, color="firebrick",
            label=mylabel,
            linestyle="--")


        axs[1].plot(self.calc.x_LinearApprox_lowE,self.calc.LinearApprox_lowE*self.StoppingFactor, color="limegreen",
            label=r"$-\frac{dE}{dx}=c v$, $c =$"+ f"{self.calc.maxsloop:.1f} eV/(Åv₀)",
            linestyle="--")    
        
        if self.calc.Approximations:
            axs[1].plot(x, self.calc.DL_StoppingEnergy*self.StoppingFactor,
                color="orchid",label=r"DL Formula from ELF\nC₁= "+ f"{self.calc.C1 :.3f}" + ", ω_p=" + f"{self.calc.MIE :.1f} eV",
                linestyle="-.")
            if np.amax(self.calc.DL_Stopping_sum_Energy) > 0.0:
                axs[1].plot(x,self.calc.DL_Stopping_sum_Energy*self.StoppingFactor,
                    color="black",label="summed DL oscillators",linestyle="dotted")
            axs[1].plot(x, self.calc.DL_Stopping_from_ELF*self.StoppingFactor,
                color="limegreen",label="DL stopping from ELF",linestyle="dotted")
        axs[2].plot(x, self.calc.StragglingEnergy, label=self.LaTex_label[2])
        if self.calc.particle == "proton":
            scalingfactor=1.0
        else:
            scalingfactor=4.0    
        axs[2].plot([x[0], x[self.calc.NStopping - 1]], [self.calc.BohrStraggling/scalingfactor, 
            self.calc.BohrStraggling/scalingfactor], color="darkgreen", label="Bohr Limit", linestyle="dashdot")
        if self.calc.particle == "proton":
            axs[2].plot( x, self.calc.Straggling_Jackson, color="firebrick",
                label="Jackson eq. 13.50", linestyle="--")
                
            
        if self.calc.Approximations:
            axs[2].plot( x, self.calc.BetheStragglingEnergy, color="firebrick",
                label="Bethe Formula\n"+ "ω_p=" + f"{self.calc.MIE :.1f}" + "eV",
                linestyle="--")
            if ( np.amax(self.calc.DL_StragglingEnergy) > 0.0) and self.calc.OscillatorsPresent:
                axs[2].plot( x, self.calc.DL_StragglingEnergy, label="DL Formula from zero width osc.\n  C₂="
                    + f"{self.calc.C2 :.3f}" + ", ω_p=" + f"{self.calc.Istraggling :.1f}"+ "eV",  color="orchid",
                    linestyle="-.")
            elif np.amax(self.calc.DL_StragglingEnergy) > 0.0:
                axs[2].plot(x, self.calc.DL_StragglingEnergy, color="orchid", linestyle="-.",
                    label="DL Formula from ELF\n,  C₂=" + 
                    f"{self.calc.C2 :.3f}" + ", ω_p ="+ f"{self.calc.Istraggling :.1f}" + " eV")
            if np.amax(self.calc.DL_Straggling_sum_Energy) > 0.0:
                axs[2].plot(x, self.calc.DL_Straggling_sum_Energy, color="black",label="summed DL oscillators", 
                    linestyle="dotted" )
            axs[2].plot(x,self.calc.DL_Straggling_from_ELF, color="limegreen",label="DL straggling from ELF", 
                linestyle="dotted")

        if self.Overplot == 1:
            if not self.calc.LogX:
                minimum= min(x)
                maximum= max(x)
                extrabit=(maximum-minimum)/50.0
                plt.xlim((minimum-extrabit,maximum+extrabit))
            axs[self.comp_option_choice].plot(self.xCompArray, self.yCompArray,
                label=self.LiteratureDescription,linestyle=self.overplotline,marker='.', color="red")
          
        plt.xlabel(self.xlabel)
        for i in range(3):
            axs[i].legend()
            if self.calc.LogX == 1:
                axs[i].set_xscale("log")
            if (self.calc.LogY == 1):
                axs[i].set_yscale("log")    
            axs[i].tick_params(axis="both", which="both", direction="in", right=1, top=1)
            axs[i].set_ylim(bottom=0)
        
        axs[1].set_ylim(top=1.1 * np.nanmax(self.calc.StoppingEnergy)*self.StoppingFactor)
        self.plotshow(plt)
   

    def plot_result10(self):
        self.my10labels=[]
        fig, axs = plt.subplots(1, 1)
        axs.ticklabel_format(axis="both", style="sci", scilimits=(-3,4))
        if(self.title_fontsize> 0): 
            axs.set_title(self.PlotDescription+" "+self.CalcDescription)
      
        axs.tick_params(direction="in", right=1, top=1)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.11)
        q_step = self.calc.UpperqLimit / 10.0
        for i in range(10):
            q_lower = i * q_step
            q_upper = (i + 1) * q_step
            if self.MyChapApp.plotchoice == "partial_DIIMFP":
                mytext=r"{:.1f} $< q <$ {:.1f}, λ_p={:6.2e} Å".format(q_lower,q_upper,self.calc.PartIntSum[i])
            else:
                mytext=r"{:.1f} $< q <$ {:.1f}, S_p={:6.2e} eV/Å".format(q_lower,q_upper,self.calc.PartIntSum[i])
            self.my10labels.append(mytext)        
            plt.plot(self.calc.x_axis, self.calc.partialresults[:, i], label=self.my10labels[i])
        if(self.calc.LogX == 1): plt.xscale("log")
        if(self.calc.LogY == 1): plt.yscale("log")        
        if(self.legendfontsize > 0): plt.legend(loc=self.legendposition)
 
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        self.plotshow(plt)

    def colorplot_result(self):
        fig, axs=plt.subplots(1, 1)
        if(self.title_fontsize> 0): 
            axs.set_title(self.PlotDescription+" "+self.CalcDescription)
        
        if (self.MyChapApp.plotchoice == "eq_plot") or (self.MyChapApp.plotchoice == "self_eq_plot"):
            ratio = self.calc.UpperqLimit / (self.calc.UpperELimit - self.calc.CenterFirstBin)
           # maximum scale and log/linear are taken care of in calc.
            print("lower lim",self.calc.LowerELimit)
            plt.imshow(
                self.calc.my_scaled_image,
                extent=[0,self.calc.UpperqLimit,self.calc.UpperELimit,self.calc.LowerELimit],
                aspect= ratio*self.figureheight/self.figurewidth
            )
           
            if (self.MyChapApp.plotchoice == "eq_plot"):        
                plt.xlabel("q (a.u.)")
            else:
                plt.xlabel(r"$q_\parallel$"+"(a.u.)")    
            plt.xlim(0, self.calc.UpperqLimit)
        else:  # self.MyChapApp.plotchoice == "dcs_omega_eq_plot" or cerenkov
            ratio = self.calc.theta_max / self.calc.UpperELimit
            plt.imshow(self.calc.my_scaled_image,
                extent=[0, self.calc.theta_max, self.calc.UpperELimit, self.calc.CenterFirstBin],
                aspect=ratio)
            
            plt.xlabel("θ(mrad)")
            plt.xlim(0, self.calc.theta_max)
        plt.ylim(self.calc.CenterFirstBin, self.calc.UpperELimit)

        plt.colorbar()
        plt.ylabel("ω (eV)")

        if self.calc.E0 < 1000:
            mylabel = "upper limit ω\n E₀= {:.2f} keV, {}".format(self.calc.E0,self.particle_LaTeX)    
        else:
            mylabel =  "upper limit ω\n E₀= {:.2f} MeV, {}".format(self.calc.E0/1000,self.particle_LaTeX)    
        if self.MyChapApp.plotchoice == "eq_plot" and self.calc.E0 > 0.0:
            self.calc.Calculate_Integration_limits()
          #  plt.plot(self.calc.xArray, self.calc.yArray, color="lime", label=mylabel + " non-relativistic")  #uncomment this if you want the non-relativistic boundary as well
           
            plt.plot(self.calc.xArray, self.calc.yArray_relativistic, color="red", label=mylabel)
            if(self.legendfontsize > 0):
                plt.legend(loc=self.legendposition,frameon=True)  
        elif  self.MyChapApp.plotchoice == "dcs_omega_eq_plot"  and self.calc.particle == "proton":  # limiting lines for protons

            if self.calc.theta_max > 1000.0/cnst.Mp:
                plt.plot([1000/cnst.Mp, 1000/cnst.Mp], [self.calc.CenterFirstBin,self.calc.UpperELimit ],
                   color="red",linestyle="dashed",  label=r"$\sin(θ_{\rm max}) =\frac{M_e}{M_p}$")
                   #https://www.physicsforums.com/threads/solving-scattering-angle-problem-when-m-m-or-m-m.28662/ for derivation
                   #or Sigmund 3.1.1
                
            maxloss= 2*(self.calc.beta_r *self.calc.gamma_r* cnst.C)**2 * cnst.HARTREE
            if maxloss < self.calc.UpperELimit :  
                plt.plot([0, self.calc.theta_max], [maxloss,maxloss],
                   color="red",linestyle="dotted", label=r"$2v^2$ * 27.211")
                   
        plt.locator_params(axis='x', nbins=3)
        plt.locator_params(axis='y', nbins=3)
        self.plotshow(plt)
        
        
    def plotshow(self, plt):
        
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(self.start_x,self.start_y+30,self.dx, self.dy)

        self.start_y += self.dy+30  #60 for toolbar and titlebar
        if self.start_y + self.dy > self.screen_height:
            self.start_y=0
            self.start_x += self.dx
            if self.start_x + self.dx > self.screen_width:
                self.start_x=40
     
        plt.show()        


    def Close_all_plots(self):
        plt.close("all")
        self.start_x=40
        self.start_y=0
        
        

   
        

# ==============================end plotting======================================      
#from matplotlib documentation
#this seems buggy but useful
        
   ## ***************************************************************************
## * INTERACTIVE KEYMAPS                                                     *
## ***************************************************************************
## Event keys to interact with figures/plots via keyboard.
## See https://matplotlib.org/stable/users/explain/interactive.html for more
## details on interactive navigation.  Customize these settings according to
## your needs. Leave the field(s) empty if you don't need a key-map. (i.e.,
## fullscreen : '')
#keymap.fullscreen: f, ctrl+f   # toggling
#keymap.home: h, r, home        # home or reset mnemonic
#keymap.back: left, c, backspace, MouseButton.BACK  # forward / backward keys
#keymap.forward: right, v, MouseButton.FORWARD      # for quick navigation
#keymap.pan: p                  # pan mnemonic
#keymap.zoom: o                 # zoom mnemonic
#keymap.save: s, ctrl+s         # saving current figure
#keymap.help: f1                # display help about active tools
#keymap.quit: ctrl+w, cmd+w, q  # close the current figure
#keymap.quit_all:               # close all figures
#keymap.grid: g                 # switching on/off major grids in current axes
#keymap.grid_minor: G           # switching on/off minor grids in current axes
#keymap.yscale: l               # toggle scaling of y-axes ('log'/'linear')
#keymap.xscale: k, L            # toggle scaling of x-axes ('log'/'linear')
#keymap.copy: ctrl+c, cmd+c     # copy figure to clipboard     
        
 
        
