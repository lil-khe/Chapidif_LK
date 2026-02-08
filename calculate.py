# generic python modules
import os
import sys
import math
import cmath
#python routines that may not be installed in the generic Python distroibution 
import numpy as np
#from scipy.integrate import cumulative_trapezoid
import time
#my stuff
import constants as cnst
print("platform:",sys.platform)
if(sys.platform == 'linux'):
    # to compile epslib issue (on Linux) in the source directory:
    # cmake -S . -B LinuxLib
    # cmake --build LinuxLib
    from LinuxLib import epslib
elif(sys.platform == 'win32'): 
    # to compile the library epslib issue (on windows) in the source directory:
    # cmake  -G "MinGW Makefiles" -S . -B WindowsLib
    # cmake --build WindowsLib 

    from WindowsLib import epslib
elif(sys.platform == 'darwin'):   
    from MacLib import epslib
      

class calculate():
    def __init__(self,parent):
        self.MyChapApp=parent
        epslib.test()
        self.init()        
    def init(self):

#########start variables associated with DF, interfaced in tab1, and written to DF file=========== 
# varialbes determining array sizes.  changing these  will require similar changes in epslib        
        self.maxOscillators: int = 250
        self.maxGOS: int =10
        self.maxKaneko:int =250
        self.maxBelkacem:int =5

        self.massunitcell = 27  #in amu
        self.specificweight = 2.7 # gr/cm^3
        self.maxEnergyDensityEffect=0.0
        self.PreciseDenityCor =False
        self.ApplySumRuleToGOS=1
        self.Merminize = 1 #0 means Lindhard, 1 means Merminize, 2 Direct approach
        self.MottCorrection = 1 # 
        self.Dispersion_choice = 0  #0 simple dispersion, 1 full dispersion
        self.Dispersion_relativistic = 1  #0 non-relativistic dispersion, 1 relativisitic dispersion (all DF)
        self.delayed_dispersion = 0  
        self.Add_Doppler_Width  = 0  # 0 constant width, 1 add Doppler width ("classical models only")
        
        self.DFmodel="Mermin"
        self.AddELF: int = 1  # 1 add elf 0 add chi
        self.Stopping_calc_quality = 0
        self.ExchangeCorrection = False
        self.Exchange_as_in_SBethe = False   #False means use method from Ashleya
        self.max_q_considered_surface: float = 200
        self.DebugMode= False
#========this concludes variables associated with DF, interfaced in tab1, and written to DF file===========

#====================== now variables for tab2=======================
   
        self.LowerELimit=0.0
        self.Stepsize = 0.5  # stepsize
        
        self.UpperELimit = 100.0
        self.MaxNPoints=1000
        self.Energy_Scale_choice=0  # 0: eV, 1: nm 2: cm⁻¹
        self.epsilon_chi_choice=0
        self.theta_max = 0.6  #mrad
        self.omega_ddcs = 10.0
        self.theta_ddcs = 100.0
        self.NThetaStep = 100
        self.particle="electron"
        self.ProjectileMass = 1
        self.E0 = 50.0  # beam energy (keV)
        self.my_updateProjectileEnergy()  
        self.q =0.05  # momentum of E plots
        self.UpperqLimit = 3.0
        self.Stepsize_qplot = 0.03
        self.Energy_qplot = 0.5
        self.q_Compton = 60.0
        self.E_Fresnel = 1.0
        self.phi_ellipsometry=70.0
        
        self.NStopping = 50
        self.IncrFactor = 1.25
        self.c_transition = 0.005
        if (sys.platform == "win32" or sys.platform == "linux"):
            self.c_transition = 0.0001  # for replacement of Mermin etc. with DL
        self.BE_for_exchange = 50.0
        self.eps_bkg=1.0  # for tauc and DL
        self.n_infty=1.123  # for Forouhi_bloomer
        self.a_TL_an=0.001  # in eV!
        self.RadiativeLosses = 0 # 0 means no radiative losses else Bethe-Heitler approach 

        self.ErrorMessage="my error"
       
#---------------------for (R)eels--------------------        
        self.Eres           = 1.0
        self.coef1          = 0.0
        self.coef2          = 0.0
        self.coef3          = 0.0
        self.EELS_thickness = 1000 # thickness for transmission EELS (angstrom)
        self.EELS           = 0 #EELS=1 is TEELS,  0 means REELS

        self.thetaIn        = 0.0
        self.thetaOut       = 45.0
        self.surf_ex_factor = 1.0
        self.fraction_DIIMFP= 1.0
        
 #---------------various       
        self.LogX = False
        self.LogY = False
        self.LogXY = False # for color plot
        self.log_range = 5
        self.max_eq = 0.0
        
        self.xArray = np.zeros(30)
        self.yArray = np.zeros(30)
        self.yArray_relativistic = np.zeros(30)
        self.velocity_projectile(self.E0)  # sets beta and gamma
        self.Compton_k_limit = 30  # was 10
        self.StoppingResultArray = np.zeros(3) # for output diimfp calc lambda stopping straggling
        self.first_electron_energy: float =round( 0.5 * 0.5**2 * cnst.HARTREE / 1000.0,6)  # E (keV) corresponding to v=0.5 a.u.
        self.first_proton_energy: float =round(0.5 * 0.5**2 * cnst.Mp * cnst.HARTREE / 1000.0,3) # E (keV) corresponding to v=0.5 a.u.
        self.w_p_TPP = 15.0   #for TPP formula

        self.U_factor=0.0  # contribution of U to the peak position
        self.Q_Kaneko_transform = 1.0
        self.l_Kaneko_transform = 1
        self.Approximations= False
        self.N_oscillator_used = self.maxOscillators
        self.First_oscillator_transform = 0
        self.Last_oscillator_transform = self.maxOscillators - 1
  
#========================== for communicating with epslib=======================

        self.NDFPAR = 5 * self.maxOscillators + 4 * self.maxGOS + 4 * self.maxBelkacem + 6 * self.maxKaneko
# length of all the oscillaotor, Gos etc arrays filled in "fill_oscillators()"
      
        self.ParArray = np.zeros(self.NDFPAR + 100)
        self.PartIntSum = np.zeros(10)
        self.OOSEnergy = np.zeros(1000)
        self.OOS = np.zeros(1000)
        # self.ELF_wide_Energy = 1300  # should match with epslib
        # self.ELF_wide = 1300   # should match with epslib
    
        
#======================================================================        
    def ZeroDF(self):
        self.InitDF()
        self.MyChapApp.update_all_tables()

    def InitDF(self):  # init variables tab1
        self.Amps = [0.0] * self.maxOscillators
        self.Omegas = [0.0] * self.maxOscillators
        self.Alphas = [1.0] * self.maxOscillators
        self.Gammas = [0.0] * self.maxOscillators
        self.Us = [0.0] * self.maxOscillators
        
        self.ConcGOS = [0.0] * self.maxGOS
        self.EdgeGOS = [100.0] * self.maxGOS
        self.nlGOS = [21] * self.maxGOS
        self.ZGOS:int = [6] * self.maxGOS
        self.N_Kaneko = [0.0]*self.maxKaneko # number of electrons in a shell
        self.Q_Kaneko = [2.0]*self.maxKaneko  # characteristic momentum of each shell
        self.width_Kaneko = [4.0]*self.maxKaneko  # characteristic width of shell
        self.l_Kaneko = [0]*self.maxKaneko # angular momentum
        self.Edge_Kaneko = [5.0]*self.maxKaneko  # minimum excitation energy of that shell (Archubi&Arista ext.)
        self.gamma_Kaneko =[1.0]*self.maxKaneko # plasmon position fudge parameter
        self.Kaneko_choice = 0 # 0 means modified Kaneko
        
        self.Conc_Belkacem = [0.0] * self.maxBelkacem
        self.w_Belkacem =  [0.0] * self.maxBelkacem
        self.gamma_Belkacem =  [0.0] * self.maxBelkacem
        self.DF_prop_text="\n"
        self.SumAiText=""
    
    def initParArray(self):
        if self.DFmodel == "Drude":
            self.DFChoice = 1
        elif self.DFmodel == "DL":
            self.DFChoice = 2
        elif self.DFmodel == "Mermin":
            self.DFChoice = 3
        elif self.DFmodel == "Vlasov":
            self.DFChoice = 4
        elif self.DFmodel == "Tauc":
            self.DFChoice = 5
        elif self.DFmodel == "TL_an":
            self.DFChoice = 6
        elif self.DFmodel == "TL_Mermin":
            self.DFChoice = 7
        elif self.DFmodel == "FB":
            self.DFChoice = 8 
        elif self.DFmodel == "BB":
            self.DFChoice = 9   
        elif self.DFmodel == "OC":
            self.DFChoice = 10        
        else:
            print("something wrong",self.DFmodel, "is not implemented") 
        self.ErrorMessage="my error 2"
        error_code = self.fill_oscillators()
        if error_code == 0:
            error_code = self.fill_remainder()
        if error_code != 0: self.calc.ErrorMessage = "initParArray retuned error code:"+ str(error_code)
        return error_code    
            

    def fill_oscillators(self):
        if self.DFmodel =="Drude" or self.DFmodel == "Tauc" or self.DFmodel == "BB"or self.DFmodel == "OC":
            self.ParArray[0] = self.eps_bkg  # not sure if this is still used
        elif self.DFmodel == "TL_an":
            self.ParArray[0] = self.a_TL_an  # this is NOT eps_bkg,array element used for  different thing
        elif self.DFmodel == "FB":
            self.ParArray[0] = self.n_infty*self.n_infty  # this is eps_bk
        try:
            for i in range(self.maxOscillators):
                self.ParArray[5 * i + 1] = self.Amps[i]
                self.ParArray[5 * i + 2] = self.Omegas[i]
                self.ParArray[5 * i + 3] = self.Gammas[i]
                self.ParArray[5 * i + 4] = self.Alphas[i]
                self.ParArray[5 * i + 5] = self.Us[i]
        except ValueError as e:
            self.ErrorMessage += "Input error Oscillator no: " + str(i+1) + "\n"  + str(e) 
            return 1   
        offset = 5 * self.maxOscillators
        try: 
            for i in range(self.maxGOS):
                self.ParArray[4 * i + offset + 1] = self.ConcGOS[i]
                self.ParArray[4 * i + offset + 2] = self.ZGOS[i]
                self.ParArray[4 * i + offset + 3] = float(self.nlGOS[i])
                self.ParArray[4 * i + offset + 4] = self.EdgeGOS[i]
        except ValueError as e:
            self.ErrorMessage += "Input error GOS no: " + str(i+1) + "\n"  + str(e)  
            return 2      
        offset += 4 * self.maxGOS
        try:
            for i in range(self.maxBelkacem):
                self.ParArray[4 * i + offset + 1] = self.Conc_Belkacem[i]
                self.ParArray[4 * i + offset + 2] = self.w_Belkacem[i]
                self.ParArray[4 * i + offset + 3] = self.gamma_Belkacem[i]
            offset += 4 * self.maxBelkacem
        except ValueError as e:
            self.ErrorMessage += "Input error Belkacem no: " + str(i+1) + "\n"  + str(e)     
            return 3             
        try:
            for i in range(self.maxKaneko):
                self.ParArray[6 * i + offset + 1] = self.N_Kaneko[i]
                self.ParArray[6 * i + offset + 2] = self.Q_Kaneko[i]
                self.ParArray[6 * i + offset + 3] = self.width_Kaneko[i]
                self.ParArray[6 * i + offset + 4] = self.Edge_Kaneko[i]
                self.ParArray[6 * i + offset + 5] = float(self.l_Kaneko[i])
                self.ParArray[6 * i + offset + 6] = self.gamma_Kaneko[i]
            
        except ValueError as e:
            self.ErrorMessage += "Input error Kaneko no: " + str(i+1) + "\n"  + str(e)   
            return 4
        return 0                

    def fill_remainder(self):
        self.CenterFirstBin = self.LowerELimit + 0.5*self.Stepsize # first energy point of plot
        self.CenterLastBin = self.UpperELimit - 0.5*self.Stepsize # last energy point
        self.NPoints = int((self.UpperELimit - self.CenterFirstBin) / self.Stepsize) + 1
        if self.NPoints >self.MaxNPoints:
            self.stepEfixed=False
            self.NPoints = self.MaxNPoints
        else:
            self.stepEfixed= True    
        self.Nqstep = int(self.UpperqLimit / self.Stepsize_qplot)
      
        # max_q_considered: the maximum momentum that is used in the integration in terms of q_min
        #  precision:          relates to quanc8 integration. Larger values more precise but take a
        # little longer
        if self.Stopping_calc_quality == 0:
            max_q_considered = 3000
            self.lin_cont_deltaE = 0.08
            self.precision = 10.0

        elif self.Stopping_calc_quality == 1:
            max_q_considered = 5000
            self.lin_cont_deltaE = 0.02
            self.precision = 100.0
        elif self.Stopping_calc_quality == 2:
            max_q_considered = 1e99
            self.lin_cont_deltaE = 0.005
            self.precision = 800.0

        self.ParArray[self.NDFPAR + 1] = self.E0 * 1000.0
        # presicion in integration routine for DIIMFP etc.
        self.ParArray[self.NDFPAR + 2] = self.precision
        # maximum of qmax considered in terms of qmin for diimfp
        self.ParArray[self.NDFPAR + 3] = max_q_considered
        # default 500    # maximum of qmax considered in terms of qmin for dsep
        self.ParArray[self.NDFPAR + 4] = self.max_q_considered_surface
       
 
        self.ParArray[self.NDFPAR + 5] = self.NPoints
        self.ParArray[self.NDFPAR + 6] = self.CenterFirstBin
        self.ParArray[self.NDFPAR + 7] = self.Stepsize
        # used to be atom  density for GOS
        self.ParArray[self.NDFPAR + 8] = self.UnitCellDensity
        self.ParArray[self.NDFPAR + 9] = self.lin_cont_deltaE
        # lower and upper limit for diimfp integration. partial diimfp and
        # stopping manually overrides these settings
        self.ParArray[self.NDFPAR + 10] = 0
        self.ParArray[self.NDFPAR + 11] = 1e99#self.highest_momentum_considered_stopping
        self.ParArray[self.NDFPAR + 12] = self.UpperqLimit
        self.ParArray[self.NDFPAR + 13] = self.Stepsize_qplot

        self.ParArray[self.NDFPAR + 14] = self.c_transition
        self.ParArray[self.NDFPAR + 15] = float(self.ExchangeCorrection)
        self.ParArray[self.NDFPAR + 16] = float(self.Kaneko_choice)
        self.ParArray[self.NDFPAR + 17] = float(self.AddELF)
        self.ParArray[self.NDFPAR + 18] = float(self.ApplySumRuleToGOS)
        self.ParArray[self.NDFPAR + 19] = float(self.Merminize)
        self.ParArray[self.NDFPAR + 20] = float(self.Dispersion_choice)
        self.ParArray[self.NDFPAR + 21] = self.maxEnergyDensityEffect
        # next parameterts for reels spectrum calculation
        self.ParArray[self.NDFPAR + 22] = self.Eres  # fwhm here
        self.ParArray[self.NDFPAR + 23] = self.coef1
        self.ParArray[self.NDFPAR + 24] = self.coef2
        self.ParArray[self.NDFPAR + 25] = self.coef3
        self.ParArray[self.NDFPAR + 26] = self.thetaIn  # degree here
        self.ParArray[self.NDFPAR + 27] = self.thetaOut
        self.ParArray[self.NDFPAR + 28] = self.surf_ex_factor
        self.ParArray[self.NDFPAR + 29] = 0.0  # empty!

        self.ParArray[self.NDFPAR + 30] = self.fraction_DIIMFP
        if self.particle== "electron":
            self.ParArray[self.NDFPAR + 31] = 0.0
        else: 
            self.ParArray[self.NDFPAR + 31] = 1.0  
        self.ParArray[self.NDFPAR + 32] = float(self.Dispersion_relativistic)
        
        self.ParArray[self.NDFPAR + 33] = self.theta_max
        self.ParArray[self.NDFPAR + 34] = float(self.NThetaStep) 
        self.ParArray[self.NDFPAR + 35] = float(self.MottCorrection)
        self.ParArray[self.NDFPAR + 36] = self.BE_for_exchange  
        self.ParArray[self.NDFPAR + 37] = float(self.Exchange_as_in_SBethe) #0 is Ashley 1(true) is SBethe
        self.ParArray[self.NDFPAR + 38] = float(self.delayed_dispersion)
        self.ParArray[self.NDFPAR + 39] = float(self.Add_Doppler_Width)
        self.ParArray[self.NDFPAR + 40] = float(self.DebugMode) #controls debugging output
        return 0

    

   
    def my_updateProjectileEnergy(self):      
       
        if self.particle == "proton":
            self.ProjectileMass = cnst.Mp
        else:
            self.ProjectileMass = 1.0
        self.velocity_projectile(self.E0)
        self.myvelocitytext="p=" + f"{self.p0_rel:.1f}"+"a.u.   v/c ="+f"{self.beta_r:.3f}"
        
    def velocity_projectile(self, Energy):  # retuns velocity, also sets self.p0_rel, self.beta_rel and self.gamma_r
       E0_au= Energy*1000.0/cnst.HARTREE 
       self.gamma_r = 1.0 + E0_au/(self.ProjectileMass*(cnst.C)**2)
       T = E0_au * (1.0 + self.gamma_r) / (2.0 * self.gamma_r * self.gamma_r)  # Egerton appendix E, T in a.u.
       v = np.sqrt(2.0 * T/ self.ProjectileMass)   # V in a.u.
       self.p0_rel = self.gamma_r * self.ProjectileMass * v
       self.beta_r= v/cnst.C
       return v   
        
    # ====================================================
    def calculate_energy_axis(self,FineMeshFactor):
        step=self.Stepsize/FineMeshFactor
        NStep=self.NPoints*FineMeshFactor
        CenterFirstBin=self.LowerELimit + 0.5*step
        CenterLastBin=self.UpperELimit - 0.5*step
        if self.stepEfixed:
            self.x_axis= np.linspace(CenterFirstBin,CenterLastBin, NStep)  
            
        else:
            linearLastE=self.LowerELimit+NStep*step
            linearXaxis= np.linspace(CenterFirstBin,linearLastE-0.5*step, NStep) 
            missing_bit=self.UpperELimit - linearLastE
            logXaxis=np.geomspace( step/20,missing_bit ,NStep)
            last_step=logXaxis[self.MaxNPoints-1]-logXaxis[self.MaxNPoints-2]
            logXaxis=np.geomspace(self.Stepsize/20,missing_bit-0.5*last_step,NStep)
            self.x_axis=linearXaxis+logXaxis
        self.xstepsize=np.gradient(self.x_axis)
        
    # ====================================================    
    
    def calculate_eps_array(self):
        NPnts=np.size(self.x_axis)
        eps = np.zeros(NPnts,dtype=np.complex128)
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        epslib.Eps1Eps2(self.ParArray,self.x_axis, eps, self.q, self.DFChoice)
        return eps
        
    def eps1eps2(self):
        FineMeshFactor=1
        self.calculate_energy_axis(FineMeshFactor)
        eps =self.calculate_eps_array()
        self.Result1=np.real(eps)
        self.Result2=np.imag(eps)

    def oneovereps1eps2(self):
        FineMeshFactor=1
        self.calculate_energy_axis(FineMeshFactor)
        eps =self.calculate_eps_array()
        self.Result1=np.real(1.0/eps)
        self.Result2=-np.imag(1.0/eps)  # so the loss function not 1/eps2
        
    def calculate_eps_array_q(self):   
        eps_q = np.zeros(self.Nqstep,dtype=np.complex128)
        qstep= self.UpperqLimit / self.Nqstep 
        self.x_axis= np.linspace(0.5*qstep,self.UpperqLimit-0.5*qstep, self.Nqstep)    
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        epslib.Eps1Eps2_q(self.ParArray,self.x_axis, eps_q, self.Energy_qplot, self.DFChoice) 
        return eps_q

    def eps1eps2_q(self):
        eps_q = self.calculate_eps_array_q()
        self.Result1 = np.real(eps_q)
        self.Result2 = np.imag(eps_q)

    def oneovereps1eps2_q(self):
        eps_q= self.calculate_eps_array_q()
        self.Result1 = np.real(1.0/eps_q)
        self.Result2 = -np.imag(1.0/eps_q)  # so the loss function not 1/eps2
            
    # ========================
    
    def eps_kk_test(self):
        if not self.stepEfixed:
            print("ERROR,KK tests requires fixed stepsize")
            self.MyChapApp.UpdateStatus("ERROR,KK tests requires fixed stepsize")
            return -1
            
        FineMeshFactor=9
        length_fine=self.NPoints*FineMeshFactor
        self.calculate_energy_axis(FineMeshFactor)
        eps_fine=self.calculate_eps_array() 
        self.Result1_fine = eps_fine.real
        self.Result2_fine = eps_fine.imag
 
        self.Result3_fine = np.zeros(length_fine)
        self.Result4_fine = np.zeros(length_fine)
        eps_last=eps_fine[length_fine-1]
       
        epslib.Kramers_Kronig_eps1_from_eps2(self.x_axis[0], self.xstepsize[0], eps_last,length_fine,
             self.Result2_fine , self.Result3_fine)

        epslib.Kramers_Kronig_eps2_from_eps1(self.x_axis[0], self.xstepsize[0], eps_last,length_fine,
             self.Result1_fine , self.Result4_fine)

        #now cast the result on the normal grid we use for plotting  
        self.recast(FineMeshFactor, only_2=False)
        return 0
        

        
        
     # #modified from https://github.com/utf/kramers-kronig/blob/master/kkr.py, attempt to teach me  to use numpy   
    # def kk_numpy_eps1_from_eps2(self,de, eps_imag, cshift=1e-3):   

        # """Calculate the Kramers-Kronig transformation on imaginary part of dielectric

        # Doesn't correct for any artefacts resulting from finite window function.

        # Args:
            # de (float): Energy grid size at which the imaginary dielectric constant
                # is given. The grid is expected to be regularly spaced.
            # eps_imag (np.array): A numpy array with dimensions (n, 3, 3), containing
                # the imaginary part of the dielectric tensor.
            # cshift (float, optional): The implemented method includes a small
                # complex shift. A larger value causes a slight smoothing of the
                # dielectric function.

        # Returns:
            # A numpy array with dimensions (n, 3, 3) containing the real part of the
            # dielectric function.
        # """
        # eps_imag = np.array(eps_imag)
        # nedos = eps_imag.shape[0]
        # cshift = complex(0, cshift)
        # w_i = np.arange(0, nedos*de, de, dtype=np.complex_)
        # w_i += 0.5*de

        # def integration_element( w_r):  # this is an inner function
            # factor = w_i / (w_i**2 - w_r**2 + cshift)
            # total = np.sum(eps_imag * factor, axis=0)
            # return total * (2/math.pi) * de + 1.0

        # return np.real([integration_element(w_r) for w_r in w_i[:]])
        
        
     # #modified from https://github.com/utf/kramers-kronig/blob/master/kkr.py   
    # def kk_numpy_eps2_from_eps1(self,de, eps_real, cshift=1e-3):   

        # eps_real = np.array(eps_real)
        # nedos = eps_real.shape[0]
        # cshift = complex(0, cshift)
        # w_i = np.arange(0, nedos*de, de, dtype=np.complex_)
        # w_i += 0.5*de

        # def integration_element( w_r):  # this is an inner function
            # factor = w_r / (w_i**2 - w_r**2 + cshift)
            # total = - np.sum((eps_real - 1.0) * factor, axis=0)
            # return total * (2/math.pi) * de 

        # return np.real([integration_element(w_r) for w_r in w_i[:]])    
    
    def recast(self, factor,only_2):   
        x_fine=self.x_axis  # first make a copy of the old x axis 
        self.x_axis =  np.zeros(self.NPoints)
        self.Result1 = np.zeros(self.NPoints)
        self.Result2 = np.zeros(self.NPoints)   
        if not only_2:
            self.Result3 = np.zeros(self.NPoints)
            self.Result4 = np.zeros(self.NPoints) 
    
        first_i=math.ceil(factor/2.0)-1  
        for i in range(self.NPoints):
            self.x_axis[i] = x_fine[first_i+i*factor]
            self.Result1[i] = self.Result1_fine[first_i+i*factor]
            self.Result2[i] = self.Result2_fine[first_i+i*factor]
            if not only_2:
                self.Result3[i] = self.Result3_fine[first_i+i*factor]
                self.Result4[i] = self.Result4_fine[first_i+i*factor]
           
        del self.Result1_fine
        del self.Result2_fine
        if not only_2:
            del self.Result3_fine
            del self.Result4_fine      
                
    
    def one_over_eps_kk_test(self):
        if not self.stepEfixed:
            print("ERROR,KK tests requires fixed stepsize")
            self.MyChapApp.UpdateStatus("ERROR,KK tests requires fixed stepsize")
            return -1
        FineMeshFactor=9
        length_fine=self.NPoints*FineMeshFactor

        self.calculate_energy_axis(FineMeshFactor)
        eps_fine=self.calculate_eps_array() 
        self.Result1_fine =  np.real(1.0/eps_fine)
        self.Result2_fine =  np.imag(1.0/eps_fine)
        self.Result3_fine = np.zeros(length_fine)
        self.Result4_fine = np.zeros(length_fine)

        one_over_eps_last=1.0/eps_fine[length_fine-1]
        epslib.Kramers_Kronig_eps1_from_eps2(self.x_axis[0], self.xstepsize[0], one_over_eps_last,
            length_fine, self.Result2_fine, self.Result3_fine)
        epslib.Kramers_Kronig_eps2_from_eps1(self.x_axis[0], self.xstepsize[0], one_over_eps_last,
            length_fine, self.Result1_fine , self.Result4_fine)    
       
        self.Result4_fine  *= -1 # now result4 contains im (-1/eps) (obtained via KK)   
        self.Result2_fine  *= -1 # now result2 contains im (-1/eps) (original)                          
        self.recast(FineMeshFactor, only_2=False)  
        return 0    

    def n_and_k_from_eps1_eps2(self):
        # check this! not sure what it means for q!=0
        #see also Wooton eq. 3.25,3.26

        self.calculate_energy_axis(1)
        eps=self.calculate_eps_array() 
        self.Result1= np.real(eps**0.5)
        self.Result2= np.imag(eps**0.5)
            
    def n_and_k_kk_test(self):
        if not self.stepEfixed:
            print("ERROR,KK tests requires fixed stepsize")
            self.MyChapApp.UpdateStatus("ERROR,KK tests requires fixed stepsize")
            return -1
        FineMeshFactor=9
        length_fine=self.NPoints*FineMeshFactor
        self.calculate_energy_axis(FineMeshFactor)
        eps_fine=self.calculate_eps_array() 
        n_complex=np.sqrt(eps_fine)
        self.Result1_fine= np.real(n_complex)  # contains n
        self.Result2_fine= np.imag(n_complex)   # contains k
        self.Result3_fine = np.zeros(length_fine)
        self.Result4_fine = np.zeros(length_fine)
        n_last= n_complex[length_fine-1]

        epslib.Kramers_Kronig_eps1_from_eps2(self.x_axis[0], self.xstepsize[0],n_last,length_fine,
            self.Result2_fine, self.Result3_fine)
        epslib.Kramers_Kronig_eps2_from_eps1(self.x_axis[0], self.xstepsize[0],n_last,length_fine,
            self.Result1_fine, self.Result4_fine)      # result4:  kktransform from n, should be k
        self.recast(FineMeshFactor,only_2=False) 
        return 0 
      
    # ========================

    def calcDIIMFP(self):
        # calculate over the energy range defined in chapidif       
        self.Result1 = np.zeros(self.NPoints)
        self.calculate_energy_axis(1)
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
       # epslib.DIIMFP(self.ParArray,self.Result1,self.DFChoice,self.StoppingResultArray)
        epslib.DIIMFP_variable_step(self.ParArray,self.x_axis/cnst.HARTREE,self.xstepsize/cnst.HARTREE, self.Result1,self.DFChoice, self.StoppingResultArray)

    def sum_rules(self):
        FineMeshFactor=9
        self.calculate_energy_axis(FineMeshFactor)
        eps_fine =self.calculate_eps_array()
    
        length_fine = self.NPoints*FineMeshFactor
        self.Result1_fine = np.zeros(length_fine)
        self.Result2_fine = np.zeros(length_fine)
        self.Result3_fine = np.zeros(length_fine)
        self.Result4_fine = np.zeros(length_fine)
        
        n_complex = eps_fine**0.5
        eps2 = np.imag(eps_fine)
        Ext = np.imag(n_complex)
        LossFunction = np.imag(-1.0/eps_fine)

        E_au=self.x_axis/ cnst.HARTREE
        stepsize_au=self.xstepsize/cnst.HARTREE
        normalisation =  2.0 * self.UnitCellDensity * cnst.BOHR**3 * cnst.PI**2 
    
        self.Result1_fine = np.cumsum(eps2*E_au/normalisation*stepsize_au)
        self.Result2_fine = np.cumsum(LossFunction*E_au/normalisation*stepsize_au)
        self.Result3_fine = np.cumsum(2*Ext*E_au/normalisation*stepsize_au)
        self.Result4_fine = np.cumsum((2.0/cnst.PI)*(LossFunction/E_au)*stepsize_au)
   
        self.recast(FineMeshFactor,only_2=False)
        
    def inertial_rules(self):  # needs some moce checking!
        FineMeshFactor=29
        self.calculate_energy_axis(FineMeshFactor)
        eps_fine =self.calculate_eps_array()
        
        length_fine = self.NPoints*FineMeshFactor
        self.Result1_fine = np.zeros(length_fine)
        self.Result2_fine = np.zeros(length_fine)
        self.Result3_fine = np.zeros(length_fine)  # not realy used, limit is plotted in slot 3
        self.Result4_fine = np.zeros(length_fine)

        Re_n = np.real(eps_fine**0.5)
        ReEps=np.real(eps_fine)
        Re_one_overeps=np.real(1.0/eps_fine) 
        eps1=np.real(eps_fine)
        prefactor = -1.0/(2.0*cnst.PI**2)
       
        self.Result1_fine = np.cumsum((Re_one_overeps-1.0)*self.xstepsize)   
        self.Result2_fine = np.cumsum((Re_n-1.0)*self.xstepsize)   
        self.Result4_fine = np.cumsum(prefactor*(eps1-1.0) *self.xstepsize)   
        #epslib.cumulative_trapezoid(Re_one_overeps-1.0         ,self.x_axis, self.Result1_fine,self.LowerELimit)
        # epslib.cumulative_trapezoid(Re_n-1.0        ,self.x_axis, self.Result2_fine,self.LowerELimit)
        # epslib.cumulative_trapezoid(prefactor*(eps1-1.0)         ,self.x_axis, self.Result4_fine,self.LowerELimit)

       
        
        self.limitingvalue =0.0
        for i in range(4):
            limit = np.imag(eps_fine[i])*self.x_axis[i] /(4*cnst.PI)
            if limit > self.limitingvalue:
                self.limitingvalue=limit

        self.recast(FineMeshFactor,only_2=False)    
           


    def Mean_Excitation_Energy(self):
        # make sure we have the latest value,  and always evaluate at zero
        # momentum
        oldq = self.q
        self.q = 0.01  # better avoid 0
        FineMeshFactor=1
        last=self.NPoints - 1
        self.calculate_energy_axis(FineMeshFactor)
        eps =self.calculate_eps_array()
        self.q = oldq
        elf=np.imag(-1/eps)
       
        # integrate Im -1/eps in  result1 and result 2 weighted by omega and
        # omega log omega
        E_au=self.x_axis /  cnst.HARTREE 
        stepsize_au = self.xstepsize/ cnst.HARTREE
        sumtop=np.cumsum((2.0 / cnst.PI)*elf*np.log(E_au) * stepsize_au)
        sumbottom=np.cumsum((2.0 / cnst.PI)*elf* stepsize_au)
        self.Result1=cnst.HARTREE * np.exp(sumtop / sumbottom)  # in eV
        self.I0 = self.Result1[last]
        self.C0 = sumbottom[last] / math.exp(sumtop[last] / sumbottom[last])
        
        
        sumtop=np.cumsum((2.0 / cnst.PI) * elf * E_au * np.log(E_au) * stepsize_au)
        sumbottom=np.cumsum((2.0 / cnst.PI) * elf * E_au * stepsize_au)
        self.Result2=cnst.HARTREE * np.exp(sumtop / sumbottom)  # in eV
        self.MIE = self.Result2[last]
        self.C1 =  sumbottom[last]  / (math.exp(sumtop[last] / sumbottom[last])) ** 2 
       

    def calccurves(self, smallqonly):
        self.CurvesEnergy = np.zeros(self.NStopping)
        self.CurvesVelocity = np.zeros(self.NStopping)
        self.IMFPEnergy = np.zeros(self.NStopping)
        self.StoppingEnergy = np.zeros(self.NStopping)
        self.CrosssectionEnergy = np.zeros(self.NStopping)
        # assuming 'Bethe dispersion' (constant  at plasmon energy then free electron dispersion)
        self.BetheIMFPEnergy = np.zeros(self.NStopping)
        self.BetheStoppingEnergy = np.zeros(self.NStopping)
        self.BetheStoppingEnergy_Salvat = np.zeros(self.NStopping)
        self.BetheStragglingEnergy = np.zeros(self.NStopping) 
        self.L_0_Bethe=np.zeros(self.NStopping) 
        self.Straggling_Jackson=np.zeros(self.NStopping) 
        self.FCOR_Salvat=np.zeros(self.NStopping) 
        self.DL_IMFPaverage_Energy = np.zeros(self.NStopping)  # DL using single average oscillator
        self.DL_IMFP_sum_Energy = np.zeros(self.NStopping)

        self.DL_StoppingEnergy = np.zeros(self.NStopping)  # DL stopping based on average oscillator
        self.DL_Stopping_sum_Energy = np.zeros(self.NStopping)  # DL stopping based on sum oscillators
        self.DL_StragglingEnergy = np.zeros(self.NStopping)  # DL straggling based on average oscillator
        self.DL_Straggling_sum_Energy = np.zeros(self.NStopping)  # DL straggling based on sum oscillators
        self.DL_Stopping_from_ELF = np.zeros(self.NStopping)
        self.DL_IMFP_from_ELF = np.zeros(self.NStopping)
        self.DL_Straggling_from_ELF = np.zeros(self.NStopping)
        self.StragglingEnergy = np.zeros(self.NStopping)
        self.TPP_IMFPEnergy = np.zeros(self.NStopping)
      
        dummy= np.zeros(2) #dummy array for DIIMFP
        self.ParArray[self.NDFPAR + 5]=0.0  # NSTep so no  DIIMFP array is stored in dummy
        
       
        startime=time.time()
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        q=0.01
        self.ELF_wide_Energy = np.zeros(1300)  # should match with epslib
        ELF_wide = np.zeros(1300)  # should match with epslib
        ELF_wideResultArray =  np.zeros(1300)# for c1,c2,c3, i1, mie and i3
        epslib.Loss_wide(self.ParArray,  self.ELF_wide_Energy,  ELF_wide,  
            q, self.DFChoice, ELF_wideResultArray )
        self.C0 = ELF_wideResultArray[0]
        self.C1 = ELF_wideResultArray[1]
        self.C2 = ELF_wideResultArray[2]
        self.I0 = ELF_wideResultArray[3]
        self.MIE = ELF_wideResultArray[4]
        self.Istraggling = ELF_wideResultArray[5]
        self.sumBethe= ELF_wideResultArray[6]
        I1_au =self.MIE / cnst.HARTREE  #mean ionization energy used to separate close and distant collision (if smallqonly is true)

        if  self.particle == "proton":
            CurrentE = self.first_proton_energy
        else:
            CurrentE = self.first_electron_energy
 

        if smallqonly:
            self.ParArray[self.NDFPAR + 11] = math.sqrt(2.0 * I1_au)
        for Ecounter in range(self.NStopping):
            text="calculating step: {} Energy (keV):{:.2e}".format(Ecounter, CurrentE)
            self.MyChapApp.UpdateStatus(text)
            self.CurvesEnergy[Ecounter] = CurrentE
            self.ParArray[self.NDFPAR + 1] = CurrentE * 1000
            self.CurvesVelocity[Ecounter] = self.velocity_projectile(CurrentE)   # V in a.u. (using relativistic kinematics)
 
            epslib.DIIMFP(self.ParArray,dummy, self.DFChoice,self.StoppingResultArray)
                
            self.IMFPEnergy[Ecounter] = self.StoppingResultArray[0]
            self.CrosssectionEnergy[Ecounter] = 1.0 / (self.UnitCellDensity * self.StoppingResultArray[0])
            self.StoppingEnergy[Ecounter] = self.StoppingResultArray[1]
            self.StragglingEnergy[Ecounter] = self.StoppingResultArray[2]
            CurrentE = CurrentE * self.IncrFactor
            
        if self.RadiativeLosses ==1 and self.particle == "electron": 
            Z=round(self.Nelec_per_UC)
            for Ecounter in range(self.NStopping):
                TotalE_MeV = (self.CurvesEnergy[Ecounter] + 511.0)/1000.0
                Factor=1.0 +  Z * TotalE_MeV/800.0  # Bethe Heitler estimate, Nikjoo book Interaction Radiation Matter pg 109
                self.StoppingEnergy[Ecounter] = self.StoppingEnergy[Ecounter] * Factor
                
      
        
        finishtime=time.time()
        text= "duration calculation {:.1f} sec".format(finishtime - startime)
        self.MyChapApp.UpdateStatus(text)
        self.Stopping_Linear_V()
        self.Bethe_L_Salvat()
        self.calculate_Straggling_Jackson()
        if self.particle == 0:
            self.TanumaPowellPenn()
       
        if self.Approximations:
            self.calculate_approximations()
        
                       

    def calculate_approximations(self):
         
        n0 = ( self.I0 * self.I0 / (4 * cnst.PI * cnst.HARTREE**2))  # electron density when interpreting I0 as a plamon energy (in e^-per au^3)
        n1 = self.MIE * self.MIE / (4 * cnst.PI * cnst.HARTREE**2)
        n2 = self.Istraggling * self.Istraggling / (4 * cnst.PI * cnst.HARTREE**2)
        self.sumBethe = n1 * self.C1 
        self.C0_from_BetheSum = ( self.sumBethe / n0)  # fraction of space filled with this oscillator so number of electrons correct.
       
        self.C1_from_BetheSum = self.sumBethe / n1
        
        self.C2_from_BetheSum = self.sumBethe / n2
        
        n1 = self.MIE * self.MIE / (4 * cnst.PI * cnst.HARTREE**2)
        self.C1_from_BetheSum = self.sumBethe / n1
 

        # if self.Projectile == 1:
            # mass = cnst.Mp
        # else:
            # mass = 1
        I1_au = self.MIE / cnst.HARTREE
        I0_au = self.I0 / cnst.HARTREE
        I2_au = self.Istraggling / cnst.HARTREE
        
        top0_zerowidth = 0.0
        bottom0_zerowidth = 0.0
        top1_zerowidth = 0.0
        bottom1_zerowidth = 0.0
        top2_zerowidth = 0.0
        bottom2_zerowidth = 0.0
        self.OscillatorsPresent = False
        if self.DFmodel  == "DL" or  self.DFmodel == 'Mermin':  
            for i in range(self.maxOscillators):
               
                Amp_i = self.Amps[i]
                if Amp_i > 0:
                    self.OscillatorsPresent = True
                    wi_p_au = self.Omegas[i] / cnst.HARTREE
                    top0_zerowidth += Amp_i * wi_p_au * math.log(wi_p_au)
                    bottom0_zerowidth += Amp_i * wi_p_au
                    top1_zerowidth += Amp_i * wi_p_au**2 * math.log(wi_p_au)
                    bottom1_zerowidth += Amp_i * wi_p_au**2
                    top2_zerowidth += Amp_i * wi_p_au**3 * math.log(wi_p_au)
                    bottom2_zerowidth += Amp_i * wi_p_au**3
            if self.OscillatorsPresent:
                I0_zerowidth = math.exp(top0_zerowidth / bottom0_zerowidth)
                C0_zerowidth = bottom0_zerowidth / I0_zerowidth
                I1_zerowidth = math.exp(top1_zerowidth / bottom1_zerowidth)
                C1_zerowidth = bottom1_zerowidth / (I1_zerowidth**2)
                I2_zerowidth = math.exp(top2_zerowidth / bottom2_zerowidth)
                C2_zerowidth = bottom2_zerowidth / (I2_zerowidth**3)
                print("for zero width oscillators:")
                print("C0=", C0_zerowidth, "I0=", I0_zerowidth * cnst.HARTREE, "eV")
                print("C1=", C1_zerowidth, "I1=", I1_zerowidth * cnst.HARTREE, "eV")
                print("C2=", C2_zerowidth, "I2=", I2_zerowidth * cnst.HARTREE, "eV")
                print("end calculate approximations")

        if  self.particle == "proton":
            ConstantA = 1.0/2.0  # notation so we are in line with imfp draft eq. 11
            ConstantB = 1
        else:
            ConstantA = 1
            ConstantB = 2
        for Ecounter in range(self.NStopping):
            velocity = self.CurvesVelocity[Ecounter]
            prefactor = self.C1 * I1_au**2 / velocity**2


    

            # start code for stopping  and IMFP  straggling in high energy limit for the case of
            # 'Bethe disperion'  using the MIE as energy    , (NOT I0 for IMFP, I1 for straggling)
            q_max = velocity / ConstantA #non-relativistic case have to think about this
            
            q_min = I1_au / velocity
            if q_min > q_max:
                q_min = q_max  # this will make logBethe equal 0
                
           
            q_c = np.sqrt(2 * I1_au)
          
            logBethe = np.log(q_max / q_min)
            self.BetheStoppingEnergy[Ecounter] = prefactor * logBethe * cnst.HARTREE / cnst.BOHR
    
            if logBethe > 0.0:
                tmp = math.log(q_c / q_min)
                tmp += 0.5 - I1_au / velocity**2

                oneoverlambda = self.C1 * I1_au / velocity**2 * tmp
                if oneoverlambda > 0.0:
                    self.BetheIMFPEnergy[Ecounter] = 1.0 / oneoverlambda  # in a.u.
                else:
                    self.BetheIMFPEnergy[Ecounter] = 1e30
                self.BetheIMFPEnergy[Ecounter] *= cnst.BOHR
            else:
                self.BetheIMFPEnergy[Ecounter] = 1e30
            if (q_c > q_min) and (q_max > q_c):
                logstraggling = np.log(q_c / q_min)
                self.BetheStragglingEnergy[Ecounter] = (
                    self.C1
                    * I1_au**2
                    / velocity**2
                    * (I1_au * logstraggling + (q_max * q_max - q_c * q_c) / 4.0)
                )
                self.BetheStragglingEnergy[Ecounter] *= (
                    cnst.HARTREE * cnst.HARTREE / cnst.BOHR
                )
            else:
                self.BetheStragglingEnergy[Ecounter] = 0
                
           
            # now IMFP and stopping and straggling based on more than one DL oscillator (ie  with gamma 0)
            oneoverlambda_total = 0.0
            stopping = 0.0
            straggling = 0.0
            for i in range(self.maxOscillators):
                Amp_i = self.Amps[i]
                if Amp_i != 0.0:
                    wi_p_au = float(self.Omegas[i]) / cnst.HARTREE
                    if velocity**2 > 2 * ConstantB * wi_p_au:
                        q_min_DL = (
                            velocity - np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                        ) / ConstantB
                        q_max_DL = (
                            velocity + np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                        ) / ConstantB
                    else:
                        q_min_DL = 1
                        q_max_DL = 1

                    logDL = math.log(q_max_DL / q_min_DL)
                    stopping += Amp_i * wi_p_au**2 / velocity**2 * logDL

                    tmp = logDL - 0.5 * math.log(
                        (2 * wi_p_au + q_max_DL**2) / (2 * wi_p_au + q_min_DL**2))
                    oneoverlambda_total += Amp_i * wi_p_au / velocity**2 * tmp

                    tmp = Amp_i * (wi_p_au**3 / velocity**2 * logDL)
                    fudge = (  1.0 ) # fudge here if equivalent oscillator has the wrong density

    
                    straggling += (fudge * Amp_i * 0.25 * wi_p_au**2 / velocity**2
                        * (q_max_DL**2 - q_min_DL**2) + tmp )

            self.DL_Stopping_sum_Energy[Ecounter] = stopping * cnst.HARTREE / cnst.BOHR
            if oneoverlambda_total > 0:
                self.DL_IMFP_sum_Energy[Ecounter] = (
                    1.0 / oneoverlambda_total
                ) * cnst.BOHR
            else:
                self.DL_IMFP_sum_Energy[Ecounter] = 1e20
            self.DL_Straggling_sum_Energy[Ecounter] = (
                straggling * cnst.HARTREE * cnst.HARTREE / cnst.BOHR
            )

            # start code for DL model single oscilator high-energy limit, first stopping
            if velocity**2 > 2 * ConstantB * I1_au:
                q_min_DL = (
                    velocity - np.sqrt(velocity**2 - 2 * ConstantB * I1_au)
                ) / ConstantB
                q_max_DL = (
                    velocity + np.sqrt(velocity**2 - 2 * ConstantB * I1_au)
                ) / ConstantB
            else:
                q_min_DL = 1
                q_max_DL = 1
            logDL = np.log(q_max_DL / q_min_DL)
            prefactor = I1_au**2 / velocity**2
            self.DL_StoppingEnergy[Ecounter] = (
                self.C1 * prefactor * logDL * cnst.HARTREE / cnst.BOHR
            )

            # now start code for IMFP from equivalent oscillator in high E limit
            if velocity**2 > 2 * ConstantB * I0_au:
                q_min_DL_IMFP = (
                    velocity - np.sqrt(velocity**2 - 2 * ConstantB * I0_au)
                ) / ConstantB
                q_max_DL_IMFP = (
                    velocity + np.sqrt(velocity**2 - 2 * ConstantB * I0_au)
                ) / ConstantB
            else:
                q_min_DL_IMFP = 1
                q_max_DL_IMFP = 1

            logDL_IMFP = np.log(q_max_DL_IMFP / q_min_DL_IMFP)

            tmp = logDL_IMFP - 0.5 * math.log(
                (2 * I0_au + q_max_DL_IMFP**2) / (2 * I0_au + q_min_DL_IMFP**2)
            )
            oneoverlambda = self.C0 * I0_au * tmp / velocity**2
            if oneoverlambda > 0.0:
                self.DL_IMFPaverage_Energy[Ecounter] = 1.0 / oneoverlambda  # in a.u.
            else:
                self.DL_IMFPaverage_Energy[Ecounter] = 1e20  # in a.u.
            self.DL_IMFPaverage_Energy[Ecounter] *= cnst.BOHR

            # now start code for straggling from equivalent oscillator in high E limit
            # this does not work for oscillators as then I2_au diverges, then use sum zero width oscillators
            if self.OscillatorsPresent:
                I2_au = I2_zerowidth
                self.C2 = C2_zerowidth
                self.Istraggling = I2_au * cnst.HARTREE

            if velocity**2 > 2 * ConstantB * I2_au:
                q_min_DL_Straggling = (
                    velocity - np.sqrt(velocity**2 - 2 * ConstantB * I2_au)
                ) / ConstantB
                q_max_DL_Straggling = (
                    velocity + np.sqrt(velocity**2 - 2 * ConstantB * I2_au)
                ) / ConstantB
            else:
                q_min_DL_Straggling = 1
                q_max_DL_Straggling = 1
            logDL_Straggling = np.log(q_max_DL_Straggling / q_min_DL_Straggling)
            self.DL_StragglingEnergy[Ecounter] = self.C2 * (
                I2_au**3 / velocity**2 * logDL_Straggling
            )
            # indeed I1_au^2 and C1 not I2_au^2 and C2 in last part, real average density not of average oscillator

            self.DL_StragglingEnergy[Ecounter] += (self.C1 * 0.25 * I1_au**2
                / velocity**2 * (q_max_DL_Straggling**2 - q_min_DL_Straggling**2))

            self.DL_StragglingEnergy[Ecounter] *= (cnst.HARTREE * cnst.HARTREE / cnst.BOHR)
        self.DL_stopping_IMFP_straggling_from_ELF()

    def DL_stopping_IMFP_straggling_from_ELF(self):
        #should work for all models
     
        if self.particle =="proton":
            ConstantB = 1
        else:
            ConstantB = 2
      
        
        for Ecounter in range(self.NStopping):
            velocity = self.CurvesVelocity[Ecounter]

            self.DL_Stopping_from_ELF[Ecounter] = 0.0
            one_over_lambda = 1e-20  # so we never divide by 0
            self.DL_Straggling_from_ELF[Ecounter] = 0.0
            for i_w in range(1300):  #should match number of steps in epslib.Loss_wide
                wi_p_au = self.ELF_wide_Energy[i_w] / cnst.HARTREE
                if i_w ==0:
                    Stepsize_au =  self.ELF_wide_Energy[0] / cnst.HARTREE
                else:    
                    Stepsize_au = (self.ELF_wide_Energy[i_w] - self.ELF_wide_Energy[i_w-1]) / cnst.HARTREE
                C_w = (2.0 / (cnst.PI * wi_p_au) * self.ELF_wide[i_w])  # G(omega) in Penn's paper for vanishing small gamma
                elec_dens = C_w * wi_p_au * wi_p_au / (4.0 * cnst.PI)
                prefactori = 4 * np.pi / velocity**2 * elec_dens
                if velocity**2 > 2 * ConstantB * wi_p_au:
                    q_min_DL = (
                        velocity - np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                    ) / ConstantB
                    q_max_DL = (
                        velocity + np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                    ) / ConstantB
                else:
                    q_min_DL = 1
                    q_max_DL = 1

                logDL = math.log(q_max_DL / q_min_DL)
                self.DL_Stopping_from_ELF[Ecounter] += (
                    prefactori * logDL * cnst.HARTREE / cnst.BOHR * Stepsize_au
                )
                tmp = logDL - 0.5 * math.log(
                    (2 * wi_p_au + q_max_DL**2) / (2 * wi_p_au + q_min_DL**2)
                )
                one_over_lambda += prefactori / wi_p_au * tmp * Stepsize_au
                self.DL_Straggling_from_ELF[Ecounter] += (
                    prefactori
                    * (wi_p_au * logDL + 0.25 * (q_max_DL**2 - q_min_DL**2))
                    * cnst.HARTREE**2
                    / cnst.BOHR
                    * Stepsize_au
                )

            self.DL_IMFP_from_ELF[Ecounter] = (1.0 / one_over_lambda) * cnst.BOHR
            
    def TanumaPowellPenn(self):
        
        bandgap = 0.0 #bandgap not implementd
        rho = self.specificweight
        if self.w_p_TPP > 0.0:
            beta_TPP_2m = (-0.1 + 0.944 / math.sqrt(self.w_p_TPP**2 + bandgap * bandgap)
                + 0.069 * rho**0.1)  # assume no gap
        # from  TPP SIA 43 689 2011
        gamma_TPP_2m = 0.191 * rho**-0.5
        U = self.w_p_TPP**2 / 829.4
        C_TPP_2m = 1.97 - 0.91 * U
        D_TPP_2m = 53.4 - 20.8 * U
        print("calculated Beta", beta_TPP_2m, "w_p^2*beta",   beta_TPP_2m * self.w_p_TPP**2,
            "gamma_TPP_2m",  gamma_TPP_2m)
                # start code tor TPP IMFP 
       
        for Ecounter in range(self.NStopping):
            if (self.particle=="electron") and (self.w_p_TPP > 0.0):  # electrons
                E_eV = self.CurvesEnergy[Ecounter] * 1000
                tmp = ( beta_TPP_2m * math.log(gamma_TPP_2m * E_eV) - C_TPP_2m / E_eV
                    + D_TPP_2m / E_eV**2)
                self.TPP_IMFPEnergy[Ecounter] = E_eV / (tmp * self.w_p_TPP**2)
            else:
                self.TPP_IMFPEnergy[Ecounter] = 0.0    
                
                
    def calculate_Straggling_Jackson(self):            
     # Jackson Classical electrodynamics eq.13.50, Salvat PRA 2022 eq. 112, only for protons  
        for Ecounter in range(self.NStopping):
            CurrentE= self.CurvesEnergy[Ecounter] # in keV
            v=self.velocity_projectile(CurrentE) #subroutine also sets self.gamma_r and self.beta_r
            gamma2=self.gamma_r*self.gamma_r
            beta2=self.beta_r*self.beta_r
            curlyB=2 * cnst.PI/v**2   
          #  straggling_per_UC= curlyB * 2 * self.Nelec_per_UC *(gamma2*beta2*cnst.C**2*(1.0-beta2/2.0))
            straggling_per_UC= curlyB * 2 * self.Nelec_per_UC *(gamma2*v**2*(1.0-beta2/2.0))
            straggling_per_UC *=  cnst.HARTREE**2 * cnst.BOHR**2# now in eV^2/Angstrom^2 per UC
            self.Straggling_Jackson[Ecounter] = straggling_per_UC*self.UnitCellDensity # now eV^2/angstrom
            
    def Stopping_Linear_V(self):
         # a linear approximation
        tmp=self.StoppingEnergy/self.CurvesVelocity   
        Imaxsloop=np.argmax(tmp)
        self.maxsloop=self.StoppingEnergy[Imaxsloop]/self.CurvesVelocity[Imaxsloop]   
        self.LinearApprox_lowE=np.zeros(Imaxsloop+2)
        self.x_LinearApprox_lowE=np.zeros(Imaxsloop+2)# filled in run_and_plot with either E or v
        for Ecounter in range(Imaxsloop+1):
            self.LinearApprox_lowE[Ecounter+1]=self.CurvesVelocity[Ecounter]* self.maxsloop
        #self.LinearApprox_lowE[0] =0.0   
                
    def get_stopping_fine(self):  # also do the same for straggling
        self.range_integration_factor=10
        self.calccurves(False) 
        self.stopping_fine = np.zeros(self.range_integration_factor*self.NStopping)
        self.straggling_fine = np.zeros(self.range_integration_factor*self.NStopping)
        self.energy_fine = np.zeros(self.range_integration_factor*self.NStopping)
        self.velocity_fine = np.zeros(self.range_integration_factor*self.NStopping)
        
        stopping_before=0.0  
        straggling_before=0.0 
        energy_before=0.0 
        counter=0    
        for i in range(self.NStopping):
            energy_after= self.CurvesEnergy[i]
            stopping_after=  self.StoppingEnergy[i]
            straggling_after = self.StragglingEnergy[i]
            Estep= (energy_after-energy_before)/self.range_integration_factor
            stopping_step= (stopping_after-stopping_before)/self.range_integration_factor
            for j in range(self.range_integration_factor):
                self.energy_fine[counter]= energy_before+Estep*(j+1)
                self.velocity_fine[counter] = self.velocity_projectile(self.energy_fine[counter])
                self.stopping_fine[counter]=stopping_before+stopping_step*(j+1)
                self.straggling_fine[counter]=straggling_before+straggling_after*(j+1)
                counter+=1
            energy_before=energy_after
            stopping_before=stopping_after
            straggling_before=straggling_after
                     
    def projectile_range(self):
        self.get_stopping_fine()
        csda_range=0.0  
        print("self.MyChapApp.runplot.x_axis_keV",self.MyChapApp.runplot.x_axis_keV)
        Energy_before=0.0 
        if self.MyChapApp.runplot.x_axis_keV: 
            self.x_axis = self.energy_fine 
        else:
            print("velocity")
            self.x_axis = self.velocity_fine 
        self.Result1 = np.zeros(self.range_integration_factor*self.NStopping)
        for i in range(self.NStopping*self.range_integration_factor):
            csda_range += (self.energy_fine[i]-Energy_before)*1000/  self.stopping_fine[i] 
            self.Result1[i]=csda_range
            Energy_before=self.energy_fine[i] 
    
    def gaussian(self,x, mu, sig):  # https://stackoverflow.com/questions/14873203/plotting-of-1-dimensional-gaussian-distribution-functio
        return (1.0 / (np.sqrt(2.0 * np.pi) * sig)) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)    
            
    def Energy_Depth_Dist(self):
        NS=20000
        # first calculate the optimum energy increment factoor based on E0 and self.first_electron_energy
        ratio_energies_involved = self.E0/self.first_electron_energy
        best_incr_factor=math.exp(math.log(ratio_energies_involved)/(self.NStopping-1))
        print("best_incr_factor",best_incr_factor)
        oldfactor=self.IncrFactor 
        self.IncrFactor = best_incr_factor*1.0001
        self.get_stopping_fine()
        self.IncrFactor = oldfactor
        csda_range=0.0 
        
        Energy_before=0.0  
        for i in range(self.NStopping*self.range_integration_factor):
            csda_range += (self.energy_fine[i]-Energy_before)*1000/  self.stopping_fine[i] 
           
            Energy_before=self.energy_fine[i] 
            if self.energy_fine[i] > self.E0:
                range_this_energy= csda_range
                break
        stepsize=range_this_energy/(0.9*NS)
        self.x_axis =np.arange(0,range_this_energy/0.9, stepsize)
  
        depth_distribution_this_energy=np.zeros(NS)
        self.Result1 = (np.zeros(NS))
        self.Result2 = (np.zeros(NS))

        current_Energy=self.E0
        straggling_integrated=0
        for Idepth in range(NS): 
            for i in range(self.NStopping*self.range_integration_factor):
                if current_Energy < 0:
                        break
                if current_Energy < self.energy_fine[i]:
                    #self.Result1[Idepth]=self.stopping_fine[i]  # this is without straggling
                    current_Energy-= self.stopping_fine[i]*stepsize/1000  # E0 in keV, stopping eV/angstrom
                    straggling_integrated +=self.straggling_fine[i]*stepsize # in eV^2
                    width_beam=math.sqrt(straggling_integrated)/self.stopping_fine[i]
                    current_depth=self.x_axis[Idepth]
                    depth_distribution_this_energy=stepsize *self.gaussian(self.x_axis,current_depth,width_beam)
                    discreet_sum=depth_distribution_this_energy.sum() # this should be 1 but this does not work is width beam not much larger than stepsize              
                    self.Result2 +=  self.stopping_fine[i]*depth_distribution_this_energy/discreet_sum
                    self.Result1[Idepth]=self.stopping_fine[i]  # this is without straggling
                    break
                 
             
    def Bethe_L_Salvat(self):   
		#subroutine aiming to reproduce high energy limit sBethe from Salvat
     
        I_au=self.MIE / cnst.HARTREE
       
        for Ecounter in range(self.NStopping):
            prefactor =  4 * np.pi / self.CurvesVelocity[Ecounter] ** 2 * self.sumBethe * cnst.HARTREE / cnst.BOHR
            CurrentE= self.CurvesEnergy[Ecounter] # in keV
            v=self.velocity_projectile(CurrentE) #subroutine also sets self.gamma_r and self.beta_r
            gamma=self.gamma_r
            gamma2=gamma*gamma
            beta2=self.beta_r*self.beta_r
            if  self.particle == "proton":
                R=1.0/(1.0+(1.0/cnst.Mp**2)+2.0*gamma/cnst.Mp)
                FCOR=math.log(R)+((gamma2-1.0)*R/(gamma*cnst.Mp))**2
                
            else:   
                FCOR=(2.0*gamma2-1.0)/gamma2+0.125*((gamma-1.0)/gamma)**2-(4.0-((gamma-1.0)/gamma)**2)*math.log(2.0)-math.log(gamma+1.0)
            
            self.L_0_Bethe[Ecounter]=math.log(2*v*v/ I_au) +  math.log(gamma2)-beta2
       
            if self.L_0_Bethe[Ecounter] < 0.0:
                self.L_0_Bethe[Ecounter] = 0.0
            self.FCOR_Salvat[Ecounter]= FCOR
            self.BetheStoppingEnergy_Salvat[Ecounter] = prefactor * (self.L_0_Bethe[Ecounter] + 0.5*FCOR)
            if self.BetheStoppingEnergy_Salvat[Ecounter] < 0.0:
                self.BetheStoppingEnergy_Salvat[Ecounter] = 0.0
            
           
               

    def shell_effect(self, smallqonly):
        self.x_axis =  np.zeros(self.NStopping)  
        self.Result1 = np.zeros(self.NStopping)  
        self.Result2 = np.zeros(self.NStopping)  
        self.Result3 = np.zeros(self.NStopping)  
        self.Result4 = np.zeros(self.NStopping)  
        self.calccurves(smallqonly)
        for Ecounter in range(self.NStopping):  
            prefactor =  4 * np.pi / self.CurvesVelocity[Ecounter] ** 2 * self.sumBethe * cnst.HARTREE / cnst.BOHR
            L_df = self.StoppingEnergy[Ecounter] / prefactor # recover L
            self.Result1[Ecounter] = L_df    # already restricted value if smallqonly == true
            if smallqonly:
                self.Result2[Ecounter] = 0.5 * (self.L_0_Bethe[Ecounter]+  0.5*self.FCOR_Salvat[Ecounter] )
            else:
                self.Result2[Ecounter] = self.L_0_Bethe[Ecounter]+  0.5*self.FCOR_Salvat[Ecounter] 
            if(self.Result2[Ecounter] < 0.0):
                self.Result2[Ecounter] = 0.0
                
            
            self.Result3[Ecounter] =  self.Result1[Ecounter] -self.Result2[Ecounter]   # "shell effect"
            self.Result4[Ecounter] = 0.5 * self.FCOR_Salvat[Ecounter]

    def CalcCompton(self):
        # actual momentum array but plotting routine expect x-axis in this array
        
         # k_limit is the largest momentum component of target electrons considered
        self.q_old=self.q
        self.q=self.q_Compton 
        if( self.Compton_k_limit > 0.5 * self.q_Compton):
            Compton_k_limit_used = 0.5 * self.q_Compton   # so elower stays >=0
        else:
            Compton_k_limit_used  = self.Compton_k_limit
        if self.Dispersion_relativistic:
            Q_recoil_au =self.calc_Qrecoil(self.q_Compton)
        else:
            Q_recoil_au=self.q_Compton**2 / 2.0
                
        E_center = Q_recoil_au * cnst.HARTREE
        width_per_au =math.sqrt(2*Q_recoil_au)
        E_lower = E_center - Compton_k_limit_used * width_per_au * cnst.HARTREE
        E_upper = E_center + Compton_k_limit_used * width_per_au * cnst.HARTREE
 
        self.x_axis= np.linspace(E_lower,E_upper,self.NPoints )  
  
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        eps=self.calculate_eps_array()
        
        self.Result1= np.imag(eps)
        self.Result1 *=  Q_recoil_au / (2*cnst.PI**2)* width_per_au   # momentum density per  a.u^3.
        self.Result1 /= (cnst.BOHR**3 * self.UnitCellDensity) # now  momentum density per atom 
        
        self.x_axis= np.linspace(-Compton_k_limit_used,Compton_k_limit_used,self.NPoints )  #x-axis in momentum a.u.
        Mom_step =  self.x_axis[1]- self.x_axis[0]
        area=np.sum(self.Result1)*Mom_step
        print("electrons per unit cell",area)
        self.q=self.q_old 

    def Fresnel_at_E(self):  # thus as a function of angle
      
        stepsize=(90.0/self.NPoints)  # in degrees
        self.x_axis= np.linspace(0.5*stepsize,90-0.5*stepsize, self.NPoints)  
        θi=self.x_axis*math.pi/180.0  # array in rad
      
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        eps= epslib.single_Eps1Eps2(self.ParArray, 0.0, self.E_Fresnel,  self.DFChoice)
     
        
#from: https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/calc/reflection.py
        n1=complex(1.0,0.0)  # number: complex index of refr. vacuum
        n2=eps**0.5    # number  complex index  of refr.  medium
        θt = np.asin(n1/n2*np.sin(θi))  #array,   θt is complex!!
        rs = (n1*np.cos(θi)-n2*np.cos(θt)) / (n1*np.cos(θi)+n2*np.cos(θt))
        rp = (n2*np.cos(θi)-n1*np.cos(θt)) / (n1*np.cos(θt)+n2*np.cos(θi))
        self.Result1 = abs(rs)
        self.Result2 = abs(rp)
        self.Result3 = np.rad2deg(np.atan(abs(rp)/abs(rs)) ) #psi (in radians
        self.Result4 = np.rad2deg(np.angle(rp/rs) )   #Delta (in radians)
        
        # for i in range(self.NPoints):
                
                   
                    # θt = cmath.asin(n1/n2*math.sin(θi))  #   θt is complex!!
                    # rs = (n1*cmath.cos(θi)-n2*cmath.cos(θt)) / (n1*cmath.cos(θi)+n2*cmath.cos(θt))
                    # rp = (n2*cmath.cos(θi)-n1*cmath.cos(θt)) / (n1*cmath.cos(θt)+n2*cmath.cos(θi))
                    # self.Result1[i]=abs(rs)
                    # self.Result2[i]=abs(rp)
                    # # self.Result1[i]=cmath.phase(rp/rs)
                    # # self.Result2[i]=math.atan(abs(rp)/abs(rs))
                    
#from: https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/calc/reflection.py

    def Fresnel_at_angle(self): # thus as a function of energy
        #from: https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/calc/reflection.py
        print("in delta psi")
        self.calculate_energy_axis(1)
        epscomplex =self.calculate_eps_array()
      
        n1=complex(1.0,0.0)  # complex index of refr. vacuum
        n2=epscomplex**0.5   #complex array n2 as function omega
       
        θi= np.deg2rad(self.phi_ellipsometry)# θi: single (real) value
        θt = np.arcsin(n1/n2*np.sin(θi))  # θt: complex array
        rs = (n1*np.cos(θi)-n2*np.cos(θt)) / (n1*np.cos(θi)+n2*np.cos(θt))
        rp = (n2*np.cos(θi)-n1*np.cos(θt)) / (n1*np.cos(θt)+n2*np.cos(θi))

        self.Result1 = np.abs(rs)#**2
        self.Result2 = np.abs(rp)#**2
        self.Result3 = np.rad2deg(np.atan(abs(rp)/abs(rs)) ) #psi (in radians
        self.Result4 = np.rad2deg(np.angle(rp/rs) )   #Delta (in radians)
        # for i in range(self.NPoints):
            # N2=epscomplex**0.5 
            # N=N2.real
            # θt = cmath.asin(N1/N*math.sin(θi))
            # rs = (N1*cmath.cos(θi)-N2*cmath.cos(θt)) / (N1*cmath.cos(θi)+N2*cmath.cos(θt))
            # rp = (N2*cmath.cos(θi)-N1*cmath.cos(θt)) / (N1*cmath.cos(θt)+N2*cmath.cos(θi))
            # self.Result1[i]=cmath.phase(rp)-cmath.phase(rs)
            # self.Result2[i]=math.atan(abs(rp)/abs(rs))
            
            
                    
    def xray_absorption(self):
       

        q_old=self.q
        self.q=0.0    
        self.n_and_k_from_eps1_eps2()
        self.q=q_old
        for i in range(self.NPoints):
            lmbda=12398/self.x_axis[i]  #lambda in angstronm
            self.Result1[i]= 4*math.pi*self.Result2[i]/lmbda  # in angstrom
            self.Result1[i]=self.Result1[i]*1e8/self.specificweight 
            
        
                                                    
            
    
    def colorplot_lossfunction(self):
        eps = np.zeros(self.NPoints,dtype=np.complex128)
        current_q = 0.5 * self.Stepsize_qplot  # initial value
        if self.NPoints > 0 and self.Nqstep > 0:
            self.my_image = np.zeros((self.NPoints, self.Nqstep))
        else:
            print("dimenson error self.NPoints", self.NPoints,"self.Nqstep", self.Nqstep)
            return
       
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        self.x_axis= np.linspace(self.CenterFirstBin,self.CenterLastBin, self.NPoints)    

        for i in range(self.Nqstep):
            epslib.Eps1Eps2(self.ParArray,self.x_axis,eps, current_q, self.DFChoice)
            if self.bulk_eq:
                oneovereps=1.0/eps
                LossFunction=-np.imag(oneovereps)
                self.my_image[:, i] = LossFunction
            else:   
                SurfLossFunction=np.imag((eps - 1.0)*(eps - 1.0) / (eps*(eps + 1.0)))
                self.my_image[:, i] = SurfLossFunction
            current_q = current_q + self.Stepsize_qplot
            
    def Calculate_Integration_limits(self):   
        MC2 = self.ProjectileMass* cnst.C*cnst.C
        for x in range(30):
            q = x * self.UpperqLimit / 29.0
            self.xArray[x] = q
            # E0 is in keV, nonreativistic calculation
            p1_min = math.sqrt(2 * self.ProjectileMass * self.E0 * 1000.0 / cnst.HARTREE) - q
 
            omega = self.E0 * 1000.0 - p1_min * p1_min * cnst.HARTREE / (2 * self.ProjectileMass)
            self.yArray[x] = omega
            # and now we try to do the same using relativistic kinematics for the same energy losses
            p1_min_rel = self.p0_rel - q
            # use E^2=(pc)^2+(M_0 C^2), fingers crossed enough precision
            E_square = (p1_min_rel*cnst.C)**2 + (MC2)**2
            E_after = np.sqrt(E_square)
            E_kin_after=E_after - MC2
            omega = self.E0 * 1000.0 - E_kin_after * cnst.HARTREE
            self.yArray_relativistic[x] = omega
 
    def scale_image(self): 
        print("scaling") 
        arraymax = np.amax(self.my_image)
        if self.max_eq > 0.0:
            arraymax=self.max_eq
     
        minimum_value=np.amin(self.my_image)
        if self.LogXY:
            minimum_value = arraymax / 10**self.log_range
            print("new minimum_value",minimum_value, "maximum",arraymax)
        self.my_scaled_image = self.my_image.clip( minimum_value,arraymax)   
        if self.LogXY:
           self.my_scaled_image=np.log10(self.my_scaled_image)
          

        

   

    def colorplot_ddcs(self):
        self.x_axis = np.zeros(self.NThetaStep)
        self.Result1 = np.zeros(self.NThetaStep)
        self.Result2 = np.zeros(self.NThetaStep)
        
        ddcs_array =np.zeros(self.NPoints * self.NThetaStep)
        self.my_image = np.zeros((self.NPoints,self.NThetaStep))
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
      
        self.stopping_longitudinal=epslib.stopping_longitudinal(self.ParArray, ddcs_array, self.DFChoice )
        print("stopping from integration over this angular-energy loss range (eV/Å, no retard.):",self.stopping_longitudinal)
        for i in range(self.NPoints):
            for j in range(self.NThetaStep):
                self.my_image[i, j] =ddcs_array[j * self.NPoints + i]
        
    def colorplot_Cerenkov(self):
        ddcs_array =np.zeros(self.NPoints * self.NThetaStep)
        self.my_image = np.zeros((self.NPoints, self.NThetaStep))
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        self.stopping_total=epslib.stopping_total(self.ParArray, ddcs_array,self.DFChoice)  # calculates the 2d ddcs array, integrates over 2pi to get stopping and this is the return value
        print("stopping from integration over this angular-energy loss range (eV/Å):",self.stopping_total)

        for i in range(self.NPoints):
            for j in range(self.NThetaStep):
                self.my_image[i, j] = ddcs_array[j * self.NPoints + i]
                
                
    def difference_due_to_Cerenkov(self):  
        self.colorplot_Cerenkov()    
        cerenkov_image= self.my_image
        self.colorplot_ddcs()
        self.my_image=cerenkov_image - self.my_image
        self.stopping_due_to_photons=self.stopping_total - self.stopping_longitudinal
        print(" stopping total", self.stopping_total, "stopping longitudinal",self.stopping_longitudinal )#, "due to photons:",self.stopping_due_to_photons )
        
          
        
        

      
    def DDCS_at_omega(self):
        self.x_axis = np.zeros(self.NThetaStep)
        self.Result1 = np.zeros(self.NThetaStep)
        self.Result2 = np.zeros(self.NThetaStep)
        self.x_axis[0]= self.theta_max
        for i in range(self.NThetaStep-1):
            self.x_axis[i+1] = self.x_axis[i]*0.95
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
 
        epslib.DDCS_at_omega(self.ParArray,self.x_axis,self.Result1,self.Result2,
            self.DFChoice,self.omega_ddcs)
            
    def DDCS_at_theta(self):
        self.x_axis = np.zeros(self.NPoints)
        self.Result1 = np.zeros(self.NPoints)
        self.Result2 = np.zeros(self.NPoints)
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        
        epslib.DDCS_at_theta(self.ParArray,self.Result1,self.Result2,
            self.DFChoice,self.theta_ddcs)            
            
    def DCS(self):
        self.x_axis = np.zeros(self.NThetaStep)
        self.Result1 = np.zeros(self.NThetaStep)
        self.Result2 = np.zeros(self.NThetaStep)
        self.Result3 = np.zeros(self.NThetaStep)
        self.x_axis[0]= self.theta_max/1000.0
        for i in range(self.NThetaStep-1):
            self.x_axis[i+1] = self.x_axis[i]*0.975
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        self.sigma_from_DCS = epslib.DCS(self.ParArray,self.x_axis,self.Result1,self.Result2,self.DFChoice ) 
        print("self.sigma_from_DCS",self.sigma_from_DCS)
        E0_au= self.E0*1000.0/cnst.HARTREE 
        gamma_r = 1.0 + E0_au/(self.ProjectileMass*(cnst.C)**2)  
        for j in range(self.NThetaStep): # put Rutherford in result3
            current_theta=self.x_axis[j] # in rad
            q=2.0*np.sin(current_theta/2.0) * self.p0_rel  ## for angles >> larger that theta_0=(|P_1|- |P_0|)/|p_0|
            self.Result3[j] = (4 * self.Nelec_per_UC * gamma_r**2 * self.ProjectileMass**2 / q**4 )/ cnst.BOHR**2 / (4.0 * cnst.PI)  # not sure about the 4 pi   
        for j in range(self.NThetaStep): # we like our x-axis in mrad 
            self.x_axis[j] *= 1000
 
    def surfaceloss(self):
        self.x_axis = np.zeros(self.NPoints)
        self.Result1 = np.zeros(self.NPoints)
        self.Result2 = np.zeros(self.NPoints)
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize

        epslib.SurfLossFunc(
            self.ParArray,
            self.Result1,
            self.q,
            self.DFChoice,
        )

        self.SurfExProb = epslib.DSEP(
            self.ParArray, self.Result2, self.DFChoice
        )

    def PartDIIMFP(self):

        self.Result1 = np.zeros(self.NPoints)
        self.x_axis= np.linspace( self.CenterFirstBin,self.CenterLastBin, self.NPoints)  
        self.partialresults =np.zeros(shape=(self.NPoints,10))
        #self.partcross = np.zeros(shape=(10,self.NPoints))

        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)

        q_step = self.UpperqLimit / 10.0
        for i in range(10):
            # here we override what was assigned in initParArray as is used by
            # standard diimfp calc, will be put back to normal values in fill_remainder
            self.ParArray[self.NDFPAR + 10] = i * q_step
            self.ParArray[self.NDFPAR + 11] = (i + 1) * q_step
            epslib.DIIMFP(self.ParArray,self.Result1,self.DFChoice,self.StoppingResultArray)

            if self.weight == 0:
                self.PartIntSum[i] = self.StoppingResultArray[0]
                for j in range(self.NPoints):
                    self.partialresults[j, i] = self.Result1[j]
            else:
                self.PartIntSum[i] = self.StoppingResultArray[1]
                for j in range(self.NPoints):
                    self.partialresults[j, i] = (self.Result1[j] * self.x_axis[j])
   

    def Penn_from_ELF(self):
        constant_stepsize=True
        self.MyChapApp.runplot.initcalc()
        current_q=self.q
        self.q=0.0
        self.oneovereps1eps2()
        self.q=current_q
        N_oscillator_used = self.Last_oscillator_transform -self.First_oscillator_transform+1 
        print("first,last",self.First_oscillator_transform,self.Last_oscillator_transform)
        if N_oscillator_used < 1:
            print(" oscillator range input error")
            return
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        
        if constant_stepsize:
            incrementfactor=1
            CurrentStepsize=(self.UpperELimit-self.LowerELimit)/N_oscillator_used
            CurrentEnergy=self.LowerELimit+0.5*CurrentStepsize
            old_Ai=self.Amps[0]  #Hack, remove
            
        else:    
            rangefactor=(self.UpperELimit)/self.CenterFirstBin
            exponent=1.0/self.N_oscillator_used
            incrementfactor=rangefactor**exponent
            CurrentStepsize=self.CenterFirstBin*(incrementfactor-1)
            if  CurrentStepsize < 0.5*self.CenterFirstBin:
                CurrentEnergy = self.CenterFirstBin
                
            else:
                CurrentEnergy = 0.5* CurrentStepsize 
            rangefactor=(self.UpperELimit)/CurrentEnergy
            exponent=1.0/self.N_oscillator_used
            incrementfactor=rangefactor**exponent
            
        sumnewAi=0.0
        for i in range(self.First_oscillator_transform,self.Last_oscillator_transform+1):
            Loss_at_E = epslib.Lossfunction_inclGOS_atE(self.ParArray, CurrentEnergy, self.DFChoice )
            g_omega_prefactor = 2.0 / (cnst.PI * CurrentEnergy)  # Penn delta function  in the case  gamma=0
            # alternative from nguyen  for gamma !=0 (J of Phys Chem 119 23627 2015) does not work Sum A_i> 1
            # g_omega_prefactor=2.0/(cnst.PI*CurrentEnergy**2*CurrentGamma) * \
            # math.sqrt(2*CurrentEnergy*(CurrentEnergy**2+CurrentGamma**2)*(math.sqrt(CurrentEnergy**2+CurrentGamma**2)-CurrentEnergy))
           
            Amp = g_omega_prefactor * Loss_at_E * CurrentStepsize
            sumnewAi+=Amp
           
            self.Amps[i] =Amp
            self.Omegas[i] =CurrentEnergy
            self.Gammas[i] = 1.*CurrentStepsize  # used to be 1.2
            self.Alphas[i]= 1.0
            self.Us[i]= 0.0
                 
            CurrentStepsize *= incrementfactor
            CurrentEnergy = CurrentEnergy + CurrentStepsize
        self.DFmodel = 'Mermin'
        self.AddELF=1
        for i in range(self.maxGOS):
            self.ConcGOS[i] = 0.0
        for i in range(self.maxKaneko):
            self.N_Kaneko[i] = 0.0
        for i in range(self.maxBelkacem):
            self.Conc_Belkacem[i] = 0.0
        print("old ai", old_Ai, "sumnewAi ",sumnewAi )   
        for i in range(self.maxOscillators):  
            self.Amps[i]=self.Amps[i]* old_Ai/ sumnewAi  #hack, remove
        self.MyChapApp.update_all_tables()    
        
        
    def convert_DL_to_Kaneko(self): #convert a set of DL (or Mermin)  oscillators to a set of Kaneko DL oscillators wilth the same optical ELF
        if not (self.DFmodel == 'DL' or self.DFmodel == 'Mermin'):
            print("requires DL or Mermin Oscillators, not Drude-Lindhard Tauc or Vlasov oscillators")
            return

        Q_used =  self.Q_Kaneko_transform
        l_used = self.l_Kaneko_transform
        if l_used == 0:
            doublefact = 1
        elif l_used == 1:
            doublefact = 3
        elif l_used == 2:
            doublefact = 5 * 3
        else:
            doublefact = 7 * 5 * 3
        w_pl_l=math.sqrt(2*(2*l_used+1)*doublefact*Q_used**3*math.exp(l_used)/((2.0*l_used)**l_used*math.sqrt(cnst.PI)))
 
     

        for i in range(self.maxKaneko):
            if self.Amps[i] > 0.0:
                A=self.Amps[i]
                W = self.Omegas[i] / cnst.HARTREE
                U = self.Us[i] /cnst.HARTREE
                elec_dens = ( A*W*W / (4.0 * cnst.PI) / (cnst.BOHR**3) )
                elec_per_unit_cell = elec_dens / self.UnitCellDensity
                
                self.N_Kaneko[i]=elec_per_unit_cell
                self.Q_Kaneko[i]=Q_used
                self.Edge_Kaneko[i]=U*cnst.HARTREE
                self.l_Kaneko[i]=l_used
                self.width_Kaneko[i]=self.Gammas[i]
                self.gamma_Kaneko[i]=W**2 /w_pl_l**2
        
        for i in range(self.maxOscillators):
            self.Amps[i]=0.0
        self.MyChapApp.update_all_tables()        
        
    def convert_DF(self):
        if self.DFmodel == 'DL' or self.DFmodel == 'Mermin':
             self.convert_DL_to_Drude()
        elif self.DFmodel == 'Drude':
            self.convert_Drude_to_DL()
        self.MyChapApp.update_all_tables()           
                
    def convert_Drude_to_DL(self):   # from  extended Drude toDrude-Lindhard
        for i in range(self.maxOscillators):
            oldA = self.Amps[i]
            if oldA != 0.0:
                old_omega = self.Omegas[i]
                U = self.Us[i] # U remains unchanged
                old_position=math.sqrt(oldA+old_omega**2 + U**2)
                self.Omegas[i] = math.sqrt(old_position**2 - U**2) 
                self.Amps[i] = oldA/self.Omegas[i]**2 
        self.DFmodel = 'DL'
        
        
    def convert_DL_to_Drude(self):    # from Drude-Lindhard to extended-Drude
         for i in range(self.maxOscillators):  
            oldA = self.Amps[i]
            if oldA != 0.0:
                old_omega = self.Omegas[i]
                U = self.Us[i]
                self.Amps[i]=old_omega**2*oldA
                old_position=math.sqrt(old_omega**2 + U**2)
                square_value= old_position**2-self.Amps[i] - U**2
                if square_value >= 0.0:
                    self.Omegas[i] = math.sqrt( old_position**2-self.Amps[i] - U**2) 
                else:     #square_value sometimes negative due to rounding errors
                    self.Omegas[i] = 0.0  
         self.DFmodel = 'Drude'
    
  


    def Change_U(self):
        self.MyChapApp.UpdateStatus("")
        Nosc = self.Last_oscillator_transform -self.First_oscillator_transform+1 
        if Nosc < 0:
            print("osc range input error")
        if Nosc > self.maxOscillators:
            Nosc=self.maxOscillators
        if self.U_factor < 0.0:
            print("U factor should be positive")
            return
        if self.U_factor < 1 or(self.U_factor==1 and self.DFmodel == "drude"):
            for i in range(self.First_oscillator_transform,self.Last_oscillator_transform+1):  # now calculate new values for new U,
                
                oldA = self.Amps[i]
                if oldA != 0.0:
                    old_omega = self.Omegas[i]
                    old_U = self.Us[i]
                    if self.DFmodel == 'Drude':
                        old_position=math.sqrt(oldA+old_omega**2 + old_U**2)
                        a=old_position**2-oldA
                        if a<0:a=0
                        shift=math.sqrt(a)
                        self.Us[i]= shift*math.sqrt(self.U_factor)
                        a=1.0-self.U_factor
                        if a<0:a=0
                        self.Omegas[i]=shift*math.sqrt(a)
                    elif self.Omegas[i] > 0.0:
                        old_position=math.sqrt(old_omega**2 + old_U**2)
                        newU=math.sqrt(self.U_factor*old_position**2)
                        self.Omegas[i] = math.sqrt(old_position**2 - newU**2)
                        self.Amps[i] = oldA * (old_omega / self.Omegas[i])**2    
                        self.Us[i]=newU 
        else:
            print(" U factor too large")
            return                 
        self.MyChapApp.update_all_tables()              

    def Calc_Os_Strength(self): 
        self.oneovereps1eps2()  # after much confusion osc. strength based on omega times loss function , see Nikjoo book eq. 19.36, not omega times eps2
        for i in range(self.NPoints):
            self.Result1[i] = (2.0 * self.x_axis[i] / (cnst.PI * self.PlasmonE**2) * self.Result2[i])  # oscillator strength based on 1/eps2, put in result1    
            
    def dyn_struct_factor(self):
        FineMeshFactor=9
        length_fine = self.NPoints*FineMeshFactor
        self.calculate_energy_axis(FineMeshFactor)
        eps_fine=self.calculate_eps_array() 
        one_over_eps=1.0/eps_fine
        self.Result1_fine = - np.imag(one_over_eps)  #this is Im[-1/eps],the loss function, NOT im[1/eps]
        self.Result2_fine = np.zeros(length_fine)
        elecdens = self.Nelec_per_UC *self.UnitCellDensity * cnst.BOHR**3  # elec/per a.u.^3
        omega_p_square=4.0 * cnst.PI * elecdens #Hartree square
        self.Result1_fine *= self.q**2 / ( cnst.PI* omega_p_square) /cnst.HARTREE
        self.Result2_fine = np.cumsum(self.Result1_fine*self.x_axis*self.xstepsize)
        self.recast( FineMeshFactor,only_2=True)    
      
    
    def DL_from_OOS(self):
        Nosc = self.Last_oscillator_transform -self.First_oscillator_transform+1 
        if Nosc < 0:
            print("osc range input error")
        if Nosc > self.maxOscillators:
            Nosc=self.maxOscillators
        step = self.Stepsize
        CenterFirstBin = self.LowerELimit + 0.5*step # first energy point of plot
        CenterLastBin  = self.UpperELimit - 0.5*step # last energy point assuming constant stepsiize
        linearLastE=self.LowerELimit+Nosc*step 
        x_axis= np.linspace(CenterFirstBin, linearLastE, Nosc) 
        missing_bit=self.UpperELimit - linearLastE
        if missing_bit > 0.0:
            logXaxis=np.geomspace( step/5,missing_bit ,Nosc)
            last_step=logXaxis[Nosc-1]-logXaxis[Nosc-2]
            logXaxis=np.geomspace(self.Stepsize/20,missing_bit-0.5*last_step,Nosc)
            x_axis=x_axis+ logXaxis
        xstepsize=np.gradient(x_axis)
        omega=x_axis[0]
        currentBin=0
        for i in range(self.First_oscillator_transform,self.Last_oscillator_transform):
            omega = x_axis[i-self.First_oscillator_transform]
            while omega >  self.OOSEnergy[currentBin]:
                currentBin +=1
            PosWithinBin =  (omega - self.OOSEnergy[currentBin-1])/ (self.OOSEnergy[currentBin] - self.OOSEnergy[currentBin-1]) 
            currentOOS = self.OOS[currentBin-1] + PosWithinBin*(self.OOS[currentBin]-self.OOS[currentBin-1])
            currentELF=currentOOS/(2.0*omega)*cnst.PI *  self.PlasmonE**2   
            g_omega_prefactor = 2.0 / (cnst.PI * omega) 
            Amp = g_omega_prefactor *currentELF * xstepsize[i-self.First_oscillator_transform]
            self.Omegas[i]=omega
            self.Amps[i]=Amp
            self.Gammas[i]=1.15*xstepsize[i-self.First_oscillator_transform]
            self.Alphas[i]=1.0
            self.Us[i]=0
             
    def PseudoChargeDensity(self):
        """calculate the distribution of the charge density within the Penn Pseudo charge picture"""
        copyq=self.q
        self.q=0.05  # make sure elf is calculated near the optical limit
        self.oneovereps1eps2()  # calculartes also re[1/eps] in result1, but later overwritten here
        self.q=copyq   # put back value defined by interface (prob. not required)
        Energy_step_au = self.Stepsize / cnst.HARTREE  # energy step in a.u.
        vol_so_far = 0.0
        charge_so_far=0.0
        i_overfull = 0.0
        for i in range(self.NPoints):

            # now energy axis in atomic units
            self.x_axis[i] = self.x_axis[i] / cnst.HARTREE
        for i in range(self.NPoints):
            j = self.NPoints - i - 1
            # G(omega from Penn eq. 9b) within the statistical interpretation
            # this is also proportional to the fraction of space filled with
            # this electron density
            GPenn = self.Result2[j] * 2 / (cnst.PI * self.x_axis[j])
            # e density of this electron loss (a.u.)^-3
            elec_dens = self.x_axis[j] ** 2 / (4 * cnst.PI)
            # electron density per angstrom^3
            self.Result1[j] = elec_dens / cnst.BOHR**3
            vol_so_far += GPenn * Energy_step_au
            charge_so_far+=elec_dens* GPenn * Energy_step_au/self.UnitCellDensity / cnst.BOHR**3
            self.x_axis[j] = vol_so_far
            if self.x_axis[j] > 1.0:
                i_overfull += 1

                self.x_axis[j] = 1.0

        self.x_axis[0] = 1.0
        self.Result1[0] = 0.0

    def ConvertToRadialPseudoChargeDensity(self):
        # should be called straight after  calc.PseudoChargeDensity(), only makes sense for pure elements
        self.r_array = np.zeros(self.NPoints)
        self.dens_array = np.zeros(self.NPoints)
        UC_Volume = 1 / self.UnitCellDensity
        for i in range(self.NPoints):
            current_volume = self.x_axis[i] * UC_Volume
            radius = (current_volume * 3.0 / (4.0 * cnst.PI)) ** (1.0 / 3.0)
            self.r_array[i] = radius
            self.dens_array[i] = self.Result1[i]
            self.x_axis[i] = self.r_array[i]  # so we can plot it by itself
            
            
    def stopping_IMFP_w_p_versus_r(self):
        self.Result3 =np.zeros(self.NPoints) # for stoppings result1 for w_p, result2 for lambda
        self.Result4 = np.zeros(self.NPoints) # for straggling result1 for w_p, result2 for lambda
        if self.particle == "proton":
            ConstantA = 1.0/2.0  # notation so we are in line with imfp draft eq. 11
            ConstantB = 1
        else:
            ConstantA = 1
            ConstantB = 2
        velocity = self.velocity_projectile(self.E0)  # in a.u.   
        
        for i in range(self.NPoints):
            current_density = self.Result1[i]
            if(current_density > 0.0):
                PlasmonE= np.sqrt(4.0 * cnst.PI * current_density * cnst.BOHR**3)
                if velocity**2 > 2 * ConstantB * PlasmonE:
                    q_min_DL = (velocity - np.sqrt(velocity**2 - 2 * ConstantB * PlasmonE)
                    ) / ConstantB
                    q_max_DL = (velocity + np.sqrt(velocity**2 - 2 * ConstantB * PlasmonE)
                    ) / ConstantB
                else:
                    q_min_DL = 1
                    q_max_DL = 1
                logDL = np.log(q_max_DL / q_min_DL)    
                tmp = logDL - 0.5 * math.log(
                    (2 * PlasmonE + q_max_DL**2) / (2 * PlasmonE + q_min_DL**2))
                oneoverlambda = PlasmonE * tmp / velocity**2    
                if oneoverlambda > 0.0:
                    Lambda = (1.0 / oneoverlambda)  # in a.u.
                else:
                    Lambda = 0.0 
                stopping =   PlasmonE**2 / velocity**2 * logDL
                tmp =(PlasmonE**3 / velocity**2 * logDL) 
                straggling = 0.25 * PlasmonE**2 / velocity**2 * (q_max_DL**2 - q_min_DL**2) + tmp 
            else:
                Lambda = 0.0 
                PlasmonE = 0.0 
                stopping=0.0
                straggling = 0.0      
            
            self.Result1[i] = PlasmonE * cnst.HARTREE
            self.Result2[i] = Lambda  * cnst.BOHR
            self.Result3[i] = stopping  * cnst.HARTREE /cnst.BOHR
            self.Result4[i] = straggling  * cnst.HARTREE * cnst.HARTREE / cnst.BOHR
   

    def stopping_versus_r(self):
        """IMFPversus r, to be called straight after  ConvertToRadialPseudoChargeDensity()
        atom in cube with volume Muffin_Tin sphere, i.e. unit cell volume"""
        Cube_length = (1.0 / self.UnitCellDensity) ** (1.0 / 3.0)
        self.MT_radius = (1.0 / self.UnitCellDensity * 3.0 / (4.0 * cnst.PI)) ** (
            1.0 / 3.0
        )
        print("MT radius", self.MT_radius, "side cube", Cube_length)
        # prepare the array for stopping power calculation with one Mermin

        # first store the calculated density versus r values safely
        currentEnergy = self.CenterFirstBin
        currentstepsize = self.Stepsize
        pointswithfixedstepsize = int(
            (self.UpperELimit - self.CenterFirstBin) / self.Stepsize
        )

        for ipoints in range( pointswithfixedstepsize):  # this should match what happens in epslib
  
            if pointswithfixedstepsize > self.max_no_lin_cont_deltaE:
                self.fixedstepsize = False
                self.ParArray[self.NDFPAR + 9] = self.lin_cont_deltaE
                currentstepsize = self.Stepsize + self.lin_cont_deltaE * currentEnergy

            else:
                self.fixedstepsize = True
                self.ParArray[self.NDFPAR + 9] = 0.0

            currentEnergy += currentstepsize
            if currentEnergy > self.UpperELimit:
                break

        diimfparray = np.zeros(self.NPoints)  # calculated but not really used

        self.stop_array = np.zeros(self.NPoints)

        for i in range(self.maxOscillators):  # make sure other oscillators in c routine have amplitude 0 
            self.ParArray[4 * i + 1] = 0.0
        for i in range(self.maxGOS):  # make sure all GOS amp are 0
            self.ParArray[4 * i + 4 * self.maxOscillators + 1] = 0.0
        for i in range(self.maxKaneko):  # make sure all GOS amp are 0
            self.ParArray[
                4 * i + 4 * self.maxOscillators + 4 * self.maxKaneko + 1
            ] = 0.0

        # make sure remainder of pararray is filled with correct data
        self.fill_remainder()

        self.ParArray[self.NDFPAR + 5] = (
            ipoints + 1
        )  # put to NPoints by fill remainder but should be ipoints+1 (variable step size)
        for i in range(self.NPoints):
            elec_dens_au = self.dens_array[i] * cnst.BOHR**3
            w_pl = cnst.HARTREE * math.sqrt(4 * cnst.PI * elec_dens_au)  # in eV
            # at distance r_array[i] we have stopping due to DF with plasmon energy w_pl
            self.ParArray[1] = 1.0
            self.ParArray[2] = w_pl
            self.ParArray[3] = (w_pl / 10)  # this is gamma, should not matter may afffect comp. time required
            epslib.DIIMFP(self.ParArray, diimfparray,self.DFChoice,self.StoppingResultArray)

            self.stop_array[i] = self.StoppingResultArray[1]

    # def impact_dep_stop(self): #currently not connected to user interface, newver called
        # """impact dependent stopping, to be called straight after  ConvertToRadialPseudoChargeDensity()
        # atom in cube with volume Muffin_Tin sphere, i.e. unit cell volume.  It is assumed that the atoms
        # are in a simple cubic lattice with nearest neighbor distance Cube_length
        # Maximum impact parameter considered is Cube_length/2. Particle impinges perpendicular to cube
        # The energy loss for crossing the cube is calculated as a function of the impact parameter.
        # Part of the cube is outside the Muffin Tin sphere.  For that part the electron density at
        # the edge of the MT shpere is taken.
        # This edge density again is determined by the lowest energy considered for the ELF. Take this not too low"""
        # # calculate charge density distribution (Penn's pseudo charge density")
        # self.PseudoChargeDensity()
        # # put it in a radial form
        # self.ConvertToRadialPseudoChargeDensity()
        # # calculate sthe local stopping for that density (this takes time)
        # self.stopping_versus_r()
        # # calculate inpact parameter dependent stopping

        # Cube_length = (1.0 / self.UnitCellDensity) ** (1.0 / 3.0)
   
        # self.Nst = 500
        # sumB2 = 0
        # sumstop = 0

        # self.b_dep_stop = np.zeros( self.Nst)
        # self.b_param = np.zeros( self.Nst)
        # steplength = Cube_length / (2 * self.Nst)
        # for i in range(self.Nst):
            # self.b_dep_stop[i] = 0.0
            # self.b_param[i] = (
                # i + 0.5
            # ) * steplength  # so start at b=Cube_length/(2*Nst) up to almost Cube_length/2
            # sumB2 += self.b_param[i] * self.b_param[i]
            # for j in range(self.Nst):
                # along = (
                    # j + 0.5
                # ) * steplength  # so from almost  central atom up to cube edge

                # r = math.sqrt(along * along + self.b_param[i] * self.b_param[i])
                # if r > self.MT_radius:
                    # r = self.MT_radius - (
                        # r - self.MT_radius
                    # )  # folding back when outside mt radius
                # for current_index in range(self.NPoints):
                    # if r > self.r_array[current_index]:
                        # break  # so if we are never outside MT sphere
                # current_index = (
                    # current_index - 1
                # )  # this means if part crystal is empty this part get stopping of zero not the next value
                # self.b_dep_stop[i] += (
                    # 2 * self.stop_array[current_index] * steplength / Cube_length
                # )  # the factor of 2 because we integrate only from edge to center
                # # not to the other edge
            # sumstop += self.b_param[i] * self.b_param[i] * self.b_dep_stop[i]
  
        # print("estimate of average stopping:", sumstop / sumB2)
        # # now we are going to calculate the contribution of bunching to straggling
        # self.straggling_from_bunching = 0.0
        # self.straggling_from_bunching_rel_vacuum = 0.0
        # for i in range(self.Nst):
            # self.straggling_from_bunching += (
                # self.b_dep_stop[i] - sumstop / sumB2
            # ) ** 2 * self.b_param[i] ** 2
            # self.straggling_from_bunching_rel_vacuum += (
                # self.b_dep_stop[i]
            # ) ** 2 * self.b_param[i] ** 2
        # self.straggling_from_bunching = self.straggling_from_bunching / sumB2
        # self.straggling_from_bunching_rel_vacuum = (
            # self.straggling_from_bunching_rel_vacuum / sumB2
        # )
        # print(self.straggling_from_bunching, "bunching per unit cell")
        # self.straggling_from_bunching = self.straggling_from_bunching / Cube_length
        # self.straggling_from_bunching_rel_vacuum = (
            # self.straggling_from_bunching_rel_vacuum / Cube_length
        # )
        # print(self.straggling_from_bunching, "bunching angstrom")
        # print(
            # self.straggling_from_bunching / self.BohrStraggling,
            # "bunching rel to Bohr straggling",
        # )
        # print(
            # self.straggling_from_bunching_rel_vacuum
            # / Cube_length
            # / self.BohrStraggling,
            # "bunching rel to Bohr straggling_rel_vacuum",
        # )
        # print("self.BohrStraggling", self.BohrStraggling)

    # def bunching_versus_E0(self):
        # # calculate charge density distribution (Penn's pseudo charge density")
        # self.PseudoChargeDensity()
        # # put it in a radial form
        # self.ConvertToRadialPseudoChargeDensity()
        # oldE0 = self.E0
       
        # if self.particle  == "proton":
            # self.E0 = 2
        # else:
            # self.E0 = self.first_electron_energy

        # StragglingResult1 = np.zeros(self.NStopping)
        # StragglingResult2 = np.zeros(self.NStopping)
        # Beam_incr = 1.25
        # NBunch = 35
        # for Ecounter in range(NBunch):
            # print("currently calculating No", Ecounter, "Energy", self.E0)
            # self.ParArray[self.NDFPAR + 1] = self.E0 * 1000
            # self.impact_dep_stop()
            # StragglingResult1[Ecounter] = (
                # self.straggling_from_bunching / self.BohrStraggling
            # )
            # StragglingResult2[Ecounter] = (
                # self.straggling_from_bunching_rel_vacuum / self.BohrStraggling
            # )
            # self.E0 = self.E0 * Beam_incr

        # self.E0 = oldE0
        # self.x_axis  = np.zeros(NBunch)
        # self.Result1 = np.zeros(NBunch)
        # self.Result2 = np.zeros(NBunch)
        # E = self.first_proton_energy
        # self.NPoints = NBunch
        # for i in range(NBunch):
            # self.x_axis[i] = E
            # self.Result1[i] = StragglingResult1[i]
            # self.Result2[i] = StragglingResult2[i]
            # E = E * Beam_incr

    def REELS_spectrum(self):
        Start_REELS = -5.0  # hard coded to start at -5
        NREELS = (int((self.UpperELimit - Start_REELS) / self.Stepsize) + 1)  # hard coded to start at -5
        NDIIMFP = int(self.UpperELimit / self.Stepsize) + 1
        self.x_axis = np.zeros(NREELS)
        self.Result1 = np.zeros(NREELS)
        self.NormDIIMFP = np.zeros(NDIIMFP)
        epslib.eps_Scaling_init(self.ParArray, self.DFChoice)
        for Ecounter in range(NREELS):
            self.x_axis[Ecounter] = -5.0 + Ecounter * self.Stepsize
            self.Result1[Ecounter] = 0.0
        epslib.calc_REELS(self.ParArray,Start_REELS, NREELS,
            self.DFChoice, self.Result1, self.EELS, self.EELS_thickness)
    def calc_Qrecoil(self, q):  #calculates (relativistically ) the kinetic energy of an initially stationary electron after it's momentum changes to q
        #see appendix J.M. Fernández-Varea et al. / Nucl. Instr. and Meth. in Phys. Res. B 229 (2005) 187–218
        # #high q approach
        # Q= math.sqrt(cnst.C**2*q**2+cnst.C**4)-cnst.C**2
        # #low q approach
        # x=0.5*(cnst.C*q/cnst.C**2)**2
        # Q_low=cnst.C**2*(x-x**2/2+x**3/2)
        # print ("q=",q,"Q high q approach",Q, "Q, low q approach", Q_low, "percentage", 100* (Q_low-Q)/Q, "q^2/2",q*q/2)
        #conclusion of above test: use q^2/2 up to 3 a.u. high-Q approach above that
        if q < 3.0:
            return q**2/2
        else:    
            return math.sqrt(cnst.C**2*q**2+cnst.C**4)-cnst.C**2
    
        
    def get_df_properties(self):  
    
        self.SumAi = 0.0
        SumGOS = 0.0
        SumKaneko = 0.0
        VolumeFractionKaneko=0.0
        SumBelkacem = 0.0
        elec_dens = 0.0
        if self.massunitcell <= 0.0:
            return
        mole_per_cm3 = self.specificweight / self.massunitcell
       
        self.UnitCellDensity = mole_per_cm3 * cnst.NAvogadro / 1.0e24
        self.UnitCellDensityText = "unit cell  (u.c.) density:\n  " + f"{self.UnitCellDensity:.4f}"+ " per Å³"
        
        self.PlasmonE= ( np.sqrt(4.0 * cnst.PI * self.UnitCellDensity * cnst.BOHR**3) * cnst.HARTREE)
        self.w_p_1e_per_uc = "ωₚ  for  1 e⁻ per u.c.: " + f"{self.PlasmonE:.2f}" + " eV"

        self.MT_radius = (1.0 /self.UnitCellDensity * 3.0 / (4.0 * cnst.PI)) ** (1.0 / 3.0)
        try:
            for i in range(self.maxGOS):
                SumGOS += float(self.ConcGOS[i])
            for i in range(self.maxBelkacem):
                SumBelkacem += float(self.Conc_Belkacem[i])
            for i in range(self.maxKaneko):
                
                if self.N_Kaneko[i] > 0.0:
                    currentN = self.N_Kaneko[i]
                    currentQ = self.Q_Kaneko[i]
                    currentG = self.gamma_Kaneko[i]
                    currentl = self.l_Kaneko[i]
                    SumKaneko += currentN
                    if self.Kaneko_choice == 0:   #modified Kaneko
                        if currentl==0:
                            doublefact=1
                        elif currentl==1: 
                            doublefact=3
                        elif currentl==2:
                            doublefact=15
                        elif (currentl==3):
                            doublefact=105

                        w_p=math.sqrt(currentG*2*(2*currentl+1)*doublefact*currentQ**3 \
                            *math.exp(currentl)/((2.0*currentl)**currentl*math.sqrt(math.pi)) )
                        w_p=w_p*cnst.HARTREE   
                        VolumeFractionKaneko += currentN*self.PlasmonE*self.PlasmonE/(w_p*w_p)  
                    else:   # original Kaneko
                        qmean= currentQ*currentN**(1.0/3.0) 
                       # tmp=currentG*qmean**3/math.sqrt(math.pi)
                        w_p = math.sqrt(currentG*qmean**3 / math.sqrt(math.pi) ) * cnst.HARTREE
                        VolumeFractionKaneko += currentN*self.PlasmonE*self.PlasmonE/(w_p*w_p)   
                         
            if self.DFmodel == "Drude":
                for i in range(self.maxOscillators):
                    elec_dens +=  self.Amps[i] / cnst.HARTREE**2 / (cnst.BOHR**3) / (4.0 * cnst.PI)
                elec_per_unit_cell = elec_dens / self.UnitCellDensity
            elif self.DFmodel == "DL" or self.DFmodel == "Mermin" or  self.DFmodel == "vlasov" :
                for i in range(self.maxOscillators):
                    if abs(self.Amps[i]) > 1e-50:
                        if self.DFmodel == "vlasov":
                            vlasov_Q=self.Omegas[i]  # for Vlasov the omega variable contains the Q value. Bit confusing, sorry
                            W = math.sqrt(vlasov_Q**3 / math.sqrt(math.pi) ) 
                        else:    
                            W = self.Omegas[i] / cnst.HARTREE
                            U = self.Us[i] /cnst.HARTREE
                        
                        if W > 0.0:
                            W_M_square = W**2 + U**2
                            A_M = W*W / W_M_square *self.Amps[i]
                        else:
                            A_M=0.0  
                            W_M_square=0.0
                        self.SumAi += A_M
                        
                        elec_dens += ( A_M * W_M_square / (4.0 * cnst.PI) / (cnst.BOHR**3) )
            elif self.DFmodel == "Tauc":
                for i in range(self.maxOscillators):
                    A = self.Amps[i]/cnst.HARTREE
                    if abs(A) > 1e-50:
                        E0 = self.Omegas[i]/ cnst.HARTREE
                        Eg = self.Us[i] /cnst.HARTREE
                        C = self.Gammas[i] / cnst.HARTREE
                        Norm_rel_to_Drude =epslib.TaucSumRule(A, C,E0,Eg)
                        TaucDens =Norm_rel_to_Drude* self.Amps[i]*self.Omegas[i] / cnst.HARTREE**2 / (cnst.BOHR**3) / (4.0 * cnst.PI)
                        elec_dens += TaucDens
                        print("elec dens from epslib:", elec_dens,"Drude limit" , self.Amps[i]*self.Omegas[i] / cnst.HARTREE**2 / (cnst.BOHR**3) / (4.0 * cnst.PI))
            elif self.DFmodel == "TL_an":
                a_TL_an_au=self.a_TL_an/cnst.HARTREE
                for i in range(self.maxOscillators):
                    A = self.Amps[i]/cnst.HARTREE
                    if abs(A) > 1e-50:
                        E0 = self.Omegas[i] / cnst.HARTREE
                        Eg = self.Us[i] /cnst.HARTREE
                        C  = self.Gammas[i] / cnst.HARTREE
                        Norm_rel_to_Drude = epslib.TL_an_SumRule(A,C,E0,Eg,a_TL_an_au)
                        print("norm rel to drude",Norm_rel_to_Drude)
                        TL_an_Dens =Norm_rel_to_Drude* self.Amps[i]*self.Omegas[i] / cnst.HARTREE**2 / (cnst.BOHR**3) / (4.0 * cnst.PI)
                        elec_dens += TL_an_Dens  
            elif self.DFmodel == "TL_Mermin":
                for i in range(self.maxOscillators):
                    A = self.Amps[i]/cnst.HARTREE
                    if abs(A) > 1e-50:
                        E0 = self.Omegas[i] / cnst.HARTREE
                        Eg = self.Us[i] /cnst.HARTREE
                        C  = self.Gammas[i] / cnst.HARTREE
                        Norm_rel_to_Drude = epslib.TaucMermin_SumRule(A, C,E0, Eg)
                        print("norm rel to drude",Norm_rel_to_Drude, "E0",E0,"Eg",Eg,"C",C,"A",A)
                        TL_Mermin_Dens =Norm_rel_to_Drude* self.Amps[i]*self.Omegas[i] / cnst.HARTREE**2 / (cnst.BOHR**3) / (4.0 * cnst.PI)
                        elec_dens += TL_Mermin_Dens                   
        except  ValueError: 
            return
        except Exception as err:
            print(f"Unexpected {err=}, {type(err)=}")
            raise     
        elec_per_unit_cell = elec_dens / self.UnitCellDensity
        self.DF_prop_text= "Oscillators:  "+  f"{elec_dens:.2f}"  + " e⁻/Å³ or\n    "  + f"{elec_per_unit_cell:.2f}" + " e⁻ per U.C."
        if self.DFmodel == "DL" or  self.DFmodel == "Mermin":
            self.SumAiText  ="    (\u03A3 c(U)Aₓ ="  + f"{self.SumAi:.3f})" 
        elif self.DFmodel == "Vlasov":
            Q=self.Omegas[0]  # omega array contains Q in the Vlasov case
            KBoltzmann=8.617e-5/cnst.HARTREE 
            temp=Q*Q/(2*KBoltzmann)
            self.SumAiText  ="Q  = "  + f"{self.Omegas[0]:.2f}" + " a.u. equiv. to "+ f"{temp:.1e}" + "K" 
                
        else:  
            self.SumAiText=""        
        density = SumGOS * self.UnitCellDensity
        self.GOStext="GOS:            "+ f"{ density:.3f}" + " e⁻/Å³ or\n     " + f"{SumGOS:.2f}"  + " e⁻ per U.C."
       
        density = SumBelkacem *self.UnitCellDensity
        self.Belkacemtext = "BelKacem: " + f"{ density:.3f}" + " e⁻/Å³ or \n    " + f"{SumBelkacem:.2f}"+ " e⁻ per U.C."
        density = SumKaneko * self.UnitCellDensity
        self.Kanekotext="Kaneko:     " + f"{ density:.3f}"+ " e⁻/Å³ or \n    " + f"{SumKaneko:.2f}" + " e⁻ per U.C.\n    " \
            + " volume fraction: " + f"{VolumeFractionKaneko:.3f}"
   
        self.Nelec_per_UC = elec_dens / self.UnitCellDensity + (SumGOS + SumKaneko + SumBelkacem) 
        self.BohrStraggling = (4 * cnst.PI * self.Nelec_per_UC * self.UnitCellDensity * cnst.BOHR**3 * cnst.HARTREE**2 / cnst.BOHR)

        self.myconversiontext1 = ("from eV/ Å to MeV/(mg/cm²):           " + f"{0.1/self.specificweight:.3f}")
        self.myconversiontext2 = ("from eV/ Å to eV/(10¹⁵ molecules/cm²): "    + f"{self.massunitcell/(self.specificweight*6.022):.3f}")
    
