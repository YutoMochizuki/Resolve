import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
import numpy as np
import scipy as sp
import scipy.optimize as so
import scipy.special
import sys



def lineDB1():
    lineDBname='Cr_Ka'
    filename="./CalLines/Cr_Ka7lor.dat"
    data=np.loadtxt(filename, comments='#')
    energyDB=data[:,0]
    lgammaDB=data[:,1]
    ampDB=data[:,2]
    return lineDBname, energyDB, lgammaDB, ampDB


def lineDB2():
    lineDBname='V_Kb'
    filename="./CalLines/V_kb3lor_Smale.dat"
    data=np.loadtxt(filename, comments='#')
    energyDB=data[:,0]
    lgammaDB=data[:,1]
    ampDB=data[:,2]
    return lineDBname, energyDB, lgammaDB, ampDB


def calcchi(params,consts,model_func,xvalues,yvalues,yerrors,species):
    model = model_func(species,xvalues,params,consts)
    chi = (yvalues - model) / (yerrors)
    return(chi)


# optimizer
def solve_leastsq(xvalues,yvalues,yerrors,param_init,consts,model_func,species):
    param_output = so.leastsq(
        calcchi,
        param_init,
        args=(consts, model_func, xvalues, yvalues, yerrors, species),
        full_output=True)
    param_result, covar_output, info, mesg, ier = param_output
    if covar_output is not None:
      error_result = np.sqrt(covar_output.diagonal())
    else :
      error_result = None
    dof = len(xvalues) - 1 - len(param_init)
    chi2 = np.sum(np.power(calcchi(param_result,consts,model_func,xvalues,yvalues,yerrors,species),2.0))
    return([param_result, error_result, chi2, dof])


def mymodel(species, x, params, consts, tailonly = False):
    #norm,gw,gain,P_tailfrac,P_tailtau,bkg1,bkg2 = params    
    norm,gw,gain,bkg1, factor = params    
    # norm : nomarlizaion
    # gw : sigma of gaussian
    # gain : gain of the spectrum 
    # P_tailfrac : fraction of tail 
    # P_tailtau : width of the low energy tail
    # bkg1 : constant of background
    # bkg2 : linearity of background    
    #initparams = [norm,gw,gain,bkg1,bkg2]
    initparams = [norm,gw,gain,bkg1, factor]
    def rawfunc(x): # local function, updated when mymodel is called 
        return linemodel(species, x, initparams, consts=consts)               
    #model_y = smear(rawfunc, x, P_tailfrac, P_tailtau, tailonly=tailonly)
    #return model_y
    return rawfunc(x)


def linemodel(species,xval,params,consts=[]):
    # norm,gw,gain,bkg1,bkg2 = params
    norm,gw,gain, bkg1, factor = params
    # norm : normalization 
    # gw : sigma of the gaussian 
    # gain : if gain changes
    # consttant facter if needed 
    species_name, energy, lgamma, amp = lineDB1()
    species_name2, energy2, lgamma2, amp2 = lineDB2()
    prob = (amp * lgamma) / np.sum(amp * lgamma) # probabilites for each lines. 
    prob2 = (amp2 * lgamma2) / np.sum(amp2 * lgamma2) # probabilites for each lines. 
    model_y = 0 
    if len(consts) == 0:
        consts = np.ones(len(energy))
        consts2 = np.ones(len(energy2))
    else:
        consts = consts
        consts2 = np.ones(len(energy2))
    for i, (ene,lg,pr,con) in enumerate(zip(energy,lgamma,prob,consts)):
        voi = voigt(xval,[ene*gain,lg*0.5,gw])
        model_y += norm * con * pr * voi
    for i, (ene,lg,pr,con) in enumerate(zip(energy2,lgamma2,prob2,consts2)):
        voi = voigt(xval,[ene*gain,lg*0.5,gw])
        model_y += norm * con * pr * voi * factor
    # background = bkg1 * np.ones(len(xval)) + (xval - np.mean(xval)) * bkg2
    #model_y = model_y #+ background
    background = bkg1 * np.ones(len(xval)) 
    model_y = model_y + background
    # print "bkg1,bkg2 = ", bkg1,bkg2, background
    return model_y


def voigt(xval,params):
    center,lw,gw = params
    # center : center of Lorentzian line
    # lw : HWFM of Lorentzian (half-width at half-maximum (HWHM))
    # gw : sigma of the gaussian 
    z = (xval - center + 1j*lw)/(gw * np.sqrt(2.0))
    w = scipy.special.wofz(z)
    model_y = (w.real)/(gw * np.sqrt(2.0*np.pi))
    return model_y



gfwhm = 6
gw = gfwhm / 2.35
norm = 1000.0
gain = 1.0
bkg1 = 1.0
#bkg2 = 0.0
#P_tailfrac = 0.25
#P_tailtau = 10
factor = 0.5
#init_params=[norm,gw,gain,P_tailfrac,P_tailtau,bkg1,bkg2]
init_params=[norm,gw,gain,bkg1, factor]
# consts = [1,1,1,1,1]
consts = []


species_use=0
try :
      species_name, energy, lgamma, amp = lineDB1()
      species_name2, energy2, lgamma2, amp2 = lineDB2()
      energymin=np.min(energy)
      energymax=np.max(energy2)
      energyrange_min=energymin-20
      energyrange_max=energymax+20
      pixel='all'
      file='spec/spec_%s.txt' % pixel
      specdata=np.loadtxt(file)
      x_raw=specdata[:,0]
      y_raw=specdata[:,1]
      x0=x_raw[(energyrange_min<x_raw) & (x_raw<energyrange_max)]
      y0=y_raw[(energyrange_min<x_raw) & (x_raw<energyrange_max)]
      yerr_raw = np.sqrt(y0)
      x = x0[y0!=0]
      y = y0[y0!=0]
      yerr = yerr_raw[y0!=0]
      model_y = mymodel(species_use,x,init_params,consts)
      
      # do fit 
      result, error, chi2, dof = solve_leastsq(x, y, yerr, init_params, consts, mymodel,species_use)
      
      # get results 
      norm = np.abs(result[0])
      gw = np.abs(result[1])
      gain = np.abs(result[2])
      #tailfrac = np.abs(result[3])
      #tailtau = np.abs(result[4])
      #bkg1 = np.abs(result[5])
      #bkg2 = np.abs(result[6])
      
      if error is not None :
        norme = np.abs(error[0])
        gwe = np.abs(error[1])
        gaine = np.abs(error[2])
      else :
        norme = 0
        gwe = 0
        gaine = 0
      
      fwhm = 2.35 * gw
      fwhme = 2.35 * gwe
      label0 = "%s (+ %s), PIXEL=%s" % (species_name,species_name2, pixel)
      label1 = " (E_obs/Emodel) -1 = " + str("%.2E(+/-%.2E)" % (gain-1,gaine)) + ", dE = " + str("%4.2f(+/-%4.2f)" % (fwhm,fwhme) + " eV (FWHM)")
      label3 = "chi/dof = " + str("%4.2f"%chi2) + "/" + str(dof) + " = " + str("%4.2f"%  (chi2/dof))
      
      fitmodel = mymodel(species_use,x,result,consts)
      plt.figure(figsize=(10,8))
      ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
      plt.title(label0 + "\n" + label1 + "\n"  + label3 )
      #plt.xlabel("EPI (eV)")
      plt.errorbar(x0, y0, yerr=yerr_raw, fmt='ko', label = "data")
      plt.plot(x, fitmodel, 'r-', label = "model")
      plt.ylabel('Counts / 0.5 eV bin')
      background = bkg1 * np.ones(len(x)) #+ (x - np.mean(x)) * bkg2
      plt.plot(x, background, 'b-', label = "background", alpha = 0.9, lw = 1)
      eye = np.eye(len(energy))
      for i, oneeye in enumerate(eye):
        plt.plot(x, mymodel(species_use,x,result,consts=oneeye), alpha = 0.7, lw = 1, linestyle="--", label = str(i+1))
      
      plt.grid(linestyle='dotted',alpha=0.5)
      plt.legend(numpoints=1, frameon=False, loc="upper left")
      ax2 = plt.subplot2grid((3,1), (2,0))    
      plt.xscale('linear')
      plt.yscale('linear')
      plt.xlabel(r'EPI (eV)')
      plt.ylabel(r'Resisual')
      resi = fitmodel - y 
      plt.errorbar(x, resi, yerr = yerr, fmt='ko')
      plt.legend(numpoints=1, frameon=False, loc="best")
      plt.grid(linestyle='dotted',alpha=0.5)    
      plt.savefig("../figures/fit4_line%s_pix%s.png" % (species_use,pixel))
      plt.clf()
      plt.close()
      #plt.show()
      f = open('linefitting_result_all.txt','a')
      f.write('%f %f %f %E %E \n' % (energy[0], fwhm, fwhme, gain-1, gaine))
      f.close()
except :
      print('Error: #' + str(species_use) + ' (line=' + species_name + ')')
